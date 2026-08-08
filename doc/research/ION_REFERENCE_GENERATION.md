# Generating the full charge-aware ion-α reference set (cations → +4, anions → −2/−3)

Companion to `METHODS_EXPLAINED.md` and `../../data/targets_ionic_solids/METHOD_LIMITS_AND_FIXES.md`.
Goal: a comprehensive reference set so XDM's charge-aware α is defined for every
main/normal oxidation state — **cations out to +4 (stern route), and as many
anions as bind (−1 stern, −2/−3 via the density/compute route).**

All tooling lives on **dev-srv** (`ssh albd@10.10.49.104`) in
`~/critic2_dev/tools/wfc_generator/`. ld1.x = `/usr/bin/ld1.x` (apt QE 6.7).
The polarizabilities enter critic2 as **ratios** α(Z,q)/α(Z,0) in `src/param.F90`
(`rstern_p1`, `rstern_p2`, … for cations; `rstern_m1`, … for anions), consumed by
`ion_alpha_stern` / `chargeaware_atpol` in `src/xdm@proc.f90`.

---

## Why the +q / −q routes differ (the asymmetry — do NOT skip)

- **Cations (+1…+4): bound, closed/near-closed shell.** The Sternheimer (CPKS)
  linear response is well-behaved. Ratios are stable and transferable. → **stern**.
- **Anions beyond −1: unbound in vacuum.** Free O²⁻/S²⁻ do not bind the extra
  electrons; uncoupled linear response **diverges** (verified: confined O²⁻ stern
  α = {982, −31, 0.99, …} garbage across boxes). So −2/−3 cannot use stern.
  → use the **density/compute (Kirkwood) route on the CONFINED ground-state
  density**, where the finite box binds the extra electron and makes ⟨r²⟩ finite.

So: **stern for neutral→any cation; compute for −2/−3 anions; −1 anions can use
either (stern is fine, table already has `rstern_m1`).**

---

## A. Cations out to +4 (stern route)

1. **Patch the charge loop** in `batch_sternheimer.py` (line ~56):
   ```python
   for q in (+1, +2, +3, +4, -1, -2):
   ```
   and add matching accumulator arrays + `emit(...)` calls for `rstern_p3`,
   `rstern_p4` next to the existing `rstern_p2` block (mirror the `rp1`/`rp2`
   pattern near line ~78). **NOTE:** the dev-srv copy is already patched through
   +2; the local mirror copy (`/root/critic2_development/tools/wfc_generator/`)
   is the pre-+2 version — sync from dev-srv before editing.

2. **Run** (ZMAX=54 covers H→Xe; bump to 86 for heavier — Cs…Rn):
   ```bash
   ssh albd@10.10.49.104
   cd ~/critic2_dev/tools/wfc_generator
   python3 batch_sternheimer.py 54 > /tmp/stern/batch.log 2>&1
   # outputs: /tmp/stern/sternheimer_fortran.txt  (rstern_pN arrays)
   #          /tmp/stern/sternheimer_alpha.dat     (absolute α, for sanity)
   ```
   A cation with fewer electrons than its charge (e.g. Li³⁺) is skipped/→1.0 —
   physically it has no valence density, α→0, harmless.

3. **Insert** each `rstern_pN(1:maxzat0)` array into `src/param.F90` right after
   `rstern_p2` (pad to `maxzat0`=123 with `1d0`, as `rstern_p2` is). Keep the
   one-line provenance comment (functional/box/ld1 version).

4. **Extend the interpolation** in `ion_alpha_stern` (`src/xdm@proc.f90`).
   Current (through +2):
   ```fortran
   qc = min(max(q,-1d0),2d0)
   if (qc >= 1d0) then
      f = qc - 1d0; ratio = (1d0-f)*rstern_p1(iz) + f*rstern_p2(iz)
   else if (qc >= 0d0) then
      f = qc;       ratio = (1d0-f) + f*rstern_p1(iz)
   else
      f = -qc;      ratio = (1d0-f) + f*rstern_m1(iz)
   end if
   ```
   To go to +4, raise the clamp to `min(max(q,-1d0),4d0)` and add ladder branches:
   `qc≥3 → p3..p4`, `qc≥2 → p2..p3`, `qc≥1 → p1..p2`. Add `rstern_p3, rstern_p4`
   to the `use param, only:` list.

5. **Validate** before trusting: the +1 ratios a fresh batch emits MUST reproduce
   the existing `rstern_p1` (Mg 0.3044, Ca 0.3209) — this proves the ld1 box +
   H=7.045 absolute miscalibration cancels in the ratio. Spot-check a couple of
   absolute α in `sternheimer_alpha.dat` against physical (Mg²⁺≈0.57 a₀³,
   Ca²⁺≈5.2 a₀³; NOT the 21 the +1-only extrapolation gave).

## B. Anions −2 / −3 (density / compute route)

Stern cannot do these. Instead generate a **confined** ld1 density and let
critic2's **compute** route take Kirkwood α = ⟨r²⟩²/N on it.

1. **Generate the confined `.rho`** (the box binds the extra electron):
   ```bash
   cd ~/critic2_dev/tools/wfc_generator
   ./gen_anion_rho_ld1.py O -2          # standardized box rmax = 3.6*R99(neutral)
   ./gen_anion_rho_ld1.py S -2
   ./gen_anion_rho_ld1.py N -3          # explicit override if needed: --rmax 5
   # batch: ./batch_anions_ld1.sh
   ```
   Box rule: `rmax = ALPHA*R99(neutral)`, ALPHA=3.6; if ld1 reports
   "convergence not achieved" it steps rmax down ×0.8 until it binds. The anion of
   charge q uses zed=Z with the isoelectronic-neutral configuration CONF[Z−q].
   Output `.rho` goes to `~/critic2_dev/dat/wfc/` (or the `wfcdir` you pass to
   critic2 via `wfcdir <dir>`).

2. **critic2 consumes it** automatically in the compute route: the confined
   density makes ⟨r²⟩ finite (the grid fix 42260ea4 already guards the vacuum
   divergence). No param.F90 edit — compute reads the density at runtime.
   Keyword: `... hirshfeld_i alpharef compute wfcdir <dir>`.

3. **Sanity**: confined S²⁻ compute α ≈ 206 a₀³ is still soft-anion-inflated vs
   physical ~35 (documented limit); stern (44.9) is closer but only exists for
   −1. This is the KNOWN soft-anion limitation — for −2/−3 the compute number is
   the best density-based estimate we have; flag it in any solid where a soft
   double-anion dominates (MgS).

---

## Coverage target (what "the full picture" means here)

| charge | route   | how                                    | status |
|--------|---------|----------------------------------------|--------|
| 0      | (free)  | `alpha_free`                           | done   |
| +1     | stern   | `rstern_p1`                            | done   |
| +2     | stern   | `rstern_p2`                            | **built this session** |
| +3,+4  | stern   | extend batch loop + `rstern_p3/p4`     | TODO (A above) |
| −1     | stern   | `rstern_m1`                            | done   |
| −2     | compute | `gen_anion_rho_ld1.py X -2`            | tooling ready (B) |
| −3     | compute | `gen_anion_rho_ld1.py X -3`            | tooling ready (B) |

Main-group priority for A/B: alkali/alkaline-earth (+1/+2, done/in-build),
group 13 (+3: Al³⁺), 14 (+4: Si⁴⁺), halogens (−1, done), chalcogens (−2),
pnictogens (−3). That spans every solid in the campaign and typical ionic crystals.
