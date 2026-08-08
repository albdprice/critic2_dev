# Charge-aware (Hirshfeld-I) XDM dispersion — complete paper trail

**Canonical publication record.** Everything done / doing / planned for the charge-aware XDM
project, with provenance, validation numbers, decisions, and references. Target: a paper with
A. Otero-de-la-Roza (aoterodelaroza) and E. Johnson. Keep this current; it is the single source
for the methods + results narrative. Companion docs (details): [`METHODS_EXPLAINED.md`],
[`REFERENCE.md`], [`METHOD_LIMITS_AND_FIXES.md`], [`ION_REFERENCE_GENERATION.md`],
repo `doc/research/xdm_hirshfeld_i_notebook.md` (chronological lab notebook), memory
`xdm_master_handoff.md` (infra/resume).

Last updated: 2026-08-08.

---

## 1. Problem statement & thesis
Standard XDM (exchange-hole dipole moment) dispersion uses **free-neutral-atom** reference
polarizabilities α, volume-scaled to the atom-in-molecule (AIM) volume:
`α_AIM = α_free · V_AIM/V_free` (Otero-de-la-Roza & Johnson, JCP 136, 174109 (2012), Eq 10).
In **ionic** solids this is wrong: a cation (Li⁺, Na⁺, Mg²⁺, Ca²⁺) is treated as a big soft neutral
atom, giving a spurious cation polarizability and C6, and **systematic over-binding**. The thesis:
replace α_free with a **charge-matched** reference α_ref(Q) tied to the Hirshfeld-I (HI) atomic
charge Q, and dispersion for ionic solids is corrected — reproducing, within XDM, the fix that
Bučko's TS/HI brought to Tkatchenko–Scheffler.

**Headline result (preliminary):** neutral XDM over-binds the ±1 alkali-halide + oxide set; the
charge-aware routes cut the cohesive-energy MAE ~3× (neutral 0.176 → stern 0.057 eV/atom vs
Csonka/Zhang). Cation C6 collapses ~2–3 orders (LiF Li 189→0.1 a.u.), anion C6 grows 3–4×
(F⁻ 7.5→31) — the exact TS→TS/HI signature (Bučko JCTC 2013 Table 1).

## 2. The six α reference routes (keyword-gated; default XDM unchanged)
Keyword surface: `xdm … hirshfeld_i [volonly] [volref] [alpharef gould|scale|compute|stern|sternws] [wfcdir <dir>]`
(molecular `xdm_wfn`); periodic `xdm grid rho … elf … hirshfeld_i alpharef … wfcdir …`.
`ialpharef` code in `xdm@proc.f90`: 0=neutral, 1=gould, 2=scale, 3=compute, 4=stern, 5=sternws.

| # | route | α_ref(Q) model | provenance | charge coverage | status |
|---|-------|----------------|-----------|-----------------|--------|
| 0 | **neutral** | α_free·V_AIM/V_free (baseline) | OdlR&J 2012 | — | shipped |
| 1 | **gould** | α_FI(Q) table, linear-N interp | Gould-Bučko JCTC 12 4644 (2016), arXiv:1604.02751 (LaTeX-source, no paywall) | −1..+1 (frozen-orbital) | commit c483a6f5 |
| 2 | **scale** | α_free·(V_AIM/V_free⁰)^{p′_Z} | Gould JCP 145 084308 (2016), arXiv:1608.04161; p′=p−0.615 | any (volume law) | commit 1457dbdd |
| 3 | **compute** | Kirkwood ⟨r²⟩²/N on confined ld1 ion | this work (Stage 2C) | any (bound density) | shipped |
| 4 | **stern** | uncoupled Sternheimer/CPKS on confined ld1 ion; ratio×α_free | this work (2C-rigorous) | −1..+2 (cations extended) | commits …/1a8dd386 |
| 5 | **sternws** | **Watson-sphere self-consistent Sternheimer** (deep-anion route) | this work (Stage 2D), 2026-08 | 0..+2 cations (via stern tables) + −1..−3 anions (raw Watson α) | commit f00e3ec9 (+ aws-fix, in progress) |

Full per-method write-up with strengths/limits: `METHODS_EXPLAINED.md`.

## 3. Reference polarizabilities — how each is generated (tooling)
All in `tools/wfc_generator/` (dev-srv `~/critic2_dev`; mirror `/root/critic2_development`).
ld1.x = `/usr/bin/ld1.x` (apt QE 6.7); gives consistent RATIOS (H abs 7.045 is a red herring —
cancels). Data enters `param.F90` as compile-time arrays.

- **compute / stern (confined ld1 ion):** `gen_anion_rho_ld1.py` (density), `gen_ion_alpha_sternheimer.py`
  (α by inverting the ld1 KS potential + radial uncoupled Sternheimer), `batch_sternheimer.py`
  → `rstern_p1,rstern_p2` (cations +1/+2), `rstern_m1` (−1). Box rmax=3.6·R99(neutral).
- **sternws (Watson-sphere self-consistent):** `gen_ion_alpha_watson_scf.py` — a **self-contained
  radial KS-DFT atom (LDA/PW92)** with a Watson sphere `V_ws(r)=−|q_net|/max(r,R_W)`,
  `R_W=⟨r³⟩^(1/3)` of the neutral (Bučko's TS/HI rule), INSIDE the SCF; then Pauli-projected
  uncoupled Sternheimer. `batch_watson.py` → `aws_m{1,2,3}` (raw α_ws anion, a₀³) +
  `vws_m{1,2,3}` (r³-moment). Cations reuse rstern_p1/p2.
  (`gen_ion_alpha_watson.py` = the FAILED post-hoc attempt, kept only for the `rmoment_neutral`
  helper — see decision log D4.)

## 4. Validation (the evidence for the paper)
### 4a. Reference α vs experiment (Tessman & Kahn, Phys Rev 92 890 (1953), in-crystal column)
Watson-sphere free-ion α (this work), converted a₀³→Å³ (÷6.7483):

| ion | route | our α (Å³) | Tessman in-crystal (Å³) |
|-----|-------|-----------|--------------------------|
| O²⁻ | sternws (free Watson) | 1.52 | 1.65 (MgO) |
| S²⁻ | sternws (free Watson) | **4.78** | 4.8–5.9 |
| S²⁻ | compute (Kirkwood) | 30 | (too soft) |
| S²⁻ | stern (−1 clamp) | 6.7 | (too stiff) |
| Mg²⁺| stern (+2, rstern_p2) | 0.045 | ~0.09 (MgO) |
| Ca²⁺| stern (+2) | 0.386 | ~0.48 (CaO) |

**sternws is the only route with a stable, first-principles deep-anion α that matches Tessman.**

### 4b. Cross-validation & solver ladder (sternws solver)
- bare-H: α = 4.5003 a₀³ (exact 4.5) — Sternheimer + grid correct.
- The independent Watson SCF reproduces the ld1-Sternheimer CATION tables: Mg²⁺ ratio 0.00475 vs
  tabulated rstern_p2 0.00461 (3%); Ca²⁺ 0.0188 vs 0.0179 (5%).
- Anion α grid-stable to 6 digits over rmax 18→30 a₀ (post-hoc approach was erratic/negative).

### 4c. C6 before→after (a.u.), ours vs Bučko TS Hirshfeld→HI (JCTC 2013 Table 1)
Ours (neutral→stern): Mg(MgO) 175→0.30, Ca(CaO) 563→5.72, O(MgO) 11.3→45.2, S(MgS) 74.5→286.
Bučko (Hirshfeld→HI): Na(NaCl) 1225→7.0, Cl 88→236, Li(LiF) 1039→0.1, F 8.3→31, K(KI) 3048→92, I 367→808.
Same signature: cation collapse ~99%, anion growth 3–4×.

### 4d. +2 divalent fix (head-node, MgO/MgS/CaO)
α(Mg²⁺) 20.2→0.31, α(Ca²⁺) 46.7→2.60 a₀³ (was unphysical +1-extrapolation). See REFERENCE.md.

### 4e. Cohesive energies — targets & status
Experimental targets (authoritative): Zhang NJP 20 063020 (2018) SI + Csonka PRB 79 155107 (2009);
KCl/CaO derived from CRC 97th. Full table in REFERENCE.md (LiF 4.40, NaF 3.93, NaCl 3.33, LiCl 3.55,
KCl 3.35, MgO 5.12, MgS 3.71, CaO 5.50 eV/atom). Bučko method overlay (JCP 2014): LiF/NaCl/MgS
TS-overbinds→TS/HI-fixes. **Full 6-option cohesive campaign on fir — see §6.**

## 5. Key design decisions (decision log)
- **D1 — charge-aware α only, not re-derived moments.** HI weights applied to volumes/α; the
  exchange-hole moments Mₗ kept on neutral weights by default (RQ; task #28).
- **D2 — a1/a2 damping NOT refit for the charge-aware routes.** Damping is trained on neutral
  covalent/vdW dimers (KB49 0.293/2.545; refdata 0.6512/1.4633) where charge-aware α ≈ neutral;
  the ionic solids are the TEST set. Mirrors Bučko TS/HI (reuses TS damping) and Gould FI (β=0.83
  fixed). Confirmed by RQ5 / task #43. **sternws uses the SAME a1/a2 as stern** — the sternws vs
  stern cohesive comparison at fixed a1/a2 IS the deliverable; any shift is anion-physics signal,
  not a damping artifact. (Optional reported robustness check: joint refit incl. ionic — not the
  headline, risks overfitting + breaks transferable-damping philosophy.)
- **D3 — +2/−2 asymmetry.** Cations (+1..+4) are bound → Sternheimer clean, extend freely (stern).
  Free anions beyond −1 are UNBOUND in vacuum → uncoupled response diverges. Resolution: the
  divergence is a vacuum artifact; the crystal Madelung field binds them. Model it with a Watson
  sphere (sternws) → self-consistent bound anion, finite response.
- **D4 — sternws must be SELF-CONSISTENT (Watson in the SCF), not post-hoc.** Adding the Watson
  potential to a hard-box ld1 solution after the fact gave erratic/negative α (O²⁻ 4.96/1.47/6.48/
  0.13/20.1; S²⁻ negative across boxes). Only putting the sphere inside the SCF works.
- **D5 — anions use the RAW Watson α, cations the α_free-renormalized ratio.** The α_free×ratio
  renormalization (which cancels the LDA-uncoupled absolute error, correct for cations) pulls the
  deep anion BELOW experiment (S²⁻ 3.2 vs 4.8–5.9 Å³) because the LDA error is not transferable
  across shell-filling. The raw Watson α already equals the in-crystal value (S²⁻ 4.78 Å³), so the
  anion branch uses aws_m{1,2,3} directly with the existing (≈1) V_AIM/V_free volume scaling.
  [Superseded the earlier "V_ws double-counting" hypothesis: V_AIM≈V_free_ld1 for the diffuse
  anion, so the volume scaling was never the problem — the α normalization was.]
- **D6 — unphysical over-fills clamped.** Anion charges that force an electron into the NEXT shell
  (F²⁻, S³⁻, Kr⁻: ratio blows to 13–15×) are set to 0 → critic2 falls back to neutral scaling.

## 6. fir campaign (the 6-option comparison)
`/scratch/albd/targets/` on fir (acct def-anatole_cpu). `run_targets.sh` = SBATCH array,
53 systems (X23 molecular crystals + 8 ionic solids) × 5 options `(neutral gould scale compute stern)`
= 265 tasks, **job 53769928** (RUNNING as of 2026-08-08). FI-relaxed via the outer-loop driver
`fixdm_val/fixdm_relax.py` (inner: pw.x B86bPBE SCF; frozen-coeff XDM per option; outer: refresh
critic2 HI coeffs at converged geometry). Collate: `collate_ecoh.py`. a1/a2 shared across options.
**sternws (6th option) queued separately** — array 1–53, option=sternws, driver pointing at
`critic2_sternws` (see §7); runs/<name>/sternws. ⚠ the divalent-STERN cells in 53769928 used the
pre-+2 binary → re-run those with the +2 binary when the campaign finishes (or swap binary first).

## 7. Infrastructure (where everything runs)
- **fir** (Alliance; `ssh fir`; acct def-anatole_cpu): all periodic compute. critic2 at
  `/scratch/albd/critic2_ca` (main binary + `critic2_sternp2` [+2] + `critic2_sternws`). **Full source
  tree now staged at `/scratch/albd/critic2_src`; rebuild = `bash -l …/fir_build.sh`** (cvmfs
  StdEnv/2023 gcc12 imkl libxc flexiblas; explicit -DLIBXC_* + -DBLA_VENDOR=FlexiBLAS). To RUN:
  `module load StdEnv/2023 gcc/12.3 imkl/2024.2.0 libxc/6.2.2 flexiblas` + `export CRITIC_HOME=
  /scratch/albd/critic2_ca`. QE 7.5 module; B86b pseudos `/scratch/albd/kb49_b86b_pp`.
- **dev-srv** (`ssh albd@10.10.49.104`): git repo `~/critic2_dev` (branch research/xdm-hirshfeld-i,
  fork albdprice/critic2_dev — COMMIT HERE ONLY, user-attributed). Build on TANK
  `/data/Iterative_hirshfeld/critic2_build`. ld1.x=/usr/bin/ld1.x. Sternheimer/Watson generators.
- **tank** (corvette-local `/tank/research/xdm_chargeaware/`): durable archive. `data/` (this doc +
  all method/reference docs + kb49 mesh/grid data + targets). `papers/solid_refs/` (all PDFs incl.
  PhysRevB.87.064110.pdf, CRC, Bučko JCTC+JCP, Tessman, Csonka, Zhang SI — mirrored 2026-08-08).
- **nibi**: old Psi4/mesh data (2FA-gated; not needed — fir self-sufficient).

## 8. Literature / references (with acquisition notes)
- **Method foundation:** Otero-de-la-Roza & Johnson, JCP 136 174109 (2012) — periodic XDM; our
  toolchain. `papers/solid_refs/oterodelaroza2012_xdm_solids.pdf`.
- **Bučko TS/HI:** JCTC 9 4293 (2013) [Letter, ion-C6 Table 1] = `bucko2013_paperI.pdf` (MISLABELED —
  it is the JCTC letter, not the PRB); JCP 141 034114 (2014) [full TS/HI; LiF/NaCl/MgS cohesive] =
  `bucko2014_TS_HI.pdf`; **PRB 87 064110 (2013)** [TS/TS+SCS on solids; ionic = NaCl+KI only] =
  `PhysRevB.87.064110.pdf`. NB: no single Bučko paper has the full alkali-halide TS/HI cohesive
  table — that was an incorrect assumption; our Zhang/Csonka set is the complete reference.
- **Ion-α databases:** Gould & Bučko JCTC 12 4644 (2016), arXiv:1604.02751 (ingested via arXiv
  LaTeX source, paywall-free); Gould JCP 145 084308 (2016), arXiv:1608.04161 (volume-scaling
  exponents); FI-MBD Gould et al. JCTC 12 5920 (2016)/arXiv:1703.08786.
- **Reference data:** Tessman-Kahn-Shockley Phys Rev 92 890 (1953) [in-crystal ionic α];
  Csonka et al. PRB 79 155107 (2009) + Zhang et al. NJP 20 063020 (2018) [cohesive E]; CRC 97th
  [KCl/CaO thermochem]. All in `papers/solid_refs/`.

## 9. Open items / plan
- [in progress] sternws aws-fix: swap param rstern_ws→aws (raw α), rebuild fir, validate
  S²⁻ AIM ≈ Tessman, commit.
- [queued] sternws fir campaign (X23 + 8 ionic solids, FI-relaxed) → 6-option cohesive comparison.
- [ ] After 53769928 finishes: swap fir main binary → +2/sternws build; re-run divalent-stern cells.
- [ ] Extend cations to +3/+4 (stern) + anions coverage (per ION_REFERENCE_GENERATION.md); task #44.
- [ ] Cohesive-energy sensitivity: sternws vs stern vs neutral at fixed a1/a2 (D2) — the headline table.
- [ ] Molecular ionic benchmark (GMTKN55 IL16/AHB21/CHB6/IONPI19); task #45.
- [ ] Optional: dynamic α(iω) for full C6 (vs the static-α approximation); task #44.

## 10. Reproducibility — exact commands
- Generate Watson anion table: `cd ~/critic2_dev/tools/wfc_generator; LD1X=/usr/bin/ld1.x python3
  batch_watson.py 54` → `/tmp/watson/rstern_ws_fortran.txt` (aws_mN, vws_mN, rstern_ws_mN).
- Single-ion check: `python3 gen_ion_alpha_watson_scf.py O -2` (add `--gridscan` for stability,
  `--check` for bare-H).
- Build fir: `rsync src → /scratch/albd/critic2_src/src/; bash -l /scratch/albd/critic2_src/fir_build.sh`.
- Run a solid one-off: `critic2_sternws c_sternws.cri out` with CRITIC_HOME set (see §7).
- Commits (fork research/xdm-hirshfeld-i): stern+2 `1a8dd386`; sternws `f00e3ec9`; aws-fix (pending).
