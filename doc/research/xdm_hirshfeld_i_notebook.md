# Charge-aware (iterative-Hirshfeld) reference handling for XDM dispersion — lab notebook

**Project:** Does replacing neutral-Hirshfeld partitioning with iterative
Hirshfeld (Hirshfeld-I) in the Becke–Johnson XDM dispersion model improve
it, and what is the correct treatment of the free-atom reference
quantities?

**Maintainers:** A. Price (critic2 fork), with A. Otero-de-la-Roza (XDM
author) consulted. Worklog kept lab-notebook style (dated entries at the
bottom); top sections are the living "paper skeleton."

---

## 0. CURRENT STATE & HANDOFF (read this first)

**Repo / infra**
- Dev machine: `dev-srv` (`ssh albd@10.10.49.104`, sudo pw `800127828`).
- critic2 fork clone: `/home/albd/critic2_dev`; build dir
  `~/critic2_dev/build` (out-of-source; rebuild `cd build && make -j$(nproc) critic2`).
- Local mirror used for editing: `/root/critic2_development` (edit here, `scp` to dev-srv, build).
- Active branch: **`research/xdm-hirshfeld-i`** (fork only; NOTHING pushed
  to upstream `aoterodelaroza/critic2`). Push: `git push origin research/xdm-hirshfeld-i`.
- This notebook is the source of truth. §1–7 = paper skeleton; §7 worklog = dated log.

**What is DONE and merged to fork master (earlier work):** HIRSHFELD_I
keyword + SCF; basindat refactor; Route-1 (Gaussian, `dat/hirshfeld_proatoms/{pbe,pbe0}`)
and Route-2 (confined ld1.x, `dat/hirshfeld_proatoms/ld1_pbe`, 137 files,
Z=1–117 q−1 + p-block q−2) anion databases; standardized `rmax=3.6·R99`;
LiH regression test `tests/009_intgrid/022_hirshfeld_i`. Generators in
`tools/wfc_generator/`. (See memories `hirshfeld-i-design`, `hirshfeld-i-anion-references`.)

**XDM/HI work — status by stage**
- **Stage 0 (grid path, `xdm grid … HIRSHFELD_I [VOLONLY] [WFCDIR]`): DONE, committed.**
  Drops `min(ratio,1)` cap, feeds HI weights to volumes + moments. Validated
  on water (compact): O V/Vfree 0.90→1.15, C6 O–O +53%. **BUT** uniform-grid
  HI SCF is unreliable for all-electron *molecular* densities (aliasing vs
  cusp-loss) — see §f. Grid path is correct for *planewave/pseudopotential*
  densities only.
- **postg/xdm_wfn audit: DONE (§g).** Key result: molecular XDM has **NO cap**
  (postg AND critic2 `calc_coefs`); cap is grid-path-only → Stage-0 removal is
  consistent, not novel. Both integrate on a Becke/Franchini **mesh** (cusp-safe).
- **Option M (mesh HI in molecular `xdm_wfn`, `xdm a1 a2 chf HIRSHFELD_I [VOLONLY] [WFCDIR dir]`):
  WORKING — mesh integration correct AND the ionic SCF now converges (commit 31c6e744, §i).**
  - Implemented: `hirshfeld` helpers `hirsh_i_prepare/refrho/qfloor/cache_clean`
    + shared loader `hi_cache_load`; keyword parse in `xdm_driver`; mesh HI SCF
    + HI volume/moment loop + neutral-moment ablation (VOLONLY) in `xdm_wfn`.
  - Mesh integration correct (NaCl `nelec,total=28.00`; §f grid problem solved).
  - **SCF STABILIZED (§i).** The blocker was a **negative reference density**
    (cubic spline of a box-confined ld1 ref undershoots <0 past the cutoff →
    `phihi=Σrefrho` tiny → weight blows up to 1e23 → limit cycle), NOT mainly
    over-diffuseness. Three fixes: (1) floor `hirsh_i_refrho` at 0 [key];
    (2) element-aware anion clamp `hirsh_i_qfloor` (Cl→−1, O→−2); (3) cation
    clamp `Z−1d-3` not `Z−1` (the latter pinned every H at Q=0) + β=0.3 mixing.
  - Converged, charge-conserving, 0 warnings: NaCl ±0.888 (48 it), LiF ±0.931
    (51), H₂O O−0.870/H+0.435 (83), CH₄ C−0.466/H+0.116 (70). NaCl V/Vfree:
    Na 0.40→0.14, Cl 1.12→1.45 (uncapped). Debug trace gated by param
    `xdm_hi_verbose` (default `.false.`).

**IMMEDIATE NEXT ACTIONS (in order)**
1. **Benchmark neutral vs HI vs HI-VOLONLY** on the test set (NaCl, LiF, H₂O,
   CH₄ — `/tmp/xdmtest/*_mhi.cri`, `*_neu.cri`; add VOLONLY variant): tabulate
   per-atom M1, V/Vfree, and total C6/energy. Goal: quantify the charge-aware
   shift and the volumes-only-vs-moments-too split (RQ2).
2. **Cross-check vs postg** (`~/projects/postg/postg`, neutral only) for the
   neutral baseline C6/volumes.
3. **Stage 1** (task #31): charge-matched volume reference V_ref(Q) in
   `calc_coefs`/`xdm_wfn` (FI volume formula; we have ion densities). Then
   **Stage 2** (#32): ion reference polarizabilities (Gould–Bučko set); decide
   whether the uncapped anion α needs an XDM a1/a2 refit.
4. Add a small ionic-molecule HI-XDM regression test (NaCl charges/C6) so the
   stabilized SCF is protected.

**Molecular HI-XDM run recipe (the working pipeline, via fchk — NOT wfx):**
```
# generate fchk:  g09 PBEPBE/def2TZVP opt (chk) -> formchk   (GAUSS_EXEDIR=~/g09 etc.)
molecule X.fchk
load X.fchk
meshtype franchini small
xdm 0.4 2.5 pbe hirshfeld_i wfcdir /home/albd/critic2_dev/dat/hirshfeld_proatoms/ld1_pbe
#   add VOLONLY after hirshfeld_i for the neutral-moments ablation; drop
#   "hirshfeld_i ..." for the neutral baseline.
```
Test fchks already on dev-srv in `/tmp/xdmtest/` (h2o, NaCl, LiF, CH4) — but
`/tmp` may be cleared; regenerate via `/tmp/xdm_bench.py` (geometries inside).
**Gotchas:** wfx-as-structure reader is buggy (`</`-tag check) → use **fchk**;
inside a chemical function reference the field by **name** (`gkin($mol)`), not `$1`;
`gkin($field)` = orbital τ for the BR hole; numpy 2.x uses `np.trapezoid`;
g09 cubegen mode-5 points are in **Ångström** (don't use — use psi4 collocation or mesh).

---

---

## 1. Question and motivation

critic2's XDM implementation partitions atomic volumes and exchange-hole
moments with **neutral** Hirshfeld weights, and sets the atomic
polarizability to

```
α_A = min( V_A^AIM / V_A^free,neutral , 1 ) · α_A^free,neutral      (current critic2 XDM)
```

i.e. the in-molecule volume is divided by the *neutral* free-atom volume,
scaled by the *neutral* free-atom polarizability, and **capped** so an
atom can only become *less* polarizable than its neutral free state.

Plain Hirshfeld is known to under-respond to oxidation state, so this
choice is expected to misrepresent **ionic / charge-transfer** systems
(oxides, halides, nitrides, layered TMDs), where the polarizable species
are precisely the anions that the cap suppresses. We have just added a
working **Hirshfeld-I** capability to critic2 (keyword `HIRSHFELD_I`,
`.rho` confined free-ion reference densities for q = −1, −2 across
Z = 1–117). The natural question: feed those charge-aware weights and
references into XDM — does it help, and how should the reference
quantities be defined?

## 2. How critic2's XDM uses Hirshfeld (code audit)

From `src/xdm@proc.f90` (volume/moment loop, ll. ~525–590):

```fortran
wei  = rhofree * rhot / max(pdens, 1d-14)     ! = w_A^Hirshfeld(r) · ρ(r)
avol = Σ wei · ri^3                            ! V_A^AIM  (3rd moment)
ml(l)= Σ wei · (ri^l − db^l)^2,  l=1,2,3        ! exchange-hole moments M_l
...
afree(A)  = free_volume(Z)                      ! ∫ ρ_neutral r^3 dr  (NEUTRAL)
alpha(A)  = min(avol/afree, 1) · alpha_free(Z)  ! capped, NEUTRAL ref
```

- `rhofree` = `agrid(Z)%interp(r)` = **neutral** free-atom density (the
  Hirshfeld numerator); `pdens` = **neutral** promolecular density
  (denominator). → standard neutral Hirshfeld weights.
- The dispersion coefficients are built from α and the moments M_l:
  `C6 = α_i α_j M1_i M1_j / (M1_i α_j + M1_j α_i)`, etc.

**Key structural observation (XDM ≠ TS):** the moments `M_l` are
integrals of the *actual molecular density* partitioned by the Hirshfeld
weight — they already adapt to the environment. Only **α** uses the
neutral volume-ratio scaling and the cap. So in XDM the charge-awareness
gap is concentrated in α (the cap + the neutral reference), unlike TS
where *all* environment dependence rides on the single volume ratio.

There is **no published charge-aware XDM**: the `min(ratio,1)` cap and
neutral reference are critic2 implementation choices, not from the XDM
papers. This work is therefore novel.

## 3. Literature synthesis — what charge-aware AIM dispersion methods do

(Deep, multi-source review; 24/25 extracted claims verified 3-0 on
primary sources. Full provenance in §7 log entry 2026-05-30a.)

**Headline:** every charge-aware method does the *opposite* of critic2's
XDM — they do **not** cap at the neutral atom; they replace the neutral
reference with a **charge-matched (fractional-ion) reference**, so anions
become **more** polarizable than the neutral atom (and than cations).

| Method | Volume-ratio denominator | Reference polarizability | Cap at neutral? | Unbound anions |
|---|---|---|---|---|
| **TS** (PRL 102, 073005, 2009) | neutral free-atom V | neutral α_free | no cap; ν stays ≲1 in practice | n/a |
| **TS+HI** (JCTC 9, 4293, 2013; JCP 141, 034114, 2014) | **fractional-ion** ref (interp. between integer charge states) | charge-dependent | **no** — anions > cations explicitly | integer-ion ref densities; *2014 confinement method = open gap (§5)* |
| **FI / MBD@rsSCS-FI** (JCTC 12, 5920, 2016; arXiv:1703.08786) | **fractional-ion** V_FI(N): α_AIM = α_FI(N)·V_eff(N)/V_FI(N) | **fractional-ion**; piecewise-linear interp of integer-ion **TDDFT** α (Gould–Bučko "minimal chemistry" DB) | **no**; Li→Li⁺(0.2 au) not Li⁰(164 au) | confining power-law `(r/r_a)^σ` + ion DB; no bound free anion required |
| **MCLF** (RSC Adv. 2019, c9ra03003d) | — (no volume ratio) | charge-dependent **scaling laws on ⟨r³⟩,⟨r⁴⟩**, neutral data only | n/a | sidestepped — no per-ion reference |

**Findings most relevant to us:**

1. **Our HI density machinery already matches the field.** FI/TS+HI build
   fractional references by *linear interpolation between integer
   charge-state densities* — exactly what `hirsh_i_eval` does — and FI
   handles unbound anions with a *confining potential*, the continuous
   analogue of our box confinement. The reference-density half is done
   and conventional.

2. **No damping refit was needed for FI** (β = 0.83 unchanged), with the
   payoff concentrated where expected: crystal-polarizability error
   157% → 43% (HI) → 23% (FI); ionic solids that *crash* under
   neutral-reference MBD (NaCl, MgO, LiF) succeed under FI; large gains
   for layered TMDs (MoS₂ ≈ RPA); only *modest* gains on
   molecular/covalent (S66/X23). So charge-aware references mainly help
   ionic/charge-transfer systems and can fix outright failures.

3. **"Self-consistent" = partitioning charge↔reference loop only.** The
   dispersion energy never feeds back into the partitioning (a claim that
   it did was refuted 3-0). Our existing HI SCF *is* the whole
   self-consistency story; XDM just consumes the converged partition.

4. **The polarizability *reference* is the hard part.** FI's correctness
   relies on a tabulated database of integer-**ion TDDFT
   polarizabilities** (first 6 rows). critic2 has α_free for **neutral
   atoms only**. A fully FI-faithful HI-XDM needs ion polarizabilities we
   do not yet have.

## 4. Research questions

- **RQ1 (cap):** how much is `min(ratio,1)` currently suppressing? Does
  dropping it (allowing anions > α_free) change ionic-system XDM
  coefficients materially? *(Stage 0)*
- **RQ2 (α only vs also moments): RESOLVED.** E. Johnson (lead PI on XDM,
  2026-05-30, pers. comm.): *"yes, you would want to compute the moment
  integrals using the HI weights as well … that should be
  straightforward, it's the reference atom bit that is not obvious."* →
  apply HI weights to **both** volumes and exchange-hole moments (default);
  the hard/open problem is the **reference atom** (Stages 1–2). We keep
  the volume-only mode as an *ablation* to quantify the moment effect.
- **RQ3 (volume reference):** does swapping V_free,neutral → charge-matched
  V_ref(Q_A) (FI volume formula, using our ion densities) improve over
  neutral-volume + dropped cap? *(Stage 1)*
- **RQ4 (polarizability reference):** is the neutral α_free anchor good
  enough, or is the FI ion-polarizability database necessary for
  accuracy? *(Stage 2)*
- **RQ5 (refit):** do XDM's a1/a2 damping parameters need refitting after
  any of the above? (FI did not need a refit; XDM ≠ MBD — verify.)

## 5. Free-ion reference densities for unbound anions — resolved (Task #27)

**Finding:** the field has no single canonical recipe; unbound-anion
reference densities for Hirshfeld-I are regularized by one of four
families, and **our box confinement is a legitimate member**:

| Approach | Mechanism | Parameter | Notes |
|---|---|---|---|
| **Watson sphere / confining potential** | external (charged shell or wall) potential binds the extra electron | sphere radius & charge | *most common in HI literature*; smooth; mimics a Madelung field |
| **Power-law confinement** `(r/r_a)^σ` | smooth external potential | r_a, σ | used by the rigorous **FI** dispersion work (Gould 2016/17) |
| **Box / hard wall (ours)** | reduce the radial box `rmax` | `rmax` (we tie it to `3.6·R99`) | simplest; crudest (sharp cutoff at wall); we showed the HI charge is flat to ~0.01 e across the gentle-confinement window, which mitigates the parameter dependence |
| **Fractional nuclear charge** (Heidar-Zadeh/Ayers, *J. Mol. Model.* 2017) | raise Z to the smallest *effective* nuclear charge that binds all electrons, then scale | effective Z* | no external cavity; arguably best-behaved density shape for deeply-reduced atoms |

Two corroborating points from the primary literature:
- The 2013/2014 TS+HI papers and the VASP docs are **silent** on the
  confinement recipe (abstracts only; the reference data is a precomputed
  database "for the first six rows except lanthanides"). So there is *no
  published canonical choice to match*; confinement is implementation-
  specific. This vindicates treating it as a free methodological knob.
- The unbound-anion problem has two distinct sub-cases worth stating in
  the paper: monoanions "physically bound but computationally unbound"
  (O⁻ at HF) vs "physically unbound but computationally bound" (N⁻ with
  diffuse DFT). Polyanions are always unbound. All confinement schemes
  exist precisely to regularize these.

**Implication for us:** our box-confinement databases are defensible and
in-family. The **fractional-nuclear-charge** scheme is the most notable
*alternative* generation route (no cavity, better density tails for
highly-reduced atoms) and is logged as a possible Stage-2+ refinement to
compare against our `rmax`-confined `.rho` set. (Sources: search-verified
abstracts of *J. Mol. Model.* 23:341 (2017), 10.1007/s00894-017-3514-6;
FI arXiv:1703.08786; VASP TS+HI wiki.)

## 6. Staged execution plan

- **Stage 0** *(cheapest, most informative)* — keyword-gated `XDM …
  HIRSHFELD_I`: feed `hirsh_i_eval` into the volume/moment loop and
  **drop the `min(ratio,1)` cap**; keep neutral α_free. Two sub-modes:
  (a) HI for α/volume only, (b) HI for volume + moments (RQ2). No refit.
  Benchmark ionic vs covalent.
- **Stage 1** — charge-matched volume reference V_ref(Q_A) from our ion
  densities (FI volume formula); neutral α_free anchor. Re-benchmark.
- **Stage 2** *(if warranted)* — integer-ion reference polarizabilities
  (Gould–Bučko set), piecewise-linear interpolation; full FI-faithful.
  Check a1/a2 refit need.

Default behaviour never changes until validated; everything behind a
keyword.

## 7. Worklog

<!-- newest entries appended below; keep dated, terse, falsifiable -->

### 2026-05-30 — project setup
- Branch `research/xdm-hirshfeld-i` off fork master (which already carries
  the merged `HIRSHFELD_I` keyword, basindat refactor, Routes 1+2 anion
  databases `dat/hirshfeld_proatoms/{pbe,pbe0,ld1_pbe}`, standardized
  `rmax = 3.6·R99` rule, and the LiH regression test).
- Audited `src/xdm@proc.f90`; §2 records the exact weight/volume/moment
  expressions and the cap. Confirmed: XDM uses neutral Hirshfeld for both
  volumes and moments; only α is capped+neutral-referenced.

### 2026-05-30a — literature review (deep multi-source, verified)
- Method: fan-out web search (5 angles) → 17 primary/secondary sources
  fetched → 77 claims → 25 adversarially verified (2/3-refute kill rule)
  → 24 confirmed, 1 killed. Synthesis in §3.
- Primary sources: Tkatchenko–Scheffler PRL **102**, 073005 (2009);
  Bučko/Lebègue/Hafner/Ángyán JCTC **9**, 4293 (2013, ct400694h) and
  JCP **141**, 034114 (2014, OSTI 22419884 / PubMed 25053308);
  Gould/Lebègue/Ángyán/Bučko JCTC **12**, 5920 (2016, acs.jctc.6b00925 /
  arXiv:1703.08786); Manz MCLF, RSC Adv. **9** (2019, c9ra03003d);
  VASP wiki TS / TS+HI / MBD-FI implementation pages.
- Killed claim (0-3): that TS+HI establishes a *dispersion→partitioning*
  self-consistent loop. It does not; self-consistency is partitioning
  charge↔reference only (Finding 3).
- Net decision input: charge-matched references + **no cap** is the
  field standard; the minimal first experiment is to drop our cap and
  feed HI weights in (Stage 0).

### 2026-05-30b — Stage 0 implementation plan (xdm@proc.f90, grid driver)
Code map of the XDM grid driver:
- Parser loop (ll. ~190–290): field selectors + a1/a2/upto/onlyc.
- Promolecular `ipdens` (neutral Σρ_B) and core densities built on the
  grid (ll. ~324–390).
- Volume/moment loop (ll. ~525–575): per nneq atom, sum
  `wei = ρ_A^free · ρ / ρ_promol`; `M_l` from the valence ρ, `V_A` from
  the all-electron ρ (irhoae) or ρ+core. `ρ_A^free` from `agrid` (neutral).
- α (l. ~588): `min(avol/afree,1)·alpha_free` (cap + neutral ref).

Stage-0 change (keyword-gated `XDM … HIRSHFELD_I [WFCDIR dir]`), done in
two sub-steps:
- **0.1 (full HI, this entry):** run the HI SCF (`hirsh_i_driver`) on the
  XDM grid; overwrite `ipdens` with the HI promolecular Σρ_B^{Q_B};
  replace the neutral `agrid` numerator with `hirsh_i_eval` (charged
  ref); **drop the cap** → `alpha = avol/afree·alpha_free`. This makes
  BOTH `V_A` and `M_l` HI-weighted (sub-mode b). Neutral `afree`,
  `alpha_free` kept (Stage 0).
  Requires `irho == sy%iref` (HI SCF integrates `s%f(s%iref)`); warn
  otherwise. Map each nneq atom → a representative cell atom for the HI
  charge state.
- **0.2 (RQ2 toggle):** add `MOMENTS` sub-option to choose whether the
  exchange-hole moments use HI (sub-mode b) or stay neutral (sub-mode a),
  to isolate the moment effect. Implemented after 0.1 builds/runs.
Default XDM path untouched.

### 2026-05-30c — Stage 0 implemented (both sub-modes), builds clean
`src/xdm@proc.f90` (`xdm_grid`):
- New keywords: `XDM … HIRSHFELD_I [MOMENTS] [WFCDIR dir] [HITOL t]`
  (alias `HI`). `dohi` gates everything; default path byte-for-byte
  unchanged.
- Runs `hirsh_i_driver` on the XDM grid (requires `rho == reference
  field`; errors otherwise). Keeps the **neutral** promolecular in
  `ipdens` and stores the **HI** promolecular Σρ_B^{Q_B} in a separate
  `phi_hi` array, so the two sub-modes are selectable per-use:
  - volume/α weight: HI when `dohi` (always, for the polarizability);
  - exchange-hole moment weight: HI only when `MOMENTS` (sub-mode b),
    else neutral (sub-mode a). → directly answers RQ2 by toggling.
- Numerator ρ_A from `hirsh_i_eval` (charged ref) vs `agrid` (neutral);
  nneq→cell-atom map for the per-atom charge state. Thread-safe (HI
  cache preloaded in the driver; eval is read-only in the OMP loop).
- α: cap dropped when `dohi` → `alpha = (avol/afree)·alpha_free`
  (neutral afree/alpha_free anchor kept; that's Stage 1/2's job).
- Cleanup via `hirsh_i_cleanup(bashi)`.
- Builds clean (gfortran 13.3.0). **Not yet run** — needs XDM-ready
  inputs (ρ + τ/ELF + BR-hole b); that's the next entry (benchmark).
- Design note for RQ2: because XDM's M_l already integrate the *real*
  density, sub-mode (a) changes M_l only through the *partition weight*
  (how the shared hole is divided), not through the reference shape —
  so we expect (a) vs (b) to differ less than neutral-vs-HI does. To be
  measured.

### 2026-05-30d — RQ2 resolved by XDM PI; code default flipped
- E. Johnson (XDM lead) confirms HI weights should partition the
  exchange-hole moments too (the "b" sub-mode), and flags the
  **reference-atom** quantities as the genuinely hard part — consistent
  with §3/§4 (charge-matched references are the crux; Stages 1–2).
- Code: `HIRSHFELD_I` now applies HI to **both** volumes and moments by
  default; new `VOLONLY` sub-keyword runs the volume-only ablation. (So
  `MOMENTS` is no longer needed; default = recommended.)
- Re-prioritization: Stage 0 stands (drop cap, HI weights incl. moments);
  the scientific weight shifts to Stage 1/2 (the reference atom). Keep
  the volume-only ablation only to *report* the moment contribution.

### 2026-05-30e — Stage 0 first result: water (PBE/def2-TZVP)
Pipeline (molecular wfx→grid; the wfx-as-structure reader has a latent
bug — `</`-tag check in wfn_private — so use **fchk**; and inside a
chemical function the field must be referenced by *name*, not `$1`):
```
molecule h2o.fchk ; load h2o.fchk id mol
load as "$mol"        80 80 80 id rho1
load as "gkin($mol)"  80 80 80 id tau1   ! orbital τ for the BR hole
reference rho1
xdm grid rho rho1 tau tau1 rhoae rho1 xa1 1 xa2 1 [hirshfeld_i wfcdir …/ld1_pbe [volonly]]
```
HI SCF converged in 13 iters. Atomic units:

| atom | V (neut) | V (HI) | V/Vfree (neut→HI) | M1 (neut) | M1 (HI) |
|---|---|---|---|---|---|
| O | 21.38 | 27.33 | 0.90 → **1.15** | 4.710 | **5.649** |
| H | 5.52 | 3.25 | 0.53 → **0.31** | 1.257 | **0.834** |

C6 (a.u.): O–O 11.46 → **17.58** (+53%); O–H 3.95 → 3.13 (−21%);
H–H 1.49 → **0.583** (−61%).

Three findings:
1. **The cap binds even for water.** HI gives O a volume ratio **1.15 > 1**;
   the old `min(ratio,1)` would clip O's polarizability (~−15%). So
   dropping the cap is not just an ionic-solid concern.
2. **Moments matter (RQ2 / Erin).** HI-weighting the exchange-hole
   moments shifts M1 by +20% (O) / −34% (H) vs the volonly ablation
   (which holds moments at neutral). So the PI's "do the moments too" is
   numerically significant, not cosmetic.
3. **Charge-transfer-consistent redistribution:** dispersion piles onto
   the electron-rich O (C6 O–O +53%) and comes off the electron-poor H
   (C6 H–H −61%) — the qualitatively correct direction.
Caveat: Stage 0 keeps the **neutral** Vfree/α_free anchor and no refit,
so these are *raw partition* changes, not validated energies. Whether
they improve accuracy needs reference C6 and the Stage-1/2 reference
atom (the hard part). Ionic set (LiF, NaCl) + CH4 next.

### 2026-05-30f — KEY FINDING: uniform-grid HI SCF is unreliable for
### all-electron *molecular* densities (ionic set blocked, reframed)
Ran the ionic set (LiF, NaCl) + CH4 through the same molecular wfn→grid
pipeline. **The Hirshfeld-I SCF gives grid/box-dependent garbage charges**
for anything but compact molecules:

NaCl, converged HI charges vs molecular-box border (bohr) / grid:
| border | grid | sum(Q) | Q(Na) | Q(Cl) |
|---|---|---|---|---|
| 6  | 100 | −4.19 | −0.80 | −3.39 |
| 10 | 120 | +0.19 | +0.47 | −0.28 |
| 16 | 160 | +2.42 | +1.70 | +0.71 |

This is **not convergence** — it is two opposing errors crossing:
- small box → **periodic aliasing** of the diffuse density tails →
  over-counted electrons → Q too negative;
- large box (coarser spacing at fixed N) → **nuclear-cusp undersampling**
  of the all-electron density → lost core electrons → Q too positive.
They cancel near border≈10, but there is no robust plateau. Water (tiny,
compact) sat in the good regime by luck (sum Q = 0.005), so the §e water
result stands; LiF/NaCl/CH4 on a uniform grid do **not**.

**Root cause & reframe.** critic2's `xdm grid` is designed for periodic
**planewave/pseudopotential** densities — smooth valence, *no
all-electron cusps* — where a uniform grid integrates cleanly. Forcing an
all-electron molecular wavefunction onto a uniform grid for the HI SCF is
the wrong tool. Two correct ways forward:
- **(P) Native substrate — planewave ionic solids** (MgO, NaCl, LiF) from
  QE/VASP: smooth densities, uniform grid is appropriate, *and* this is
  the regime where charge-aware XDM matters most (TS+HI/FI evidence). The
  current `XDM GRID HIRSHFELD_I` code applies directly.
- **(M) Mesh-based molecular HI** — port the HI SCF + weight evaluation to
  the atom-centred mesh used by the analytic `xdm_wfn` path (`meshmod`,
  Franchini grids), which handles cusps by construction. Bigger code
  change, but gives clean molecular benchmarks (LiF, NaCl, organics).

Decision pending (user/Alberto/Erin): pursue (P), (M), or both. The
Stage-0 grid code is correct and stays; it is the molecular *uniform-grid*
input that is inadequate, not the partitioning. This belongs in the paper
as a methods caveat: HI-XDM volumes require an integration grid that
resolves the reference density used in the SCF.

### 2026-05-30g — postg + xdm_wfn audit → choose M; the cap is grid-only
Read **postg** (Erin/Alberto's reference molecular-XDM code,
`~/projects/postg`) and critic2's own molecular path `xdm_wfn`:
- **No cap anywhere in molecular XDM.** postg sets the atomic
  polarizability as `mol%v(i)*frepol(z)/frevol(z)` and critic2's
  `calc_coefs` as `atpol = v(i)*alpha_free(z)/frevol(z,chf)` — both
  `V_AIM/V_free·α_free` with **no `min(ratio,1)`**. The cap exists *only*
  in critic2's `xdm_grid`. → Stage-0's cap removal aligns the grid path
  with postg *and* critic2's molecular path; it is not a new assumption.
  (Notebook §3 said "no cap" is the TS+HI/FI field standard; postg shows
  XDM itself already does this for molecules.)
- Both integrate volumes/moments on an **atom-centred Becke/Franchini
  mesh** (`m%gen(sy%c,mesh_type,mesh_level)` in `xdm_wfn`), which is
  cusp-safe — exactly the fix for the §f uniform-grid failure.
- `frevol(z,chf)` is functional-dependent (per-XC free-volume tables);
  `alpha_free`/`frepol` are neutral free-atom polarizabilities. Both are
  **neutral** references — the charge-matched versions are the open
  Stage-1/2 problem Erin flagged.

**Decision: implement M** (charge-aware HI in critic2's mesh `xdm_wfn`),
keep P (grid) for periodic/planewave solids. postg is not needed as a
separate engine but is the perfect external **cross-check** (and an
option to port HI into later if the community wants it).

**M design.** `xdm_wfn` currently: mesh → `m%f(:,1)=ρ`, `m%f(:,4)=b`;
per-atom neutral ref `agrid(iz)%interp(r,…)` → promol `m%f(:,2)` →
Hirshfeld weight `ρ_A/promol` → `v(i)`, `mm(l,i)`. To add HI:
1. extract a reusable helper `hirsh_i_refrho(iz, Q, r) → ρ` (charge-aware
   reference density at distance r, from the cached integer-(Z,q) `.rho`
   grids + linear interpolation in Q) — refactor out of `hirsh_i_eval`;
2. **mesh HI SCF** (cusp-safe): iterate Q_A using mesh integration of the
   HI Hirshfeld populations until converged (reuses the mesh `m` and the
   helper);
3. swap the neutral `agrid` ref for `hirsh_i_refrho(iz,Q_A,r)` in the
   promol and per-atom weight; HI weights drive **both** v and mm (per
   Erin); neutral `frevol/α_free` anchor kept (Stage 0 semantics).
Keyword: extend the molecular `XDM … chf` form (and/or the `xdm` driver)
with `HIRSHFELD_I [WFCDIR …]`. Default unchanged.

### 2026-05-30h — M implemented; mesh fixes integration but SCF diverges (ionic)
Built M in `src/xdm@proc.f90`: helpers `hirsh_i_prepare/refrho/cache_clean`
(+ shared `hi_cache_load`) in `hirshfeld`; `xdm a1 a2 chf HIRSHFELD_I
[VOLONLY] [WFCDIR dir]` parsed in `xdm_driver`; mesh HI SCF + HI
volume/moment loop in `xdm_wfn`. Compiles clean; default path untouched.
- **Win:** NaCl mesh integrates to `nelec,total = 28.00` (vs the uniform
  grid's 35.6 in §f) — the cusp/aliasing integration problem is SOLVED by
  the Becke/Franchini mesh.
- **Blocker:** the HI SCF **diverges for NaCl**. Bare iteration: max|dQ|
  ~4.7, wrong-sign charges. With linear mixing β=0.5: WORSE (max|dQ|→1e13,
  Q_Na→−14). Water converges fine. Diagnosis: HI instability for ionic
  systems — Q runs toward the −5 cache clamp → Cl⁻⁵-like ultra-diffuse
  reference → w_Cl≈1 everywhere → runaway. The "reference-atom" hard part
  (Erin). Next: smaller/Anderson mixing, physical per-iteration Q clamp,
  and check whether the confined `ld1_pbe` anion refs are over-diffuse.
  Committed as-is (builds; ionic SCF known-unstable) so state is preserved.

### 2026-05-30i — Option M SCF STABILIZED; ionic molecules converge (commit 31c6e744)
The ionic divergence is fixed. The earlier "diffuse-reference runaway"
diagnosis was only half right — the dominant cause was a **negative
reference density**, not over-diffuseness. Three fixes in `xdm_wfn` +
`hirshfeld@proc`:
1. **`hirsh_i_refrho` floored at zero.** The cubic spline of a box-confined
   ld1 density undershoots **below zero** just past the confinement cutoff.
   Since `ni = Σ w·(refrho/phihi)·ρ` is bounded by N only when every
   `refrho ≥ 0`, a single negative spline value made `phihi = Σ refrho`
   tiny/negative at a mesh point and blew the weight up to `max|dQ| ~ 1e23`,
   which kicked a stable ~15-iteration limit cycle (Cl→+16, Na→−1, then a
   slow march back, repeat). Flooring `refrho` at 0 removed it entirely; the
   SCF then converged **monotonically**. *This is the key fix.*
2. **Element-aware anion floor** (`hirsh_i_qfloor`, new public fn): clamp Q
   to the most-negative *available* reference per element. Probes the WFCDIR
   for `sym_q-N.{rho,wfc}`. Cl→−1 (halogens have no q−2 — chemically correct,
   they don't bind a 2nd e⁻), O/N/C/group13–15→−2. Stops Q escaping to the
   best-effort `read_db` extrapolations near the −5 cache clamp.
3. **Cation clamp `Z−1d-3`** (was `Z−1`): `Z−1` pinned **all hydrogens at
   exactly Q=0** (H has Z=1 ⇒ Z−1=0), wrongly forcing H to be ≤0 when it is a
   cation in H₂O/CH₄. `Z−1d-3` matches `hi_clampq` and lets the q∈(0,1)
   interpolation toward the zero-density bare-proton reference work.
   Plus linear charge mixing β=0.3.

**Converged HI charges** (mesh `franchini small`, `ld1_pbe` refs, all
charge-conserving, 0 warnings):

| molecule | charges | iters |
|---|---|---|
| NaCl | Na +0.888 / Cl −0.888 | 48 |
| LiF  | Li +0.931 / F −0.931  | 51 |
| H₂O  | O −0.870 / H +0.435   | 83 |
| CH₄  | C −0.466 / H +0.116   | 70 |

These are textbook Hirshfeld-I magnitudes (HI > plain Hirshfeld). **Charge-
aware volumes**, NaCl, V/Vfree: Na 0.40→**0.14** (cation contracts), Cl
1.12→**1.45** (anion expands past 1 — exactly what the dropped `min(ratio,1)`
cap used to clip). C6 contribution NaCl −3.892e-4 (neutral) → −3.681e-4 (HI);
H₂O −1.207e-4 → −8.39e-5. Default XDM path + grid regression test unchanged.
- **Stage 0/M deliverable is now functional end-to-end.** A debug per-
  iteration trace is gated behind `xdm_hi_verbose` (param, default `.false.`).
- **Next:** (1) benchmark neutral vs HI vs HI-VOLONLY on the test set
  (volumes, M1, C6) + cross-check the neutral numbers vs postg; (2) decide
  whether the larger anion volumes/uncapped α need an XDM a1/a2 refit
  (Stage 1 #31, Stage 2 #32). The SCF-stability blocker (#34) is closed.

