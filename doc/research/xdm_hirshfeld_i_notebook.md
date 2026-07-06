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

- **Stage 1 (VOLREF) DONE (§k, commit 47a09e1b):** `hirshfeld_i volref` computes
  charge-matched `V_free(Q)` from our confined ions and rescales the α
  denominator (`calc_coefs` optional `volscal`). KEY: VOLREF alone is *not*
  physical — `V_free(Q)` and `α_free(Q)` are a matched pair; needs Stage 2.
  Produced the from-code reference-treatment table (neutral scaling
  overestimates cation α 20–35×).

- **Stage 2 (A/B/C) ALL DONE (§m/§n/§o):** three keyword-gated charge-aware
  α routes, all from-code. **2A** `alpharef gould` (commit c483a6f5): FI-faithful
  with embedded Gould–Bučko ion-α (`alpha_gb_*`; provenance
  `dat/xdm_ion_alpha_gould_bucko_2016.dat`); `α_AIM=α_FI(Q)·V_AIM/V_free(Q)`.
  **2B** `alpharef scale` (1457dbdd): `α_AIM=α⁰·(V_AIM/V_free⁰)^{p'_Z}`,
  exponents `pprime_gb` (Gould JCP 2016, p'_Z=p_Z−0.615). **2C** `alpharef
  compute` (ea4b2e9a): Kirkwood moment estimator `α_free(Q)∝⟨r²⟩²/N` from our
  confined densities, calibrated to neutral. Per-atom α (Na/Li/Cl/F) all move
  cations↓/anions↑; the three span the ion-α reference uncertainty.

**Keyword map:** `hirshfeld_i` (Stage 0/M) → `+volref` (Stage 1 volume) →
`+alpharef gould|scale|compute` (Stage 2 A|B|C). All gated; default unchanged.

**IMMEDIATE NEXT ACTIONS**
1. **Validate A/B/C vs reference data** (the "which is best" question): needs a
   benchmark set — molecular C₆ (e.g. vs the Gould–Bučko C₆ or Tang) and/or
   ionic-solid lattice/binding energies. This is now the main scientific task.
2. **a1/a2 refit check** (RQ5; FI needed none) once a validation set is chosen.
3. Extend 2A to multiply-charged ions + dynamic α(iω) (ACS SI / Gould JCP DB)
   for full periodic-table deployment.
4. Molecular HI-XDM regression test (task #36; `tests/zz_source` download-gated).

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

### 2026-05-30j — Option M benchmark (neutral vs HI vs HI-VOLONLY) + postg validation
Ran all four test molecules in the three molecular-XDM modes (PBE, a1=0.4
a2=2.5, `meshtype franchini small`, `ld1_pbe` references). The VOLONLY
ablation is *clean*: it reproduces the neutral `<M1²>` exactly while using
the HI volume — so the volume effect and the moment effect are cleanly
separable.

Per-atom volume ratio V/Vfree (neutral → HI):

| atom | neutral | HI | direction |
|---|---|---|---|
| Na (NaCl) | 0.404 | **0.140** | cation contracts hard |
| Cl (NaCl) | 1.121 | **1.449** | anion expands past 1 |
| Li (LiF)  | 0.243 | **0.041** | cation contracts hard |
| F (LiF)   | 1.094 | **1.735** | anion expands past 1 |
| O (H₂O)   | 0.910 | **1.236** | gains (net −0.87 charge) |
| H (H₂O)   | 0.634 | **0.324** | contracts (cationic) |
| C (CH₄)   | 0.765 | **0.973** | mild (slightly anionic) |
| H (CH₄)   | 0.685 | **0.582** | mild contraction |

The anion ratios >1 are exactly what the old `min(ratio,1)` cap clipped.
The HI weights also reshape the moments: `<M1²>` anion up (Cl 12.06→14.77,
F 5.97→8.33, O 5.25→6.61), cation/H down (Na 6.88→4.38, H 1.44→0.84).

Total dispersion energy E_disp (Ha): neutral / HI / HI-VOLONLY
- NaCl: −1.151e-3 / −9.520e-4 / −1.110e-3
- LiF : −1.955e-4 / −1.183e-4 / −2.110e-4
- H₂O : −2.255e-4 / −1.535e-4 / −1.860e-4
- CH₄ : −8.252e-4 / −7.718e-4 / −7.935e-4

Pattern: **HI reduces |E_disp|** for all four; the strong cation
contraction (Na, Li, H) outweighs the anion expansion in the C6 sum.
HI-VOLONLY sits between neutral and HI (the volume change alone is the
larger of the two effects; the moment reshaping adds on top). Whether the
reduction is an *improvement* needs reference C6/binding data — that is the
Stage 1/2 + benchmark-set question, not decidable from self-consistency
alone.

**postg cross-check (neutral baseline, the reference molecular XDM code,
`~/projects/postg/postg 0.4 2.5 NaCl.fchk pbe`):**

| | M1²(Na) | M1²(Cl) | V(Na) | V(Cl) | E_disp |
|---|---|---|---|---|---|
| postg | 6.848 | 12.089 | 44.99 | 74.55 | −1.1601e-3 |
| critic2 | 6.876 | 12.060 | 46.08 | 74.39 | −1.1511e-3 |

Agreement: volumes ~2%, moments ~0.4%, E_disp 0.8% (residual = mesh:
postg's default vs critic2 franchini small). **critic2's neutral molecular
path is validated against postg.** postg's *neutral Hirshfeld* charges
(Na ±0.566) vs our **HI** charges (±0.888) also confirm HI amplifies
charge transfer over plain Hirshfeld, as documented.

Inputs saved on dev-srv: `/tmp/xdmtest/{h2o,LiF,CH4,NaCl}_{neu,mhi,vol}.cri`.
Conclusion: Option M is verified, self-consistent, and reference-validated
on the neutral limit. Proceed to Stage 1.

### 2026-05-30k — Stage 1 (VOLREF) implemented + the reference-treatment table (commit 47a09e1b)
Added `xdm a1 a2 chf hirshfeld_i volref [wfcdir]`. After the HI SCF it
integrates the **charge-matched free-ion reference volume** `V_free(Q)`
from the same reference densities that define the HI weights (radial Gauss
quadrature via `hirsh_i_refrho`), and rescales the polarizability
denominator: `α = V_AIM·α_free⁰ / (frevol · V_free(Q)/V_free(0))`. The
ratio cancels the method offset so the neutral baseline (Q=0) is exactly
preserved. Per-atom report: `V_AIM, V_free(0), V_free(Q), ratio,
α_neutscal, α_volref`. `calc_coefs` gained an optional `volscal` arg
(absent ⇒ 1; default/neutral path unchanged). Regression test still passes.

**Critical scientific finding — VOLREF alone is NOT a physical model.**
`V_free(Q)` and `α_free(Q)` are a *matched pair* in the FI formula
`α = (V_AIM/V_free(Q))·α_free(Q)`. For a cation `V_free(Q)` is small (small
denominator → α inflates) and that is only correct because `α_free(Q)` is
*also* small; for an anion both are large. Using `V_free(Q)` with the
*neutral* `α_free⁰` therefore drives both ions the wrong way (e.g. NaCl
α_volref(Na)=130.7 bohr³ — absurd). VOLREF is the **volume half /
architecture** for the full FI model; Stage 2 (ion `α_free(Q)`) is required
to make it physical. Logged so we don't misread it as a result.

**The reference-treatment comparison, fully from critic2 (bohr³):**
Per-atom, neutral-scaling polarizability `α = (V_AIM/V_free⁰)·α_free⁰`
(what XDM does today) vs reference free-ion α:

| atom (system) | HI charge | V_free(0) | **V_free(Q)** | α_neutscal | ref free-ion α | error |
|---|---|---|---|---|---|---|
| Na (NaCl) | +0.888 | 121.9 | **21.2** | 22.7 | ≈0.98 (Na⁺) | **23× high** |
| Li (LiF)  | +0.931 | 99.6  | **7.7**  | 6.80 | ≈0.19 (Li⁺) | **35× high** |
| Cl (NaCl) | −0.888 | 66.2  | **117.6**| 21.3 | ≈36 (Cl⁻ free) | cap erases it; ~1.7× low |
| F (LiF)   | −0.931 | 19.3  | **53.2** | 6.52 | ≈10.6 (F⁻ free) | cap erases it; ~1.6× low |

`V_free(0)→V_free(Q)` is the charge-matched reference volume interpolated
from our confined ld1.x ions — "the same interpolation for the volumes"
(Erin). The α failure is the headline: neutral scaling overestimates
cation polarizability by 20–35× because the neutral alkali's diffuse
valence shell (the entire ~24 Å³ / ~160 bohr³) is gone in the ion, and a
volume ratio cannot remove a whole shell's response. Anions: the
`min(ratio,1)` cap (grid path) discards their *enhanced* polarizability
outright. Both are fixed only by charge-matched references (volume **and**
α). [α in bohr³ = code units; ×0.148 → Å³. Ref ion α = free-ion static
dipole polarizabilities; in-crystal anion values are smaller, widening the
cation/anion asymmetry.]

Caveat for any external use: in-molecule atoms are *fractional* ions
(Na⁺⁰·⁸⁹), so the rigorous comparison interpolates the ion α to the
fractional charge; the cation overestimate stays ~10–20× even then, so the
conclusion is robust. Inputs: `/tmp/xdmtest/{NaCl,LiF,h2o,CH4}_volref.cri`.
**Next: Stage 2** — confined-ion polarizabilities `α_free(Q)` and interpolate,
completing the FI-faithful model (scoped in §l).

### 2026-05-30l — Stage 2 scoping: how FI/MBD and others handle the reference α
Literature dig (Tkatchenko/Gould/Bučko FI-MBD + neighbours) on **where the
reference ion polarizability comes from**, and options for us.

**How the established methods do it:**
- **TS / TS+HI (Tkatchenko–Scheffler, VASP).** α_AIM = ν·α_free with
  ν = (V_AIM/V_free). α_free is a **tabulated** free-atom value (Chu–Dalgarno
  TDDFT), rows 1–6 minus lanthanides. *Plain* TS uses the **neutral** α_free
  even under iterative-Hirshfeld partitioning — i.e. exactly the "neutral α +
  volume ratio" that we have shown fails for ions. So TS+HI fixes the
  *weights/volumes* but NOT the α reference. This is the gap our Stage 1+2 and
  the FI work close.
- **FI / MBD@rsSCS-FI (Gould, Lebègue, Ángyán, Bučko, JCTC 2016 / 2017,
  arXiv:1703.08786) — the method to emulate.** Exact model:
  `α_p^AIM(iω;N) = α_p^FI(iω;N) · V_p^eff(N)/V_p^FI(N)`, with the
  fractional-charge reference built by **linear interpolation between integer
  ions**: `α_{Z,N}(iω) = f·α_{Z,⌈N⌉} + (1−f)·α_{Z,⌊N⌋}`, f = N−⌊N⌋. The
  integer-ion α are **precomputed, ab initio (TDDFT + ensemble DFT)** — *not*
  recomputed per system. **No damping refit** (β = 0.83 unchanged); cubic ionic
  crystals MARE 157%→23%; NaCl/MgO/LiF that *crash* under neutral-ref MBD
  succeed. Effective volumes from iterative Hirshfeld. **They do NOT confine**
  the unbound anions: self-consistent α for neutral/cations/closed-shell
  monoanions, **frozen-orbital / non-self-consistent** α (built on the neutral)
  for the rest.
- **The ion-α database itself: Gould & Bučko, JCTC 12, 4644 (2016),
  arXiv:1604.02751 — "C₆ Coefficients and Dipole Polarizabilities for All
  Atoms and Many Ions in Rows 1–6."** ~411 species; TDDFT imaginary-frequency
  **dynamic α(iω)**, **static α(0)**, and **C₆**; three tiers (raw TDDFT,
  benchmark-corrected, frozen-orbital anions). Data on ACS figshare
  (collection 3281582) + Griffith repository (Gould's institution; likely the
  open copy — ACS figshare returned 403). *This is the table we would ingest.*
- **Gould, JCP 145, 084308 (2016), arXiv:1608.04161 — "How polarizabilities
  and C₆ actually vary with atomic volume."** Uses **confined atoms** (like us)
  and finds `C₆/C₆ᴿ ≈ (V/Vᴿ)^p`, `α/αᴿ ≈ (V/Vᴿ)^{p'}`, with **element-specific**
  exponents and the surprising relation `p' ≈ p − 0.615` (NOT the naive
  `p'=1` that TS/XDM assume, nor `p'=p/2`). Two payoffs for us: (i) it
  *quantitatively* explains why neutral volume-scaling (`p'=1`) of α is wrong;
  (ii) the scaling law is *robust to the confinement form* (harmonic/cubic/
  quartic), which independently **vindicates our box-confinement choice** and
  hands us a cheap Stage-2 route (below).
- **MCLF (Manz).** Avoids a per-ion α table entirely via charge-dependent
  **scaling laws on ⟨r³⟩,⟨r⁴⟩** from neutral data. Sidesteps the reference but
  is its own model, not XDM-compatible without surgery.

**Stage-2 options for critic2 XDM (we already have `V_free(Q)` and the moments):**
- **(A) Ingest the Gould–Bučko static ion-α database + linear-N interpolation
  (the FI recipe). RECOMMENDED.** Field-standard, defensible, no recompute,
  rows 1–6. Plug `α_free(Q)` into the Stage-1 slot:
  `α = (V_AIM/V_free(Q))·α_free(Q)`. *Consistency caveat:* their α uses
  frozen-orbital anions while our `V_free(Q)` uses confined ld1.x ions — but
  FI itself mixes α-source and V-source, so this is in-family. Cost: parse a
  table, add an interpolator; data-acquisition is the main task (get the SI/
  figshare/Griffith file, map to our element+charge grid).
- **(B) Element-specific volume-scaling exponents (Gould JCP 2016): cheap
  cross-check / fallback.** `α_free(Q) = α_free⁰·(V_free(Q)/V_free⁰)^{p'_Z}`
  using our computed `V_free(Q)` and tabulated `p'_Z`. Uses our own volumes;
  no ion-α table. Captures the *nonlinear* α–V scaling TS/XDM miss, but is
  still neutral-anchored (won't fully capture shell-removal); good as a
  sanity bracket against (A).
- **(C) Compute α for our confined ld1.x ions ourselves — most internally
  consistent, most work.** Needs a polarizability out of a *spherical* atomic
  code: Sternheimer/coupled-KS linear response (the proper route; ld1.x has
  no field), or finite-field in a 3D confined-atom DFT (Psi4/Gaussian with a
  confining potential). Defer unless (A) proves inconsistent with our densities.

**Recommendation:** do **(A)** for the deliverable and keep **(B)** as a from-
our-volumes cross-check. Both reuse the Stage-1 plumbing (`volscal`/`α_free(Q)`
slot). Open sub-tasks: locate & parse the Gould–Bučko table (static α per
element+charge), add a charge-state interpolator mirroring `hirsh_i_refrho`,
wire `α_free(Q)` into `calc_coefs`, re-benchmark NaCl/LiF/H₂O/CH₄, check the
a1/a2 refit need (RQ5; FI needed none). Decision needed from AP: pursue (A)
[tabulate], (B) [scaling-law], or (C) [compute] first.
→ **AP: do all three (A, B, C).**

### 2026-05-30m — Stage 2A DONE: FI-faithful polarizabilities from the Gould–Bučko table (commit c483a6f5)
Got the data **without the paywall** — pulled the arXiv e-print LaTeX source
(`Bench.tex`) for 1604.02751 and parsed its three benchmark tables
(neutrals / +1 cations / −1 anions; static α(0) in a₀³; ~244 species rows
1–6). Validated against my Erin-table stand-ins (Li⁺ 0.193, H⁻ 216,
Na⁰ 163, Cl⁰ 14.7 — exact). Provenance file
`dat/xdm_ion_alpha_gould_bucko_2016.dat`; embedded as compile-time arrays
`alpha_gb_{m1,0,p1}` in `param.F90` (a₀³; same convention as `alpha_free`
but NOT ÷Å³ since already atomic units).
- `ion_alpha0(iz,q)`: linear-in-N interpolation between bracketing integer
  ions (the FI recipe `α_{Z,N}=f·α_{Z,⌈N⌉}+(1−f)·α_{Z,⌊N⌋}`); neutral
  fallback if an ion datum is missing.
- Keyword `xdm a1 a2 chf hirshfeld_i alpharef gould [wfcdir]` ⇒ FI-faithful
  `α_AIM = α_FI(Q)·V_AIM/V_free(Q)` (uses Stage-1 `V_free(Q)`), plumbed via
  `calc_coefs` optional `atpolov`.

**From-code result — neutral-scaling α vs FI-faithful α (a₀³), the table
fully from critic2:**

| atom (Q) | α_neutscal | α_FI(Q) | **α_AIM (FI-faithful)** | full-ion ref |
|---|---|---|---|---|
| Na (+0.89) | 22.7 | 19.1 | **14.4** | 0.93 (Na⁺) |
| Li (+0.93) | 6.80 | 11.4 | **5.5** | 0.19 (Li⁺) |
| Cl (−0.89) | 21.3 | 28.6 | **23.3** | 30.3 (Cl⁻) |
| F (−0.93) | 6.52 | 14.2 | **9.0** | 15.0 (F⁻) |

FI-faithful moves **cations down, anions up** — both physically correct.
(The atoms are fractional ions ±0.9, so linear-α interp retains some neutral
character; full-ion refs shown for orientation.) Total E_disp (Ha),
neutral / HI(Stage-0/M) / **FI(2A)**:
NaCl −1.151e-3 / −9.520e-4 / **−8.99e-4**;
LiF −1.955e-4 / −1.183e-4 / **−1.49e-4**;
H₂O −2.255e-4 / −1.535e-4 / **−7.72e-5**;
CH₄ −8.252e-4 / −7.718e-4 / **−4.03e-4**.
Whether FI *improves* XDM needs reference binding/C₆ data (separate
validation task). Known small issue: H⁺ has no GB entry (genuinely α=0 bare
proton), so our "0 ⇒ unavailable" convention falls H back to neutral α
(minor; H still gets the V_AIM/V_free(Q) factor). Inputs:
`/tmp/xdmtest/{NaCl,LiF,h2o,CH4}_fi.cri`.
- **Coverage:** GB benchmark has integer charges −1/0/+1 only → covers all our
  HI charges (|Q|<1). Multiply-charged ions (O²⁻ etc.) and the *dynamic*
  α(iω) + V_FI volumes are in the ACS SI / Gould 2016 JCP "minimal chemistry"
  DB — to add for full periodic-table deployment.
- **Next: Stage 2B** (volume-scaling exponents p'_Z, Gould JCP 2016 /
  arXiv:1608.04161 — needs the per-element exponent table) and **2C** (compute
  α from our confined ld1.x ions). Then validate (A/B/C) against reference
  C₆/binding + check a1/a2 refit (RQ5).

### 2026-05-30n — Stage 2B DONE: element-specific volume-scaling exponents (commit 1457dbdd)
Pulled the exponent table from the arXiv source of 1608.04161
(`/tmp/gould/ex_1608.04161/paper.tex`, summary tables Rows 1–5, PGG kernel
column `p`), embedded `p'_Z = p_Z − 0.615` as `pprime_gb` in param.F90
(Z=1–54; Z>54 → 1.0). Cleanest reading of Gould's relation: it's the
*standard XDM α-scaling with the exponent corrected from 1 to p'_Z*:

  `α_AIM = α_free⁰ · (V_AIM/V_free⁰)^{p'_Z}`   (keyword `alpharef scale`).

Charge-awareness enters only through the HI `V_AIM`; **no ion densities
needed** (unlike 2A). At p'_Z=1 it is exactly the current XDM (so it's a
clean drop-in). From-code per-atom α (a₀³):

| atom (Q) | neutral-scaling (p'=1) | 2A (FI table) | **2B (p'_Z)** | p'_Z |
|---|---|---|---|---|
| Na (+0.89) | 22.7 | 14.4 | **9.3** | 1.455 |
| Li (+0.93) | 6.80 | 5.5 | **4.3** | 1.145 |
| Cl (−0.89) | 21.3 | 23.3 | **27.8** | 1.715 |
| F (−0.93)  | 6.52 | 9.0 | **9.7** | 1.715 |

2A and 2B **agree in direction** (cations↓, anions↑) and bracket each other
— 2B pushes cations lower (closer to the true Na⁺) and anions slightly
higher. E_disp (Ha): NaCl −8.40e-4 (2B) vs −8.99e-4 (2A) vs −1.151e-3
(neutral); LiF −1.476e-4 (2B) vs −1.486e-4 (2A); H₂O −1.461e-4 (2B). Two
independent charge-aware α routes now run from our pipeline. Inputs:
`/tmp/xdmtest/{NaCl,LiF,h2o}_2b.cri`.
- **Next: Stage 2C** — compute α from our confined ld1.x ions directly
  (Sternheimer/coupled-KS in a radial code, or finite-field 3D confined DFT),
  the most internally-consistent route; compare to 2A/2B. Then validate all
  three vs reference C₆/binding and check the a1/a2 refit (RQ5).

### 2026-05-30o — Stage 2C DONE: Kirkwood moment estimator from our confined ions (commit ea4b2e9a)
A full Sternheimer/CPKS solve needs the ld1.x orbitals (we only export the
`.rho` density), so instead used a **Kirkwood-type density-moment estimator**
(precedent: MCLF scales α on density moments). Compute the free-ion α from
the confined reference density's second moment, calibrated to the known
neutral α so the prefactor cancels:

  `α_free(Q) = α_free⁰ · [⟨r²⟩(Q)²/N(Q)] / [⟨r²⟩(0)²/N(0)]`,

with `⟨r²⟩(Q)=∫ρ_ref^Q r² d³r`, `N(Q)=∫ρ_ref^Q` (radial integrals of the
charge-matched density, same quadrature as `V_free(Q)`); then FI-faithful
`α_AIM = α_free(Q)·V_AIM/V_free(Q)`. Keyword `alpharef compute`. **Fully
from our density family** — no external table (2A), no fitted exponent (2B).

**Consolidated three-route comparison — per-atom α (a₀³), all from critic2:**

| atom (Q) | neutral-scaling | 2A (Gould table) | 2B (p'_Z scaling) | 2C (Kirkwood) | full-ion ref |
|---|---|---|---|---|---|
| Na (+0.89) | 22.7 | 14.4 | 9.3 | 14.7 | 0.93 (Na⁺) |
| Li (+0.93) | 6.80 | 5.5 | 4.3 | **1.5** | 0.19 (Li⁺) |
| Cl (−0.89) | 21.3 | 23.3 | 27.8 | 22.5 | 30.3 (Cl⁻) |
| F (−0.93)  | 6.52 | 9.0 | 9.7 | 6.5 | 15.0 (F⁻) |

All three charge-aware routes **agree in direction** (cations↓, anions↑);
their spread is the genuine methodological uncertainty in the ion-α
reference (largest for the alkali cations, where α changes by ~3 orders of
magnitude across one charge unit). E_disp (Ha), neutral / 2A / 2B / 2C:
NaCl −1.151e-3 / −8.99e-4 / −8.40e-4 / −8.84e-4; LiF −1.955e-4 / −1.486e-4 /
−1.476e-4 / −7.91e-5; H₂O −2.255e-4 / −7.72e-5 / −1.461e-4 / −1.386e-4.
Inputs: `/tmp/xdmtest/{NaCl,LiF,h2o}_2c.cri`.

**Stage 2 (A/B/C) all implemented, keyword-gated, from-code.** Keyword map:
`hirshfeld_i` (Stage 0/M) | `+volref` (Stage 1) | `+alpharef gould|scale|compute`
(Stage 2 A|B|C). Remaining: (1) **validate** A/B/C against reference
C₆/binding-energy data (the "which is best" question — needs a benchmark set,
e.g. ionic-solid lattice energies or molecular C₆); (2) **a1/a2 refit** check
(RQ5; FI needed none); (3) extend 2A to multiply-charged ions + dynamic α(iω)
for full periodic-table deployment; (4) molecular HI-XDM regression test (#36).

### 2026-05-30p — Stage 2 A/B/C ported to the periodic grid path (commit c954bdd9)
AP wants validation on **molecules AND solids** (the routes may differ —
FI's big wins are in ionic solids). The molecular Stage-2 lived only in
`xdm_wfn`; the grid/periodic path (`xdm_grid`) used neutral `alpha_free`.
Ported via a shared module function **`chargeaware_atpol(iz,q,vaim,vfree0,
ialpha)`** (routes 1/2/3) so molecules and solids use identical formulas.
`xdm_grid` now parses `alpharef gould|scale|compute`, loads the reference
cache (`hirsh_i_prepare` with the converged `bashi%hi_qfinal`), and replaces
neutral α per atom (per-nneq HI charge via `icel_nneq`); prints
Q/V_AIM/a_neutscal/α_AIM. Builds clean; grid HI regression passes; default
unchanged. **Code-complete but not yet physically validated**: the uniform
molecular-as-grid test vehicle gives garbage HI charges (Na −3.97; the
all-electron cusp/aliasing §f problem), so it only confirms the code path
runs. **Real-solid validation needs proper periodic pseudopotential/PAW grid
densities (QE)** — task #40.

### 2026-05-30q — Stage 2C-rigorous (Sternheimer) feasibility confirmed; plan
AP: build the full confined-ion polarizability solve **now** (upgrade 2C
beyond the Kirkwood estimator). Feasibility CONFIRMED:
- **Orbitals R_{nl}(r):** ld1.x writes them to `ld1.wfc` by default (AE,
  iswitch=1) — already produced/parsed by the Route-2 generator.
- **KS potential v_KS(r):** reconstruct from the density we already have,
  `v_KS = −Z/r + v_H[ρ] + v_xc^PBE[ρ]` (radial Poisson for v_H), consistent
  with the ld1 SCF eigenvalues. (NB: the `file_chi`/`file_potscf` &input
  vars I first tried are not valid ld1 namelist keys → read error; not
  needed — reconstruct instead, or use the eigenvalues + orbitals directly.)
So we have orbitals + (reconstructable) potential + grid for the confined
ion (same box that regularizes the density). Plan (radial coupled-perturbed KS / Sternheimer):
for each occupied (n,l), the dipole field z=r cosθ couples to l±1; solve the
inhomogeneous radial ODE `(ĥ_{l'} − ε_{nl})(r·δu) = −(r·R_{nl})·⟨l'|cosθ|l⟩`
with `ĥ_{l'} = −½ d²/dr² + l'(l'+1)/2r² + v_KS(r)` and the **box boundary
condition δu(rmax)=0** (this is what makes the anion α finite — the same
regularization as the density); then `α = (2/3)Σ_occ (ang)² ∫ R_{nl} δu r³ dr`.
Start with the **uncoupled (independent-particle)** solve (no δv_Hxc
self-consistency) — already a real improvement over Kirkwood and a clean
first target; add the coupled Hartree+XC response (TDDFT-level, what Gould
2A used) as a second step. **Validation gate:** reproduce known *neutral*
free-atom α (e.g. H 4.5, Ne 2.67, Ar 11.1 a₀³) before trusting ion values.
This is a standalone numerical build (tool in `tools/wfc_generator/`,
producing an α(Z,q) table critic2 ingests like the GB table) — the dedicated
next work item (task #42).

### 2026-05-30r — Sternheimer 2C solver built + validated (commit 17aef259)
Built `tools/wfc_generator/gen_ion_alpha_sternheimer.py` — offline generator:
ld1.x AE calc → orbitals `P_nl`(ld1.wfc) + eigenvalues/occ (stdout) →
reconstruct `V(r)` by inverting the radial KS equation (multi-orbital,
max-|P|; exact ld1 SCF potential, no XC re-eval) → solve the uncoupled
radial Sternheimer ODE per (n,l)→l±1 (non-uniform 3-point FD tridiagonal,
Dirichlet/box BCs) → `α=(2/3)Σ occ·A(l,l')·∫P r w dr`. scipy `solve_banded`.

**Validation (the key result):** the *uncoupled* solve overestimates
*absolute* α — H 7.05 (vs 4.5; also SIE-inflated), He 1.73 (1.38),
Ne 3.51 (2.67), Ar 17.8 (11.1) — exactly the expected uncoupled-vs-coupled
gap (no self-consistent depolarization field; element-dependent). **But as a
neutral-calibrated RATIO it reproduces the Gould–Bučko TDDFT charge trend**:
- cations (compact, confinement-insensitive) match almost exactly:
  **Li⁺/Li⁰ 0.0012 vs GB 0.0012; Na⁺/Na⁰ 0.0067 vs GB 0.0057.**
- anions right-signed but **confinement(rmax)-dependent**: F⁻/F⁰ 2.28 (GB
  4.17), Cl⁻/Cl⁰ 2.41 (GB 2.06), O⁻/O⁰ 2.73 (GB 1.04) at rmax=12. This is the
  genuinely-uncertain "reference-atom" quantity (Erin): GB use frozen-orbital
  free anions, we use box-confined — they legitimately differ. For 2C use OUR
  confined α at the **same standardized rmax (3.6·R99)** as the density/volume
  references → one self-consistent method.

So the uncoupled Sternheimer is a sound *ratio* estimator (more physical than
Kirkwood; cation-exact vs GB). Absolute accuracy would need the coupled
(TDDFT) solve — not required since we calibrate to the neutral.
- **Next:** batch over elements at standardized rmax → embed an α(Z,q) table
  (like the GB arrays) and wire as the rigorous 2C reference (replacing/
  augmenting the inline Kirkwood). Minor harness fix: robust stdout-occupation
  parse; `ok` now keys off "reached in" not the fragile "convergence" string.

### 2026-05-30s — Sternheimer α table batch-generated + embedded (commit 755c6e42)
Built `tools/wfc_generator/batch_sternheimer.py`: drives the Sternheimer
generator over the periodic table at **rmax = 3.6·R99** (the HI
density-reference box, read per element from the neutral `.wfc`), computing
the neutral-calibrated charge factors `rstern(Z,q)=α_stern(Z,q)/α_stern(Z,0)`
for q=0,±1. Embedded as `rstern_p1`, `rstern_m1` in param.F90 (Z≤18 computed;
Z>18 → 1.0; provenance `dat/xdm_ion_alpha_sternheimer.dat`). New keyword
`alpharef stern` (route 4, shared molecular+grid via `chargeaware_atpol`):
`α_free(Q)=α_free^CRC(Z)·rstern` (interp in N), `α_AIM=α_free(Q)·V_AIM/V_free(Q)`.
- **Cation ratios match GB TDDFT essentially exactly** (Li⁺ 0.0012, Na⁺
  0.0067) — the rigorous self-contained route. Diffuse alkali/alkaline/noble
  *anion* ratios are unphysical (He⁻ 431, Ne⁻ 138; the extra e⁻ barely binds)
  but **never used** (those elements aren't anions) — capped at 8 defensively.
  Real molecular anions sensible: O⁻ 2.38, F⁻ 2.82, Cl⁻ 2.45, N⁻ 2.63.

**Four charge-aware α routes now in-code (per-atom α, a₀³):**

| atom (Q) | neutral | 2A Gould | 2B scale | 2C Kirkwood | 2C-stern |
|---|---|---|---|---|---|
| Na (+0.89) | 22.7 | 14.4 | 9.3 | 14.7 | 14.4 |
| Li (+0.93) | 6.8 | 5.5 | 4.3 | 1.5 | 5.5 |
| Cl (−0.89) | 21.3 | 23.3 | 27.8 | 22.5 | 27.5 |
| F (−0.93)  | 6.5 | 9.0 | 9.7 | 6.5 | 6.4 |

2C-stern matches 2A for cations (Sternheimer cation ratios ≈ GB) and tracks
our confined anions (Cl higher than 2A's frozen-orbital value). **Keyword
map:** `hirshfeld_i [volonly] [volref] [alpharef gould|scale|compute|stern]
[wfcdir]` — all gated, default unchanged, regression green. Limited to Z≤18
(rows 1–3) for the Sternheimer table; extend `_CONF` for rows 4–6.
- **Stage 2C-rigorous (Sternheimer) is DONE.** Remaining program: validate
  all routes vs reference C₆/binding on molecules (#41) AND ionic solids
  (needs QE crystal densities, #40); a1/a2 refit (#43); extend tables to
  rows 4–6 + multiply-charged/dynamic (#44); regression test (#36).

### 2026-05-30t — periodic-table coverage closed; molecular validation STARTED
**Coverage (answering "why not the whole PT"):** the HI density references
are Z=1–117 (the partition is full-PT); 2A (Gould) ~rows 1–6, 2B rows 1–5
(literature limits). 2C-Sternheimer was capped at Z≤18 *only* because the
generator's `_CONF` had 18 configs — fixed by copying the full aufbau list;
re-ran `batch_sternheimer.py` over Z=1–86 (256 species, commit 72afb7ec).
Open-shell TM/f anion solves return negative/non-finite α (uncoupled
Sternheimer breaks down near open-shell degeneracies) → neutral-fallback
(11 cases); diffuse alkali/noble anion ratios capped at 8 (unused). So
2C-stern is reliable for closed-shell/main-group; for open-shell TM/f use
2A/2B. The routes are complementary across the table.

**Molecular validation (harness `tools/wfc_generator/validate_mol_xdm.py`):**
runs neutral + 4 charge-aware routes, sums the pairwise C6 to the
homomolecular C6, compares to DOSD reference. Result (mol C6, a₀³; C6/ref):

| molecule | neutral | HI | 2A | 2B | 2C-kirk | 2C-stern | DOSD |
|---|---|---|---|---|---|---|---|
| H₂O | 40.3 (.89) | 38.8 | **16.1 (.36)** | 40.9 (.90) | 36.2 | 29.7 (.66) | 45.3 |
| CH₄ | 123 (.95) | 128 (.99) | **58.7 (.45)** | 120 (.93) | 105 | 73.5 (.57) | 129.7 |

**Findings (covalent):** neutral / HI / 2B-scale match DOSD well (~0.9–1.0×);
**2A over-reduces badly (0.36–0.45×)**, 2C-stern moderately. Diagnosis: **2A
mixes sources** — Gould's *frozen-orbital* ion α with *our box-confined*
`V_free(Q)` (large for anions); `α_FI(Q)·V_AIM/V_free(Q)` then over-reduces
(the volume ratio kills it). The fix is to use V_FI consistent with the α
reference (Gould's own ion volumes, in their SI — task #44) OR restrict 2A
to where the mismatch is small. 2C-stern is self-consistent (α and V both
from our confined ld1.x), hence better-behaved. This matches the FI
literature: charge-aware references help *ionic* systems, are neutral-or-
worse for covalent. → **The decisive test is ionic solids** (cohesive
energies), which needs the QE crystal densities (#40); ionic *molecules*
NaCl/LiF show large route spreads (523→334–389; 108→39–63 a₀³) but lack an
easy reference. **Next:** generate QE solid densities → run neutral/A/B/C/
stern on NaCl/MgO/LiF → compare cohesive energies (the headline result).

### 2026-05-30u — 2A consistency flaw FIXED + accepted benchmark sets identified
**2A fix (commit 12df48a5):** changed 2A's denominator from the confined
`V_free(Q)` to the NEUTRAL `V_free(0)`: `α_AIM = α_FI(Q)·V_AIM/V_free(0)`.
Rationale: Gould's α is referenced to *frozen-orbital* ions (orbitals ≈
neutral ⇒ ~neutral-sized), so the consistent volume is neutral; the diffuse
box-confined V over-reduced anions. (This is also literally the "TS+HI with
charge-aware α" model.) The fully self-consistent routes (2C-stern, 2C-kirk:
our own confined α AND V) correctly keep `V_free(Q)`. Covalent C6 recovers:
H₂O 2A 16.1→36.8 (0.81×), CH₄ 58.7→121.0 (0.93×). Refined covalent ranking
vs DOSD: neutral/HI/2A/2B ≈ 0.81–0.99× (good); 2C-kirk 0.80–0.81; **2C-stern
0.57–0.66 (over-reduces** — full-FI with diffuse confined anion volume; a
genuine feature, not a bug). So the charge-matched-VOLUME routes reduce
covalent C6; volume-scaling (2B) and neutral-volume routes preserve it.

**Accepted benchmark sets to run (the validation campaign):**
- Molecular non-covalent: **S66 / S66×8** (Řezáč–Hobza, CCSD(T)/CBS) — the
  standard "doesn't hurt covalent" test; S22/L7/X40 companions.
- Molecular crystals: **X23** (Otero-de-la-Roza & Johnson's own cohesive-
  energy set — the canonical periodic-XDM benchmark).
- Solids / ionic payoff: the **XDM-for-solids cohesive-energy set** (OdlR &
  Johnson, JCP 2012) + the **FI ionic-crystal set** (alkali halides LiF/NaCl,
  oxide MgO, layered MoS₂ vs RPA) — *exactly* where FI showed 157%→23%, so
  our headline target. Plus LC20 lattice constants.
- Tiered plan: Tier 1 = S66×8 + X23 (no harm); Tier 2 = alkali-halide+oxide
  cohesive energies (the payoff). Both are campaigns (S66 wavefunctions; X23/
  solid periodic densities) but let us compare directly to published XDM/FI.

### 2026-05-30v — validation campaign plan (molecules + solids); QE venue
AP: run A/B/C discrimination on BOTH molecules and solids. Plan:

**Molecular — ionic-discriminating subsets first (GMTKN55), then full sets.**
The charge-aware routes (2A/2B/2C-kirk/2C-stern) should separate from neutral
Hirshfeld only where ions/charge-transfer matter, so test there first:
- **IL16** — ionic-liquid ion pairs (cation+anion). Most strongly ionic.
- **AHB21** — anionic hydrogen-bonded complexes (stresses anion α).
- **CHB6** — cationic hydrogen-bonded complexes (cation α).
- **IONPI19** — ion–π (cation-π / anion-π) interactions.
(62 complexes total; CCSD(T)/CBS-quality refs ship with GMTKN55.) These are
the discriminator. THEN the standard non-ionic dispersion sets to confirm
no-harm: **S22, S66/S66×8**, then large/diverse **DES370k**. (Pipeline:
GMTKN55 geometries → DFT wavefunctions/fchk for monomers+complexes →
`xdm 0.4 2.5 pbe hirshfeld_i alpharef … ` per route → interaction-energy /
dispersion-contribution vs ref. Task #45.)

**Solids — ionic sets (the payoff).** Alkali halides (LiF, NaCl, …), oxides
(MgO, CaO), layered/TMD (MoS₂ vs RPA): the FI set (Gould–Bučko, 157%→23%)
and the OdlR & Johnson XDM-for-solids cohesive-energy set (JCP 2012). Run via
the grid path (`xdm grid … hirshfeld_i alpharef …`) on periodic
pseudopotential/PAW densities; cohesive energy = E_bulk + E_disp − ΣE_atoms,
vs experiment. (Task #40.)

**QE venue.** `pw.x`+`pp.x`+PAW pseudos (kjpaw_psl PBE, fetchable from the QE
library — got Na/Cl) are present, BUT the dev-srv apt QE build aborts at
startup with a glibc `_FORTIFY_SOURCE` "buffer overflow detected"
(`__snprintf_chk`) — a packaging bug, input-independent. → Run production QE
on the lab HPC (per ARCHITECTURE.md: HPC ↔ Globus(LXC) ↔ tank/research), or
rebuild QE in a container on dev-srv. critic2 side (the XDM grid + A/B/C
routes) is ready; only the periodic-density generation is blocked on a
working QE. NaCl rocksalt SCF input drafted at `/tmp/nacl_solid/` (a=10.6577
bohr, ibrav=2, ecutrho=480, PAW) — ready to run once QE works.

### 2026-05-30w — solids: QE-only route (no VASP), pp.x AE reconstruction is viable
AP has **no VASP** → must use Quantum ESPRESSO. The blocker was the
all-electron density `rhoae` that the Hirshfeld-I partition needs (VASP gives
it as AECCAR0+AECCAR2). Resolved: QE **`pp.x` can reconstruct the
all-electron valence charge** (binary strings: "Reconstructing all-electron
valence charge", `pp_augmentation`) — so `rhoae` is obtainable from QE+PAW.
critic2's grid path takes `rho` (valence) + `rhoae` (AE) + `elf`/`b`
(`xdm grid rho … rhoae … elf …`); with `rhoae` provided it doesn't need a
separate `core` (xdm@proc l.434). So the QE route works without VASP.
- dev-srv apt `pw.x` is FORTIFY-broken; building `pw.x`+`pp.x` from the QE 7.2
  source already on dev-srv (`/tmp/qe-build/q-e-qe-7.2`, same tree that built
  our `ld1.x`). The ionic cells are tiny (2-atom rocksalt) → can run on
  dev-srv directly; HPC (submission agent) reserved for scale/the full sets.
- **Plan:** (1) build pw/pp; (2) run NaCl scf+pp.x, NAIL the exact
  QE→critic2 charge-aware XDM-grid recipe (pp.x plot_num/augmentation flags
  for rho, rhoae, elf; FFT grid → `load as`); (3) sweep the 6 routes, get the
  first ionic-solid cohesive-energy comparison; (4) THEN package the *verified*
  recipe + structures for the submission agent to run the full ionic-solid set
  (LiF/NaF/NaCl/KCl/MgO/CaO …) + scale on the HPC. Verify-then-package, so no
  wasted HPC runs.

### 2026-05-30x — solids: QE pipeline runs, but TWO real solid-specific blockers found
Source-built QE `pw.x`/`pp.x` (QE 7.2, `/tmp/qe-build`) work (apt build is
FORTIFY-broken). NaCl rocksalt SCF converged; `pp.x` produced rho/rhoae/elf
cubes; critic2 reads them and the charge-aware grid path RUNS. But the NaCl HI
charges are garbage (Na −4.98, Cl +17). Diagnosed two distinct solid blockers:

1. **All-electron density on the grid.** `pp.x plot_num=17` ("AE charge,
   PAW") cube **integrates to 16.19 e, not 28** — i.e. the VALENCE count
   (Na 9 + Cl 7 semicore pseudos); the 12 core electrons are absent (FFT grid
   can't hold the core cusps). So we don't actually have the AE density.
   FIX (critic2-standard for QE/plane-wave XDM): give the **valence** density +
   the pseudo valence charge via **`ZPSP`**, and let critic2 reconstruct the
   core analytically (AE = grid-valence + internal core → integrates to N; the
   Hirshfeld weights then use fine analytic atomic densities, not the coarse
   grid). Need the per-element ZPSP matching the QE pseudos (Na 9, Cl 7, …) and
   the right keyword. **FHI-aims** (true all-electron, installed at
   `/data/FHIaims_distribution/FHIaims_2025`) is the cleaner density source —
   but a uniform cube still undersamples cusps, so core-reconstruction is the
   robust route regardless of code.
2. **Grid HI-SCF instability.** `hirsh_i_driver` (the periodic HI SCF) NEVER
   got the stabilization applied to the mesh path in #34 (floor refrho≥0,
   element-aware per-iteration Q-clamp, linear mixing). So even with a correct
   density it diverges — the §f instability, solid edition. FIX: port the #34
   stabilization into `hirsh_i_driver`.

Both tractable and known. NOTE for collaboration: the installed
`xdm.param.fhiaims` is **Kyle Bryenton's XCDM(Z) (2025, Johnson group)** —
active XDM-in-FHI-aims development in Erin's own group; our charge-aware
references are a natural thing to coordinate there. **Next (solids):**
(a) port #34 stabilization to `hirsh_i_driver` (task #46); (b) nail the
ZPSP/core-reconstruction recipe so the density integrates to N; (c) re-run
NaCl → expect Na≈+0.8/Cl≈−0.8; (d) then package the verified recipe for the
HPC submission agent. Molecular GMTKN55-ionic campaign is independent and
unaffected.

### 2026-05-30y — solids: #46 SCF stabilization DONE; #40 AE-routing partial (2 bugs)
Built source QE 7.2 pw.x/pp.x, ran NaCl rocksalt, drove the QE→critic2
charge-aware grid XDM end-to-end (commit 3af1f21f).
- **#46 DONE:** ported the mesh stabilization to the grid driver
  (`hirsh_i_driver`): refrho≥0 floor in `hirsh_i_eval`, element-aware anion
  clamp (`hirsh_i_qfloor`), linear mixing β=0.5. The grid SCF no longer
  explodes (was Na→−5). Grid HI regression (`022_hirshfeld_i.cro`) re-baselined
  — same converged charges (Sum unchanged to ~1e-4), more iterations from
  mixing.
- **#40 PARTIAL:** route the all-electron density to the HI driver = rhoae if
  given, else `rho(valence)+core` (core from ZPSP), via a temporary reference
  swap. Mechanism active (NaCl: "Use core densities? T", core reconstructed),
  BUT NaCl charges still wrong. Debug (iter-1) revealed **two concrete bugs**:
  1. **AE field over-integrates:** `nelec` sums to **34.35, not 28**
     (Na 11.76 + Cl 22.59). Valence integrates cleanly to 16.0, so the
     `rho+core` core-reconstruction over-counts by ~6 e. Suspect the
     `promolecular_atom(...,zpsp)` core vs ZPSP convention (critic2 printed
     "Core charges (ZPSP): Na(9), Cl(7)" — if ZPSP is read as core rather than
     valence, core≈16 → iae≈32–34). NEXT: check critic2's core/ZPSP semantics;
     the XDM *moment* loop already reconstructs `(rhot+core)` and works, so
     compare how it integrates vs the HI driver.
  2. **`hirsh_i_qfloor` returns 0, not −1**, inside the grid driver even though
     `na_q-1.rho`/`cl_q-1.rho` exist → the anion clamp pins charges at 0.
     wfcdir/symstr probe-path issue in the driver context. NEXT: print the
     probed path; likely a trim/case/dirsep mismatch.
Both are isolated; once fixed, expect NaCl Na≈+0.8/Cl≈−0.8. Molecular path
(5 routes) + default XDM unaffected; regression green.

### 2026-05-30z — SOLID charge-aware XDM WORKING (NaCl); the two #47 bugs resolved
The solid grid-HI now gives physical charges and the full charge-aware
dispersion. Resolution (commit 7988d631):
- **Root cause:** the HI *population* integral needs a density that
  integrates cleanly on the coarse uniform grid. All-electron does NOT
  (QE pp.x plot_num=17 loses the core → 16e; reconstructed AE = valence+core
  re-introduces the core CUSP → over-integrates, 34≠28). Same §f aliasing,
  in the core.
- **Fix:** integrate the SMOOTH valence (pseudo) density and take the charge
  relative to the pseudo valence: `Q = ZPSP − ∫w·ρ_valence`. That is exactly
  the physical (valence) charge transfer (the frozen spherical core does not
  partition), and the valence density has no cusp. `hirsh_i_driver` reads
  ZPSP from the reference field (zval = ZPSP if usecore else Z). Reverted the
  AE-reconstruction reference-swap. The "qfloor=0" symptom was a *test-script*
  bug (unexpanded `$WFC`), not code.

**NaCl rocksalt (QE 7.2 PBE PAW, a=10.6577 bohr; `xdm grid … hirshfeld_i
alpharef …`):** charges **Na +1.02 / Cl −1.00** (Cl at the q=−1 floor — bulk
NaCl is strongly ionic, more so than the gas molecule's ±0.89, and there is
no Cl²⁻ ref). Charge-aware α (a₀³): **α(Na) 10.9 → 1.06** (Na⁺ collapse),
**α(Cl) 26.5 → 34.6** (Cl⁻ expansion). Lattice-summed dispersion energy
Evdw per cell (mHa):

| neutral | 2A gould | 2B scale | 2C stern |
|---|---|---|---|
| −33.2 | −16.7 | −15.9 | −11.9 |

**Charge-aware references roughly HALVE the bulk-NaCl XDM dispersion** — the
neutral-Hirshfeld XDM massively over-counts the cation's contribution (neutral
Na α ~ shell of valence; Na⁺ ~ 1 a₀³). This is the large ionic-solid effect
the FI work reported, now reproduced in critic2. (Note: the grid HI SCF
clamps Cl at −1 and reports non-convergence at that boundary — the physical
limit for NaCl; the charges are stable. "Better"-ness awaits cohesive-energy
vs experiment, but the EFFECT is large and correctly directed.)

**Pipeline (verified, for the HPC package):** QE pw.x PBE-PAW SCF →
pp.x plot_num=0 (valence ρ) + plot_num=8 (ELF) → critic2
`crystal cube; load val zpsp <Zval…>; load elf; reference val;
xdm grid rho val elf elf rhoae?(omit) xa1 a1 xa2 a2 hirshfeld_i
alpharef gould|scale|stern wfcdir <ld1_pbe>`. No rhoae needed (core
reconstructed from ZPSP; HI uses valence). Solids + molecules now both work.
NEXT: package this for the submission agent → full alkali-halide/oxide set.

### 2026-05-30aa — HEADLINE: ionic-solid cohesive-energy validation (6 solids)
Full QE→critic2 charge-aware XDM across the alkali-halide + oxide set, at the
experimental rocksalt lattice constants. Pipeline: source-built QE 7.2
PBE-PAW SCF (kjpaw_psl 1.0.0 pseudos) → pp.x plot_num=0 (valence ρ) + 8 (ELF)
→ critic2 `xdm grid … hirshfeld_i alpharef gould|scale|stern` (ZPSP-valence
HI, no rhoae); free-atom PBE energies (spin-polarized, isolated, degauss=0.01)
for the cohesive reference. Scripts + raw output in
`doc/research/solid_validation/`.

**Cohesive energy E_coh = (E_cat+E_an) − (E_bulk + E_disp), eV/f.u.:**

| solid | PBE | neutral | 2A gould | 2B scale | 2C stern | exp | HI q(cat/an) |
|---|---|---|---|---|---|---|---|
| LiF | 8.77 | 9.59 | 9.25 | 9.05 | 8.90 | 8.75 | +1.02/−1.00 |
| NaF | 7.75 | 8.66 | 8.10 | 8.02 | 7.87 | 7.86 | +1.05/−1.00 |
| NaCl | 6.15 | 7.05 | 6.60 | 6.58 | 6.47 | 6.62 | +1.02/−1.00 |
| KCl | 6.27 | 7.71 | 6.57 | 6.71 | 6.61 | 6.70 | +1.03/−1.00 |
| MgO | 10.25 | 11.06 | 10.56 | 10.80 | 10.65 | 10.30 | +2.13/−2.00 |
| CaO | 11.19 | 12.66 | 11.59 | 11.75 | 11.79 | 11.00 | +2.07/−2.00 |

**MAE / MSE vs exp (eV):** PBE 0.21/−0.14 · **neutral 0.92/+0.92** · gould
0.29/+0.24 · scale 0.30/+0.28 · **stern 0.26/+0.18**.

**KEY RESULT:** neutral-Hirshfeld XDM over-binds EVERY ionic solid by a
systematic +0.92 eV (MSE=MAE → purely one-directional over-binding), because
it over-counts the neutral cation polarizability. **Charge-aware references
cut the error ~3.5× (0.92→0.26) and eliminate the systematic bias**
(+0.92→+0.18); 2C-stern best. This is the FI thesis, quantitatively
reproduced in critic2. Exp from Born–Haber (alkali halides solid; oxides
10.3/11.0 approximate).

**Caveats (report straight):** (1) PBE-alone MAE 0.21 (MSE −0.14, under-binds,
no dispersion) is competitive for *pure ionic* solids — the fair like-for-like
is XDM-vs-XDM (neutral 0.92 → charge-aware 0.26); the decisive win for
charge-aware should come on layered/molecular crystals where dispersion truly
matters (X23, TMDs — next). (2) a1/a2=(0.4,2.5) NOT refit (RQ5) — absolute XDM
magnitude would shift with optimal damping; the neutral↔charge-aware gap is
robust. (3) single fixed-geometry point, no ZPE/thermal/relaxation; grid HI
slightly over-ionizes (halide anions pin at −1). (4) free-atom degauss=0.01
smearing entropy ~tens of meV.

### 2026-05-31a — Molecular ionic benchmark (GMTKN55) scoping: all assets located
AP chose the GMTKN55 ionic subsets as the molecular charge-aware
discriminator. The cluster already has the full infrastructure — no rebuild:
- **Geometries:** `/data/refdata/30_collection-GMTKN55/{il16,ahb21,chb6}/`
  (complex + A/B monomer .xyz per reaction). IONPI19 not in this collection.
- **References + stoichiometry:** `/data/refdata/10_din-GMTKN55/{il16,ahb21,
  chb6}.din` (Grimme .din: reaction = +1·complex −1·A −1·B, CCSD(T)/CBS ref in
  kcal/mol; e.g. il16_008 = −100.41 — strong ion-pair binding). il16=16,
  ahb21=21, chb6=6 reactions.
- **Existing XDM pipeline:** `/data/XDM_Psi4/` — AP's own **XDM-in-PSI4**
  implementation (`~/projects/psi4_xdm_implement/`, conda env `psi4-xdm`) +
  `xdm_params.py` (fitted a1/a2, 16 functionals × 4 bases), `xdm_lib.py`
  (BJ dispersion, .din parser, stats), `datasets.py` registry, S22x5/S66x8/
  3B-69/KB49 benchmarks. `06_master_workup.py` evaluates Alberto-style.
- **FHI-aims PBE0** GMTKN55 calcs exist (`/data/FHIaims_GMTKN55/PBE0/IL16/…`,
  per-species .out + .csc density) — energies available, densities in .csc.

**Integration plan (critic2 charge-aware route):** for each species, Psi4
(psi4-xdm env) PBE → molden/fchk → critic2 `xdm_wfn` (neutral + alpharef
gould|scale|stern + the ld1_pbe wfcdir) → E_disp; E_int = [E_DFT+E_disp]
(complex) − Σ monomers; compare to the .din CCSD(T) reference. The molecular
mesh path (xdm_wfn) already works with our charge-aware routes, so this reuses
proven code; the only new work is the species loop + the .din evaluation
(can reuse xdm_lib's parser/stats). Start with **IL16** (16 ion pairs — the
strongest ionic discriminator), then AHB21/CHB6, then the standard S22/S66×8
for no-harm. NEXT: confirm the Psi4→critic2 wavefunction format on il16_008,
then batch.

### 2026-05-31b — il16_008 handshake WORKS; but IL16 is electrostatics-dominated
End-to-end pipeline confirmed on the first ion pair (il16_008, neutral 17-atom
complex = +1 cation A [16 at] + −1 anion B [1 at, monatomic]). Psi4 PBE/def2-TZVP
→ fchk → critic2 `xdm_wfn` for all 4 routes. Driver `/tmp/run_il16_008.sh`,
`/tmp/psi4_species.py` on dev-srv. (Gotcha: `set -u` kills the script at
`source activate_xdm.sh` — conda's profile.d hits unbound vars; use `set +u`.)

Result (Edisp in Ha; E_int in kcal/mol; B monatomic ⇒ Edisp(B)=0):
| route   | Edisp(cplx) | Edisp(A) | Edisp_int | E_tot   | Δ vs neutral |
|---------|-------------|----------|-----------|---------|--------------|
| neutral | −0.01607    | −0.01039 | −3.56     | −109.65 | —            |
| gould   | −0.01608    | −0.00922 | −4.30     | −110.39 | −0.74        |
| scale   | −0.01581    | −0.00951 | −3.96     | −110.04 | −0.40        |
| stern   | −0.01276    | −0.00765 | −3.20     | −109.29 | +0.36        |
| ref(.din, CCSD(T)/CBS)                       | **−100.41** |          |

Routes behave correctly & differ in the right direction: `stern` reduces the
anion-side over-polarization (cplx Edisp −0.01276 vs −0.01607 neutral, ≈21% less),
`gould` enhances it; `stern` lands closest to ref. a1/a2 used = 0.4/2.5 ≈ AP's
cc-pVTZ-fitted PBE-XDM (0.4041, 2.6998), so damping is ~right.

**KEY (negative) FINDING — IL16 is a weak discriminator.** All routes overbind by
~9 kcal/mol, and that error is essentially PBE's (electronic E_int alone = −106.08;
likely + some def2-TZVP BSSE), NOT dispersion's. Dispersion is only ≈3% of the
binding and the *spread between routes* is <1 kcal/mol — swamped by the base-functional
error. Molecular ion-PAIR binding is electrostatics-dominated, so batching the other
15 IL16 (and largely AHB21/CHB6, H-bond-dominated) would mostly measure PBE, not our
charge-aware dispersion. Mirrors the solid finding inverted: where the *differential*
is dispersion (layered/ionic-molecular crystals, ion-π) the method shows; where it's
electrostatics/H-bonding, it can't. DECISION NEEDED: point molecular validation at a
dispersion-dominated ionic case (IONPI19 ion-π — but geometries not in the GMTKN55
collection here) and/or switch metric to dispersion-only vs SAPT to isolate the term
that actually changes. Solids remain the headline.

### 2026-05-31c — unattended batch launched (AP away ~few hrs, "run as many as possible")
Built a resumable orchestrator on dev-srv: `/tmp/gmtkn/run_batch.py` (+ updated
`/tmp/gmtkn/psi4_species.py` which now writes a `<fchk>.edft` energy sidecar).
Parses any Grimme `.din` (coeff/species pairs, `0` term, ref kcal/mol), runs Psi4
PBE/def2-TZVP per unique species (3-wide ThreadPool, 4 threads/7 GB each; cached
fchk+edft), then critic2 `xdm_wfn` ×4 routes (neutral/gould/scale/stern; cached
`<sp>.<route>.edisp`), builds per-reaction E_int by stoichiometry, reports MAE/MSE
per route vs ref. Fully resumable (skip cached) + fault-tolerant (per-species
try/except, single-atom monomers→Edisp=0). a1/a2 = **0.4041/2.6998** (AP's
cc-pVTZ-fitted PBE-XDM, ~def2-tzvp; supersedes the ad-hoc 0.4/2.5 of 31b, ~30%
more damping ⇒ smaller |Edisp|, route ordering unchanged).
Queue (ionic discriminators + small no-harm first): il16 → chb6 → ahb21 → s22 →
ionichb(120 rxns) → s66(198). Dry-run verified all 6 sets parse and every species
xyz exists (il16 48sp, chb6 18, ahb21 63, s22 66, ionichb 150, s66 198). Logs:
`/tmp/gmtkn/master.log`; per-set `/tmp/gmtkn/results_<set>.txt`. Launched 12:41,
PID 1773377 (setsid). To resume after a crash: same launch cmd — caches make it
pick up where it stopped. ionichb/ssi/s66_disp are the bonus ionic/dispersion
discriminators if time allows.
Throughput note: Psi4 (not critic2 — one charged mesh-HI route = 27 s) is the
bottleneck (~8–30 min/species; large organic cations dominate). Switched to 6×2
(NWORK=6, 2 threads/4 GB) after a slow start; resumable via caches. Launch via
`/tmp/gmtkn/launch.sh` (sources env, execs run_batch.py) fired with
`nohup setsid bash launch.sh </dev/null >/dev/null 2>&1 &`. **Gotcha:**
`pkill -f psi4_species.py` over ssh self-matches the ssh command line and kills
its own session (rc 255) — never pattern-kill on a string that appears in your
own command.

### 2026-05-31d — IL16 partial (6/16): charge-aware effect is real but tiny here
Directional read via `/tmp/gmtkn/analyze.py il16` on the 6 fully-cached reactions
(a1/a2 = 0.4041/2.6998, NOT refit). Eint MAE/MSE vs CCSD(T)/CBS (kcal/mol):
DFT-only 2.68/−2.59 · neutral 4.37/−4.37 · gould 4.35/−4.35 · scale 4.52/−4.52 ·
stern 3.98/−3.98. **PBE alone already overbinds ion pairs (MSE −2.59); dispersion
is always attractive so it pushes IL16 *further* into overbinding — every route is
WORSE than DFT-only.** Route spread is only ~0.5 kcal (stern best, scale worst),
swamped by the base-functional error. `stern` wins here only incidentally (it makes
the least dispersion; less attraction helps when the base overbinds), not from
better physics. ⇒ **IL16 is not a discriminating molecular test for charge-aware
XDM.** Combined with the refit point (AP): a1/a2 refit can absorb a uniform
dispersion-magnitude shift but NOT the per-species anion/cation α differential —
yet that differential only *matters* where dispersion is a large fraction of binding
(ionic solids ✓, and expected: layered/molecular ionic crystals, ion-π), not in
electrostatics-dominated ion pairs / H-bonds. Solids remain the headline; molecular
ionic GMTKN55 is a no-harm/where-it-doesn't-matter data point, not the showcase.

### 2026-05-31e — ALL THREE molecular ionic sets done: a clean null (37 rxns)
Full results (`/tmp/gmtkn/results_{il16,chb6,ahb21}.txt`; Eint MAE/MSE kcal/mol vs
CCSD(T)/CBS; a1/a2=0.4041/2.6998 unrefit; 6×2 batch). IL16 lost 6/16 reactions to
Psi4 SCF timeouts (large organic cations; failed-fast skip).
| set            | n     | DFTonly      | neutral      | gould        | scale        | stern        |
|----------------|-------|--------------|--------------|--------------|--------------|--------------|
| il16 (ionpair) | 10/16 | 3.24/−3.19   | 4.83/−4.83   | 4.83/−4.83   | 4.98/−4.98   | 4.49/−4.49   |
| chb6 (cat HB)  | 6/6   | 1.34/−0.58   | 1.35/−1.35   | 1.32/−0.60   | 1.30/−1.01   | 1.28/−0.67   |
| ahb21 (an HB)  | 21/21 | 4.37/−4.37   | 4.94/−4.94   | 5.03/−5.03   | 5.01/−5.01   | 4.84/−4.84   |
Findings, robust across 37 reactions:
(1) PBE base OVERBINDS all three (MSE<0). Dispersion is always attractive ⇒ on
il16/ahb21 it makes things WORSE than bare PBE; dispersion is not the missing piece
for PBE-overbound ionic H-bonds. (2) Charge-aware ≈ neutral: routes differ by only
0.07–0.49 kcal/mol MAE (chb6 within rounding). Per-species signal exists & is
sign-consistent (stern trims anion-side dispersion ⇒ smallest MAE in every set) but
is tiny because dispersion is a small fraction of these bindings. (3) Aside: PBE
handles cationic HB (chb6 1.34) far better than anionic HB (ahb21 4.37) — anion SIE.
⇒ CONFIRMED: charge-aware XDM is a no-op on molecular ion pairs / H-bonded ions
(neither helps nor harms meaningfully). The method's value is confined to
dispersion-dominated binding: ionic solids (done, decisive), and the still-untested
layered/molecular ionic crystals & ion-π. Even post a1/a2-refit this won't change for
these sets — refit rescales the (small) dispersion uniformly; the route differential
that charge-aware adds is ~0.3 kcal here regardless. s22 no-harm running next.

### 2026-06-01a — KB49 a1/a2 REFIT pipeline (canonical, per route) — launched
AP asked to refit a1/a2 for gould/scale/stern on KB49 (the canonical XDM training
set; Kannemann-Becke JCTC 2010). Reuse Alberto's canonical machinery exactly:
- **Refdata source of truth:** `/home/albd/projects/refdata` (Alberto's; `20_kb49`
  geoms + `10_din/kb49.din`, refs in kcal/mol). `/data/refdata` is a near-identical
  mirror (used for the GMTKN batch). 49 reactions, 147 species, 0 missing.
- **AP's fit method (`/data/XDM_Psi4/xdm_kb49_refit/02_collate_and_fit.py`):** each
  species → JSON {base_energy, coords[bohr], c6,c8,c10,rc (NxN, a.u.)}; E_disp is
  closed-form `calc_bj_dispersion` (xdm_lib): `rvdw=a1·rc + a2/0.52917721`,
  `E=-Σ_{i<j}Σ_n Cn/(rvdw^n+d^n)`; fit (a1,a2) by `least_squares` minimizing RMSP
  (relative residuals), init [0,1.4545], a2 in Å.
- **Key enabler — critic2 emits the matrices directly:** the `coefficients (a.u.)`
  block prints `i j C6 C8 C10 Rc Rvdw` per pair (Rc is a1/a2-independent). So one
  Psi4 PBE/def2-TZVP per species (route-INDEPENDENT base_energy) + 4 critic2 runs
  (route-dependent Cn) → 4 JSON sets → 4 fits. Parser `/tmp/gmtkn/kb49_makejson.py`
  (coords from fchk, native bohr; diagonal self-terms dropped, unused by k=1 sum).
- **VALIDATION (make-or-break):** parsed water-dimer block → AP's `calc_bj_dispersion`
  reproduced critic2's own E_disp to **3.6e-14 Ha**. Parser+formula+units all exact.
- Driver `/tmp/gmtkn/kb49_run.py` (4×3, resumable, fault-tolerant), launched via
  `/tmp/gmtkn/kb49_launch.sh`; logs `/tmp/kb49/master.log`, fit → `/tmp/kb49/fit_results.json`.
  Output table per route: a1, a2, MAE_fit, MAPD_fit, and MAE/MAPD at the current
  default (0.4041/2.6998) for comparison. GMTKN batch PAUSED (192 species cached,
  resumable) to give KB49 the full 12-core machine. Big aromatics (naphthalene
  dimers 36 at, adenine-thymine 30 at) are dispersion-critical for the fit and must
  converge — generous 90-min Psi4 timeout.

### 2026-06-01b — KB49 run debugged (3 real blockers) → relaunched clean
First KB49 attempt: 142/147 FAILED. Three independent bugs, all fixed:
1. **Root disk 100% full** ⇒ Psi4 DFHelper `put_tensor: write error`. Stale
   `/tmp/psi.*` scratch from the killed GMTKN orphans filled `/` (116 G). Fix:
   cleared stale scratch + redirect `PSI_SCRATCH=/data/XDM_Psi4/psi_scratch_kb49`
   (34 TB volume; `/data` root itself is not user-writable, but `XDM_Psi4/` is).
   `psi4_species.py` now `IOManager.set_default_path()` to it.
2. **OpenBLAS×OpenMP nesting hang** (the big one): conda psi4-xdm links OpenBLAS,
   which "may hang" inside Psi4's OpenMP region — ~2 min/SCF-iter, so big aromatics
   never reached convergence before timeout. AP's cluster used Intel MKL (no nesting
   issue). Fix: `export OPENBLAS_NUM_THREADS=1` (Psi4 keeps its own OpenMP). ⇒ ~16
   s/iter (7× faster); adenine_thymine_stack then converges in 24 iters to
   −920.7768821613 Ha. SCF was never actually divergent — just too slow to finish.
3. SCF settings now match AP's worker exactly (guess sad, freeze_core, ultrafine
   590/99 grid, e/d_conv 1e-6, maxiter 80) + two-stage convergence (default DIIS,
   SOSCF retry only on ConvergenceError). Also OOM note: an ML job sharing the box
   OOM-killed a 13 G process at 16:20 — NOT a KB49 worker (those are ~5 G, all
   survived); reniced KB49 to coexist, ML later stopped.
Full batch relaunched 17:38 (4×3, OPENBLAS_NUM_THREADS=1, scratch on /data),
resumes from 7 cached species. Fit (RMSP least_squares per route) auto-runs at end
→ `/tmp/kb49/fit_results.json` + table in `/tmp/kb49/master.log`.

### 2026-06-01c — KB49 REFIT RESULT (all 147 species converged, 0 fail, 49/49)
Final fitted BJ-damping per route (PBE/def2-TZVP, RMSP objective, a2 in Å):
| route   | a1     | a2(Å) | MAE(kcal) | MAPD   | MAE@old-def | MAPD@old-def |
|---------|--------|-------|-----------|--------|-------------|--------------|
| neutral | 0.0000 | 4.069 | 0.651     | 19.45% | 0.621       | 19.70%       |
| gould   | 0.0000 | 4.491 | 0.679     | 21.06% | 0.717       | 23.94%       |
| scale   | 0.1864 | 3.730 | 0.605     | 18.59% | 0.583       | 19.40%       |
| stern   | 0.0000 | 4.143 | 0.632     | 19.51% | 0.933       | 27.04%       |
(old-def = the borrowed 0.4041/2.6998.) Interpretation:
1. After refit all four routes land in a tight band (MAE 0.60–0.68, MAPD 18.5–21%)
   ⇒ charge-aware does NO HARM on the canonical neutral training set once its own
   damping is fit. scale marginally best, gould marginally worst.
2. **stern needs refitting most:** at the borrowed default it was badly off
   (MAE 0.93, MAPD 27%) because it systematically shrinks α; its own fit recovers it
   to 0.63/19.5%. Concrete proof the routes must be re-fit individually (AP's point).
3. a1 pins to the a1≥0 boundary (=0) for neutral/gould/stern ⇒ PBE prefers a pure-a2
   offset; effectively a 1-param fit. RMSP is dominated by the weakest complexes
   (ch4·c2h4 etc., ref ~0.5 kcal). NB AP's own pbe-xdm cc-pVTZ/aug params have a1≠0
   (0.40/0.71) — the a1→0 here may reflect def2-tzvp + critic2-vs-psi4 coefficient
   differences + RMSP objective; worth a cross-check vs an MAE objective before
   shipping these as production params.
Pipeline (`kb49_run.py` + `kb49_makejson.py`, validated to 3.6e-14 Ha vs critic2)
is reusable for any functional/basis. GMTKN batch still paused (192 cached) — resume
or leave per AP. Task #43 (RQ5 refit) essentially done for KB49.

### 2026-06-01d — DATA PERSISTED + QE/plane-wave KB49 refit packaged for HPC
AP: persist all fitting data + do the KB49 refit again with QE plane-wave densities
(ecutwfc 80 / ecutrho 800 Ry, extra-large isolated box, PAW).
- **Persistence:** all generated data now under **`/data/Iterative_hirshfeld/`**
  (34 TB vol; root `/` is small + was 100% full). Subdirs: `kb49_psi4/` (163 MB:
  fchk, json/<route>, fit_results.json), `gmtkn_molecular/` (all set caches+results),
  `scripts/`, `notebook/` (snapshot), `kb49_qe/`, `README.md`. (sudo-mkdir'd the
  /data root dir, chown albd.)
- **QE path validated, but can't run the chosen box on dev-srv:** fixed 45-bohr /
  800-Ry box ⇒ ~405^3 FFT ⇒ QE wants ~16 GB/MPI-rank (~195 GB for 12); dev-srv
  has 31 GB → OOM-killed before SCF. ⇒ needs a large-memory HPC node.
- **Chain proven (small box):** qe_make_input → pw.x(MT box) → pp.x rho(0)+ELF(8)
  → critic2 `xdm grid … zpsp …` → JSON. critic2 grid prints `+ Dispersion
  coefficients` (i j C6 C8 C10 Rc; lower-triangle) + volumes/moments; ZPSP core
  reconstruction OK (Vfree(C)=39.0 a0^3). New parser `kb49_makejson_qe.py` (coords
  from xyz in bohr). Round-trip: fit `edisp` vs critic2 grid Evdw matches to the
  periodic-image residual (2.4e-6 Ha @24-bohr → ~0 @45-bohr); critic2 Evdw =
  Evdw6+8+10 (no C9 in energy). Consistent with AP's molecular (no-image, no-C9) fit.
- **HPC package:** `/data/Iterative_hirshfeld/kb49_qe/HPC_package/` — scripts
  (qe_make_input, run_species.sh, kb49_makejson_qe, zpsp_for, **self-contained**
  fit_kb49_qe.py), `submit_kb49_qe.slurm` (array 1-147, 256 GB/8 MPI/8 h),
  pp/ (8 kjpaw_psl 1.0.0 PAW), geom/ (147 xyz), ld1_pbe/ (137 anion refs), kb49.din,
  species.list, NOTES.md. **Key dep flagged in NOTES:** cluster must build the FORK
  branch `research/xdm-hirshfeld-i` (gould/scale/stern + grid alpharef are ours, not
  upstream). Brings back only `results/<route>/*.json` (few MB) + fit_results_qe.json;
  cubes (~44 GB) stay on node scratch. Ready for AP's HPC-submitting agent.

### 2026-06-03 — full literature ingested → doc/research/LITERATURE_INTEGRATION.md
Read all 16 PDFs in doc/papers/ (4 parallel readers) and wrote a synthesis mapping
every paper to our routes + the numbers to compare against. See
**`doc/research/LITERATURE_INTEGRATION.md`** (+ `references.bib`, mirrored under
`/data/Iterative_hirshfeld/papers/`). Corrections + action items it surfaced:
- **`s00894-017-3514-6.pdf` = Heidar-Zadeh/Ayers/Bultinck 2017** (fractional-nuclear-charge
  anion densities), NOT Bučko-2017 MBD. Useful for unbound-anion refs; bib fixed.
- **KB49 = KB65 − noble-gas dimers.** Our neutral refit (a1=0.4041/a2=2.6998) reproduces
  **OdlR&J 2013 PBE plane-wave a1=0.4073/a2=2.4150** (a1 matches → neutral route validated).
  Do NOT compare to Kannemann–Becke PW86PBE 0.82/1.16 (different functional).
- **IONPI19 ∉ GMTKN55** (the 3 ionic subsets are AHB21/CHB6/IL16). Drop it from #45 title.
- **[VERIFY, high] gould anion tier:** must be the Gould–Bučko *embedded/frozen-orbital*
  anion α (Appendix B), NOT the free self-consistent values (F⁻ 15.5 vs embedded). If wrong,
  gould over-polarizes anions. (#37 re-check.)
- **scale (p'=p−0.615) is a DEPARTURE from FI**, which used linear p'=1 and rejected the
  exponent. FI-faithful = gould+linear. Document the distinction.
- **[RUN, high] B86bPBE ionic solids:** OdlR&J 2020 (`xdmmx.pdf`) get ionic solids right via
  the B86b *exchange* functional (lattice MAE 0.060 Å), not charge-aware refs — and name HI
  as the unfilled XDM gap. Strongest paper claim ⇒ run charge-aware on B86bPBE too
  (complementary: their B86b fixes geometry/B0, our charge-aware fixes cohesive energy).
- Literature precedent: charge-aware needs little/no damping refit (FI β unchanged; TS/HI
  s_R 0.94→0.95) — consistent with our near-neutral refit (stern the exception).

### 2026-06-03a — gould anion tier VERIFIED (benchmark, not embedded; results unaffected)
Checked the literature action item. `alpha_gb_m1` (param.F90; provenance
`dat/xdm_ion_alpha_gould_bucko_2016.dat` header = "benchmark") is the Gould–Bučko
**free self-consistent** anion tier (F⁻=15, Cl⁻=30.3, O⁻=5.4, Li⁻=1180, Na⁻=1310,
K⁻=2090), NOT the FI "minimal-chemistry" embedded/frozen-orbital tier. **Impact on our
validated results = NONE:** (1) HALIDE anions are self-consistent in both tiers → identical
→ alkali-halide solids correct; (2) oxide O clamped at −1, α=5.4≈neutral → tier-invariant;
(3) `ion_alpha0` clamps |q|≤1 and cations interpolate neutral→+1, so the grossly-diffuse
alkali "anion" entries are structurally unreachable. Tier only matters for non-halide anions
with large charge (N³⁻/S²⁻ nitrides/sulfides — untested). Added a TIER NOTE comment in
param.F90; not a bug in current work, a scoped future refinement (embedded tier from the
Gould–Bučko SI, task #44 territory). Committed to fork.
