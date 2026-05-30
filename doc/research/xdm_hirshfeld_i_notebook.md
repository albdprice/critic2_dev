# Charge-aware (iterative-Hirshfeld) reference handling for XDM dispersion — lab notebook

**Project:** Does replacing neutral-Hirshfeld partitioning with iterative
Hirshfeld (Hirshfeld-I) in the Becke–Johnson XDM dispersion model improve
it, and what is the correct treatment of the free-atom reference
quantities?

**Maintainers:** A. Price (critic2 fork), with A. Otero-de-la-Roza (XDM
author) consulted. Worklog kept lab-notebook style (dated entries at the
bottom); top sections are the living "paper skeleton."

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

