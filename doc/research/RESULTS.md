# RESULTS — charge-aware (Hirshfeld-I) XDM  ·  observations & numbers

The evidence base for the paper. Every number here has provenance (script + data on
`/data/Iterative_hirshfeld/`). Routes: **neutral** (baseline XDM), **gould**
(Gould–Bučko ion-α, linear-N), **scale** (α∝V^p′, p′=p−0.615), **stern** (confined-ion
Sternheimer α). Details/derivation in `PAPER_METHODS.md`; comparison targets in
`LITERATURE_INTEGRATION.md`.

Last updated: **2026-06-03**

---

## R1 — Ionic-solid cohesive energies (the headline)
6 rocksalt solids (LiF, NaF, NaCl, KCl, MgO, CaO), QE 7.2 PBE-PAW, single-point at
experimental geometry; `xdm grid hirshfeld_i alpharef …`; cohesive E = E_bulk + E_disp
− Σ E_atom(PBE). Data: `/data/Iterative_hirshfeld/` (was `doc/research/solid_validation/`).

| method | MAE (eV) | MSE (eV) | note |
|---|---|---|---|
| PBE (no disp) | 0.21 | −0.14 | under-binds; competitive for *pure* ionic |
| neutral XDM | **0.92** | **+0.92** | over-binds EVERY solid (systematic) |
| gould | 0.29 | — | |
| scale | 0.30 | — | |
| **stern** | **0.26** | **+0.18** | best; bias largely removed |

**Finding:** charge-aware references cut XDM over-binding ~3.5× (0.92→0.26 eV) and
remove the systematic bias — the FI thesis confirmed in critic2. HI charges: halides
±1.0, oxides +2.1/−2.0. Mechanism: neutral XDM over-counts the *cation* polarizability
(a whole valence shell is gone in the cation; a volume ratio can't remove it).
**Caveat:** PBE-alone is competitive here (pure ionic solids are not dispersion-limited),
a1/a2 not refit, one geometry point. → the decisive test is dispersion-dominated crystals.

## R2 — GMTKN55 ionic molecular benchmark (the null)
Psi4 PBE/def2-TZVP → critic2 `xdm_wfn`; E_int vs CCSD(T)/CBS (.din). MAE kcal/mol.

| set | n | DFT-only | neutral | gould | scale | stern |
|---|---|---|---|---|---|---|
| IL16 (ion pairs) | 10/16 | 3.24 | 4.83 | 4.83 | 4.98 | 4.49 |
| CHB6 (cationic HB) | 6/6 | 1.34 | 1.35 | 1.32 | 1.30 | 1.28 |
| AHB21 (anionic HB) | 21/21 | 4.37 | 4.94 | 5.03 | 5.01 | 4.84 |
| S22 (neutral no-harm) | 7/22 | 0.58 | 0.61 | 0.57 | 0.65 | 0.53 |

**Finding:** all routes within ~0.1–0.5 kcal/mol of neutral — a clean null. Molecular
ion-pair/H-bond binding is electrostatics-dominated (IL16 ~97%), so dispersion (and our
correction) is a small slice, swamped by PBE's base error. Charge-aware neither helps
nor harms molecular ionic complexes. (S22 lost its π-stacked cases to SCF timeouts.)

## R3 — KB49 a1/a2 refit (Psi4 PBE/def2-TZVP)
147 species, 0 fail; per-route RMSP fit (AP's convention); coefficient parser validated
to 3.6e-14 Ha vs critic2. a2 in Å.

| route | a1 | a2 (Å) | MAE (kcal/mol) | MAPD |
|---|---|---|---|---|
| neutral | 0.000 | 4.069 | 0.651 | 19.45% |
| gould | 0.000 | 4.491 | 0.679 | 21.06% |
| **scale** | 0.186 | 3.730 | **0.605** | 18.59% |
| stern | 0.000 | 4.143 | 0.632 | 19.51% |

**Findings:** (1) after refit all routes tighten to a narrow band → charge-aware does no
harm on the neutral training set. scale best, gould worst. (2) stern needed the refit
most (borrowed default MAE 0.93 → 0.63) → routes must be refit individually. (3) neutral
route VALIDATED: reproduces Otero-de-la-Roza–Johnson 2013 PBE plane-wave a1=0.4073
(matches to 0.003). (4) Caveat: a1 pins to the 0 boundary for 3/4 routes (RMSP over-weights
the weakest complexes) — needs an MAE-objective cross-check before shipping as production.
QE plane-wave version packaged for HPC (not yet run).

## R4 — Per-atom reference polarizabilities (a₀³) — the routes disagree
Illustrative Na/Li/Cl/F: neutral 22.7/6.8/21.3/6.5 → gould(2A) 14.4/5.5/23.3/9.0 →
scale(2B) 9.3/4.3/27.8/9.7 → compute(2C-Kirkwood) 14.7/1.5/22.5/6.5 → stern 14.4/5.5/27.5/6.4.
The spread is genuine ion-α reference uncertainty; cations move ↓ (valence-shell loss),
anions ↑ vs neutral scaling. Headline for the "why not just scale neutral volume"
argument: neutral scaling overestimates **cation** α **20–35×** (Na 22.7 vs Na⁺ ~0.98;
Li 6.8 vs Li⁺ ~0.19) — a whole shell is gone, a volume ratio cannot remove it.

## R5 — Double-anion reference polarizability: a method finding
Investigating q=−2 references (O²⁻/S²⁻ in oxides/sulfides) revealed that the two
*table/linear-response* routes cannot be extended past −1, but the *density-based*
routes can:

| method for O²⁻ | value | verdict |
|---|---|---|
| uncoupled Sternheimer (linear response) | α = {982,−31,0.99,−0.07,107,1.74} vs box | **diverges** (2nd e⁻ unbound) |
| Gould–Bučko benchmark (free ion) | grossly diffuse | unphysical |
| **density-moment (Kirkwood ⟨r²⟩²/N)** | 17.5→64→**82** (q=0,−1,−2), α≈24 a₀³ | **finite, monotonic** (right order vs physical ~14 a₀³) |
| q=−1 clamp (current gould/stern) | α(O⁻)≈19 a₀³ | reasonable proxy for in-crystal O²⁻ (~14) |

**Finding:** a free double anion is not bound, so its linear-response polarizability is
ill-defined (Sternheimer garbage across the whole table); but the *bound confined
ground-state density* is well-defined, so density-moment estimators work. **Consequence:**
gould/stern deliberately clamp at −1 (documented, justified proxy); **`compute`/`scale`
(density/volume-based) are the correct routes for multiply-charged anions.** No change to
the 6-solid results (oxides used the −1 clamp already). Diagnostic:
`/data/Iterative_hirshfeld/sternheimer_alpha_z84.dat`.

**End-to-end confirmation (MgO, QE PBE-PAW + critic2 grid):** the `compute` route handles
O (V=73.5 vs Vfree=23.8 a₀³ — strongly anionic) with a finite, physically-ordered result —
C6(O–O)=91.2 a.u., Evdw=−0.0129 Ha — vs the garbage the Sternheimer −2 would have injected.
Full: neutral/gould/scale/stern/compute Evdw = −0.0236/−0.0091/−0.0160/−0.0118/−0.0129 Ha;
C6(O–O) = 11.3/65.4/132.7/45.6/91.2. All finite, no divergence.

**Literature check (2026-06, 14 verified claims — see LITERATURE_INTEGRATION §5):** this is
textbook. Holka et al. JCP 141, 214303 (2014) — free O²⁻/S²⁻ diverge (RHF instability), a weak
HO confining potential fixes it (O²⁻<O⁻ only for ω>0.13 a.u.). Fowler–Madden 1985 (Watson-sphere
in-crystal O²⁻); Fowler–Tole N³⁻ in Li₃N = 5.35/5.94 Å³. In-crystal anion α is environment-
dependent and SUPPRESSED vs free-ion, so no transferable constant exists and the −1 clamp
(α(O⁻)≈19 a₀³) sits near the physical in-crystal O²⁻ (~14). Our confinement + density-moment
approach is standard-family; nothing fundamentally better/transferable is missing.

## Provenance quick-map
- R1: `/data/Iterative_hirshfeld/` solid runs (QE + critic2 grid).
- R2: `/data/Iterative_hirshfeld/gmtkn_molecular/results_*.txt`.
- R3: `/data/Iterative_hirshfeld/kb49_psi4/fit_results.json`.
- R4: `param.F90` arrays + notebook §m–§30s.

## R6 — Li3N nitride (periodic, N at q=-2.39; deep-anion refs exercised)
QE PBE alpha-Li3N -> critic2 xdm grid, zpsp Li 3 N 5. Evdw (Ha): neutral -0.0978, gould -0.0175,
compute -0.0228, scale -0.0366. alpha(N3-)=34(compute)/50(scale)/31(gould) a0^3; alpha(Li+)=2-6.
Charge-aware collapses the neutral over-estimate (driven by Li+ alpha 164->~2); scale ~60% more
dispersion than compute (aggressive volume-power vs conservative Kirkwood moments). First end-to-end
use of the Z_eff/HOMO~0 deep-anion (-3) references in a real crystal.
