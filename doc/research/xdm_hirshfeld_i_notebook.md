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
- **RQ2 (α only vs also moments):** XDM-specific, no literature guidance
  (TS/FI have no exchange-hole moments). Should HI weights partition only
  the polarizability/volume, or **also** the exchange-hole moments M_l?
  The moments already integrate the real density — quantify how much the
  *weight* vs the *reference* matters. *(Stage 0, two sub-modes)*
- **RQ3 (volume reference):** does swapping V_free,neutral → charge-matched
  V_ref(Q_A) (FI volume formula, using our ion densities) improve over
  neutral-volume + dropped cap? *(Stage 1)*
- **RQ4 (polarizability reference):** is the neutral α_free anchor good
  enough, or is the FI ion-polarizability database necessary for
  accuracy? *(Stage 2)*
- **RQ5 (refit):** do XDM's a1/a2 damping parameters need refitting after
  any of the above? (FI did not need a refit; XDM ≠ MBD — verify.)

## 5. Open literature gap (being closed — Task #27)

The 2014 TS+HI paper's method for obtaining its **integer free-ion
reference densities for unbound anions** was not pinned down by the
automated review (only FI's 2017 power-law confinement is documented).
This is the most load-bearing comparison for our box-confinement choice.
→ reading Bučko 2013/2014 methods directly; result logged below.

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

