# Literature integration — charge-aware (Hirshfeld-I) XDM

How the 16 papers in `doc/papers/` map onto our project (routes neutral/gould/scale/stern,
the KB49 a1/a2 refit, and the molecular + ionic-solid validation), the **numbers we must
compare against**, and the **action items / corrections** the close reading surfaced.
Companion to `xdm_hirshfeld_i_notebook.md` (the worklog). Built from a full read of all
PDFs (Jun 2026).

PDFs are at `doc/papers/` and mirrored at `/data/Iterative_hirshfeld/papers/pdf/`.

---

## 0. Headline corrections the reading forced (READ FIRST)

1. **`s00894-017-3514-6.pdf` is mislabeled in our bib.** It is **Heidar-Zadeh, Ayers &
   Bultinck, *J. Mol. Model.* 23:348 (2017), "Fractional nuclear charge approach to
   isolated anion densities for Hirshfeld partitioning"** — NOT Bučko 2017 periodic MBD.
   It's actually very useful (unbound-anion reference densities), but the real Bučko-2017
   MBD paper is **not** in the folder. Bib updated.

2. **"KB49" vs "KB65".** Kannemann–Becke 2010 fit on **65** complexes (S22 + 10 rare-gas +
   12 NC31/05 + 21 JB-vdW), PW86PBE, BR-hole, a1=0.82/a2=1.16 Å, MAE 0.33 kcal/mol.
   **KB49 = KB65 minus the noble-gas dimers** and is the set Otero-de-la-Roza & Johnson
   2013 (and we) refit on. So the correct neutral-route comparison is **OdlR&J 2013 PBE
   *plane-wave* a1=0.4073, a2=2.4150 Å** — which our refit (a1=0.4041, a2=2.6998)
   reproduces on a1 to ~0.003 (a2 ~0.28 Å higher). **Do NOT compare our a1/a2 to KB's
   PW86PBE 0.82/1.16** — different functional + different hole. Our neutral refit is
   validated; the small a2 offset is the charge-α-shift + our exact KB49 subset.

3. **IONPI19 is NOT in GMTKN55.** The three GMTKN55 ionic molecular subsets are
   **AHB21 (21), CHB6 (6), IL16 (16)** — which is exactly what we ran. IONPI19 (ion-π) is
   a separate later set; treat it as optional/external, not part of GMTKN55. (Task #45
   title should drop IONPI19.)

4. **gould-route anion tier — VERIFIED 2026-06 (benchmark tier; results unaffected).**
   Confirmed: `alpha_gb_m1` in `param.F90` (provenance `dat/xdm_ion_alpha_gould_bucko_2016.dat`,
   header "benchmark") is the **free self-consistent** anion tier (F⁻=15, Cl⁻=30.3, Li⁻=1180,
   Na⁻=1310), NOT the FI embedded/frozen-orbital tier. **Impact: none for anything validated.**
   (a) **Halide anions are self-consistent in BOTH tiers ⇒ identical**, so the alkali-halide
   solids (F⁻/Cl⁻) are correct. (b) Oxide O is clamped at −1 with α=5.4≈neutral ⇒ tier-invariant.
   (c) `ion_alpha0` clamps |q|≤1 and cations interpolate neutral→+1, so the grossly-diffuse
   electropositive-metal "anion" entries are **structurally never reached**. The tier only
   differs for **non-halide anions with large charge (N³⁻, S²⁻ — not in our test set)**; strict
   FI-faithfulness there = swap `alpha_gb_m1` for the embedded tier (Gould–Bučko SI). Documented
   in the `param.F90` comment.

5. **scale-route is a departure from FI, by design.** FI-MBD uses **linear** volume scaling
   (p'=1) and Gould explicitly tested the p'≈p−0.615 exponent (his 2016 JCP) and reported it
   "did not noticeably improve results" for MBD. So our **scale** route (p'=p−0.615) is a
   *more aggressive, distinct* hypothesis, not a reproduction of FI. Expect scale ≠ gould.
   The "FI-faithful" emulation is the **gould route with linear (p'=1) volume scaling**.

6. **ionic-solid functional — we use PBE; the canonical XDM ionic-solid paper uses B86bPBE.**
   OdlR&J 2020 (`xdmmx.pdf`) argue the **exchange functional (B86b)**, not charge-aware
   references, is the lever for ionic-solid non-bonded repulsion; PBE-XDM over-binds. Our
   0.92→0.26 eV cohesive-MAE win is on **PBE**. Strongest paper-ready claim ⇒ also run
   **B86bPBE** charge-aware and show the effects are complementary (their B86b fixes
   geometry/B0; our charge-aware fixes cohesive energy). Action item below.

---

## 1. Per-paper digest, grouped by role

### A. XDM foundations (what our code computes)
- **Johnson & Becke 2005, JCP 123, 024101** (`024101`) — the XDM seed: C6 from the
  exchange-hole dipole moment ⟨d²ₓ⟩ and α; Hirshfeld partition of the moment; α apportioned
  by hole-moment fraction. Uses an *energy-cutoff* damping (κ), **not** BJ a1/a2. Role:
  definitional ancestor of our "neutral" α and moment partition.
- **Becke & Johnson 2007, JCP 127, 124108** (`124108`) — the canonical C6/C8/C10 +
  ⟨M_ℓ²⟩ moments + **BJ damping R_vdw=a1·R_c+a2** + the **volume-scaling α_i=(⟨r³⟩_i/⟨r³⟩_free)·α_free**
  that our "scale" route generalizes. **No cap in the original equations** — the min(ratio,1)
  is a later safeguard, so dropping it is a return to the original form (justify for anions
  where V_AIM>V_free legitimately). Flags a C10 λ(1,3) erratum worth checking in `calc_coefs`.
- **Otero-de-la-Roza & Johnson 2012, JCP 136, 174109** (`174109`) — first **periodic /
  plane-wave XDM** (the lineage of our `xdm_grid`): uniform grids, BR hole, **PAW
  all-electron ρ,τ reconstruction** (= our ZPSP core-cusp gotcha), lattice sums, forces/stress.
  Neutral, no charge-awareness — the baseline we improve. PBE fit there: a1=0, a2=3.879 Å.
- **Otero-de-la-Roza & Johnson 2013, JCP 138, 204109** (`204109`) — **the source of our KB49
  per-functional a1/a2** and the functional-matched V_free convention. **PBE plane-wave
  a1=0.4073, a2=2.4150** (our refit target). Recommends BLYP/B3LYP/LC-ωPBE-XDM.
- **Otero-de-la-Roza & Johnson 2020, JCP 153, 054121** (`xdmmx`) — **XDM for IONIC SOLIDS**
  (our headline domain): 20 alkali halides + CsX polymorphs. Uses **B86bPBE / B86bPBE-25X**.
  Lattice-const MAE 0.060 / 0.039 Å; bulk-modulus MAE ~4.9 / 4.7 GPa; dispersion decisive for
  CsCl B1↔B2 ranking. **Neutral references, no HI** — explicitly names HI as the fix for TS
  ionic volumes but doesn't implement it for XDM. This is the gap our work fills.

### B. The TS / iterative-Hirshfeld lineage (closest prior art)
- **Tkatchenko–Scheffler 2009, PRL 102, 073005** (`PhysRevLett...`) — TS: α^eff/α^free =
  V^eff/V^free (neutral Hirshfeld), C6 ∝ ratio². Neutral Chu–Dalgarno free-atom refs, no ions,
  no cap. C6 MARE 5.5% / 1225 pairs; one damping param s_R=0.94 (PBE), S22 MAE 13 meV. = TS
  analogue of our neutral baseline; α-scaling exponent is **1** (vs our scale's p−0.615).
- **Bučko et al. 2013, JCTC 9, 4293** (`improved-density-dependent...`) — **TS+HI**: swap in
  iterative-Hirshfeld weights so ionic AIM volumes become physical (LiF: Li ratio 0.009 / F
  1.817; charges ≈ ±1.0). **Reference α stays NEUTRAL** — charge enters *only* via the HI
  volume ratio. Unbound anions handled by **Watson spheres**. No real damping refit. Helps
  ionic/adsorption decisively (MgO/H₂O −54→−46.5 vs −46 DMC); **slightly hurts** covalent
  (S22 9.3→10.6%). = our "scale with p'=1" essentially; our α-replacement routes go beyond it.
- **Bučko et al. 2014, JCP 141, 034114** (`034114`) — full TS+HI validation: octet ionicity
  series, alkali halides, layered crystals (TMD overbinding 74–120%→43–58%), hydrides
  (LiBH₄ vol −48%→−2.6%). HI charges vs Born: NaCl ±1.02 vs 1.10, MgO ±2.11 vs 1.98, LiF ±1.01
  vs 1.05. s_R only 0.94→0.95. **"TS/HI−TS difference grows with ionicity"** = our exact
  ionic-decisive / molecular-null narrative. Richest solid-benchmark template to mirror.
- **Heidar-Zadeh, Ayers & Bultinck 2017, JMM 23:348** (`s00894...`) — **fractional-nuclear-charge**
  bound anion densities: raise Z to the smallest Z_eff that binds N electrons (zero-EA),
  then coordinate-rescale to restore the true cusp. Alternative to Watson spheres / our box
  confinement & `hirsh_i_qfloor`. Confirms even O⁻/N⁻ are ~unbound at HF ⇒ a clamp/confinement
  is necessary, not optional.

### C. The Gould charge-aware reference data (our gould & scale routes)
- **Gould, Lebègue, Ángyán, Bučko 2016, JCTC 12, 5920** (`a-fractionally-ionic...`) — **FI-MBD,
  the method we emulate.** α_AIM(iω;N)=α_FI(iω;N)·V_eff(N)/V_FI(N); α_FI by **linear-N interp
  between integer ions**; HI populations; **β=0.83 unchanged** from neutral. Crystal-α MARE
  **TS 157% → HI 43% → FI 23%**. Li⁰=164→embedded ~5≈Li⁺ 0.2 a₀³. NaCl/MgO/LiF that crash
  under neutral MBD succeed under FI. Uses **linear** volume scaling (rejected p'-exponent).
- **Gould & Bučko 2016, JCTC 12, 3603** (`c6-coefficients...`) — **the ion-α/C6 database** our
  gould route reads. TDDFT α(iω) for ~411 species rows 1–6; static α(0) Tables 2–4
  (Li⁺=0.193, Na⁺=0.930, K⁺=5.05, F⁻=15.5*, Cl⁻=30.3*, O=5.20, F=3.60 a₀³ — *free anion tier*);
  **two-Lorentzian** dynamic model in SI; **embedded/frozen-orbital anion tier (Appendix B)**.
  Empirical C6≈Ξ(α_Xα_Y)^0.73. **Use the embedded anion tier for HI references** (see §0.4).
- **Gould 2016, JCP 145, 084308** (`084308`) — **volume scaling, our scale route.**
  α/αᴿ=(V/Vᴿ)^p', C6/C6ᴿ=(V/Vᴿ)^p, with **element-specific p' = p − 0.615** (not the
  rule-of-thumb p'=1). Tabulated p (PGG): C 2.00, N 2.12, O 2.24, F 2.33, Cl 2.33, Na 2.07,
  Li 1.76 ⇒ p'≈ C 1.39, O 1.63, F 1.72, Cl 1.72, Li 1.15. Confined neutral atoms, **rows 1–5
  (Z≤54) only**; ion-charge transfer assumed same-Z exponent (the main caveat).

### D. Partitioning & benchmarks (our pipeline scaffolding)
- **Bultinck et al. 2007, JCP 126, 144111** (`144111`) — **Hirshfeld-I definition**: w^i=ρ^{i-1}/ρ_mol^{i-1},
  fractional refs by **linear two-integer interpolation** (= our `hirsh_i_eval`), Δ<0.0005
  convergence, unique convex minimum (start-independent — sanity check for our SCF). Anions
  "quite problematic" ⇒ motivates `hirsh_i_qfloor`.
- **Kannemann & Becke 2010, JCTC 6, 1081** (`van-der-waals...intermolecular`) — KB65/KB49 set
  + BJ damping form; BR≫XX for molecules. **No charged complexes** (the gap we fill); their
  MAE 0.33 (BR) is the neutral floor we approach on covalent subsets, not beat on ionic.
- **Goerigk et al. 2017, PCCP 19, 32184** (`c7cp04913g`) — **GMTKN55**: AHB21 (avg |ΔE| 22.5),
  CHB6 (26.8), IL16 (109.0 kcal/mol) — all CCSD(T)/CBS; **WTMAD-2 normalizer 56.84 kcal/mol**.
  IL16 is electrostatics-dominated (dispersion a small fraction) ⇒ tests that charge-aware
  doesn't *break* things; AHB21/CHB6 are the dispersion-sensitive ionic tests. Compare our
  routes to GGA-D3(BJ) (BLYP/revPBE), same ladder rung — not the double-hybrid winners.
- **Manz et al. 2019, RSC Adv 9, 19297** (`c9ra03003d`) — **MCLF**, the *reference-free* foil:
  charge-aware α/C6 from in-material ⟨r³⟩,⟨r⁴⟩,N power-laws (no per-ion table), DDEC6
  partition. C6 MARE 4.45% / 1225 pairs (vs TS/IH 8.6%). Its killer argument — free O²⁻ etc.
  are unbound — is exactly **why our stern route confines the ion** (we keep a physical
  per-ion reference instead of avoiding it). Note MCLF prefers **DDEC6 > HI > Hirshfeld**.

---

## 2. Comparison-target table (numbers to reproduce / beat / cite)

| Quantity | Literature target | Source | Our status / action |
|---|---|---|---|
| PBE plane-wave KB49 a1/a2 (neutral) | **0.4073 / 2.4150 Å** | OdlR&J 2013 | ✓ ours 0.4041/2.6998 — a1 matches; check a2 offset |
| KB-set MAE (neutral, well-behaved) | 0.33 (BR) / 0.53 (XX) kcal/mol | KB 2010 | ours 0.65 (KB49, charge-α) — approach on covalent subset |
| Crystal-α MARE TS→HI→FI | **157% → 43% → 23%** | Gould FI 2016 | our analogue = cohesive MAE 0.92→0.26 eV (different metric) |
| Alkali-halide lattice-const MAE | B86bPBE-XDM **0.060 Å**, 25X 0.039 | OdlR&J 2020 | **run B86bPBE; confirm geometry doesn't regress** |
| Bulk-modulus MAE | ~4.9 / 4.7 GPa | OdlR&J 2020 | not yet computed (optional) |
| HI charges (partition sanity) | NaCl ±1.02, MgO ±2.11, LiF ±1.01 | Bučko 2014 | spot-check our HI-SCF charges match |
| Integer-ion static α (gould route) | Li⁺ 0.193, Na⁺ 0.930, K⁺ 5.05, F⁻/Cl⁻ embedded-tier | Gould–Bučko 2016 | **verify we ingested embedded tier (§0.4)** |
| scale exponents p' (=p−0.615) | C 1.39, O 1.63, F 1.72, Cl 1.72, Li 1.15 | Gould 2016 JCP | verify Stage-2B table matches |
| FI volume scaling | **linear, p'=1** (NOT p−0.615) | Gould FI 2016 | gould=linear is the FI-faithful route |
| Damping refit need (charge-aware) | β unchanged (FI); s_R 0.94→0.95 (TS/HI) | Gould/Bučko | consistent w/ our near-neutral refit; stern is the exception |
| GMTKN55 ionic subset sizes / refs | AHB21 21/22.5, CHB6 6/26.8, IL16 16/109.0 kcal/mol | Goerigk 2017 | ✓ matches our runs; report WTMAD-2 (norm 56.84) |
| C6 MARE yardstick (charge-aware) | MCLF 4.45%, TS/IH 8.6%, D4 3.8% | Manz 2019 | optional: compute our HI-XDM C6 MARE on DOSD 1225 |

---

## 3. How the literature pins down our four routes

- **neutral** = Becke–Johnson 2007 α=(V/V_free)·α_free with NEUTRAL ref + cap. Baseline.
  Damping = OdlR&J 2013. ✓ validated by reproducing their PBE a1.
- **gould** = Gould–Bučko 2016 integer-ion α, linear-N interp (Bultinck/FI form). With
  **linear volume scaling this is the FI-faithful route** (Gould FI 2016 eq 1). **Must use the
  embedded-anion tier.** Closest to the published FI method.
- **scale** = Gould 2016 JCP element-specific p'=p−0.615. **More aggressive than FI** (which
  rejected this). A distinct hypothesis we are testing for XDM specifically.
- **stern** = our confined-ion Sternheimer/CPKS α. Motivated by the unbound-anion problem
  (Bultinck 2007, Manz 2019, Heidar-Zadeh 2017): keep a per-ion reference but make it
  physical via confinement (cf. Bučko's Watson spheres, Heidar-Zadeh's Z_eff). No exact
  literature equivalent — our methodological contribution. Needed the a1/a2 refit most.
- **partition (all routes)** = Hirshfeld-I (Bultinck 2007), charge-matched volumes
  (V_free(Q), our Stage 1). Unbound-anion references: our box confinement / qfloor, with
  Watson-sphere (Bučko) and Z_eff (Heidar-Zadeh) as documented alternatives.

---

## 4. Action items the reading generated

1. **[verify, high] gould anion tier** — confirm Stage-2A (#37) ingested the
   embedded/frozen-orbital anion α from Gould–Bučko, not the free self-consistent values
   (F⁻ 15.5 vs embedded). If wrong, gould over-polarizes anions. (§0.4)
2. **[run, high] B86bPBE ionic solids** — re-run the alkali-halide/oxide validation with
   B86bPBE (the OdlR&J-2020 recommended solid functional) to show charge-aware is
   complementary to the exchange-functional lever, and confirm geometries match their
   0.060 Å lattice MAE. (§0.6, table row 4)
3. **[verify, med] scale p' table** — cross-check Stage-2B p' against Gould's tabulated p
   (C 2.00→p'1.39, etc.), and note rows 6–7 fall back (Z>54 has no tabulated p).
4. **[frame, med] FI-faithful = gould+linear** — document that the *reproduction* of FI is
   gould with p'=1, and scale (p−0.615) is a separate, more aggressive test.
5. **[partition sanity, med] HI charges** — spot-check our HI-SCF charges vs Bučko 2014
   (NaCl ±1.02, MgO ±2.11, LiF ±1.01).
6. **[report, low] WTMAD-2** — express the GMTKN55 ionic results with the published
   WTMAD-2 normalizer (56.84 kcal/mol) and compare to GGA-D3(BJ), not double hybrids.
7. **[bib, done] corrections** — s00894 = Heidar-Zadeh (not Bučko 2017); KB49 vs KB65;
   IONPI19 ∉ GMTKN55. Fix task #45 title.
8. **[optional] coverage gaps (#44)** — multiply-charged ions (>+3/<−2), dynamic α(iω)
   two-Lorentzian params, row 6–7 exponents are all absent from the Gould data; gould/scale
   must fall back there (to stern or neutral).

---
## 5. Double-anion / multiply-charged-anion references (2026-06 literature check)

Verified (deep-research, 14 claims survived 3-0 adversarial verification; synthesis by hand).
Answers "are we missing a better method for O2-/S2- references?" — **No; our approach is standard-family.**

- **The free double-anion divergence is textbook.** Holka, Urban, Neogrady & Paldus, JCP 141,
  214303 (2014): free O2-/S2- give "unrealistically high" polarizabilities from RHF instability;
  a weak harmonic confining potential stabilizes them. O2- is lower-energy than O- only for
  confinement omega > 0.13 a.u. == exactly our Sternheimer divergence + confinement fix.
- **Confinement is THE accepted workaround, in three established flavors:** (a) external HO
  confining potential (Holka 2014, CCSD(T) — the gold-standard O2-/S2- values); (b) Watson-sphere
  / stabilizing-shell + coupled HF (Fowler & Madden 1985, the canonical in-crystal O2-);
  (c) frozen neutral-atom potential (Gould-Bucko minimal-chemistry). Our box (3.6*R99) is in this family.
- **No single transferable O2-/S2- value exists — it is environment-dependent.** In-crystal anion
  alpha is SUPPRESSED vs free-ion and rises with decreasing coordination / increasing bond length
  (Fowler-Madden; jp068257s 2007). Reference values to compare against: Holka 2014 (CCSD(T) confined
  O2-/S2-); Fowler-Tole N3- in Li3N = 36.1/40.1 a.u. ~ 5.35/5.94 A^3.
- **Implication for us:** (1) the -1 clamp is defensible because in-crystal O2- is pulled DOWN toward
  the O- value (alpha(O-)~19 a0^3 ~ suppressed in-crystal O2- ~14). (2) The density-moment route
  (compute / MCLF-style <r3>,<r4>) is a recognized reference-free family; our Kirkwood <r2>^2/N on the
  confined density belongs to it. (3) Nothing fundamentally better/transferable exists — the accurate
  values are environment-specific and use the same confinement idea. Method is sound.
- **For the strictly-best O2-/S2- reference:** cite/compare Holka 2014 (CCSD(T) confined). A future
  refinement would use the actual IN-CRYSTAL AIM density moments (MCLF-style) rather than a confined
  free-ion reference, which would capture the environment-dependence automatically.

Provenance: /tmp/claude-.../tasks/wjbtw6evi.output (14 verified claims + sources).
