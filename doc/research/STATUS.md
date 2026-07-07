# STATUS — charge-aware (Hirshfeld-I) XDM  ·  living board

Quick-scan status. Update at the end of each session. Detail lives in the lab
notebook (`xdm_hirshfeld_i_notebook.md`, chronological), numbers in `RESULTS.md`,
method prose in `PAPER_METHODS.md`, related work + comparison targets in
`LITERATURE_INTEGRATION.md`.

Last updated: **2026-06-03**

---

## ✅ DONE (implemented + validated)
- **Code — four charge-aware α routes** in molecular (`xdm_wfn`) and periodic
  (`xdm_grid`) paths, keyword-gated: `hirshfeld_i [volonly] [volref]
  [alpharef gould|scale|compute|stern] [wfcdir]`. Cap dropped; HI weights on
  volumes + moments. Shared `chargeaware_atpol` helper (molecule = solid formula).
- **Hirshfeld-I SCF** stabilized in both mesh and grid paths (anion clamp, refrho≥0
  floor, mixing). Anion `.rho` reference densities: whole periodic table (Route 2
  confined ld1.x, `dat/hirshfeld_proatoms/ld1_pbe`).
- **Reference data** in `param.F90`: Gould–Bučko ion α (2A gould), volume exponents
  p′ (2B scale), Sternheimer ratios (2C stern). gould anion tier verified = benchmark
  (halides tier-invariant → validated results correct).
- **Ionic-solid cohesive-energy validation (6 solids)** — the headline win.
- **GMTKN55 ionic molecular benchmark (IL16/CHB6/AHB21)** — clean null.
- **KB49 a1/a2 refit (Psi4 PBE/def2-TZVP)** — done, 147 species; scale best.
- **Literature** — all 16 method papers read + integrated.
- All data persisted `/data/Iterative_hirshfeld/`; committed to fork
  `research/xdm-hirshfeld-i`.

## 🔄 RUNNING / EXTERNAL
- **QE plane-wave KB49 refit** — packaged for HPC (`/data/Iterative_hirshfeld/
  kb49_qe/HPC_package/`); waiting on the HPC agent (can't run on dev-srv, 31 GB).
- **GMTKN batch** — PAUSED (192 species cached; ionichb/s66 unfinished, low value).

## 🔲 TODO (priority order)
1. **[SCIENCE, decisive] Dispersion-dominated benchmark** — X23 molecular crystals
   and/or layered/TMD solids. This is where FI/TS+HI report the biggest charge-aware
   wins AND where PBE-alone fails. **Our current central claim is under-supported
   without it** (see gap analysis). Needs periodic densities (QE/FHI-aims).
2. **[SCIENCE, high] B86bPBE ionic solids** — OdlR&J-2020 fix ionic solids via B86b
   *exchange*, not references. Show charge-aware helps on top (complementarity) or a
   reviewer says "just use B86b." Needs B86bPBE-XDM a1/a2 + B86bPBE densities.
3. **[SCIENCE, med] Geometry/properties** — we did single-point cohesive energy at
   fixed geometry; add relaxed lattice constants + bulk moduli (vs OdlR&J-2020
   0.06 Å / ~5 GPa) to strengthen the solid case.
4. **[METHOD, med] a1/a2 refit robustness** — MAE-vs-RMSP cross-check (a1→0 boundary);
   validate refit params on an independent set; finish the QE plane-wave refit.
5. **[CODE, med] Molecular HI-XDM regression test** (#36) — needs shipped wfx/fchk.
5b. **[CODE, low] Anderson/DIIS mixing for the grid HI-SCF** — linear mixing (β=0.2) now
   converges soft multi-site charge-transfer modes (e.g. Li3N Li1b↔Li2c) but needs ~215
   iters; Anderson/Pulay on the per-atom charge vector would do it in ~20–30 and remove the
   β sensitivity. Do if the periodic multi-site ionic set grows. See notebook 2026-06-03i.
6. **[RESOLVED 2026-06] Multiply-charged anion references.** O²⁻/S²⁻ done (embedded α + density routes). **N³⁻/P³⁻/As³⁻ and C⁴⁻/Si⁴⁻ now covered** via the Z_eff generator (Heidar-Zadeh, HOMO≈0 zero-EA criterion) → density routes. Deep refs are compact frozen-orbital constructs (env-dependent α; small effect on final energies via vfq cancellation). Investigated:
   the linear-response (Sternheimer) and free-benchmark (Gould) routes cannot give a
   double-anion α (unbound 2nd electron → diverges/diffuse). The **density-based routes
   (`compute`/`scale`) handle any charge via the bound confined density** and are the
   correct choice; gould/stern clamp at −1 (documented, justified proxy ≈ physical
   in-crystal O²⁻). See RESULTS.md R5. Remaining option (low priority): strictly
   FI-faithful gould at −2 via the frozen-orbital embedded tier (external data).
7. **[low] Multiply-charged ions + dynamic α(iω)** (#44); resume/drop GMTKN ionichb/s66.

## ⛔ BLOCKED / DEPENDENCIES
- QE plane-wave refit → HPC access (large-memory node).
- Dispersion-dominated solids → periodic densities pipeline (have QE small-cell; big
  cells / molecular crystals need HPC or careful cell setup).

## ❓ OPEN QUESTIONS (scientific)
- Does charge-aware XDM beat neutral XDM on a set where dispersion *dominates* and
  PBE-alone fails? (the paper's crux — untested)
- Is the ionic-solid win from charge-aware references *additive* to the B86b-exchange
  fix, or do they overlap?
- Do the refit a1/a2 (esp. scale) generalize off KB49, or is scale's edge KB49-specific?
