# Charge-aware XDM — complete implementation & paper handoff

**READ THIS FIRST to resume.** Single source for continuing the code and writing the paper.
Covers: what's done, where everything lives (all machines), how to access the clusters, the data
needed, how to build/run/validate each code, the a1/a2 refit, all validation numbers, in-flight
jobs, and the paper narrative. Companion docs: `PAPER_METHODS.md` (repo, canonical methods),
`PAPER_TRAIL.md`, `REFERENCE.md`, `METHODS_EXPLAINED.md`, `ION_REFERENCE_GENERATION.md`,
`DATA_INDEX.md` (tank data map), memory `xdm_master_handoff.md`. Last updated 2026-08-24.

---

## 0. THESIS / WHAT THIS IS
Charge-aware (Hirshfeld-I) reference polarizabilities for XDM dispersion. Standard XDM scales the
FREE-NEUTRAL-atom α by the AIM volume; for ions this over-polarizes cations and over-binds. We
replace α with a CHARGE-MATCHED reference tied to the iterative-Hirshfeld (HI) charge Q_i, via
selectable routes. Same idea as Bučko's TS/HI, and it dovetails with Alberto's parallel **NEXDM**
(non-empirical Kirkwood-in-XDM). Target paper with A. Otero-de-la-Roza (aoterodelaroza) + E. Johnson.

Routes (ialpha code): 0=**hineutral** (HI partitioning + neutral α — the baseline that isolates the
partitioning change) · 1=gould · 2=scale · 3=**compute** (Kirkwood ⟨r²⟩²/N) · 4=**stern**
(Sternheimer ratios) · 5=**sternws** (Watson-sphere self-consistent Sternheimer; the deep-anion route).
The three we actively push into postg/QE: **compute, stern, sternws** (+ hineutral baseline).

---

## 1. CLUSTER / MACHINE ACCESS
- **corvette** = the local box these sessions run on (cwd `/root/critic2_development`, the editing
  MIRROR — NOT a git repo). Tank mounted at `/tank`. Reaches all others; `ssh fir` works from here.
- **fir** (Alliance cluster): `ssh fir`. Acct **def-anatole_cpu** (NOT plain def-anatole). SLURM.
  Everything under `/scratch/albd/`. Modules: `StdEnv/2023 gcc/12.3`, `quantumespresso/7.5`,
  `gaussian/g16.c01` (g16!), `imkl/2024.2.0 libxc/6.2.2 flexiblas`, `scipy-stack`. NO 2FA snag.
- **dev-srv**: `ssh albd@10.10.49.104` (sudo pw **800127828**). Git repo `~/critic2_dev`
  (branch research/xdm-hirshfeld-i, fork **albdprice/critic2_dev** — COMMIT/PUSH HERE ONLY,
  user-attributed only). Build critic2 on TANK `/data/Iterative_hirshfeld/critic2_build`
  (root disk full). ld1.x = `/usr/bin/ld1.x` (apt QE6.7). **Gaussian: g09 at `/home/albd/g09/g09`**
  (g09root=/home/albd; `. /home/albd/g09/bsd/g09.profile`); the `/home/albd/g16` dir is INCOMPLETE.
- **tank** (corvette-local `/tank/research/xdm_chargeaware/`): DURABLE archive (fir/nibi scratch is
  Lustre-flaky + purge-prone → tank is truth). Also `/tank/research/refdata/` = the canonical fit set.
- **nibi**: old Psi4/mesh data, 2FA-gated, not needed.
- Heredoc/`(`/`$` over ssh mangles → write scripts locally + scp/rsync.

---

## 2. GIT / REPOS
- **critic2 fork** albdprice/critic2_dev, branch research/xdm-hirshfeld-i (dev-srv ~/critic2_dev;
  mirror /root/critic2_development). Charge-aware XDM in `src/xdm@proc.f90` + `src/param.F90`.
  Key commits: 42260ea4 (vacuum-α fix), 1a8dd386 (stern+2), f00e3ec9 (sternws), **c89eee14**
  (sternws raw-aws fix), cdc9e9fa (scrub), 89c6eb92 (PAPER_METHODS), f3342c24 (gen README).
- **postg fork** = `/root/critic2_development/postg_alastair/` (git repo, Erin/Alberto/Kyle
  toolchain + our `chargeaware.f90`). My working copy: `/root/critic2_development/postg/`
  (identical chargeaware.f90). postg is NOT the critic2 fork — separate codebase.
- **QE 7.6** = `/root/critic2_development/q-e/` (staged to fir? no — built on dev-srv `~/qe_ca`).
- **RULES:** push only to the fork, NEVER upstream. Commits user-attributed only, no
  co-author/attribution lines. No AI-assistant or scratchpad-path references in tracked files;
  grep tracked files for such references before any commit. No hosted-webpage publishing for this project.

---

## 3. THE CHARGE-AWARE MODULE (portable core — validated bit-exact vs critic2)
`chargeaware.f90` (postg) / `ca_xdm.f90` (QE, same code, module renamed) + `ca_tables.inc`
(the α tables extracted verbatim from critic2 param.F90 via a script; alpha_free is in Å³ ⇒ /bohr³).
Self-contained (numpy-free Fortran): `ca_load_refs(dir)`, `ca_refrho(Z,Q,r)` (log-grid 4-pt Lagrange
+ linear-in-Q blend of bracketing integer refs), `ca_hirshfeld_i(...)` (Picard HI SCF, β=0.2,
tol=1e-4), `ca_alpha(Z,Q,vaim,ialpha)` (the routes), `free_moments` (V_free(Q) + Kirkwood moments),
`ca_run` (HI SCF + repartition volumes/moments), `ca_setup(ialpha,dir)`. **Reference `.rho` format:**
`<sym>_q<±N>.rho`, columns `r[bohr]  rho[e/bohr³]`, log grid; q signed (`o_q-2.rho`, `na_q+1.rho`).

Generators (`tools/wfc_generator/` on dev-srv ~/critic2_dev): `gen_anion_rho_ld1.py <Sym> <q>`
(confined ld1.x density; q=0/+/−), `gen_ion_alpha_sternheimer.py`+`batch_sternheimer.py`
(rstern_p1/p2/m1), `gen_ion_alpha_watson_scf.py`+`batch_watson.py` (aws_m1/m2/m3, the Watson-sphere
self-consistent LDA atom). See ION_REFERENCE_GENERATION.md.

---

## 4. postg — DONE + VALIDATED
- Files (in `/root/critic2_development/postg/` and postg_alastair/): `chargeaware.f90`,
  `ca_tables.inc`; wired into `postg.f90` (keyword parse + `ca_run` after evalwfn + α printout) and
  `wfnmod.f90` (edisp uses `ca_alpha`). Makefile: added `chargeaware.o` + `-ffree-line-length-none`.
- **Run:** `postg <a1> <a2> <wfx> <func> ca {hineutral|compute|stern|sternws} <refsdir>`
  (func e.g. pbe0/b3lyp). Neutral baseline: omit the `ca …`.
- **Forces/stress: automatic** — edisp computes f/q from the (now charge-aware) c6/c8/c10; frozen-
  coefficient, identical to neutral XDM (no dC_n/dR). So geometry opts work directly.
- **Validation vs critic2 (NaCl molecule, same ld1 refs):** HI charge Na +0.8940 vs critic2 0.89473
  (0.08%); anion α ~0.05%, cation α ~1% (ld1-vs-read_db cation ref). α ROUTES bit-exact (O²⁻ sternws
  9.23832=9.23832, S²⁻ 31.45507). Total Evdw differs ~10% CODE-BASELINE (postg vs critic2 exchange-hole
  M3 differs 40% — present in neutral too; NOT a port error). Build on fir: `module load StdEnv/2023
  gcc/12.3; make` in postg_src → `/scratch/albd/kb49_ca/postg_src/postg`.

---

## 5. QE 7.6 — code complete, building (dev-srv ~/qe_ca)
- Files: `PW/src/ca_xdm.f90` + `PW/src/ca_tables.inc` (staged), edits in `PW/src/xdm_dispersion.f90`:
  (a) HI SCF on the FFT grid in `energy_xdm` (grid × periodic-image `lvec` loop, `mp_sum`, Q=Z−N on
  rhoae, Picard); (b) volume/moment repartition with HI weights (`promol_hi`, `ca_refrho`);
  (c) α line → `ca_alpha`; (d) **hineutral** = ca_ialpha=0 → neutral α on HI volume, NO cap
  (matches postg). Makefile: `ca_xdm.o` added to PWOBJS; make.depend has `xdm_dispersion.o : ca_xdm.o`.
  nspin=1 only for now (campaign is closed-shell PAW).
- **Activation via ENV VARS:** `XDM_CA_ROUTE={hineutral|compute|stern|sternws}` + `XDM_CA_REFS=<dir>`
  (no QE namelist change). Standard `vdw_corr='xdm'` + a1/a2 still apply.
- **Forces/stress: automatic** — `force_xdm`/`stress_xdm` return fsave/ssave computed in energy_xdm
  from the charge-aware coefficients (frozen-coefficient). Opts work once built.
- **BUILD (dev-srv ~/qe_ca):** `./configure` (openblas+MPI+internal FFT; make.inc MUST persist —
  earlier backgrounded configures got clobbered, run it FOREGROUND) → `make depend` (generates the
  43 make.depend, incl. the ca_xdm dep) → external libs are `-j`-RACE-PRONE (devxlib/MBD `mkdir`
  collide): build them SERIALLY `for t in libfox libmbd libdevx libbeef; do make -j1 $t; done`
  (if devxlib fails: `rm -f devxlib/make_devx devxlib/src/*.{o,mod,a}; make libdevx`), THEN
  `make -j$(nproc) pw`. Binary → `bin/pw.x`. **STATUS 2026-08-24: devxlib fixed, PW compiling;
  `ca_xdm.o`+`xdm_dispersion.o` compile is the code-validity checkpoint (not yet confirmed).**
- **VALIDATE (when built):** run pw.x scf on MgO (PAW, B86b pseudos) with `vdw_corr='xdm'` +
  `XDM_CA_ROUTE=sternws XDM_CA_REFS=<refs>`; compare HI charges (Mg +2/O −2) + O²⁻ α + Evdw to
  critic2 `xdm grid` on MgO (numbers in headnode_cmp/MgO/ + campaign). MgO pw.x setups on fir.

---

## 6. a1/a2 REFIT — postg + Gaussian, on fir (in flight)
**Goal:** fit BJ damping a1/a2 for each functional/basis × route, so we can then run geometry OPTS.
KB49 = single-point interaction energies (NO opts, NO forces for the fit). Downstream opts use forces.
- **Dir:** `/scratch/albd/kb49_ca/` on fir. `geom/` (147 .xyz = 49 dimers + 98 monomers, from
  tank kb49_mesh_psi4/b3lyp/geom), `refs/` (191 .rho, combined ld1 neutral+cation+anion; backed to
  tank `reference_densities/combined_ld1_qall/`), `kb49.din` (references), `postg_src/postg` (built).
- **Jobs:** g16 array **56582586** (2 func PBE0[pbe1pbe]/B3LYP × 147 geoms, **def2-TZVPD**,
  output=wfx) → postg array **56582587** (afterany dep; 4 routes hineutral/compute/stern/sternws →
  `.cro`). Scripts: `run_g16.sh`, `run_postg_kb49.sh`.
- **FIT:** `python3 /scratch/albd/kb49_ca/fit_kb49_ca.py` (module load scipy-stack). Matches the
  CANONICAL Octave convention in `/tank/research/refdata/50_fit/` (dev-srv ~/projects/refdata):
  reads postg's "coefficients and distances" block (like `reader_g09.m`), `energy_bj` (rvdw=a1·rc+a2)
  and **`energy_bj0`** (rvdw=a2, ZERO-damping, 1-param), metrics **MAD / RMS / MAPD** (=MAPE) as in
  `fit_report_full.m`. VALIDATED: energy_bj reproduces postg's own Evdw exactly. **If a1 fits ≤0 for
  a route → auto-refit that route with bj0** (a1 collapse = BJ wants negative = use zero-damping).
  Output table: `func route damp a1 a2(Ang) RMS MAD MAPD%`.
- **When done:** `sbatch`-less `python3 fit_kb49_ca.py` → the a1/a2 table. Then repeat the same
  workflow for QE (env-var route) + the def2-QZVP basis if wanted; and X23/solid OPT tests.

---

## 7. IONIC-SOLID CAMPAIGN (the headline result — DONE)
- fir jobs **53769928** (5 opts) + **53790775** (sternws) DRAINED. Results pulled to tank
  `campaign_runs/` (198M; QE wfc scratch left on fir). Redo of 49 failed cells + 3 stale
  divalent-stern = job **55975629**. Collated ionic table `campaign_runs/RESULTS_ionic_cohesive.txt`.
- **Result:** neutral over-binds (MAE ~0.37 eV/atom); sternws best (MAE ~0.09), esp. divalent
  (MgO 5.51→5.20, MgS 4.21→3.81, CaO 6.19→5.72 vs expt 5.12/3.71/5.50). α vs Tessman in-crystal:
  S²⁻ 4.66 Å³ (Tessman 4.8–5.9), O²⁻ 1.37. Full table + provenance in REFERENCE.md / PAPER_TRAIL.md.

---

## 8. DATA — what the paper needs, all on tank (`/tank/research/xdm_chargeaware/`)
- `data/reference_densities/` — 218 anion .rho + 118 neutral wfc + generators + SUMMARY.md;
  `combined_ld1_qall/` (191, neutral+cation+anion, the postg/QE refs). `_meta/pseudopotentials/`.
- `data/campaign_runs/` — ionic + X23 cohesive results (result.json + scf.out + xdm.out + cubes 239G).
- `data/headnode_cmp/{MgO,MgS,CaO}/` — the α-validation reference outputs.
- `data/alpha_generation/{watson,stern}/` — raw α-table generation data.
- `papers/solid_refs/` — all cited PDFs incl. PhysRevB.87.064110.pdf, Bučko JCTC+JCP, Tessman,
  Csonka, Zhang SI, CRC. `/tank/research/refdata/` — the canonical KB49/GMTKN55 .din + Octave fit.
- postg + QE source: NOT on tank (in git/local); back up if desired.

---

## 9. IN-FLIGHT / NEXT STEPS (priority order)
1. **QE build finish** (dev-srv ~/qe_ca) → confirm ca_xdm.o+xdm_dispersion.o compile → MgO validate.
2. **KB49 refit** (fir 56582586→56582587) → `fit_kb49_ca.py` → a1/a2 table (PBE0/B3LYP × 4 routes,
   def2-TZVPD; watch for bj0 switches).
3. Repeat refit for **QE** (same din+geoms, QE energies via env-var route) + optionally def2-QZVP.
4. With a1/a2: **geometry-opt tests** (X23 crystals, ionic solids) in both codes — forces wired.
5. Molecular ionic benchmark (GMTKN55 IL16/AHB21/CHB6/IONPI19). Cations +3/+4 references.
6. Paper: methods (PAPER_METHODS.md), results (ionic cohesive + α-vs-Tessman + KB49 damping +
   the postg/QE cross-validation vs critic2), the NEXDM alignment (compute=Kirkwood, sternws best α).

---

## 10. PAPER-READY VALIDATION NUMBERS (quick reference)
- α routes bit-exact postg↔critic2 (anions), ~1% cations (ref-source).
- HI charges postg↔critic2: 0.08%.
- Ionic cohesive MAE: neutral 0.37 → sternws 0.09 eV/atom.
- α vs Tessman in-crystal: S²⁻ sternws 4.66 (compute 30 too soft, stern-clamp 6.7 too stiff), O²⁻ 1.37.
- Watson SCF cross-check: reproduces cation stern tables (Mg²⁺ 0.00475 vs 0.00461); bare-H α=4.5003.
- Bučko C6 signature reproduced (cation collapse ~99%, anion grow 3–4×).
