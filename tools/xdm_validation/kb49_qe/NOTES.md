# KB49 charge-aware XDM refit — QE plane-wave package (HPC handoff)

## Goal
Refit the XDM Becke-Johnson damping parameters **(a1, a2)** for each of the four
reference routes — **neutral, gould, scale, stern** — on the **KB49** benchmark
(49 intermolecular complexes, 147 species), using **QE PBE plane-wave** densities
instead of Gaussian/Psi4. Companion to the already-done Psi4/def2-TZVP refit
(see `../kb49_psi4/fit_results.json`).

## Protocol (as specified)
- **DFT:** PBE, plane-wave (QE pw.x v7.x).
- **Cutoffs:** `ecutwfc=80 Ry`, `ecutrho=800 Ry`.
- **Box:** single fixed **cubic cell, 45 bohr (~24 A)** for every species,
  molecule centered, `assume_isolated='mt'` (Martyna-Tuckerman) + gamma point —
  no mirror-image interaction. (Edit `BOX` in the SLURM script to change.)
- **Pseudopotentials:** PAW `kjpaw_psl 1.0.0` PBE, shipped in `pp/`
  (H,C,N,O,F,S,Cl,Si; valences 1,4,5,6,7,6,7,4).
- **Dispersion coefficients:** critic2 `xdm grid` on the QE valence density
  (`pp.x plot_num=0`) + ELF (`plot_num=8`), with **ZPSP core reconstruction**
  (valence charges above). Routes select the reference polarizability:
  neutral (default), or charge-aware `hirshfeld_i alpharef {gould|scale|stern}`.
- **Fit:** minimize **RMSP** (relative % residuals) via `least_squares`, exactly
  AP's `02_collate_and_fit` convention (replicated self-contained in
  `scripts/fit_kb49_qe.py`). E_int = sum_k coeff_k (E_base_k + E_disp_k(a1,a2)).

## CRITICAL dependency — build OUR critic2 fork
The `gould/scale/stern` routes and the `alpharef`/`hirshfeld_i` grid keywords are
**additions in the fork `albdprice/critic2_dev`, branch `research/xdm-hirshfeld-i`**
— they are NOT in upstream critic2. Build it on the cluster:
```
git clone https://github.com/albdprice/critic2_dev.git
cd critic2_dev && git checkout research/xdm-hirshfeld-i
mkdir build && cd build && cmake .. && make -j critic2     # needs gfortran, LAPACK/BLAS, libxc, OpenMP
```
Point `CRITIC2` (in the SLURM script) at `build/src/critic2`.
The anion reference densities it needs for the charge-aware routes are shipped
here as `ld1_pbe/` — `WFC` already points to it. (They also live in the repo at
`dat/hirshfeld_proatoms/ld1_pbe`.)

## Resources / why HPC
The 45-bohr / 800-Ry box gives a ~405^3 dense FFT grid ⇒ QE reports **~16 GB per
MPI rank** for the grid alone (more for the big aromatic dimers). This OOM-killed a
31 GB workstation. Use a **large-memory node (>=256 GB)**. Per-species walltime
ranges from ~minutes (HF, CH4) to a few hours (naphthalene dimers `c10h8_c10h8_*`
36 atoms, `adenine_thymine_*` 30 atoms). The SLURM array (1 task/species) is sized
8 MPI ranks / 256 GB / 8 h — tune to your cluster.

## Run
1. Build critic2 fork (above); confirm QE pw.x/pp.x v7.x available.
2. Edit `submit_kb49_qe.slurm`: `--account`, `module load` lines, `PW`/`PP`/`CRITIC2`
   paths. (`species.list` already lists all 147; `kb49.din` is the reference set.)
3. `sbatch submit_kb49_qe.slurm`  → writes `results/<route>/<species>.json` (tiny).
4. When all 147 finish: `python3 scripts/fit_kb49_qe.py results kb49.din`
   → prints per-route a1/a2/MAE/MAPD and writes `fit_results_qe.json`.

## What to bring back  —  THE ENTIRE PACKAGE DIRECTORY
Bring back the **whole job directory** (it self-contains inputs, scripts, and all
generated outputs), so the complete run is preserved exactly as executed and all
workup can be redone locally. After the array finishes it contains:
- `scripts/`, `pp/`, `geom/`, `ld1_pbe/`, `kb49.din`, `species.list`,
  `submit_kb49_qe.slurm`  — the exact inputs/config used.
- `results/<route>/<species>.json` + `fit_results_qe.json`  — distilled fit data.
- `runs/<species>/`  — the **complete per-species record**: `scf.in`, `scf.out`,
  `pp_val.in/out`, `pp_elf.in/out`, `xdm_<route>.cri/out` (all 4 routes),
  `base_energy_ha.txt`, and the **gzipped rho + ELF cubes** (`val.cube.gz`,
  `elf.cube.gz`). Cubes ~150 MB raw → ~10-20 MB gz; ~6-8 GB total. Kept so future
  critic2 experiments (new routes / a1,a2 / re-partitioning) reuse the densities
  without re-running QE (`gunzip` → feed to critic2 `xdm grid`).
- `logs/` — SLURM stdout.
Stage the returned tree on dev-srv at `/data/Iterative_hirshfeld/kb49_qe/HPC_run/`.
(Node-scratch `out/` — QE `.save` wavefunctions, hundreds of GB — is NOT staged;
it's a regenerable intermediate and the cubes carry the density we need.)

## Validation provenance (done on dev-srv, small box)
- The QE→pp.x→critic2-grid→JSON chain was validated end-to-end on ch4 in a small
  (feasible) box: critic2's grid `Evdw` = Evdw6+8+10 exactly, and the fit formula
  (`scripts/fit_kb49_qe.py` `edisp`) reproduces it to the periodic-image residual
  (~2e-6 Ha at 24 bohr → negligible at 45 bohr).
- The equivalent Gaussian-path parser was validated to 3.6e-14 Ha vs critic2.
- ZPSP core reconstruction verified (free-atom volumes sane, e.g. Vfree(C)=39.0 a0^3).

## Files
- `scripts/qe_make_input.py`   xyz → pw.x input (45-bohr MT box, 80/800, PAW)
- `scripts/run_species.sh`     one species: pw.x → pp.x(rho,elf) → critic2 x4 → JSON
- `scripts/kb49_makejson_qe.py` critic2 grid stdout → AP-format JSON (coords in bohr from xyz)
- `scripts/zpsp_for.py`        ZPSP token list per species
- `scripts/fit_kb49_qe.py`     RMSP fit per route → fit_results_qe.json
- `submit_kb49_qe.slurm`       SLURM array driver
- `species.list` (147)  `kb49.din`  `pp/` (8 PAW)  `geom/` (147 xyz)  `ld1_pbe/` (anion refs)
