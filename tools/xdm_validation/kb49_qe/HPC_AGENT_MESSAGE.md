# Task: run a KB49 charge-aware XDM refit with Quantum ESPRESSO (plane-wave)

Hi — please run the job below on the cluster and bring back a small results set.
It's a parameter-fitting run: 147 small/medium molecules, PBE plane-wave SCF, then a
lightweight post-process. Everything is pre-built and self-contained in one tarball.

## 1. Get the package
On dev-srv (`albd@10.10.49.104`):
`/data/Iterative_hirshfeld/kb49_qe/KB49_QE_HPC_package.tgz`  (3.1 MB)
Copy it to the cluster and unpack into a job dir:
```
scp albd@10.10.49.104:/data/Iterative_hirshfeld/kb49_qe/KB49_QE_HPC_package.tgz .
mkdir kb49_qe && tar xzf KB49_QE_HPC_package.tgz -C kb49_qe && cd kb49_qe
```
Read `NOTES.md` first — it has the full protocol. This message is the short version.

## 2. Prerequisites on the cluster
- **Quantum ESPRESSO 7.x** (`pw.x`, `pp.x`) — module or build.
- **Python 3** with **numpy + scipy** (for the final fit only; not on compute nodes necessarily).
- **critic2 — our fork, specific branch** (REQUIRED; the dispersion routes are our
  additions, not in upstream):
  ```
  git clone https://github.com/albdprice/critic2_dev.git
  cd critic2_dev && git checkout research/xdm-hirshfeld-i
  mkdir build && cd build && cmake .. && make -j critic2
  ```
  Needs gfortran, LAPACK/BLAS, libxc, OpenMP. Binary → `build/src/critic2`.
  (The anion reference data it reads is shipped in the package as `ld1_pbe/`.)

## 3. Configure + submit
Edit the top of `submit_kb49_qe.slurm` for your cluster:
- `--account`, and the `module load` line (QE + compilers/MPI),
- `PW`, `PP` (QE binaries), `CRITIC2` (the fork build from step 2).
Then:
```
sbatch submit_kb49_qe.slurm        # array 1-147, one molecule per task
```
Each task writes `results/<route>/<species>.json` (tiny). 4 routes: neutral, gould,
scale, stern. `species.list` (147) and `kb49.din` are already in place.

## 4. Resources (why this needs a big-memory node)
The protocol is a **fixed 45-bohr cubic box at ecutrho=800 Ry** (Martyna-Tuckerman
isolated). That's a ~405³ FFT grid ⇒ QE needs **~16 GB per MPI rank** (more for the
biggest molecules). The SLURM script requests **256 GB / 8 MPI ranks / 8 h** per task
— please adjust to your partition. Small molecules finish in minutes; the naphthalene
dimers (`c10h8_c10h8_*`, 36 atoms) and `adenine_thymine_*` (30 atoms) can take a few
hours. (For reference: this OOM-killed a 31 GB workstation, hence the cluster.)

## 5. Final step (login node, after the array finishes)
```
python3 scripts/fit_kb49_qe.py results kb49.din
```
Prints a 4-row table (route → a1, a2, MAE, MAPD) and writes `fit_results_qe.json`.
Expected: 49/49 reactions per route, MAPD roughly ~15-20%, a1/a2 in O(0.1-1)/O(1-5 Å).

## 6. What to send back — THE WHOLE JOB DIRECTORY (full provenance)
Please return the **entire job directory** so we keep the complete run exactly as
executed and do the workup locally. After the array + fit it contains:
- inputs/config: `scripts/`, `pp/`, `geom/`, `ld1_pbe/`, `kb49.din`, `species.list`,
  `submit_kb49_qe.slurm`;
- `results/<route>/*.json` + `fit_results_qe.json`;
- `runs/<species>/` — full per-species record: QE `scf.in/out`, `pp_*.in/out`,
  critic2 `xdm_<route>.cri/out` (4 routes), `base_energy_ha.txt`, and the **gzipped
  rho + ELF cubes** (`val.cube.gz`, `elf.cube.gz`) — ~6-8 GB total for all 147;
- `logs/`.
Easiest: `tar czf kb49_qe_run.tgz <jobdir>` and send that; we'll stage it at
`/data/Iterative_hirshfeld/kb49_qe/HPC_run/`. Only thing to skip: QE's node-scratch
`out/`/`.save` wavefunctions (hundreds of GB, regenerable — the cubes carry the
density we need).

## 7. Likely gotchas (already handled, just FYI)
- `pp.x` is run **serially** in `run_species.sh` (MPI cube-writing was unreliable) —
  leave it that way.
- If a task fails, it's almost always pw.x SCF non-convergence or OOM on a big
  molecule — bump `--mem` / walltime and re-run just those array indices
  (`sbatch --array=<failed ids> submit_kb49_qe.slurm`); the pipeline is per-species
  independent and re-running is safe.
- One missing species only drops its one reaction from that route's fit; not fatal.

Thanks! Ping back the `results/` + `fit_results_qe.json` and I'll fold them into the
analysis. Questions → see `NOTES.md` in the package.
```
```
