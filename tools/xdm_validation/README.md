# tools/xdm_validation — charge-aware XDM benchmark & refit pipelines

Scripts that generated the charge-aware XDM validation (see
`doc/research/xdm_hirshfeld_i_notebook.md` §31b–06-03 and
`doc/research/LITERATURE_INTEGRATION.md`). Heavy inputs/outputs (fchks, cubes,
JSONs, results) are NOT in the repo — they live on dev-srv at
`/data/Iterative_hirshfeld/` (34 TB volume). This directory is the reproducible
*code* record.

## gmtkn/ — GMTKN55 ionic molecular benchmark (result: clean null)
- `run_batch.py`  — per-species Psi4 PBE/def2-TZVP → fchk → critic2 `xdm_wfn`
  ×4 routes (neutral/gould/scale/stern), per-reaction E_int vs .din, MAE/MSE.
  Resumable, fault-tolerant, 3-wide (later 6×2). Sets: il16/chb6/ahb21/s22.
- `analyze.py`    — partial-cache directional table for a set.
- `psi4_species.py` — Psi4 PBE/def2-TZVP worker (fchk + energy sidecar).
  Gotchas baked in: `OPENBLAS_NUM_THREADS=1` (OpenBLAS×OpenMP nesting hangs),
  two-stage SCF (default → SOSCF retry), PSI_SCRATCH on a large volume.

## kb49_psi4/ — KB49 a1/a2 refit via Psi4 (DONE)
- `kb49_run.py`    — 147-species driver + per-route RMSP fit (AP's convention).
- `kb49_makejson.py` — critic2 coefficient block (mesh + grid) → AP-format JSON
  ({base_energy, coords[bohr], c6/c8/c10/rc}); validated to 3.6e-14 Ha vs critic2.
- Result: scale a1=0.186/a2=3.730 Å/MAE 0.605 (best); neutral validated vs
  Otero-de-la-Roza–Johnson 2013 PBE plane-wave (a1=0.4073).

## kb49_qe/ — KB49 refit via QE plane-wave (PACKAGED for HPC)
Fixed 45-bohr/ecutrho=800 box OOMs a 31 GB node → needs ≥256 GB HPC node.
Self-contained HPC package (this dir + pseudos/geoms staged on /data):
- `qe_make_input.py` — xyz → pw.x input (MT-isolated, 80/800, PAW kjpaw_psl 1.0.0).
- `run_species.sh`   — pw.x → pp.x rho+ELF → critic2 grid-XDM ×4 → JSON; stages
  the full per-species record (I/O + gzipped cubes) for return.
- `kb49_makejson_qe.py`, `zpsp_for.py`, `fit_kb49_qe.py` (self-contained RMSP fit).
- `submit_kb49_qe.slurm`, `NOTES.md`, `HPC_AGENT_MESSAGE.md` — the handoff.
- **Cluster must build THIS fork branch** (gould/scale/stern grid keywords are ours).

Reference bibliography of the methods: `doc/papers/references.bib`; per-paper
integration with our plan: `doc/research/LITERATURE_INTEGRATION.md`.
