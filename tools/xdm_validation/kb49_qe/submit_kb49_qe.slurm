#!/bin/bash
#SBATCH --job-name=kb49_qe_xdm
#SBATCH --array=1-147
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8          # MPI ranks for pw.x
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G                   # 45-bohr / 800-Ry box: ~16 GB/rank for the FFT grid
#SBATCH --time=08:00:00
#SBATCH --output=logs/%A_%a.out
# ---- EDIT THESE FOR YOUR CLUSTER ----
#SBATCH --account=def-anatole
module purge
module load StdEnv/2023 gcc/12.3 openmpi quantumespresso/7.2   # or your QE module
export PKG=$SLURM_SUBMIT_DIR
export PW=$(command -v pw.x)            # or absolute path to your QE 7.x pw.x
export PP=$(command -v pp.x)
export CRITIC2=$HOME/critic2_dev/build/src/critic2   # research/xdm-hirshfeld-i build (see NOTES)
export WFC=$PKG/ld1_pbe                 # anion reference densities (shipped in this package)
export SCRATCH=${SLURM_TMPDIR:-/tmp}/kb49_$SLURM_ARRAY_TASK_ID   # fast node-local work
# Persistent full-record staging: every input+output per species (gzipped cubes;
# ~6-8 GB total over all 147). MUST survive the job (project space, NOT SLURM_TMPDIR).
export RUNS=$PKG/runs
export NP=$SLURM_NTASKS_PER_NODE
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BOX=45
# --------------------------------------
mkdir -p logs
sp=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PKG/species.list")
echo "[$SLURM_ARRAY_TASK_ID] $sp on $(hostname)"
bash "$PKG/scripts/run_species.sh" "$sp"
