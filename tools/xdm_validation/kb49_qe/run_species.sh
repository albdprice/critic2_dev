#!/bin/bash
# Run ONE KB49 species end-to-end on an HPC compute node:
#   pw.x SCF (45-bohr MT box, 80/800 PAW)  ->  pp.x rho + ELF cubes
#   ->  critic2 grid-XDM x4 routes (neutral/gould/scale/stern) -> AP-format JSON.
# Writes only small results (JSON + base energy) to $PKG/results; cubes stay in scratch.
#
# Required env (set in the SLURM script):
#   PKG       package root (this dir's parent)
#   PW, PP    full paths to QE pw.x / pp.x
#   CRITIC2   full path to the research/xdm-hirshfeld-i critic2 binary
#   WFC       path to anion reference densities (dat/hirshfeld_proatoms/ld1_pbe)
#   SCRATCH   node-local scratch with space for ~2x150MB cubes
#   NP        MPI ranks for pw.x ; OMP threads/rank via OMP_NUM_THREADS
#   BOX       cubic box side in bohr (default 45)
set -e
sp=$1
: "${PKG:?}" "${PW:?}" "${PP:?}" "${CRITIC2:?}" "${WFC:?}" "${SCRATCH:?}"
BOX=${BOX:-45}; NP=${NP:-$SLURM_CPUS_PER_TASK}; export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
GEOM=$PKG/geom; W=$SCRATCH/$sp; mkdir -p "$W/out"; cd "$W"

# 1. pw.x SCF
python3 "$PKG/scripts/qe_make_input.py" "$GEOM/$sp.xyz" "$BOX" "$sp" "$W/out" scf.in
mpirun -np "$NP" "$PW" -in scf.in > scf.out 2>&1
grep -q "JOB DONE" scf.out || { echo "$sp: pw.x did not finish"; exit 1; }
E_RY=$(grep '!.*total energy' scf.out | tail -1 | awk '{print $5}')
E_HA=$(python3 -c "print($E_RY*0.5)")

# 2. pp.x rho (plot_num=0) + ELF (plot_num=8)  -- run SERIAL (MPI cube write unreliable)
for tag in val:0 elf:8; do nm=${tag%%:*}; pn=${tag##*:}
  printf '&inputpp\n prefix="%s", outdir="%s/out", plot_num=%s, filplot="tmp_%s"\n/\n&plot\n iflag=3, output_format=6, fileout="%s.cube"\n/\n' \
    "$sp" "$W" "$pn" "$nm" "$nm" > pp_$nm.in
  "$PP" -in pp_$nm.in > pp_$nm.out 2>&1
done

# 3. critic2 grid-XDM x4 routes -> JSON
ZP=$(python3 "$PKG/scripts/zpsp_for.py" "$GEOM/$sp.xyz")
for route in neutral gould scale stern; do
  if [ "$route" = neutral ]; then kw=""; else kw="hirshfeld_i alpharef $route wfcdir $WFC"; fi
  printf 'crystal %s/val.cube\nload %s/val.cube id rho1 zpsp %s\nload %s/elf.cube id elf1\nreference rho1\nxdm grid rho rho1 elf elf1 xa1 0.4041 xa2 2.6998 %s\n' \
    "$W" "$W" "$ZP" "$W" "$kw" > xdm_$route.cri
  "$CRITIC2" xdm_$route.cri > xdm_$route.out 2>&1
  mkdir -p "$PKG/results/$route"
  python3 "$PKG/scripts/kb49_makejson_qe.py" xdm_$route.out "$GEOM/$sp.xyz" "$E_HA" "$PKG/results/$route/$sp.json"
done

# 4. STAGE the COMPLETE per-species run record (every input + output, exactly as run)
#    to PERSISTENT storage ($RUNS, default $PKG/runs; NOT node scratch). This is the
#    full provenance for local workup + cube reuse. Cubes gzipped (~150MB -> ~15MB).
RUNS=${RUNS:-$PKG/runs}; R="$RUNS/$sp"; mkdir -p "$R"
# inputs + text outputs (small): QE scf, pp.x rho/elf, critic2 cri+out for all 4 routes
cp scf.in scf.out pp_val.in pp_val.out pp_elf.in pp_elf.out "$R"/ 2>/dev/null || true
cp xdm_*.cri xdm_*.out "$R"/ 2>/dev/null || true
echo "$E_HA" > "$R/base_energy_ha.txt"
# the density (rho) + ELF cubes, gzipped
gzip -c val.cube > "$R/val.cube.gz"
gzip -c elf.cube > "$R/elf.cube.gz"
# the JSONs produced for this species (all routes), alongside the central results/
for route in neutral gould scale stern; do
  [ -f "$PKG/results/$route/$sp.json" ] && cp "$PKG/results/$route/$sp.json" "$R/$sp.$route.json"
done
echo "$sp DONE  E_base=$E_HA Ha  record->$R"
