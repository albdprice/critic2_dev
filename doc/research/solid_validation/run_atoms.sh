#!/bin/bash
# Free-atom PBE total energies (spin-polarized, isolated, same PAW pseudos)
# for cohesive energies: E_coh = (E_bulk + Evdw)/Nfu - (E_cat + E_an).
set -u
PW=/tmp/qe-build/q-e-qe-7.2/bin/pw.x
PPDIR=/tmp/solids/pp
export OMP_NUM_THREADS=8
cd /tmp/solids; mkdir -p atoms; cd atoms

# sym  pseudo  tot_magnetization(unpaired e-)
ATOMS=(
"Li Li.pbe-s-kjpaw_psl.1.0.0.UPF 1"
"F  F.pbe-n-kjpaw_psl.1.0.0.UPF 1"
"Na Na.pbe-spn-kjpaw_psl.1.0.0.UPF 1"
"Cl Cl.pbe-n-kjpaw_psl.1.0.0.UPF 1"
"K  K.pbe-spn-kjpaw_psl.1.0.0.UPF 1"
"Mg Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF 0"
"Ca Ca.pbe-spn-kjpaw_psl.1.0.0.UPF 0"
"O  O.pbe-n-kjpaw_psl.1.0.0.UPF 2"
)
echo "# free-atom PBE energies (Ry)"
for row in "${ATOMS[@]}"; do
  read s p m <<< "$row"
  d=/tmp/solids/atoms/$s; mkdir -p $d/out; cd $d
  cat > a.in <<EOF
&control
 calculation='scf', prefix='$s', outdir='$d/out', pseudo_dir='$PPDIR'
/
&system
 ibrav=1, celldm(1)=24.0, nat=1, ntyp=1, ecutwfc=60, ecutrho=480,
 occupations='smearing', smearing='gaussian', degauss=0.01,
 nspin=2, tot_magnetization=$m, assume_isolated='mt'
/
&electrons
 conv_thr=1d-6, mixing_beta=0.3, electron_maxstep=300
/
ATOMIC_SPECIES
 $s 1.0 $p
ATOMIC_POSITIONS crystal
 $s 0.0 0.0 0.0
K_POINTS gamma
EOF
  $PW -in a.in > a.out 2>&1
  e=$(grep "! *total energy" a.out | tail -1 | awk '{print $5}')
  conv=$(grep -q "convergence has been achieved" a.out && echo ok || echo NOTCONV)
  printf "%-3s E= %-18s %s\n" "$s" "${e:-FAIL}" "$conv"
  cd /tmp/solids/atoms
done
echo "ATOMS DONE"
