#!/bin/bash
# Charge-aware XDM on the alkali-halide/oxide set (rocksalt B1).
# QE PBE-PAW SCF -> pp.x (valence rho, ELF) -> critic2 grid HI-XDM (neutral + A/B/C routes).
set -u
PW=/tmp/qe-build/q-e-qe-7.2/bin/pw.x
PP=/tmp/qe-build/q-e-qe-7.2/bin/pp.x
CR=/home/albd/critic2_dev/build/src/critic2
WFC=/home/albd/critic2_dev/dat/hirshfeld_proatoms/ld1_pbe
PPDIR=/tmp/solids/pp
export OMP_NUM_THREADS=6
cd /tmp/solids

# name a_Ang  catSym catZval catPP   anSym anZval anPP
SOLIDS=(
"LiF 4.0270 Li 3 Li.pbe-s-kjpaw_psl.1.0.0.UPF   F  7 F.pbe-n-kjpaw_psl.1.0.0.UPF"
"NaF 4.6342 Na 9 Na.pbe-spn-kjpaw_psl.1.0.0.UPF F  7 F.pbe-n-kjpaw_psl.1.0.0.UPF"
"NaCl 5.6402 Na 9 Na.pbe-spn-kjpaw_psl.1.0.0.UPF Cl 7 Cl.pbe-n-kjpaw_psl.1.0.0.UPF"
"KCl 6.2931 K 9 K.pbe-spn-kjpaw_psl.1.0.0.UPF   Cl 7 Cl.pbe-n-kjpaw_psl.1.0.0.UPF"
"MgO 4.2113 Mg 10 Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF O 6 O.pbe-n-kjpaw_psl.1.0.0.UPF"
"CaO 4.8105 Ca 10 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF O 6 O.pbe-n-kjpaw_psl.1.0.0.UPF"
)

printf "%-5s %-10s %-22s %-22s\n" SOLID ROUTE "charges(cat/an)" "Evdw(Ha)"
for row in "${SOLIDS[@]}"; do
  read name aA cS cZ cP aS aZ aP <<< "$row"
  d=/tmp/solids/$name; mkdir -p $d/out; cd $d
  abohr=$(python3 -c "print(f'{$aA/0.52917720859:.5f}')")
  cat > scf.in <<EOF
&control
 calculation='scf', prefix='$name', outdir='$d/out', pseudo_dir='$PPDIR', tprnfor=.true.
/
&system
 ibrav=2, celldm(1)=$abohr, nat=2, ntyp=2, ecutwfc=60, ecutrho=480, occupations='fixed'
/
&electrons
 conv_thr=1d-8, mixing_beta=0.3
/
ATOMIC_SPECIES
 $cS 1.0 $cP
 $aS 1.0 $aP
ATOMIC_POSITIONS crystal
 $cS 0.0 0.0 0.0
 $aS 0.5 0.5 0.5
K_POINTS automatic
 6 6 6 1 1 1
EOF
  if [ ! -f $name.val.cube ] || [ ! -f $name.elf.cube ]; then
    $PW -in scf.in > scf.out 2>&1
    if ! grep -q "JOB DONE" scf.out; then echo "$name SCF FAILED"; cd /tmp/solids; continue; fi
    for pn in 0 8; do tag=$([ $pn = 0 ] && echo val || echo elf)
      cat > pp_$tag.in <<EOF
&inputpp
 prefix='$name', outdir='$d/out', plot_num=$pn
/
&plot
 iflag=3, output_format=6, fileout='$name.$tag.cube'
/
EOF
      $PP -in pp_$tag.in > pp_$tag.out 2>&1
    done
  fi
  etot=$(grep "! *total energy" scf.out | tail -1 | awk '{print $5}')
  echo "# $name a=$aA Ang  E_scf=$etot Ry"
  for route in "neutral|" "gould|hirshfeld_i alpharef gould wfcdir $WFC" "scale|hirshfeld_i alpharef scale wfcdir $WFC" "stern|hirshfeld_i alpharef stern wfcdir $WFC"; do
    rname=${route%%|*}; rkw=${route#*|}
    cat > xdm_$rname.cri <<EOF
crystal $d/$name.val.cube
load $d/$name.val.cube id rho1 zpsp $cS $cZ $aS $aZ
load $d/$name.elf.cube id elf1
reference rho1
xdm grid rho rho1 elf elf1 xa1 0.4 xa2 2.5 $rkw
EOF
    $CR xdm_$rname.cri > xdm_$rname.out 2>&1
    ev=$(awk '/Evdw =/{print $3; exit}' xdm_$rname.out)
    ch=$(awk -v c="$cS" -v a="$aS" '/Stage 2|ALPHAREF/{f=1} f&&($2==c||$2==a)&&$3~/E/{printf "%s ",$3}' xdm_$rname.out | head -c 40)
    printf "%-5s %-8s %-22s %s\n" "$name" "$rname" "$ch" "$ev"
  done
  cd /tmp/solids
done
echo "ALL DONE"
