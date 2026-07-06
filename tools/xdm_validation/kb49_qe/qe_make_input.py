#!/usr/bin/env python3
# xyz -> QE pw.x scf input: fixed cubic box (celldm bohr), MT-isolated, 80/800,
# PAW kjpaw_psl 1.0.0, gamma. Molecule centered. usage: qe_make_input.py <xyz> <a_bohr> <prefix> <outdir> <out.in>
import sys
xyz, a, prefix, outdir, out = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
PP = "/data/Iterative_hirshfeld/kb49_qe/pp"
pp = {"H":"H.pbe-kjpaw_psl.1.0.0.UPF","C":"C.pbe-n-kjpaw_psl.1.0.0.UPF",
      "N":"N.pbe-n-kjpaw_psl.1.0.0.UPF","O":"O.pbe-n-kjpaw_psl.1.0.0.UPF",
      "F":"F.pbe-n-kjpaw_psl.1.0.0.UPF","S":"S.pbe-n-kjpaw_psl.1.0.0.UPF",
      "Cl":"Cl.pbe-n-kjpaw_psl.1.0.0.UPF","Si":"Si.pbe-n-kjpaw_psl.1.0.0.UPF"}
mass = {"H":1.008,"C":12.011,"N":14.007,"O":15.999,"F":18.998,"S":32.06,"Cl":35.45,"Si":28.085}
L = open(xyz).read().split("\n"); n = int(L[0])
ats = []
for ln in L[2:2+n]:
    t = ln.split()
    if len(t) >= 4:
        ats.append((t[0], float(t[1]), float(t[2]), float(t[3])))
sym = sorted({s for s, *_ in ats})
cen = a * 0.52917720859 / 2.0
cx = sum(x for _,x,_,_ in ats)/len(ats); cy = sum(y for _,_,y,_ in ats)/len(ats); cz = sum(z for _,_,z,_ in ats)/len(ats)
with open(out, "w") as o:
    o.write(f'&control\n calculation="scf", prefix="{prefix}", outdir="{outdir}",\n')
    o.write(f' pseudo_dir="{PP}", tprnfor=.false., disk_io="low"\n/\n')
    o.write(f'&system\n ibrav=1, celldm(1)={a}, nat={len(ats)}, ntyp={len(sym)},\n')
    o.write(f' ecutwfc=80, ecutrho=800, assume_isolated="mt"\n/\n')
    o.write('&electrons\n conv_thr=1d-8, mixing_beta=0.3, electron_maxstep=200\n/\n')
    o.write("ATOMIC_SPECIES\n")
    for s in sym: o.write(f" {s} {mass[s]} {pp[s]}\n")
    o.write("ATOMIC_POSITIONS angstrom\n")
    for s,x,y,z in ats: o.write(f" {s} {x-cx+cen:.6f} {y-cy+cen:.6f} {z-cz+cen:.6f}\n")
    o.write("K_POINTS gamma\n")
print(f"{out}: nat={len(ats)} ntyp={len(sym)} a={a}bohr")
