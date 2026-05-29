#!/usr/bin/env python3
# Robustness / demonstration test: plain Hirshfeld vs Hirshfeld-I across a set
# of molecules, focusing on ionic systems where plain Hirshfeld badly
# under-polarizes. PBE/def2-TZVP density (g09), Route-2 ld1_pbe anion refs.
import os, re, subprocess, sys

G = os.path.expanduser("~/g09")
ENV = dict(os.environ, GAUSS_EXEDIR=G, GAUSS_BSDDIR=G+"/bsd",
           LD_LIBRARY_PATH=G+":"+G+"/bsd", GAUSS_SCRDIR="/tmp",
           PATH=G+":"+os.environ["PATH"])
CRIT = os.path.expanduser("~/critic2_dev/build/src/critic2")
WFCDIR = os.path.expanduser("~/critic2_dev/dat/hirshfeld_proatoms/ld1_pbe")
WD = "/tmp/moltest"; os.makedirs(WD, exist_ok=True)
ZOF = {"H":1,"Li":3,"C":6,"N":7,"O":8,"F":9,"Na":11,"Cl":17}

# name: (charge, mult, geometry block in zmatrix or cartesian)
MOLS = {
 "LiF": (0,1,"Li\nF 1 1.56"),
 "LiH": (0,1,"Li\nH 1 1.60"),
 "NaF": (0,1,"Na\nF 1 1.93"),
 "NaCl":(0,1,"Na\nCl 1 2.36"),
 "HF":  (0,1,"H\nF 1 0.92"),
 "CO":  (0,1,"C\nO 1 1.13"),
 "H2O": (0,1,"O\nH 1 0.96\nH 1 0.96 2 104.5"),
 "NH3": (0,1,"N\nH 1 1.01\nH 1 1.01 2 107.\nH 1 1.01 2 107. 3 120."),
 "CH4": (0,1,"C\nH 1 1.09\nH 1 1.09 2 109.47\nH 1 1.09 2 109.47 3 120.\nH 1 1.09 2 109.47 3 -120."),
}

def run(cmd, **kw): return subprocess.run(cmd, shell=True, env=ENV, capture_output=True, text=True, **kw)

def make_density(name, chg, mult, geom):
    gjf = f"{WD}/{name}.gjf"
    open(gjf,"w").write(
        f"%chk={WD}/{name}.chk\n%mem=2GB\n%nproc=4\n"
        f"#p PBEPBE/def2TZVP opt density=current\n\n{name}\n\n{chg} {mult}\n{geom}\n\n")
    r = run(f"{G}/g09 < {gjf} > {WD}/{name}.log 2>&1")
    if "Normal termination" not in open(f"{WD}/{name}.log").read():
        return None
    run(f"{G}/formchk {WD}/{name}.chk {WD}/{name}.fchk")
    run(f"{G}/cubegen 4 fdensity=scf {WD}/{name}.fchk {WD}/{name}.cube 200 h")
    return f"{WD}/{name}.cube"

def parse_pop(out):
    """Return list of (Z, Pop) from the integrated-atomic-properties table."""
    rows=[]
    for l in out.splitlines():
        m = re.match(r"\s*\d+\s+\d+\s+\d+\s+([A-Za-z]+)\s+(\d+)\s+\S+\s+([-0-9.E+]+)", l)
        if m:
            rows.append((m.group(1), int(m.group(2)), float(m.group(3))))
    return rows

def charges(cube, method):
    cri=f"{WD}/run_{method}.cri"
    body = f"molecule {cube}\nload {cube}\nintegrable 1\n"
    body += "hirshfeld\n" if method=="h" else f"hirshfeld_i wfcdir {WFCDIR}\n"
    open(cri,"w").write(body)
    out = run(f"{CRIT} {cri}").stdout
    # take the LAST integrated-atomic-properties table
    rows = parse_pop(out)
    return rows  # [(name,Z,Pop),...]

print(f"{'mol':5} {'atom':4} {'q_Hirsh':>9} {'q_HirshI':>9}   note")
for name,(chg,mult,geom) in MOLS.items():
    cube = make_density(name,chg,mult,geom)
    if not cube:
        print(f"{name:5} --- SCF/opt failed"); continue
    h  = charges(cube,"h")
    hi = charges(cube,"hi")
    # aggregate per unique atom name keeping order; here just per-atom rows
    for (nm,Z,popH),(nm2,Z2,popI) in zip(h,hi):
        qH=Z-popH; qI=Z2-popI
        print(f"{name:5} {nm:4} {qH:9.3f} {qI:9.3f}")
    print()
