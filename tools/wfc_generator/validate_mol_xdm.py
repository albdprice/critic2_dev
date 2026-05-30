#!/usr/bin/env python3
"""Molecular A/B/C validation harness: run neutral + the four charge-aware
routes for each molecule, extract the molecular C6 (full double sum over
atom pairs, the homomolecular dispersion coefficient) and E_disp, and
compare the molecular C6 to reference DOSD values where available."""
import subprocess, re, os, sys

CRITIC = os.path.expanduser("~/critic2_dev/build/src/critic2")
WFC = os.path.expanduser("~/critic2_dev/dat/hirshfeld_proatoms/ld1_pbe")
FDIR = "/tmp/xdmtest"

# reference homomolecular C6 (a.u.), dipole-oscillator-strength distributions
REF_C6 = {"h2o": 45.3, "CH4": 129.7}   # Meath et al. / standard XDM validation set

ROUTES = [("neutral", "pbe"),
          ("HI",       "pbe hirshfeld_i wfcdir %s" % WFC),
          ("2A-gould", "pbe hirshfeld_i alpharef gould wfcdir %s" % WFC),
          ("2B-scale", "pbe hirshfeld_i alpharef scale wfcdir %s" % WFC),
          ("2C-kirk",  "pbe hirshfeld_i alpharef compute wfcdir %s" % WFC),
          ("2C-stern", "pbe hirshfeld_i alpharef stern wfcdir %s" % WFC)]

def fchk(mol):
    for c in (mol, mol.upper(), mol.capitalize()):
        p = os.path.join(FDIR, c + ".fchk")
        if os.path.exists(p): return p
    raise FileNotFoundError(mol)

def run(mol, tail):
    f = fchk(mol)
    inp = (f"molecule {f}\nload {f} id mol\nreference mol\n"
           f"meshtype franchini small\nxdm 0.4 2.5 {tail}\n")
    cri = f"/tmp/val_{mol}.cri"
    open(cri, "w").write(inp)
    out = subprocess.run([CRITIC, cri], capture_output=True, text=True).stdout
    # molecular C6 = sum over all ordered atom pairs of C6_ij (off-diag x2)
    c6 = {}
    inblk = False
    for ln in out.splitlines():
        if ln.strip().startswith("# i  j"): inblk = True; continue
        if inblk:
            m = re.match(r"\s*(\d+)\s+(\d+)\s+([-0-9.E+]+)\s+[-0-9.E+]+", ln)
            if m:
                i, j, v = int(m.group(1)), int(m.group(2)), float(m.group(3))
                c6[(i, j)] = v
            elif ln.strip() == "#" or "contribution" in ln:
                inblk = False
    molc6 = 0.0
    for (i, j), v in c6.items():
        molc6 += v if i == j else 2.0 * v
    ed = re.search(r"dispersion energy \(Ha\)\s+([-0-9.E+]+)", out)
    edisp = float(ed.group(1)) if ed else float("nan")
    return molc6, edisp

mols = sys.argv[1:] or ["h2o", "CH4", "NaCl", "LiF"]
print(f"{'molecule':9s} {'route':9s} {'mol_C6(au)':>12s} {'E_disp(Ha)':>14s} {'C6/ref':>8s}")
for mol in mols:
    ref = REF_C6.get(mol)
    for name, tail in ROUTES:
        try:
            c6, ed = run(mol, tail)
            rr = f"{c6/ref:6.3f}" if ref else "   -  "
            print(f"{mol:9s} {name:9s} {c6:12.2f} {ed:14.4e} {rr:>8s}")
        except Exception as e:
            print(f"{mol:9s} {name:9s}  ERROR {e}")
    if ref: print(f"{'':9s} {'REF-DOSD':9s} {ref:12.2f}")
    print()
