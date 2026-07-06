#!/usr/bin/env python3
"""
batch_sternheimer.py -- batch-drive gen_ion_alpha_sternheimer over the
periodic table to build the Stage-2C-rigorous reference.

For each element Z it uses rmax = 3.6 * R99(neutral) -- the SAME standardized
confinement box used to generate the Hirshfeld-I density references -- and
computes the (uncoupled) Sternheimer static polarizability for charges
q = 0, +1, -1 (and -2 as a DIAGNOSTIC only -- see below).

*** FINDING (2026-06): q=-2 is NOT ingested into param.F90. ***
The uncoupled Sternheimer diverges for double anions: the second extra electron
of a free A2- is unbound, so its linear-response polarizability is ill-defined.
Verified numerically -- confined O2- alpha = {982,-31,0.99,-0.07,107,1.7} across
rmax={6,8,10,12,14,16}; the whole -2 column of this batch is garbage (huge/negative).
rstern_m2 is emitted for the record but critic2 clamps ion_alpha_stern at q=-1.
Multiply-charged anions are instead handled by the DENSITY-based routes
(alpharef compute / scale), which use the bound confined ground-state density
(well-defined; e.g. Kirkwood <r2>^2/N for O is 17.5 -> 64 -> 82 for q=0,-1,-2,
monotonic and finite). See doc/research/xdm_hirshfeld_i_notebook.md.

It stores the neutral-CALIBRATED ratios
    rstern(Z,q) = alpha_stern(Z,q) / alpha_stern(Z,0)
so the charge-aware free-ion polarizability used by critic2 is
    alpha_free(Z,q) = alpha_free^CRC(Z) * rstern(Z,q)
(the uncoupled absolute overestimate cancels in the ratio; validated against
the Gould-Bucko TDDFT trend -- cations near-exact). Writes:
  - a human/provenance table (element charge alpha rmax ratio)
  - Fortran data arrays rstern_m1, rstern_p1 for param.F90.
"""
import os, sys, math
import numpy as np
import gen_ion_alpha_sternheimer as g

WFCDIR = os.path.expanduser("~/critic2_dev/dat/wfc")
ZMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 54
OUT = "/tmp/stern/sternheimer_alpha.dat"

rows = []   # (Z, sym, q, alpha, rmax, ratio)
a0 = {}     # Z -> alpha_stern(neutral)
rmax_of = {}
for Z in range(1, ZMAX + 1):
    sym = g.ELEM[Z] if Z < len(g.ELEM) else None
    if sym is None:
        continue
    try:
        rmax = 3.6 * g.r99_neutral(sym, WFCDIR)
    except Exception as e:
        print(f"# {sym}: no neutral wfc ({e}); skip", file=sys.stderr); continue
    rmax_of[Z] = rmax
    a, ok, msg = g.compute_alpha(sym, 0, rmax)
    if not ok:
        print(f"# {sym} q0 FAILED: {msg}", file=sys.stderr); continue
    a0[Z] = a
    rows.append((Z, sym, 0, a, rmax, 1.0))
    for q in (+1, -1, -2):
        aq, ok, msg = g.compute_alpha(sym, q, rmax)
        if not ok:
            print(f"# {sym} q{q:+d} FAILED: {msg}", file=sys.stderr); continue
        rows.append((Z, sym, q, aq, rmax, aq / a))
    gr = lambda qq: next((rt for (zz,_,q,_,_,rt) in rows if zz==Z and q==qq), float('nan'))
    print(f"{sym:>2} Z={Z:>2} rmax={rmax:6.2f} a0={a0.get(Z,float('nan')):10.3f} "
          f"r(+1)={gr(1):.4f} r(-1)={gr(-1):.4f} r(-2)={gr(-2):.4f}", flush=True)

# provenance table
with open(OUT, "w") as f:
    f.write("# Sternheimer (uncoupled CPKS) confined-ion static polarizabilities\n")
    f.write("# rmax = 3.6*R99(neutral) (the HI density-reference box). alpha in a0^3.\n")
    f.write("# ratio = alpha(Z,q)/alpha(Z,0) -- the neutral-calibrated charge factor.\n")
    f.write("# Z sym q alpha_a0^3 rmax ratio\n")
    for (Z, sym, q, a, rmax, rt) in sorted(rows):
        f.write(f"{Z:3d} {sym:>2} {q:+d} {a:12.4f} {rmax:8.3f} {rt:10.5f}\n")
print(f"\nwrote {OUT}", file=sys.stderr)

# Fortran arrays (default ratio 1.0 = no charge correction).
# rstern_m2 falls back to the -1 ratio where no bound -2 datum exists, so that
# ion_alpha_stern interpolating -1 -> -2 degrades to the (correct) -1 clamp there.
rp1 = [1.0]*(ZMAX+1); rm1 = [1.0]*(ZMAX+1); m2 = {}
for (Z, sym, q, a, rmax, rt) in rows:
    if q == 1: rp1[Z] = rt
    elif q == -1: rm1[Z] = rt
    elif q == -2: m2[Z] = rt
rm2 = list(rm1)
for Z, rt in m2.items():
    rm2[Z] = rt
def emit(name, arr):
    body = []
    vals = [f"{arr[z]:.5g}d0" for z in range(1, ZMAX+1)]
    for i in range(0, len(vals), 6):
        body.append("      " + ", ".join(vals[i:i+6]) + ",&")
    s = "\n".join(body).rstrip(",&") + " /)"
    return f"  real*8, parameter :: {name}(1:maxzat0) = (/&\n{s}\n"
with open("/tmp/stern/sternheimer_fortran.txt", "w") as f:
    f.write(emit("rstern_p1", rp1))
    f.write("\n")
    f.write(emit("rstern_m1", rm1))
    f.write("\n")
    f.write(emit("rstern_m2", rm2))
print("wrote /tmp/stern/sternheimer_fortran.txt", file=sys.stderr)
