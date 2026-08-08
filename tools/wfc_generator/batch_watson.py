#!/usr/bin/env python3
"""
batch_watson.py -- generate the Watson-sphere anion reference tables for critic2.

For each element Z, run the self-consistent Watson-sphere SCF (gen_ion_alpha_watson_scf)
for the neutral and anions q=-1,-2,-3, and emit the CONSISTENT REFERENCE PAIR:
    aws_mN(Z) = alpha_ws(Z,-N)            [a0^3, raw Watson-confined polarizability]
    vws_mN(Z) = <r^3>_ws(Z,-N)            [a0^3, r^3-moment of the same Watson density]
so critic2 can form  alpha_AIM = aws_mN * V_AIM / vws_mN.  Using the Watson-confined
V_ref (vws), NOT the diffuse ld1-box V_free(Q), avoids double-counting the crystal
compression the Watson radius already encodes (the S2- AIM->3.2 A^3 bug).

Also emits rstern_ws_mN = alpha_ws(Z,-N)/alpha_ws(Z,0) (kept for reference/back-compat).

Usage:  python3 batch_watson.py [ZMAX]   (default 54; H..Xe)
Unphysical over-fills (electron forced into the next shell: F2-, S3-, Kr-) are set to
0 => critic2 falls back to the plain ratio route / neutral scaling.
"""
import sys, os
import numpy as np
import gen_ion_alpha_watson_scf as W
import gen_ion_alpha_sternheimer as G

ZMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 54
WFCDIR = os.path.expanduser("~/critic2_dev/dat/wfc")
MAXZ = 123

def a_ws(sym, q):
    """Return (alpha_a0^3, V_ws_a0^3) or (None,None)."""
    try:
        a, conf, Rw, vws = W.compute(sym, q, WFCDIR, rmax=26.0, npts=3640)
        if a is None or not np.isfinite(a) or a <= 0:
            return None, None
        return a, vws
    except Exception:
        return None, None

aws = {1:[0.0]*(MAXZ+1), 2:[0.0]*(MAXZ+1), 3:[0.0]*(MAXZ+1)}
vws = {1:[0.0]*(MAXZ+1), 2:[0.0]*(MAXZ+1), 3:[0.0]*(MAXZ+1)}
rm  = {1:[1.0]*(MAXZ+1), 2:[1.0]*(MAXZ+1), 3:[1.0]*(MAXZ+1)}

for Z in range(1, ZMAX+1):
    sym = G.ELEM[Z]
    a0, v0 = a_ws(sym, 0)
    if a0 is None:
        print(f"{sym:2} Z={Z:2}  neutral FAILED"); continue
    line = f"{sym:2} Z={Z:2}  a0={a0:8.3f} v0={v0:8.2f}"
    for N in (1, 2, 3):
        aq, vq = a_ws(sym, -N)
        # physical band: a real anion completes the p/s shell; over-fills blow up.
        if aq is not None and 0.25 < aq/a0 < 3.0:
            aws[N][Z] = aq; vws[N][Z] = vq; rm[N][Z] = aq/a0
            line += f"  m{N}: a={aq:7.2f} V={vq:7.1f} r={aq/a0:.3f}"
        elif aq is not None:
            line += f"  m{N}:[skip {aq/a0:.1f}]"
    print(line, flush=True)

def emit(f, name, arr, fmt="{:.5g}"):
    f.write(f"\n  real*8, parameter :: {name}(1:maxzat0) = (/&\n")
    vals = [f"{fmt.format(arr[z])}d0" for z in range(1, MAXZ+1)]
    for i in range(0, len(vals), 6):
        end = ",&" if i+6 < len(vals) else " /)"
        f.write("      " + ", ".join(vals[i:i+6]) + end + "\n")

with open("/tmp/watson/rstern_ws_fortran.txt", "w") as f:
    for N in (1,2,3): emit(f, f"aws_m{N}", aws[N])
    for N in (1,2,3): emit(f, f"vws_m{N}", vws[N])
    for N in (1,2,3): emit(f, f"rstern_ws_m{N}", rm[N])
print("\nwrote /tmp/watson/rstern_ws_fortran.txt  (aws_mN, vws_mN, rstern_ws_mN)")
