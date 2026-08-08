#!/usr/bin/env python3
"""
batch_watson.py -- generate the rstern_ws anion ratio tables for critic2.

For each element Z, run the self-consistent Watson-sphere SCF polarizability
(gen_ion_alpha_watson_scf) for the neutral and for anions q=-1,-2,-3, and emit
    rstern_ws_mN(Z) = alpha_ws(Z,-N) / alpha_ws(Z,0)
as Fortran parameter arrays (padded to maxzat0 with 1d0).  These feed the
`alpharef sternws` route: cations keep rstern_p1/p2 (bound, fine); anions of
charge <= -1 use the Watson-stabilized response, which is stable and matches
Tessman in-crystal alpha for the deep anions (O2-, S2-).

Usage:  python3 batch_watson.py [ZMAX]   (default 54; H..Xe)
Only p-block anion formers give meaningful ratios; closed/ill-defined cases -> 1.
"""
import sys, math
import numpy as np
import gen_ion_alpha_watson_scf as W
import gen_ion_alpha_sternheimer as G

ZMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 54
WFCDIR = __import__("os").path.expanduser("~/critic2_dev/dat/wfc")
MAXZ = 123

def a_ws(sym, q):
    try:
        a, conf, Rw = W.compute(sym, q, WFCDIR, rmax=26.0, npts=3640)
        if a is None or not np.isfinite(a) or a <= 0:
            return None
        return a
    except Exception:
        return None

rm = {1: [1.0]*(MAXZ+1), 2: [1.0]*(MAXZ+1), 3: [1.0]*(MAXZ+1)}
alpha0 = {}
for Z in range(1, ZMAX+1):
    sym = G.ELEM[Z]
    a0 = a_ws(sym, 0)
    alpha0[Z] = a0
    if a0 is None:
        print(f"{sym:2} Z={Z:2}  neutral FAILED"); continue
    line = f"{sym:2} Z={Z:2}  a0={a0:8.3f}"
    for N in (1, 2, 3):
        aq = a_ws(sym, -N)
        # A physical anion completes the current p (or s) shell; forcing an
        # electron into the NEXT shell (e.g. Br2-, Kr-, Se3-) is unbound even in
        # the Watson sphere and the response blows up. Guard: keep only ratios in
        # a physical band; otherwise leave 1.0 (critic2 falls back to neutral).
        if aq is not None and 0.25 < aq/a0 < 3.0:
            rm[N][Z] = aq/a0
            line += f"  m{N}={aq/a0:6.3f}(a={aq:7.2f})"
        elif aq is not None:
            line += f"  m{N}=[skip {aq/a0:.1f}]"
    print(line, flush=True)

def emit(name, arr):
    print(f"\n  real*8, parameter :: {name}(1:maxzat0) = (/&")
    vals = [f"{arr[z]:.5g}d0" for z in range(1, MAXZ+1)]
    for i in range(0, len(vals), 6):
        chunk = ", ".join(vals[i:i+6])
        end = ",&" if i+6 < len(vals) else " /)"
        print(f"      {chunk}{end}")

with open("/tmp/watson/rstern_ws_fortran.txt", "w") as f:
    import contextlib
    with contextlib.redirect_stdout(f):
        emit("rstern_ws_m1", rm[1]); emit("rstern_ws_m2", rm[2]); emit("rstern_ws_m3", rm[3])
print("\nwrote /tmp/watson/rstern_ws_fortran.txt")
