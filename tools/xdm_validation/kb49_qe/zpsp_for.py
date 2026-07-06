#!/usr/bin/env python3
# Emit critic2 ZPSP token list ("C 4 H 1 ...") for the elements in an xyz,
# matching the kjpaw_psl 1.0.0 PBE valence charges.
import sys
ZV = {"H":1,"C":4,"N":5,"O":6,"F":7,"S":6,"Cl":7,"Si":4}
L = open(sys.argv[1]).read().split("\n"); n = int(L[0])
els = []
for ln in L[2:2+n]:
    t = ln.split()
    if len(t) >= 4 and t[0] not in els:
        els.append(t[0])
print(" ".join(f"{e} {ZV[e]}" for e in els))
