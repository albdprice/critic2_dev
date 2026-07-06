#!/usr/bin/env python3
"""Fit XDM BJ-damping (a1,a2) per route on KB49 from the per-species JSONs in
results/<route>/, minimizing RMSP via least_squares -- self-contained replica of
Alberto/AP's 02_collate_and_fit convention (no external imports beyond numpy/scipy).
usage: python fit_kb49_qe.py [results_dir] [kb49.din]"""
import sys, os, json, glob
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.distance import pdist, squareform

HARTREE_TO_KCAL = 627.50947
BOHR = 0.52917720859
ROUTES = ["neutral", "gould", "scale", "stern"]

def parse_din(path):
    rx = []; L = [l.strip() for l in open(path) if l.strip() and not l.startswith('#')]
    i = 0
    while i < len(L):
        if L[i] == "-111": break
        comp = [];
        while L[i] != "0":
            comp.append((float(L[i]), L[i+1])); i += 2
        i += 1; ref = float(L[i]); i += 1
        rx.append({"components": comp, "ref_e": ref, "name": comp[0][1]})
    return rx

def edisp(d, a1, a2):
    if not d["c6"]: return 0.0
    coords = np.array(d["coords"]); c6=np.array(d["c6"]); c8=np.array(d["c8"]); c10=np.array(d["c10"]); rc=np.array(d["rc"])
    a2b = a2 / BOHR
    dd = squareform(pdist(coords)); rv = a1*rc + a2b
    i, j = np.triu_indices(len(coords), k=1)
    return -np.sum(c6[i,j]/(rv[i,j]**6+dd[i,j]**6) + c8[i,j]/(rv[i,j]**8+dd[i,j]**8) + c10[i,j]/(rv[i,j]**10+dd[i,j]**10))

def fit_route(route, rdir, din):
    td = {}
    for r in din:
        for _, n in r["components"]:
            p = f"{rdir}/{route}/{n}.json"
            if n not in td and os.path.exists(p) and os.path.getsize(p) > 0:
                try: td[n] = json.load(open(p))
                except Exception: pass
    valid = [r for r in din if all(n in td for _, n in r["components"])]
    if not valid: return None
    def resid(p):
        a1, a2 = p; out = []
        for r in valid:
            cs = r["components"]
            ib = sum(c*td[n]["base_energy"] for c, n in cs)
            idsp = sum(c*edisp(td[n], a1, a2) for c, n in cs)
            out.append(((ib+idsp)*HARTREE_TO_KCAL - r["ref_e"]) / r["ref_e"])
        return np.array(out)
    res = least_squares(resid, [0.0, 1.4545], bounds=([0,0],[np.inf,np.inf]), method="trf")
    a1, a2 = res.x
    if a1 <= 1e-5: a1 = 0.0
    ae = []; pe = []
    for r in valid:
        cs = r["components"]
        calc = (sum(c*td[n]["base_energy"] for c,n in cs)+sum(c*edisp(td[n],a1,a2) for c,n in cs))*HARTREE_TO_KCAL
        ae.append(abs(calc-r["ref_e"])); pe.append(abs((calc-r["ref_e"])/r["ref_e"])*100)
    return dict(route=route, n=len(valid), a1=a1, a2=a2, mae=float(np.mean(ae)), mapd=float(np.mean(pe)))

if __name__ == "__main__":
    rdir = sys.argv[1] if len(sys.argv) > 1 else "results"
    din  = sys.argv[2] if len(sys.argv) > 2 else "kb49.din"
    D = parse_din(din)
    print("# KB49 charge-aware XDM refit (QE PBE plane-wave, ecutwfc80/ecutrho800, MT box)")
    print(f"# {'route':<8} {'n':>4} {'a1':>8} {'a2(A)':>8} {'MAE(kcal)':>10} {'MAPD%':>8}")
    res = {}
    for rt in ROUTES:
        s = fit_route(rt, rdir, D)
        if s:
            res[rt] = s
            print(f"  {s['route']:<8} {s['n']:>4} {s['a1']:>8.4f} {s['a2']:>8.4f} {s['mae']:>10.4f} {s['mapd']:>7.2f}%")
    json.dump(res, open(f"{rdir}/../fit_results_qe.json", "w"), indent=2)
