#!/usr/bin/env python3
# Partial directional analysis: for a set, build E_int per route from whatever
# species/routes are already cached, over the reactions that are FULLY cached.
import os, sys
SET = sys.argv[1]
DINMAP = {
  "il16": "/data/refdata/10_din-GMTKN55/il16.din",
  "chb6": "/data/refdata/10_din-GMTKN55/chb6.din",
  "ahb21": "/data/refdata/10_din-GMTKN55/ahb21.din",
  "s22": "/data/refdata/10_din/s22.din",
  "ionichb": "/data/refdata/10_din/ionichb.din",
  "s66": "/data/refdata/10_din/s66.din",
}
din = DINMAP[SET]; sd = f"/tmp/gmtkn/{SET}"
ROUTES = ["neutral", "gould", "scale", "stern"]
K = 627.5094740631
L = [l.strip() for l in open(din) if l.strip() and not l.strip().startswith('#')]
rx = []; i = 0
while i < len(L):
    t = []
    while i < len(L) and L[i] != '0':
        t.append((float(L[i]), L[i+1])); i += 2
    i += 1; rx.append((float(L[i]), t)); i += 1

def edft(s):
    p = f"{sd}/{s}.fchk.edft"
    return float(open(p).read()) if os.path.exists(p) else None

def edsp(s, r):
    p = f"{sd}/{s}.{r}.edisp"
    return float(open(p).read()) if os.path.exists(p) else None

hdr = f"{'rxn':<14}{'ref':>9}{'DFT':>9}" + "".join(f"{r:>9}" for r in ROUTES)
print(hdr)
errs = {r: [] for r in ROUTES}; ed_e = []; done = 0
for ref, t in rx:
    es = {s: edft(s) for _, s in t}
    if any(v is None for v in es.values()):
        continue
    dd = {r: {s: edsp(s, r) for _, s in t} for r in ROUTES}
    if any(dd[r][s] is None for r in ROUTES for _, s in t):
        continue
    done += 1
    ei_dft = sum(c * es[s] for c, s in t) * K
    row = {r: sum(c * (es[s] + dd[r][s]) for c, s in t) * K for r in ROUTES}
    print(f"{t[0][1]:<14}{ref:9.2f}{ei_dft:9.2f}" + "".join(f"{row[r]:9.2f}" for r in ROUTES))
    ed_e.append(ei_dft - ref)
    for r in ROUTES:
        errs[r].append(row[r] - ref)
print(f"--- {SET}: {done}/{len(rx)} reactions fully cached ---")
if done:
    mae = lambda e: sum(abs(x) for x in e) / len(e)
    mse = lambda e: sum(e) / len(e)
    print(f"{'DFT-only':<10} MAE={mae(ed_e):6.2f}  MSE={mse(ed_e):+6.2f}")
    for r in ROUTES:
        print(f"{r:<10} MAE={mae(errs[r]):6.2f}  MSE={mse(errs[r]):+6.2f}")
