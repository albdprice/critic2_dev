#!/usr/bin/env python3
"""KB49 charge-aware XDM refit pipeline.
Per species: Psi4 PBE/def2-tzvp -> fchk + base_energy (cached). Per route:
critic2 xdm -> parse coefficient block -> AP-format JSON (cached). Then fit
a1/a2 per route (neutral/gould/scale/stern) by minimizing RMSP, exactly as
AP's 02_collate_and_fit.py. Resumable + fault-tolerant."""
import os, sys, subprocess, time, json, traceback
from concurrent.futures import ThreadPoolExecutor, as_completed

WORK   = "/tmp/kb49"
GEOM   = "/home/albd/projects/refdata/20_kb49"
DIN    = "/home/albd/projects/refdata/10_din/kb49.din"
PSI4PY = "/tmp/gmtkn/psi4_species.py"
MKJSON = "/tmp/gmtkn/kb49_makejson.py"
CR     = "/home/albd/critic2_dev/build/src/critic2"
WFC    = "/home/albd/critic2_dev/dat/hirshfeld_proatoms/ld1_pbe"
sys.path.insert(0, "/data/XDM_Psi4")
from xdm_lib import parse_universal_din, calc_bj_dispersion, HARTREE_TO_KCAL
from scipy.optimize import least_squares
import numpy as np

ROUTES = [("neutral", ""),
          ("gould",  f"hirshfeld_i alpharef gould wfcdir {WFC}"),
          ("scale",  f"hirshfeld_i alpharef scale wfcdir {WFC}"),
          ("stern",  f"hirshfeld_i alpharef stern wfcdir {WFC}")]
NWORK, MEMGB, NTHR = 4, "7", 3
A1, A2 = "0.4041", "2.6998"          # arbitrary; coefficients are a1/a2-independent

def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def run_psi4(sp):
    fchk = f"{WORK}/{sp}.fchk"; edft = fchk + ".edft"
    if os.path.exists(edft):
        return
    xyz = f"{GEOM}/{sp}.xyz"
    subprocess.run([sys.executable, PSI4PY, xyz, fchk, MEMGB, str(NTHR)],
                   capture_output=True, text=True, timeout=5400)
    if not os.path.exists(edft):
        raise RuntimeError(f"psi4 failed {sp}")

def run_route(sp, route, kw):
    js = f"{WORK}/json/{route}/{sp}.json"
    if os.path.exists(js):
        return
    fchk = f"{WORK}/{sp}.fchk"; cri = f"{WORK}/{sp}_{route}.cri"; out = f"{WORK}/{sp}_{route}.out"
    with open(cri, 'w') as f:
        f.write(f"molecule {fchk}\nload {fchk} id mol\nreference mol\n"
                f"meshtype franchini small\nxdm {A1} {A2} pbe {kw}\n")
    with open(out, 'w') as fo:
        subprocess.run([CR, cri], stdout=fo, stderr=subprocess.STDOUT, timeout=3600)
    r = subprocess.run([sys.executable, MKJSON, out, fchk, fchk + ".edft", js],
                       capture_output=True, text=True, timeout=120)
    if not os.path.exists(js):
        raise RuntimeError(f"makejson failed {sp}/{route}: {r.stderr[-200:]}")

def do_species(sp):
    try:
        run_psi4(sp)
        for rn, kw in ROUTES:
            run_route(sp, rn, kw)
        log(f"  ok  {sp}")
        return sp, True
    except Exception as e:
        log(f"  FAIL {sp}: {e}")
        return sp, False

def fit_route(route, din_data):
    """Load JSONs for a route, fit a1/a2 (RMSP), return summary dict."""
    td = {}
    for r in din_data:
        for _, n in r['components']:
            if n not in td:
                p = f"{WORK}/json/{route}/{n}.json"
                if os.path.exists(p) and os.path.getsize(p) > 0:
                    try:
                        td[n] = json.load(open(p))
                    except Exception as e:
                        log(f"  bad JSON {route}/{n}: {e}")
    valid = [r for r in din_data if all(n in td for _, n in r['components'])]
    if not valid:
        return None
    def resid(params):
        a1, a2 = params; out = []
        for r in valid:
            cs = r['components']
            ib = sum(c * td[n]['base_energy'] for c, n in cs)
            idsp = sum(c * calc_bj_dispersion(td[n], a1, a2) for c, n in cs)
            calc = (ib + idsp) * HARTREE_TO_KCAL
            out.append((calc - r['ref_e']) / r['ref_e'])
        return np.array(out)
    res = least_squares(resid, [0.0, 1.4545], bounds=([0.0, 0.0], [np.inf, np.inf]), method='trf')
    a1o, a2o = res.x
    if a1o <= 1e-5: a1o = 0.0
    # stats at fitted params + at the current default (0.4041/2.6998)
    def stats(a1, a2):
        ae, pe = [], []
        for r in valid:
            cs = r['components']
            ib = sum(c * td[n]['base_energy'] for c, n in cs)
            idsp = sum(c * calc_bj_dispersion(td[n], a1, a2) for c, n in cs)
            calc = (ib + idsp) * HARTREE_TO_KCAL
            ae.append(abs(calc - r['ref_e'])); pe.append(abs((calc - r['ref_e']) / r['ref_e']) * 100)
        return np.mean(ae), np.mean(pe)
    mae_f, mapd_f = stats(a1o, a2o)
    mae_d, mapd_d = stats(0.4041, 2.6998)
    return dict(route=route, n=len(valid), a1=a1o, a2=a2o,
                mae_fit=mae_f, mapd_fit=mapd_f, mae_def=mae_d, mapd_def=mapd_d)

def main():
    os.makedirs(WORK, exist_ok=True)
    for rn, _ in ROUTES:
        os.makedirs(f"{WORK}/json/{rn}", exist_ok=True)
    din = parse_universal_din(DIN)
    species = sorted({n for r in din for _, n in r['components']})
    log(f"=== KB49: {len(din)} reactions, {len(species)} species, {NWORK}x{NTHR} ===")
    with ThreadPoolExecutor(max_workers=NWORK) as ex:
        futs = {ex.submit(do_species, s): s for s in species}
        for f in as_completed(futs):
            f.result()
    log("=== all species done; fitting ===")
    print("\n# KB49 charge-aware XDM refit (PBE/def2-TZVP)")
    print(f"# {'route':<8} {'n':>4} {'a1':>8} {'a2(A)':>8}  {'MAE_fit':>8} {'MAPD_fit':>9}  {'MAE@def':>8} {'MAPD@def':>9}")
    res = {}
    for rn, _ in ROUTES:
        s = fit_route(rn, din)
        if s:
            res[rn] = s
            print(f"  {s['route']:<8} {s['n']:>4} {s['a1']:>8.4f} {s['a2']:>8.4f}  "
                  f"{s['mae_fit']:>8.4f} {s['mapd_fit']:>8.2f}%  {s['mae_def']:>8.4f} {s['mapd_def']:>8.2f}%")
    json.dump(res, open(f"{WORK}/fit_results.json", "w"), indent=2)
    log("=== KB49 DONE ===")

if __name__ == "__main__":
    main()
