#!/usr/bin/env python3
# Charge-aware XDM batch over GMTKN55-style .din sets.
# For each set: parse reactions, run Psi4 PBE/def2-TZVP per unique species
# (cached fchk+energy), run critic2 xdm_wfn for 4 routes (cached Edisp),
# then build interaction energies per the .din stoichiometry and report
# MAE per route vs the high-level reference. Resumable + fault-tolerant.
import os, sys, subprocess, time, traceback
from concurrent.futures import ThreadPoolExecutor, as_completed

WORK   = "/tmp/gmtkn"
PSI4PY = f"{WORK}/psi4_species.py"
CR     = "/home/albd/critic2_dev/build/src/critic2"
WFC    = "/home/albd/critic2_dev/dat/hirshfeld_proatoms/ld1_pbe"
A1, A2 = "0.4041", "2.6998"          # AP cc-pVTZ-fitted PBE-XDM (~def2-tzvp)
K      = 627.5094740631             # Ha -> kcal/mol
NWORK  = 6                          # concurrent species (2 threads / 4 GB each)
MEMGB, NTHR = "4", 2
ROUTES = [("neutral", ""),
          ("gould",  f"hirshfeld_i alpharef gould wfcdir {WFC}"),
          ("scale",  f"hirshfeld_i alpharef scale wfcdir {WFC}"),
          ("stern",  f"hirshfeld_i alpharef stern wfcdir {WFC}")]

def log(msg):
    t = time.strftime("%H:%M:%S")
    print(f"[{t}] {msg}", flush=True)

def parse_din(path):
    """Return list of reactions: [(ref_kcal, [(coeff, species), ...]), ...]."""
    lines = [l.strip() for l in open(path) if l.strip() and not l.strip().startswith('#')]
    rxns, i = [], 0
    while i < len(lines):
        terms = []
        while i < len(lines) and lines[i] != '0':
            coeff = float(lines[i]); name = lines[i+1]; i += 2
            terms.append((coeff, name))
        i += 1                       # skip the '0'
        ref = float(lines[i]); i += 1
        rxns.append((ref, terms))
    return rxns

def run_psi4(sp, geomdir, setdir):
    fchk = f"{setdir}/{sp}.fchk"; edft = fchk + ".edft"
    if os.path.exists(edft):
        return float(open(edft).read().strip())
    xyz = f"{geomdir}/{sp}.xyz"
    if not os.path.exists(xyz):
        raise FileNotFoundError(xyz)
    r = subprocess.run([sys.executable, PSI4PY, xyz, fchk, MEMGB, str(NTHR)],
                       capture_output=True, text=True, timeout=5400)
    if not os.path.exists(edft):
        raise RuntimeError(f"psi4 failed {sp}: {r.stderr[-400:]}")
    return float(open(edft).read().strip())

def run_critic2(sp, route, kw, setdir):
    cache = f"{setdir}/{sp}.{route}.edisp"
    if os.path.exists(cache):
        return float(open(cache).read().strip())
    fchk = f"{setdir}/{sp}.fchk"; cri = f"{setdir}/{sp}_{route}.cri"
    with open(cri, 'w') as f:
        f.write(f"molecule {fchk}\nload {fchk} id mol\nreference mol\n"
                f"meshtype franchini small\nxdm {A1} {A2} pbe {kw}\n")
    r = subprocess.run([CR, cri], capture_output=True, text=True, timeout=1800)
    ed = 0.0
    for ln in r.stdout.splitlines():
        if "dispersion energy (Ha)" in ln:
            ed = float(ln.split()[3]); break
    else:
        # single-atom species legitimately have no dispersion energy line
        if "is a single atom" not in r.stdout and nat_of(fchk) > 1:
            raise RuntimeError(f"critic2 no Edisp {sp}/{route}: {r.stdout[-300:]}")
    open(cache, 'w').write(f"{ed:.12e}\n")
    return ed

def nat_of(fchk):
    try:
        for ln in open(fchk):
            if ln.startswith("Number of atoms"):
                return int(ln.split()[-1])
    except Exception:
        pass
    return 2

def do_species(sp, geomdir, setdir):
    """Psi4 then all routes for one species. Returns (sp, edft, {route:edisp}) or (sp, None, err)."""
    try:
        edft = run_psi4(sp, geomdir, setdir)
        ed = {rn: run_critic2(sp, rn, kw, setdir) for rn, kw in ROUTES}
        log(f"  ok  {sp}  EDFT={edft:.6f}  Edisp={ {k:round(v,5) for k,v in ed.items()} }")
        return sp, edft, ed
    except Exception as e:
        log(f"  FAIL {sp}: {e}")
        return sp, None, str(e)

def run_set(name, din, geomdir):
    setdir = f"{WORK}/{name}"; os.makedirs(setdir, exist_ok=True)
    rxns = parse_din(din)
    species = sorted({s for _, terms in rxns for _, s in terms})
    log(f"=== SET {name}: {len(rxns)} reactions, {len(species)} species ===")
    edft, edisp = {}, {}
    with ThreadPoolExecutor(max_workers=NWORK) as ex:
        futs = {ex.submit(do_species, s, geomdir, setdir): s for s in species}
        for fut in as_completed(futs):
            sp, e, d = fut.result()
            if e is not None:
                edft[sp] = e; edisp[sp] = d
    # build per-reaction interaction energies
    out = [f"# SET {name}  (Eint kcal/mol; ref = high-level)",
           f"# {'rxn':<16} {'ref':>9} {'DFT':>9} " +
           " ".join(f"{rn:>9}" for rn, _ in ROUTES) + "   " +
           " ".join(f"d{rn:>8}" for rn, _ in ROUTES)]
    errs = {rn: [] for rn, _ in ROUTES}; errdft = []; nok = 0
    for ridx, (ref, terms) in enumerate(rxns):
        if any(s not in edft for _, s in terms):
            out.append(f"# rxn{ridx} SKIPPED (missing species)"); continue
        edft_int = sum(c * edft[s] for c, s in terms) * K
        row = {rn: sum(c * (edft[s] + edisp[s][rn]) for c, s in terms) * K
               for rn, _ in ROUTES}
        cplx = terms[0][1]
        out.append(f"  {cplx:<16} {ref:9.2f} {edft_int:9.2f} " +
                   " ".join(f"{row[rn]:9.2f}" for rn, _ in ROUTES) + "   " +
                   " ".join(f"{row[rn]-ref:9.2f}" for rn, _ in ROUTES))
        errdft.append(edft_int - ref)
        for rn, _ in ROUTES: errs[rn].append(row[rn] - ref)
        nok += 1
    def stats(es):
        if not es: return "n/a"
        mae = sum(abs(x) for x in es)/len(es); mse = sum(es)/len(es)
        return f"MAE={mae:6.2f} MSE={mse:+6.2f}"
    out.append(f"# --- {name}: {nok}/{len(rxns)} reactions ---")
    out.append(f"#   DFT-only   {stats(errdft)}")
    for rn, _ in ROUTES:
        out.append(f"#   {rn:<10} {stats(errs[rn])}")
    res = "\n".join(out)
    open(f"{WORK}/results_{name}.txt", 'w').write(res + "\n")
    log(f"=== {name} DONE: {nok}/{len(rxns)} ===\n{res}\n")
    return res

# --- queue: ionic discriminators first, then no-harm ---
DIN  = "/data/refdata/10_din-GMTKN55"
GDIN = "/data/refdata/10_din"
COL  = "/data/refdata/30_collection-GMTKN55"
REF  = "/data/refdata"
QUEUE = [
    ("il16",  f"{DIN}/il16.din",     f"{COL}/il16"),      # ionic (ion pairs)
    ("chb6",  f"{DIN}/chb6.din",     f"{COL}/chb6"),       # ionic (cationic HB)
    ("ahb21", f"{DIN}/ahb21.din",    f"{COL}/ahb21"),      # ionic (anionic HB)
    ("s22",   f"{GDIN}/s22.din",     f"{REF}/20_s22"),     # neutral no-harm (small)
    ("ionichb", f"{GDIN}/ionichb.din", f"{REF}/20_ionichb"), # ionic (big: 120 rxns)
    ("s66",   f"{GDIN}/s66.din",     f"{REF}/20_s66"),     # neutral no-harm (big)
]
sets_to_run = sys.argv[1:] if len(sys.argv) > 1 else [q[0] for q in QUEUE]
for nm, din, gd in QUEUE:
    if nm not in sets_to_run:
        continue
    if not (os.path.exists(din) and os.path.isdir(gd)):
        log(f"=== SKIP {nm}: missing din ({os.path.exists(din)}) or geomdir ({os.path.isdir(gd)}) ===")
        continue
    try:
        run_set(nm, din, gd)
    except Exception:
        log(f"=== SET {nm} CRASHED ===\n{traceback.format_exc()}")
log("ALL SETS COMPLETE")
