#!/usr/bin/env python3
"""
gen_anion_rho_ld1.py -- Route 2 generator for Hirshfeld-I anion reference
densities: confined numerical all-electron atomic DFT via ld1.x (Quantum
ESPRESSO), written in the ".rho" radial-density format read by critic2's
grid1_read_rho.

Same numerical scalar-relativistic PBE atomic method that produced
critic2's shipped neutral/cation dat/wfc set, so neutral, cation and
anion references are all consistent.

A free anion is unbound for semilocal functionals (ld1.x: "convergence
not achieved" at the default large box). The extra electron is bound by
CONFINEMENT: shrinking the radial box rmax. The box is standardized to

    rmax = alpha * R99(neutral)

where R99 is the radius enclosing 99% of the *neutral* atom's density
(from the same ld1.x dat/wfc set) and alpha is a single global constant.
If that box does not bind, rmax is stepped down by 0.8x until ld1.x
converges (those cases are flagged in the .rho header).

The anion of element Z, charge q (N = Z - q electrons) uses the true
nuclear charge zed=Z and the configuration of the isoelectronic neutral
atom (N electrons). Per-orbital occupations for building the density are
read from ld1.x's own output (authoritative for all blocks).

Usage (ld1.x in PATH):
    gen_anion_rho_ld1.py O -1                 # standardized rmax
    gen_anion_rho_ld1.py O -2 --rmax 5        # explicit rmax override
"""
import os, sys, math, argparse, subprocess, re
import numpy as np

ALPHA_DEFAULT = 3.6   # rmax = ALPHA * R99(neutral); see doc/HIRSHFELD_I.md

SYMS = ["h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p",
 "s","cl","ar","k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga",
 "ge","as","se","br","kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag",
 "cd","in","sn","sb","te","i","xe","cs","ba","la","ce","pr","nd","pm","sm","eu",
 "gd","tb","dy","ho","er","tm","yb","lu","hf","ta","w","re","os","ir","pt","au",
 "hg","tl","pb","bi","po","at","rn","fr","ra","ac","th","pa","u","np","pu","am",
 "cm","bk","cf","es","fm","md","no","lr","rf","db","sg","bh","hs","mt","ds","rg",
 "cn","nh","fl","mc","lv","ts","og"]
ZOF = {s:i+1 for i,s in enumerate(SYMS)}

# neutral ground-state configurations indexed by electron count N (=Z for the
# neutral). Mirrors tools/wfc_generator/gen.m. Anion of charge q uses CONF[N].
_CONF = [
 "1s1","1s2","[He] 2s1","[He] 2s2","[He] 2s2 2p1","[He] 2s2 2p2","[He] 2s2 2p3",
 "[He] 2s2 2p4","[He] 2s2 2p5","[He] 2s2 2p6","[Ne] 3s1","[Ne] 3s2","[Ne] 3s2 3p1",
 "[Ne] 3s2 3p2","[Ne] 3s2 3p3","[Ne] 3s2 3p4","[Ne] 3s2 3p5","[Ne] 3s2 3p6",
 "[Ar] 4s1","[Ar] 4s2","[Ar] 3d1 4s2","[Ar] 3d2 4s2","[Ar] 3d3 4s2","[Ar] 3d5 4s1",
 "[Ar] 3d5 4s2","[Ar] 3d6 4s2","[Ar] 3d7 4s2","[Ar] 3d8 4s2","[Ar] 3d10 4s1",
 "[Ar] 3d10 4s2","[Ar] 3d10 4s2 4p1","[Ar] 3d10 4s2 4p2","[Ar] 3d10 4s2 4p3",
 "[Ar] 3d10 4s2 4p4","[Ar] 3d10 4s2 4p5","[Ar] 3d10 4s2 4p6","[Kr] 5s1","[Kr] 5s2",
 "[Kr] 4d1 5s2","[Kr] 4d2 5s2","[Kr] 4d4 5s1","[Kr] 4d5 5s1","[Kr] 4d5 5s2",
 "[Kr] 4d7 5s1","[Kr] 4d8 5s1","[Kr] 4d10","[Kr] 4d10 5s1","[Kr] 4d10 5s2",
 "[Kr] 4d10 5s2 5p1","[Kr] 4d10 5s2 5p2","[Kr] 4d10 5s2 5p3","[Kr] 4d10 5s2 5p4",
 "[Kr] 4d10 5s2 5p5","[Kr] 4d10 5s2 5p6","[Xe] 6s1","[Xe] 6s2","[Xe] 5d1 6s2",
 "[Xe] 4f1 5d1 6s2","[Xe] 4f3 6s2","[Xe] 4f4 6s2","[Xe] 4f5 6s2","[Xe] 4f6 6s2",
 "[Xe] 4f7 6s2","[Xe] 4f7 5d1 6s2","[Xe] 4f9 6s2","[Xe] 4f10 6s2","[Xe] 4f11 6s2",
 "[Xe] 4f12 6s2","[Xe] 4f13 6s2","[Xe] 4f14 6s2","[Xe] 4f14 5d1 6s2",
 "[Xe] 4f14 5d2 6s2","[Xe] 4f14 5d3 6s2","[Xe] 4f14 5d4 6s2","[Xe] 4f14 5d5 6s2",
 "[Xe] 4f14 5d6 6s2","[Xe] 4f14 5d7 6s2","[Xe] 4f14 5d9 6s1","[Xe] 4f14 5d10 6s1",
 "[Xe] 4f14 5d10 6s2","[Xe] 4f14 5d10 6s2 6p1","[Xe] 4f14 5d10 6s2 6p2",
 "[Xe] 4f14 5d10 6s2 6p3","[Xe] 4f14 5d10 6s2 6p4","[Xe] 4f14 5d10 6s2 6p5",
 "[Xe] 4f14 5d10 6s2 6p6","[Rn] 7s1","[Rn] 7s2","[Rn] 6d1 7s2","[Rn] 6d2 7s2",
 "[Rn] 5f2 6d1 7s2","[Rn] 5f3 6d1 7s2","[Rn] 5f4 6d1 7s2","[Rn] 5f6 7s2",
 "[Rn] 5f7 7s2","[Rn] 5f7 6d1 7s2","[Rn] 5f9 7s2","[Rn] 5f10 7s2","[Rn] 5f11 7s2",
 "[Rn] 5f12 7s2","[Rn] 5f13 7s2","[Rn] 5f14 7s2","[Rn] 5f14 7s2 7p1",
 "[Rn] 5f14 6d2 7s2","[Rn] 5f14 6d3 7s2","[Rn] 5f14 6d4 7s2","[Rn] 5f14 6d5 7s2",
 "[Rn] 5f14 6d6 7s2","[Rn] 5f14 6d7 7s2","[Rn] 5f14 6d8 7s2","[Rn] 5f14 6d9 7s2",
 "[Rn] 5f14 6d10 7s2","[Rn] 5f14 6d10 7s2 7p1","[Rn] 5f14 6d10 7s2 7p2",
 "[Rn] 5f14 6d10 7s2 7p3","[Rn] 5f14 6d10 7s2 7p4","[Rn] 5f14 6d10 7s2 7p5",
 "[Rn] 5f14 6d10 7s2 7p6",
 # N=119 (Og- anion, isoelectronic with hypothetical element 119): add 8s1
 "[Rn] 5f14 6d10 7s2 7p6 8s1"]
CONF_BY_N = {i+1: c for i, c in enumerate(_CONF)}


def r99_neutral(sym, wfcdir):
    """99%-enclosure radius of the neutral atom from its dat/wfc density."""
    f = os.path.join(wfcdir, f"{sym if len(sym)==2 else sym+'_'}_pbe.wfc")
    if not os.path.exists(f):
        f = os.path.join(wfcdir, f"{sym}__pbe.wfc")
    ls = open(f).readlines()
    norb = int(ls[0]); occ = [float(x) for x in ls[2].split()]; ng = int(ls[4])
    d = np.array([[float(x) for x in ls[5+i].split()] for i in range(ng)])
    r = d[:,0]; psi = d[:,1:1+norb]
    rho = (np.array(occ)[None,:]*psi**2).sum(1)/(4*math.pi*r**2)
    P = 4*math.pi*r**2*rho
    cum = np.concatenate([[0], np.cumsum(0.5*(P[1:]+P[:-1])*np.diff(r))]); cum /= cum[-1]
    return float(r[np.searchsorted(cum, 0.99)])


def run_ld1(outdir, sym, Z, conf, rmax, xmin, dx, ld1):
    """Run ld1.x at a given box; return (ok, stdout)."""
    inp = os.path.join(outdir, f"{sym}.ld1in")
    with open(inp,"w") as f:
        f.write("&input\n")
        f.write(f"  title='{sym}',\n  zed={Z}.,\n  rel=1,\n")
        f.write(f"  config='{conf}',\n  iswitch=1,\n  dft='PBE'\n")
        f.write(f"  max_out_wfc=99,\n")   # output ALL orbitals (needed for heavy atoms)
        f.write(f"  xmin={xmin},\n  dx={dx},\n  rmax={rmax},\n/\n")
    wfc = os.path.join(outdir, "ld1.wfc")
    if os.path.exists(wfc): os.remove(wfc)
    res = subprocess.run(f"{ld1} < {os.path.basename(inp)}", shell=True, cwd=outdir,
                         capture_output=True, text=True)
    out = res.stdout + res.stderr
    ok = ("convergence not achieved" not in out.lower()) and os.path.exists(wfc)
    return ok, out


def parse_occ(stdout):
    """Per-orbital occupations from ld1.x stdout lines like '2P 1( 5.00)'."""
    occ = {}
    for m in re.finditer(r"\b(\d[spdfSPDF])\s+\d+\(\s*([0-9.]+)\)", stdout):
        occ[m.group(1).upper()] = float(m.group(2))
    return occ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("symbol"); ap.add_argument("charge", type=int)
    ap.add_argument("--alpha", type=float, default=ALPHA_DEFAULT)
    ap.add_argument("--rmax", type=float, default=0.0, help="explicit rmax (overrides alpha rule)")
    ap.add_argument("--dx", type=float, default=0.005)
    ap.add_argument("--xmin", type=float, default=-7.0)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--wfcdir", default="/home/albd/critic2_dev/dat/wfc")
    ap.add_argument("--ld1", default="ld1.x")
    ap.add_argument("--nout", type=int, default=1500)
    args = ap.parse_args()

    sym = args.symbol.lower()
    if sym not in ZOF: sys.exit(f"unknown element {sym}")
    Z = ZOF[sym]; q = args.charge; N = Z - q
    if N not in CONF_BY_N: sys.exit(f"no configuration for N={N} ({sym} q{q:+d})")
    conf = CONF_BY_N[N]
    os.makedirs(args.outdir, exist_ok=True)
    tag = f"{sym}_q{q:+d}"

    # standardized confinement box, with step-down fallback
    if args.rmax > 0:
        rtry = [args.rmax]
    else:
        r99 = r99_neutral(sym, args.wfcdir)
        rmax0 = round(args.alpha * r99, 1)
        rtry = [round(rmax0*0.8**k, 1) for k in range(8)]   # step down until it binds
    rmax = None; stdout = ""
    for rm in rtry:
        if rm < 1.5: break
        ok, stdout = run_ld1(args.outdir, sym, Z, conf, rm, args.xmin, args.dx, args.ld1)
        if ok: rmax = rm; break
    if rmax is None:
        sys.exit(f"{tag}: ld1.x did not converge at any rmax in {rtry}")
    standardized = (args.rmax <= 0 and abs(rmax - rtry[0]) < 1e-6)

    # orbitals from ld1.wfc, occupations from ld1.x stdout
    lines = open(os.path.join(args.outdir,"ld1.wfc")).readlines()
    orblabels = lines[0].split()[2:]
    data = np.array([[float(x) for x in ln.split()] for ln in lines[1:] if ln.strip()])
    r = data[:,0]; psi = data[:,1:]
    occd = parse_occ(stdout)
    occ = np.array([occd.get(lab.upper(), 0.0) for lab in orblabels])
    rho = (occ[None,:]*psi**2).sum(axis=1) / (4*math.pi*r**2)

    # resample to clean log grid, rescale to N
    rgrid = np.exp(np.linspace(math.log(r[1]), math.log(r[-1]), args.nout))
    rhog = np.clip(np.interp(rgrid, r, rho), 0.0, None)
    integ = np.trapezoid(4*math.pi*rgrid**2*rhog, rgrid)
    if integ > 1e-6: rhog *= N/integ

    outf = os.path.join(args.outdir, f"{tag}.rho")
    flag = "standardized" if standardized else "stepped-down(unbound at standard box)"
    with open(outf,"w") as f:
        f.write("# critic2 radial atomic density (Route 2: confined ld1.x)\n")
        f.write(f"# element {sym} Z {Z} q {q} nelec {N} config '{conf}' "
                f"method PBE/ld1.x rmax {rmax} ({flag})\n")
        f.write(f"# integrated_electrons {integ:.6f} rescaled_to {N}\n")
        f.write(f"{len(rgrid)}\n")
        for a,b in zip(rgrid, rhog):
            f.write(f"{a:.14e} {b:.14e}\n")
    print(f"{outf}: N={N} integ={integ:.4f} rmax={rmax} {flag}")

if __name__ == "__main__":
    main()
