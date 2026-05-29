#!/usr/bin/env python3
"""
gen_anion_rho_ld1.py -- Route 2 generator for Hirshfeld-I anion reference
densities: confined numerical all-electron atomic DFT via ld1.x (Quantum
ESPRESSO), written in the same ".rho" radial-density format read by
critic2's grid1_read_rho.

Unlike Route 1 (Gaussian basis confinement), this uses the SAME numerical
scalar-relativistic atomic method (ld1.x, PBE) that produced critic2's
shipped neutral/cation dat/wfc set -- so neutral, cation, and anion
references are all consistent. The otherwise-unbound extra electron of a
semilocal-functional anion is confined by reducing the radial box (rmax):
a free anion (rmax=100) fails to converge, while a moderate box
(rmax ~ 8-15 bohr) binds it.

The anion of element Z with charge q (N = Z - q electrons) is computed
with the true nuclear charge zed=Z and the electron configuration of the
isoelectronic neutral atom (N electrons).

Usage (needs ld1.x in PATH):
    gen_anion_rho_ld1.py O -1 --rmax 12 --outdir rho_ld1
"""
import os, sys, math, argparse, subprocess
import numpy as np

# isoelectronic neutral ground-state configurations, indexed by electron
# count N (1..38 covered; mirrors tools/wfc_generator/gen.m).
CONF_BY_N = {
    1:"1s1", 2:"1s2", 3:"[He] 2s1", 4:"[He] 2s2", 5:"[He] 2s2 2p1",
    6:"[He] 2s2 2p2", 7:"[He] 2s2 2p3", 8:"[He] 2s2 2p4", 9:"[He] 2s2 2p5",
    10:"[He] 2s2 2p6", 11:"[Ne] 3s1", 12:"[Ne] 3s2", 13:"[Ne] 3s2 3p1",
    14:"[Ne] 3s2 3p2", 15:"[Ne] 3s2 3p3", 16:"[Ne] 3s2 3p4", 17:"[Ne] 3s2 3p5",
    18:"[Ne] 3s2 3p6", 19:"[Ar] 4s1", 20:"[Ar] 4s2", 21:"[Ar] 3d1 4s2",
    22:"[Ar] 3d2 4s2", 23:"[Ar] 3d3 4s2", 24:"[Ar] 3d5 4s1", 25:"[Ar] 3d5 4s2",
    26:"[Ar] 3d6 4s2", 27:"[Ar] 3d7 4s2", 28:"[Ar] 3d8 4s2", 29:"[Ar] 3d10 4s1",
    30:"[Ar] 3d10 4s2", 31:"[Ar] 3d10 4s2 4p1", 32:"[Ar] 3d10 4s2 4p2",
    33:"[Ar] 3d10 4s2 4p3", 34:"[Ar] 3d10 4s2 4p4", 35:"[Ar] 3d10 4s2 4p5",
    36:"[Ar] 3d10 4s2 4p6", 37:"[Kr] 5s1", 38:"[Kr] 5s2",
}
SYMS = ["h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p",
        "s","cl","ar","k","ca","sc","ti","v","cr","mn","fe","co","ni","cu",
        "zn","ga","ge","as","se","br","kr"]
ZOF = {s:i+1 for i,s in enumerate(SYMS)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("symbol"); ap.add_argument("charge", type=int)
    ap.add_argument("--rmax", type=float, default=0.0,
                    help="confinement box radius (bohr); 0 = auto (12 for q=-1, 9 for q<=-2)")
    ap.add_argument("--dx", type=float, default=0.005)
    ap.add_argument("--xmin", type=float, default=-7.0)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--ld1", default="ld1.x")
    ap.add_argument("--nout", type=int, default=1500, help="points in output .rho (resampled log grid)")
    args = ap.parse_args()

    sym = args.symbol.lower()
    if sym not in ZOF: sys.exit(f"unknown element {sym}")
    Z = ZOF[sym]; q = args.charge; N = Z - q
    if N not in CONF_BY_N: sys.exit(f"no configuration for N={N}")
    conf = CONF_BY_N[N]
    os.makedirs(args.outdir, exist_ok=True)
    tag = f"{sym}_q{q:+d}"
    wfc = os.path.join(args.outdir, "ld1.wfc")

    # Confinement box. A free anion (large rmax) is unbound and will not
    # converge; shrink the box until it binds. Start from the requested
    # (or charge-dependent default) rmax and step down.
    if args.rmax > 0:
        rtry = [args.rmax]
    else:
        rtry = [12.0, 10.0, 8.0, 6.0, 5.0, 4.0] if q == -1 else \
               [9.0, 7.0, 5.0, 4.0]
    rmax = None
    for rm in rtry:
        inp = os.path.join(args.outdir, f"{tag}.ld1in")
        with open(inp,"w") as f:
            f.write("&input\n")
            f.write(f"  title='{sym}',\n  zed={Z}.,\n  rel=1,\n")
            f.write(f"  config='{conf}',\n  iswitch=1,\n  dft='PBE'\n")
            f.write(f"  xmin={args.xmin},\n  dx={args.dx},\n  rmax={rm},\n/\n")
        if os.path.exists(wfc): os.remove(wfc)
        res = subprocess.run(f"{args.ld1} < {os.path.basename(inp)}",
                             shell=True, cwd=args.outdir, capture_output=True, text=True)
        out = (res.stdout + res.stderr).lower()
        if "convergence not achieved" not in out and os.path.exists(wfc):
            rmax = rm
            break
    if rmax is None:
        sys.exit(f"{tag}: ld1.x did not converge at any rmax in {rtry}")

    # parse ld1.wfc: header '# r <ORB> <ORB> ...', then rows r, psi_orb...
    lines = open(wfc).readlines()
    labels = lines[0].split()[1:]          # e.g. ['r','2P','2S','1S'] -> drop 'r'
    orblabels = lines[0].split()[2:]
    data = np.array([[float(x) for x in ln.split()] for ln in lines[1:] if ln.strip()])
    r = data[:,0]; psi = data[:,1:]        # (npt, norb), columns match orblabels

    # occupations per orbital from the configuration
    occ = occupations(conf, orblabels)
    # radial density rho(r) = sum_i occ_i * psi_i(r)^2 / (4 pi r^2)
    rho = (occ[None,:]*psi**2).sum(axis=1) / (4*math.pi*r**2)

    # resample onto a clean log grid and rescale to N
    rgrid = np.exp(np.linspace(math.log(r[1]), math.log(r[-1]), args.nout))
    rhog = np.interp(rgrid, r, rho)
    integ = np.trapezoid(4*math.pi*rgrid**2*rhog, rgrid)
    if integ > 1e-6: rhog *= N/integ

    outf = os.path.join(args.outdir, f"{tag}.rho")
    with open(outf,"w") as f:
        f.write("# critic2 radial atomic density (Route 2: confined ld1.x)\n")
        f.write(f"# element {sym} Z {Z} q {q} nelec {N} config '{conf}' "
                f"method PBE/ld1.x rmax {rmax}\n")
        f.write(f"# integrated_electrons {integ:.6f} rescaled_to {N}\n")
        f.write(f"{len(rgrid)}\n")
        for a,b in zip(rgrid, rhog):
            f.write(f"{a:.14e} {b:.14e}\n")
    print(f"{outf}: N_target={N} N_integ={integ:.4f} (err {integ-N:+.4f}) rmax={rmax} conf='{conf}'")

def occupations(conf, orblabels):
    """Occupation per output orbital label (e.g. '2P') from a config string."""
    core = {"[He]":{"1S":2}, "[Ne]":{"1S":2,"2S":2,"2P":6},
            "[Ar]":{"1S":2,"2S":2,"2P":6,"3S":2,"3P":6},
            "[Kr]":{"1S":2,"2S":2,"2P":6,"3S":2,"3P":6,"3D":10,"4S":2,"4P":6}}
    occ = {}
    for tok in conf.split():
        if tok in core:
            occ.update(core[tok])
        else:
            n = tok[0]; l = tok[1].upper(); o = int(tok[2:])
            occ[f"{n}{l}"] = occ.get(f"{n}{l}",0) + o
    return np.array([float(occ.get(lab.upper(),0)) for lab in orblabels])

if __name__ == "__main__":
    main()
