#!/usr/bin/env python3
"""
gen_anion_rho_zeff.py -- Heidar-Zadeh effective-nuclear-charge generator for
DEEP anion (q <= -3) Hirshfeld-I reference densities, where the confined true-Z
ld1.x SCF simply diverges (a free triple/quadruple anion is too unbound to
converge at ANY box; verified for N3-,P3-,As3-,C4-,Si4- down to rmax=1.0 a0).

Method (Heidar-Zadeh, Ayers & Bultinck, J. Mol. Model. 23:348, 2017):
  1. find the smallest EFFECTIVE nuclear charge Z_eff >= Z whose SCF converges AND
     whose HOMO is genuinely bound (e_HOMO <= 0, the zero-EA / marginally-bound
     point) in the confinement box -- NOT merely the first that converges, which
     can be an unbound "box state" (HOMO>0, electron pinned at the wall);
  2. take the bound density rho(r; Z_eff);
  3. coordinate-rescale to restore the TRUE-Z nuclear cusp and spatial scale
     [their Eq. 10]:  rho(r; Z) = (Z/Z_eff)^3 rho((Z/Z_eff) r; Z_eff)
     (exactly preserves the integral to N; gives the true cusp -2Z; and, since
     Z<Z_eff, EXPANDS the density -- the true anion is more diffuse than the
     more-bound Z_eff auxiliary, which is the physically-correct trend);
  4. renormalize to N, write the .rho (same format as gen_anion_rho_ld1.py).

Same ld1.x/PBE numerics as the q=-1,-2 refs, so the whole -1..-4 set is
consistent. Use this for q=-3 (N,P,As,Sb) and q=-4 (C,Si,Ge).

Usage: gen_anion_rho_zeff.py N -3 --outdir <ld1_pbe>
"""
import os, sys, math, argparse, re
import numpy as np
import gen_anion_rho_ld1 as g


def homo_energy(stdout):
    """Highest-energy OCCUPIED orbital eigenvalue e(Ry) from ld1.x stdout (the HOMO).
    Lines look like '  2 1   2P 1( 6.00)   0.1088   0.0544   1.4807' (n l NL occ eRy eHa eeV).
    Returns +inf if none parsed (treat as unbound)."""
    es = []
    for m in re.finditer(r"\d+\s+\d+\s+\d[SPDFspdf]\s+\d+\(\s*([0-9.]+)\)\s+(-?[0-9.]+)", stdout):
        if float(m.group(1)) > 0:      # occupied
            es.append(float(m.group(2)))
    return max(es) if es else float("inf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("symbol"); ap.add_argument("charge", type=int)
    ap.add_argument("--alpha", type=float, default=3.0)  # deep anions: slightly tighter box
    # than the -1/-2 refs (3.6) keeps Z_eff low -> a more diffuse (physical) auxiliary.
    ap.add_argument("--dx", type=float, default=0.005)
    ap.add_argument("--xmin", type=float, default=-7.0)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--wfcdir", default="/home/albd/critic2_dev/dat/wfc")
    ap.add_argument("--ld1", default="ld1.x")
    ap.add_argument("--nout", type=int, default=1500)
    ap.add_argument("--zstep", type=float, default=0.25)
    ap.add_argument("--zmax_off", type=float, default=5.0, help="max Z_eff-Z to scan")
    ap.add_argument("--homo_tol", type=float, default=0.02, help="HOMO<=tol Ry = bound (zero-EA)")
    args = ap.parse_args()

    sym = args.symbol.lower()
    if sym not in g.ZOF: sys.exit(f"unknown element {sym}")
    Z = g.ZOF[sym]; q = args.charge; N = Z - q
    if N not in g.CONF_BY_N: sys.exit(f"no config for N={N} ({sym} q{q:+d})")
    conf = g.CONF_BY_N[N]
    os.makedirs(args.outdir, exist_ok=True)
    tag = f"{sym}_q{q:+d}"

    r99 = g.r99_neutral(sym, args.wfcdir)
    rmax = round(args.alpha * r99, 1)          # same standard box as q=-1,-2

    # HOMO~0 (zero-EA) criterion [Heidar-Zadeh]: scan Z_eff upward and take the SMALLEST
    # Z_eff whose SCF converges AND whose HOMO is genuinely bound (e_HOMO <= homo_tol Ry).
    # This rejects converged-but-unbound "box states" (HOMO>0, electron pinned at the wall)
    # that a pure convergence test would wrongly accept, and pins the reference to the
    # marginally-bound (most diffuse) frozen-orbital anion -- the correct free-ion proxy.
    zeff = None; stdout = ""; ehomo = None
    for k in range(int(args.zmax_off/args.zstep) + 1):
        ztry = Z + k*args.zstep
        ok, out = g.run_ld1(args.outdir, sym, ztry, conf, rmax, args.xmin, args.dx, args.ld1)
        if ok and homo_energy(out) <= args.homo_tol:
            zeff = ztry; stdout = out; ehomo = homo_energy(out); break
    if zeff is None:
        sys.exit(f"{tag}: no Z_eff in [{Z},{Z+args.zmax_off}] gave a BOUND HOMO at rmax={rmax}")

    # density at Z_eff
    lines = open(os.path.join(args.outdir, "ld1.wfc")).readlines()
    orblabels = lines[0].split()[2:]
    data = np.array([[float(x) for x in ln.split()] for ln in lines[1:] if ln.strip()])
    r = data[:, 0]; psi = data[:, 1:]
    occd = g.parse_occ(stdout)
    occ = np.array([occd.get(lab.upper(), 0.0) for lab in orblabels])
    rho = (occ[None, :]*psi**2).sum(axis=1) / (4*math.pi*r**2)

    # Heidar-Zadeh coordinate rescale: rho_Z(r) = (Z/Zeff)^3 rho_Zeff((Z/Zeff) r)
    s = Z / zeff                               # < 1 -> expands the density
    rgrid = np.exp(np.linspace(math.log(r[1]), math.log(r[-1]/s), args.nout))
    rho_resc = (s**3) * np.clip(np.interp(s*rgrid, r, rho, left=float(rho[0]), right=0.0), 0.0, None)
    integ = np.trapezoid(4*math.pi*rgrid**2*rho_resc, rgrid)
    if integ > 1e-6: rho_resc *= N/integ

    outf = os.path.join(args.outdir, f"{tag}.rho")
    with open(outf, "w") as f:
        f.write("# critic2 radial atomic density (Route 2b: Heidar-Zadeh Z_eff, deep anion)\n")
        f.write(f"# element {sym} Z {Z} q {q} nelec {N} config '{conf}' "
                f"method PBE/ld1.x Z_eff {zeff:.2f} e_HOMO {ehomo:.4f}Ry rmax {rmax} (cusp-rescaled)\n")
        f.write(f"# integrated_electrons {integ:.6f} rescaled_to {N}\n")
        f.write(f"{len(rgrid)}\n")
        for a, b in zip(rgrid, rho_resc):
            f.write(f"{a:.14e} {b:.14e}\n")
    print(f"{outf}: N={N} Z_eff={zeff:.2f} e_HOMO={ehomo:+.4f}Ry integ={integ:.4f} rmax={rmax}")


if __name__ == "__main__":
    main()
