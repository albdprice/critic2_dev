#!/usr/bin/env python3
"""
gen_anion_rho.py -- generate spherically-averaged atomic (an)ion densities
for use as Hirshfeld-I reference pro-atoms in critic2.

Route 1 (Gaussian-basis confinement): runs an isolated-atom DFT SCF in
Psi4 at PBE/def2-TZVP (or PBE0), evaluates the all-electron density on a
Lebedev angular grid x logarithmic radial grid via basis-function
collocation (rho = phi^T D phi), angularly averages to a spherical
rho(r), and writes a critic2 ".rho" radial-density file.

The finite Gaussian basis confines the (otherwise unbound, for semilocal
functionals) extra electron of an anion -- the standard HORTON ProAtomDB
approach.

Output (".rho"), one file per (element, charge), rho normalized so that
4*pi*int r^2 rho dr = N exactly:
    # critic2 radial atomic density
    # element <sym> Z <Z> q <q> nelec <N> mult <m> method <lot>
    # integrated_electrons <before-rescale> rescaled_to <N>
    <ngrid>
    <r_1>   <rho_1>
    ...
r is on a log grid r_i = a*exp(b*(i-1)) (bohr); rho in electrons/bohr^3.

Run inside the psi4 environment:
    source ~/projects/psi4_xdm_implement/activate_xdm.sh
    python gen_anion_rho.py O -1 --lot PBE --outdir rho_pbe
"""
import os, sys, math, argparse
import numpy as np
import psi4

SYMS = ["h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p",
        "s","cl","ar","k","ca","sc","ti","v","cr","mn","fe","co","ni","cu",
        "zn","ga","ge","as","se","br","kr"]
ZOF = {s:i+1 for i,s in enumerate(SYMS)}

# ground-state spin multiplicity vs electron count N (isoelectronic neutral)
MULT_BY_N = {
    1:2, 2:1, 3:2, 4:1, 5:2, 6:3, 7:4, 8:3, 9:2, 10:1,
    11:2, 12:1, 13:2, 14:3, 15:4, 16:3, 17:2, 18:1, 19:2, 20:1,
    21:2, 22:3, 23:4, 24:7, 25:6, 26:5, 27:4, 28:3, 29:2, 30:1,
    31:2, 32:3, 33:4, 34:3, 35:2, 36:1, 37:2, 38:1,
}

def lebedev26():
    pts, wts = [], []
    a1, a2, a3 = 1.0/21.0, 4.0/105.0, 9.0/280.0
    for s in (+1.0,-1.0):
        pts += [(s,0,0),(0,s,0),(0,0,s)]; wts += [a1]*3
    r2 = 1.0/math.sqrt(2.0)
    for sx in (+1,-1):
        for sy in (+1,-1):
            pts += [(sx*r2,sy*r2,0),(sx*r2,0,sy*r2),(0,sx*r2,sy*r2)]; wts += [a2]*3
    r3 = 1.0/math.sqrt(3.0)
    for sx in (+1,-1):
        for sy in (+1,-1):
            for sz in (+1,-1):
                pts.append((sx*r3,sy*r3,sz*r3)); wts.append(a3)
    pts = np.array(pts); wts = np.array(wts)
    assert abs(wts.sum()-1.0) < 1e-12
    return pts, wts

def build_density_evaluator(wfn):
    basis = wfn.basisset()
    Dtot = np.array(wfn.Da()) + np.array(wfn.Db())
    nbf = basis.nbf()
    extents = psi4.core.BasisExtents(basis, 1e-10)
    def rho_at(pts):
        pts = np.ascontiguousarray(pts, float)
        n = len(pts)
        xv = psi4.core.Vector.from_array(pts[:,0].copy())
        yv = psi4.core.Vector.from_array(pts[:,1].copy())
        zv = psi4.core.Vector.from_array(pts[:,2].copy())
        wv = psi4.core.Vector.from_array(np.ones(n))
        block = psi4.core.BlockOPoints(xv, yv, zv, wv, extents)
        bf = psi4.core.BasisFunctions(basis, n, nbf)
        bf.compute_functions(block)
        phi = np.array(bf.basis_values()["PHI"])   # (n, nbf)
        return np.einsum("pm,mn,pn->p", phi, Dtot, phi)
    return rho_at

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("symbol"); ap.add_argument("charge", type=int)
    ap.add_argument("--lot", default="PBE", choices=["PBE","PBE0"])
    ap.add_argument("--basis", default="def2-tzvp")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--nrad", type=int, default=1200)
    ap.add_argument("--rmin", type=float, default=1e-3)
    ap.add_argument("--rmax", type=float, default=15.0)
    args = ap.parse_args()

    sym = args.symbol.lower()
    if sym not in ZOF: sys.exit(f"unknown element {sym}")
    Z = ZOF[sym]; q = args.charge; N = Z - q
    if N < 1: sys.exit("non-positive electron count")
    if N not in MULT_BY_N: sys.exit(f"no multiplicity for N={N}")
    mult = MULT_BY_N[N]
    func = {"PBE":"PBE","PBE0":"PBE0"}[args.lot]
    ref = "rks" if mult == 1 else "uks"

    os.makedirs(args.outdir, exist_ok=True)
    psi4.set_memory("3 GB"); psi4.core.be_quiet()
    psi4.core.set_output_file(os.path.join(args.outdir, f"{sym}_q{q:+d}_{args.lot}.psi4out"), False)

    mol = psi4.geometry(f"""
{q} {mult}
{sym.capitalize()} 0.0 0.0 0.0
units bohr
symmetry c1
""")
    psi4.set_options({"basis":args.basis, "reference":ref, "scf_type":"df",
                      "e_convergence":1e-8, "d_convergence":1e-8,
                      "maxiter":300, "dft_spherical_points":590,
                      "dft_radial_points":120})
    try:
        e, wfn = psi4.energy(func, molecule=mol, return_wfn=True)
    except Exception as ex:
        sys.exit(f"SCF failed for {sym} q{q:+d} {args.lot}: {ex}")

    rho_at = build_density_evaluator(wfn)

    nrad = args.nrad
    b = math.log(args.rmax/args.rmin)/(nrad-1)
    rgrid = args.rmin*np.exp(b*np.arange(nrad))
    angpts, angw = lebedev26()
    nang = len(angw)
    allpts = (rgrid[:,None,None]*angpts[None,:,:]).reshape(-1,3)  # (nrad*nang,3)
    rho_pts = rho_at(allpts).reshape(nrad, nang)
    rho = np.clip((rho_pts*angw[None,:]).sum(axis=1), 0.0, None)

    integ = np.trapezoid(4*math.pi*rgrid**2*rho, rgrid)
    rho_scaled = rho * (N/integ) if integ > 1e-6 else rho

    outf = os.path.join(args.outdir, f"{sym}_q{q:+d}.rho")
    with open(outf,"w") as f:
        f.write("# critic2 radial atomic density\n")
        f.write(f"# element {sym} Z {Z} q {q} nelec {N} mult {mult} method {args.lot}/{args.basis}\n")
        f.write(f"# integrated_electrons {integ:.6f} rescaled_to {N}\n")
        f.write(f"{nrad}\n")
        for r,d in zip(rgrid, rho_scaled):
            f.write(f"{r:.14e} {d:.14e}\n")
    print(f"{outf}: E={e:.6f} N_target={N} N_integ={integ:.4f} (err {integ-N:+.4f}) mult={mult} ref={ref}")

if __name__ == "__main__":
    main()
