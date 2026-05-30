#!/usr/bin/env python3
"""
gen_ion_alpha_sternheimer.py -- Stage 2C-rigorous generator.

Compute the static dipole polarizability alpha(0) of a *confined* atom/ion
from a Quantum ESPRESSO ld1.x all-electron calculation, by solving the
radial coupled-perturbed Kohn-Sham (Sternheimer) equations in the
independent-particle (uncoupled) approximation.

It is an OFFLINE dataset generator (like gen_anion_rho_ld1.py): for each
(element, charge) it prints alpha(0) in a0^3, to be tabulated and ingested
by critic2's charge-aware XDM (the `alpharef compute` route), analogous to
the embedded Gould-Bucko table but computed from our own confined ld1.x
ions -- so densities, volumes and polarizabilities all come from the same
box-confined atomic method.

Method (length gauge, uncoupled response):
  For each occupied orbital (n,l) [reduced radial fn P_nl(r) = r R_nl(r),
  eigenvalue eps_nl] and each dipole-allowed l' = l +/- 1, solve
      [ -1/2 d2/dr2 + l'(l'+1)/(2 r^2) + V(r) - eps_nl ] w(r) = - r P_nl(r)
  with Dirichlet boundary conditions w=0 at both ends (the box at rmax is
  exactly what makes the anion polarizability finite -- the same
  regularization as the density). Then
      alpha = (2/3) * Sum_nl occ_nl * Sum_l' A(l,l') * Int P_nl(r) r w(r) dr
      A(l, l+1) = (l+1)/(2l+1),  A(l, l-1) = l/(2l+1).
  The overall constant is fixed/checked against hydrogen (alpha = 4.5).

The KS potential V(r) is reconstructed by inverting the radial KS equation
from the ld1.x orbitals + eigenvalues (multi-orbital, max-amplitude
selection), so it is exactly the ld1 self-consistent potential -- no XC
functional re-evaluation needed.
"""
import os, sys, math, argparse, subprocess, re
import numpy as np
from scipy.linalg import solve_banded

LDX = os.environ.get("LD1X", "ld1.x")
LCHAR = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}

# ground-state configuration by electron number N (N=Z-q), aufbau; mirrors
# tools/wfc_generator/gen_anion_rho_ld1.py
_CONF = [
 "1s1", "1s2", "[He] 2s1", "[He] 2s2", "[He] 2s2 2p1", "[He] 2s2 2p2",
 "[He] 2s2 2p3", "[He] 2s2 2p4", "[He] 2s2 2p5", "[He] 2s2 2p6",
 "[Ne] 3s1", "[Ne] 3s2", "[Ne] 3s2 3p1", "[Ne] 3s2 3p2", "[Ne] 3s2 3p3",
 "[Ne] 3s2 3p4", "[Ne] 3s2 3p5", "[Ne] 3s2 3p6",
]

ELEM = ["", "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",
        "P","S","Cl","Ar","K","Ca"]


def run_ld1(sym, Z, conf, rmax, outdir, xmin=-7.0, dx=0.008, rel=0):
    inp = os.path.join(outdir, "ld1.in")
    with open(inp, "w") as f:
        f.write("&input\n")
        f.write(f"  title='{sym}',\n  zed={Z}.,\n  rel={rel},\n")
        f.write(f"  config='{conf}',\n  iswitch=1,\n  dft='PBE'\n")
        f.write(f"  max_out_wfc=99,\n")
        f.write(f"  xmin={xmin},\n  dx={dx},\n  rmax={rmax},\n/\n")
    wfc = os.path.join(outdir, "ld1.wfc")
    if os.path.exists(wfc): os.remove(wfc)
    res = subprocess.run(f"{LDX} < {os.path.basename(inp)}", shell=True, cwd=outdir,
                         capture_output=True, text=True)
    out = res.stdout + res.stderr
    # ld1 success = the wavefunction file was (re)written and the SCF reached
    # its error threshold ("final scf error ... reached"). The bare word
    # "convergence" appears in benign contexts, so don't key off it.
    ok = os.path.exists(wfc) and ("reached in" in out.lower())
    return ok, out, wfc


def parse_eigen(stdout):
    """Return list of (n, l, occ, eps_Ha) from the ld1.x eigenvalue table."""
    orbs = []
    for m in re.finditer(
        r"^\s*(\d+)\s+(\d+)\s+\d?[SPDFG]\s+\d+\(\s*([0-9.]+)\)\s+"
        r"([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)", stdout, re.M):
        n = int(m.group(1)); l = int(m.group(2)); occ = float(m.group(3))
        eHa = float(m.group(5))
        orbs.append((n, l, occ, eHa))
    return orbs


def parse_wfc(wfc):
    """Return r[ngrid] and dict label->P_nl(r) from ld1.wfc."""
    with open(wfc) as f:
        header = f.readline()
    labels = header.replace("#", "").split()[1:]   # skip the 'r' column label
    data = np.loadtxt(wfc, comments="#")
    r = data[:, 0]
    P = {labels[i]: data[:, i + 1] for i in range(len(labels))}
    return r, P, labels


def deriv2(f, r):
    """Second derivative on a (non-uniform) grid, 3-point formula."""
    d2 = np.zeros_like(f)
    hm = r[1:-1] - r[:-2]
    hp = r[2:] - r[1:-1]
    d2[1:-1] = 2.0 * (f[2:] / (hp * (hp + hm)) - f[1:-1] / (hp * hm)
                      + f[:-2] / (hm * (hp + hm)))
    d2[0] = d2[1]; d2[-1] = d2[-2]
    return d2


def reconstruct_V(r, orbs_data, eps_l):
    """Reconstruct the KS potential V(r) by inverting the radial KS equation
    from each orbital, picking at each point the orbital with the largest
    |P| (farthest from a node) for robustness.
    orbs_data: list of (P_array, l, eps). Returns V(r)."""
    npt = len(r)
    bestamp = np.zeros(npt)
    V = np.full(npt, np.nan)
    for (P, l, eps) in orbs_data:
        d2P = deriv2(P, r)
        with np.errstate(divide="ignore", invalid="ignore"):
            Vi = eps + (0.5 * d2P - l * (l + 1) / (2.0 * r**2) * P) / P
        amp = np.abs(P)
        take = amp > bestamp
        V[take] = Vi[take]
        bestamp[take] = amp[take]
    # clean NaN/inf by interpolation; smooth lightly
    good = np.isfinite(V)
    V = np.interp(r, r[good], V[good])
    # light 3-point smoothing to tame FD noise from d2P
    Vs = V.copy()
    Vs[1:-1] = 0.25 * V[:-2] + 0.5 * V[1:-1] + 0.25 * V[2:]
    return Vs


def solve_sternheimer(r, P, l, lp, V, eps):
    """Solve [-1/2 d2/dr2 + lp(lp+1)/(2r^2) + V - eps] w = -r P, w=0 at ends.
    Non-uniform 3-point FD -> tridiagonal solve. Returns w(r)."""
    n = len(r)
    hm = np.empty(n); hp = np.empty(n)
    hm[1:] = r[1:] - r[:-1]; hp[:-1] = r[1:] - r[:-1]
    hm[0] = hm[1]; hp[-1] = hp[-2]
    lo = np.zeros(n); di = np.zeros(n); up = np.zeros(n)
    cen = lp * (lp + 1) / (2.0 * r**2) + V - eps
    # -1/2 w'' coefficients (3-point non-uniform)
    a = 2.0 / (hm * (hm + hp))      # w_{i-1}
    c = 2.0 / (hp * (hm + hp))      # w_{i+1}
    b = -(a + c)                    # w_i
    lo[1:-1] = -0.5 * a[1:-1]
    di[1:-1] = -0.5 * b[1:-1] + cen[1:-1]
    up[1:-1] = -0.5 * c[1:-1]
    rhs = -r * P
    # Dirichlet BCs
    di[0] = 1.0; up[0] = 0.0; rhs[0] = 0.0
    di[-1] = 1.0; lo[-1] = 0.0; rhs[-1] = 0.0
    ab = np.zeros((3, n))
    ab[0, 1:] = up[:-1]   # super
    ab[1, :] = di         # diag
    ab[2, :-1] = lo[1:]   # sub
    w = solve_banded((1, 1), ab, rhs)
    return w


def alpha_sternheimer(r, P, labels, orbs, prefac=2.0/3.0):
    # map label -> l; build per-orbital arrays aligned to the eigen table
    # eigen table order may differ from wfc column order; match by (l) and
    # radial node count is complex -> instead match by label using n,l.
    lab_by_nl = {}
    # ld1.wfc labels are like '2P','1S'; build l from label
    Pcols = {lab: P[lab] for lab in labels}
    # reconstruct V from all orbitals present
    orbs_data = []
    # need eps per label: match eigen (n,l) to label via n and lchar
    nl_to_eps = {(n, l): eHa for (n, l, occ, eHa) in orbs}
    nl_to_occ = {(n, l): occ for (n, l, occ, eHa) in orbs}
    lab_nl = {}
    for lab in labels:
        n = int(lab[0]); l = LCHAR[lab[1].upper()]
        lab_nl[lab] = (n, l)
        if (n, l) in nl_to_eps:
            orbs_data.append((Pcols[lab], l, nl_to_eps[(n, l)]))
    V = reconstruct_V(r, orbs_data, None)
    alpha = 0.0
    contribs = []
    for lab in labels:
        n, l = lab_nl[lab]
        if (n, l) not in nl_to_eps: continue
        eps = nl_to_eps[(n, l)]; occ = nl_to_occ[(n, l)]
        if occ <= 0: continue
        Pnl = Pcols[lab]
        for lp, A in [(l + 1, (l + 1) / (2.0 * l + 1.0)),
                      (l - 1, l / (2.0 * l + 1.0))]:
            if lp < 0: continue
            w = solve_sternheimer(r, Pnl, l, lp, V, eps)
            I = np.trapezoid(Pnl * r * w, r)
            term = prefac * occ * A * (-I)   # -I so alpha>0 (source has -)
            alpha += term
            contribs.append((lab, lp, term))
    return alpha, contribs, V


def r99_neutral(sym, wfcdir):
    """99%-enclosure radius of the neutral atom from its dat/wfc density."""
    sym = sym.lower()
    f = os.path.join(wfcdir, f"{sym if len(sym)==2 else sym+'_'}_pbe.wfc")
    if not os.path.exists(f):
        f = os.path.join(wfcdir, f"{sym}__pbe.wfc")
    ls = open(f).readlines()
    norb = int(ls[0]); occ = [float(x) for x in ls[2].split()]; ng = int(ls[4])
    d = np.array([[float(x) for x in ls[5+i].split()] for i in range(ng)])
    r = d[:, 0]; psi = d[:, 1:1+norb]
    rho = (np.array(occ)[None, :]*psi**2).sum(1)/(4*math.pi*r**2)
    P = 4*math.pi*r**2*rho
    cum = np.concatenate([[0], np.cumsum(0.5*(P[1:]+P[:-1])*np.diff(r))]); cum /= cum[-1]
    return float(r[np.searchsorted(cum, 0.99)])


def compute_alpha(sym, q, rmax, rel=0, outdir=None):
    """Return (alpha_a0^3, ok, msg) for element sym at charge q, confined at
    rmax. ok=False on ld1/parse failure."""
    Z = ELEM.index(sym)
    N = Z - q
    if N < 1 or N > len(_CONF):
        return None, False, f"unsupported N={N}"
    conf = _CONF[N - 1]
    if outdir is None:
        outdir = f"/tmp/stern/{sym}_q{q}"
    os.makedirs(outdir, exist_ok=True)
    ok, out, wfc = run_ld1(sym, Z, conf, rmax, outdir, rel=rel)
    if not ok:
        return None, False, "ld1 failed"
    try:
        orbs = parse_eigen(out)
        r, P, labels = parse_wfc(wfc)
        alpha, contribs, V = alpha_sternheimer(r, P, labels, orbs)
    except Exception as e:
        return None, False, f"solve error: {e}"
    return alpha, True, conf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("symbol"); ap.add_argument("charge", type=int)
    ap.add_argument("--rmax", type=float, default=None)
    ap.add_argument("--wfcdir", default=os.path.expanduser("~/critic2_dev/dat/wfc"))
    ap.add_argument("--auto-rmax", action="store_true",
                    help="rmax = 3.6 * R99(neutral) (the density-reference box)")
    ap.add_argument("--rel", type=int, default=0)
    args = ap.parse_args()
    sym = args.symbol; q = args.charge
    rmax = args.rmax
    if rmax is None:
        rmax = 3.6 * r99_neutral(sym, args.wfcdir) if args.auto_rmax else 20.0
    alpha, ok, msg = compute_alpha(sym, q, rmax, rel=args.rel)
    if not ok:
        sys.exit(f"FAILED {sym} q{q}: {msg}")
    print(f"{sym} q{q:+d} rmax={rmax:.3f}  alpha(0) = {alpha:.4f} a0^3")


if __name__ == "__main__":
    main()
