#!/usr/bin/env python3
"""
gen_ion_alpha_watson.py -- Watson-sphere Sternheimer polarizability generator.

Same machinery as gen_ion_alpha_sternheimer.py (confined ld1.x + reconstructed
KS potential + radial uncoupled Sternheimer response), but the ion is stabilized
by a WATSON SPHERE instead of a bare hard box:

  A shell of compensating charge q_comp = -q_net (the ion's net charge) sits at
  radius R_W and contributes the electron potential energy
      V_ws(r) = -|q_net| / max(r, R_W)      (attractive: -|q_net|/R_W inside).
  For an anion (q<0, net charge q<0) this cancels the anion's outward monopole
  so the extra electrons become genuinely bound => the linear response converges
  where the bare free anion diverges. For cations/neutrals q_net>=0 the sphere is
  empty (V_ws=0) and this reduces to the plain confined-Sternheimer route.

  R_W = <r^3>^(1/3) of the NEUTRAL atom (Bučko's TS/HI Watson radius rule,
  PRB 87 064110 / JCP 141 034114), computed from the dat/wfc neutral density.

Implementation reuses the validated helpers from gen_ion_alpha_sternheimer:
  run_ld1, parse_eigen, parse_wfc, reconstruct_V, solve_sternheimer.
The Watson potential is added to the reconstructed V; each orbital's eigenvalue
is shifted to first order by <P|V_ws|P> so (H - eps) stays consistent for the
diffuse anion HOMO.  Output: alpha(0) in a0^3, tabulated as ratios for critic2.

Usage:
  gen_ion_alpha_watson.py O -2                 # R_W from neutral, box = k*R_W
  gen_ion_alpha_watson.py O -2 --boxscan       # print alpha vs box to check stability
"""
import os, sys, math, argparse
import numpy as np
import gen_ion_alpha_sternheimer as G


def rmoment_neutral(sym, wfcdir, p=3):
    """<r^p>^(1/p) of the neutral atom from its dat/wfc density (a0)."""
    sym_l = sym.lower()
    f = os.path.join(wfcdir, f"{sym_l if len(sym_l)==2 else sym_l+'_'}_pbe.wfc")
    if not os.path.exists(f):
        f = os.path.join(wfcdir, f"{sym_l}__pbe.wfc")
    ls = open(f).readlines()
    norb = int(ls[0]); occ = [float(x) for x in ls[2].split()]; ng = int(ls[4])
    d = np.array([[float(x) for x in ls[5+i].split()] for i in range(ng)])
    r = d[:, 0]; psi = d[:, 1:1+norb]
    rho = (np.array(occ)[None, :] * psi**2).sum(1) / (4*math.pi*r**2)
    Pr = 4*math.pi*r**2*rho                       # radial density, integrates to N
    N = np.trapezoid(Pr, r)
    rp = np.trapezoid(r**p * Pr, r) / N
    return rp**(1.0/p)


def alpha_watson(r, P, labels, orbs, Rw, qnet):
    """Sternheimer alpha with a Watson-sphere potential added to the KS V.
    qnet<0 for an anion => attractive shell of |qnet|. Returns (alpha, V, Vws)."""
    # reconstruct KS V (same as the plain generator)
    nl_to_eps = {(n, l): eHa for (n, l, occ, eHa) in orbs}
    nl_to_occ = {(n, l): occ for (n, l, occ, eHa) in orbs}
    lab_nl, orbs_data = {}, []
    for lab in labels:
        n = int(lab[0]); l = G.LCHAR[lab[1].upper()]
        lab_nl[lab] = (n, l)
        if (n, l) in nl_to_eps:
            orbs_data.append((P[lab], l, nl_to_eps[(n, l)]))
    V = G.reconstruct_V(r, orbs_data, None)
    # Watson sphere: attractive shell of |qnet| at R_W (only for net-negative ions)
    Vws = np.where(r <= Rw, -abs(qnet)/Rw, -abs(qnet)/r) if qnet < 0 else np.zeros_like(r)
    Vtot = V + Vws
    alpha = 0.0
    for lab in labels:
        n, l = lab_nl[lab]
        if (n, l) not in nl_to_eps:
            continue
        eps = nl_to_eps[(n, l)]; occ = nl_to_occ[(n, l)]
        if occ <= 0:
            continue
        Pnl = P[lab]
        # first-order eigenvalue shift from the Watson potential
        nrm = np.trapezoid(Pnl*Pnl, r)
        eps_w = eps + np.trapezoid(Pnl*Vws*Pnl, r)/nrm
        for lp, A in [(l+1, (l+1)/(2.0*l+1.0)), (l-1, l/(2.0*l+1.0))]:
            if lp < 0:
                continue
            w = G.solve_sternheimer(r, Pnl, l, lp, Vtot, eps_w)
            I = np.trapezoid(Pnl * r * w, r)
            alpha += (2.0/3.0) * occ * A * (-I)
    return alpha, V, Vws


def compute(sym, q, wfcdir, kbox=2.5, rel=0, outdir=None):
    """alpha(0) [a0^3] for sym^q via Watson-sphere Sternheimer.
    Box rmax = kbox * R_W so the soft Watson tail (not the wall) does the binding."""
    Z = G.ELEM.index(sym); N = Z - q
    if N < 1 or N > len(G._CONF):
        return None, False, f"unsupported N={N}", None
    Rw = rmoment_neutral(sym, wfcdir, p=3)
    conf = G._CONF[N-1]
    rmax = kbox * Rw
    if outdir is None:
        outdir = f"/tmp/stern_ws/{sym}_q{q}"
    os.makedirs(outdir, exist_ok=True)
    ok, out, wfc = G.run_ld1(sym, Z, conf, rmax, outdir, rel=rel)
    if not ok:
        return None, False, "ld1 failed", Rw
    try:
        orbs = G.parse_eigen(out)
        r, P, labels = G.parse_wfc(wfc)
        alpha, V, Vws = alpha_watson(r, P, labels, orbs, Rw, qnet=q)
    except Exception as e:
        return None, False, f"solve error: {e}", Rw
    return alpha, True, conf, Rw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("symbol"); ap.add_argument("charge", type=int)
    ap.add_argument("--wfcdir", default=os.path.expanduser("~/critic2_dev/dat/wfc"))
    ap.add_argument("--kbox", type=float, default=2.5)
    ap.add_argument("--rel", type=int, default=0)
    ap.add_argument("--boxscan", action="store_true",
                    help="print alpha over kbox in {1.5,2,2.5,3,4} to check stability")
    args = ap.parse_args()
    if args.boxscan:
        Rw = rmoment_neutral(args.symbol, args.wfcdir, p=3)
        print(f"{args.symbol} q{args.charge:+d}  R_W(<r^3>^1/3, neutral) = {Rw:.3f} a0")
        for k in (1.5, 2.0, 2.5, 3.0, 4.0):
            a, ok, msg, _ = compute(args.symbol, args.charge, args.wfcdir, kbox=k, rel=args.rel)
            print(f"  kbox={k:>4}  rmax={k*Rw:6.2f}  alpha = "
                  + (f"{a:9.4f} a0^3" if ok else f"FAIL ({msg})"))
        return
    a, ok, msg, Rw = compute(args.symbol, args.charge, args.wfcdir, kbox=args.kbox, rel=args.rel)
    if not ok:
        sys.exit(f"FAILED {args.symbol} q{args.charge}: {msg}")
    print(f"{args.symbol} q{args.charge:+d}  R_W={Rw:.3f}  rmax={args.kbox*Rw:.2f}  "
          f"alpha(0) = {a:.4f} a0^3")


if __name__ == "__main__":
    main()
