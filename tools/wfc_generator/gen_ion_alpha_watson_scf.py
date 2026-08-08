#!/usr/bin/env python3
"""
gen_ion_alpha_watson_scf.py -- self-consistent Watson-sphere polarizability.

A self-contained radial Kohn-Sham DFT atom (LDA, PW92 correlation) on a uniform
radial grid, with an optional WATSON SPHERE inside the SCF, followed by the
uncoupled radial Sternheimer response (reused from gen_ion_alpha_sternheimer).

Unlike the post-hoc approach (gen_ion_alpha_watson.py), the Watson potential
    V_ws(r) = -|q_net| / max(r, R_W)          [q_net<0 for an anion]
is included in EVERY SCF iteration, so the orbitals AND eigenvalues are those of
the genuinely bound, sphere-stabilized ion -- which is what makes the deep-anion
(-2,-3) linear response converge instead of diverging.  R_W = <r^3>^(1/3) of the
neutral atom (Bučko TS/HI Watson-radius rule).

Validation ladder (run --check):
  bare-H  -> alpha = 4.50 a0^3 (Sternheimer + grid)
  neutrals-> alpha in line with experiment / the ld1 stern generator
  F-,Cl-  -> compare to Bučko free-ion / Tessman in-crystal
  O2-,S2- -> STABLE across grid + physical (Tessman ~1.7/5 A^3)

Output: alpha(0) in a0^3.  Ratios alpha_ws(Z,q)/alpha_ws(Z,0) feed critic2 as
`rstern_ws` (the 6th alpharef option, `alpharef sternws`).
"""
import os, sys, math, argparse
import numpy as np
from scipy.linalg import eigh_tridiagonal
import gen_ion_alpha_sternheimer as G

# ---- LDA XC (Slater exchange + PW92 correlation, spin-unpolarized) ----
_A, _a1 = 0.031091, 0.21370
_b1, _b2, _b3, _b4 = 7.5957, 3.5876, 1.6382, 0.49294

def _eps_c(rs):
    s = np.sqrt(rs)
    Q1 = 2*_A*(_b1*s + _b2*rs + _b3*rs*s + _b4*rs*rs)
    return -2*_A*(1+_a1*rs)*np.log(1+1.0/Q1)

def vxc_lda(n):
    n = np.maximum(n, 1e-12)
    # exchange potential
    vx = -(3.0*n/math.pi)**(1.0/3.0)
    # correlation potential via V_c = eps_c - (rs/3) d eps_c/d rs (finite diff)
    rs = (3.0/(4*math.pi*n))**(1.0/3.0)
    d = 1e-4*rs
    dec = (_eps_c(rs+d) - _eps_c(rs-d))/(2*d)
    vc = _eps_c(rs) - (rs/3.0)*dec
    return vx + vc

# ---- electron configuration parser (n,l)->occ ----
_CORE = {"He":"1s2","Ne":"1s2 2s2 2p6","Ar":"1s2 2s2 2p6 3s2 3p6",
         "Kr":"1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6",
         "Xe":"1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6"}
_LMAP = {"s":0,"p":1,"d":2,"f":3,"g":4}

def parse_config(conf):
    toks = conf.replace("[","").replace("]"," ").split()
    occ = {}
    expanded = []
    for t in toks:
        if t in _CORE: expanded += _CORE[t].split()
        else: expanded.append(t)
    for t in expanded:
        n = int(t[0]); l = _LMAP[t[1]]; o = float(t[2:])
        occ[(n,l)] = occ.get((n,l),0.0) + o
    return occ

# ---- radial KS solve on a uniform grid ----
def solve_ks(r, h, V, lmax, nmax_per_l):
    """Return dict (n,l)->(eps, u[r]) for the lowest orbitals, u normalized
    so int u^2 dr = 1 (reduced radial fn u=r*R)."""
    orb = {}
    for l in range(lmax+1):
        diag = 1.0/h**2 + l*(l+1)/(2.0*r**2) + V
        off  = np.full(len(r)-1, -0.5/h**2)
        k = nmax_per_l[l]
        if k < 1: continue
        w, v = eigh_tridiagonal(diag, off, select='i', select_range=(0, k-1))
        for idx in range(k):
            n = l + 1 + idx
            u = v[:, idx]
            u = u/np.sqrt(np.trapezoid(u*u, r))
            orb[(n,l)] = (w[idx], u)
    return orb

def hartree(r, sigma):
    """V_H(r) from sigma(r)=sum occ u^2 (=4 pi r^2 n).  V_H=(1/r)int_0^r sigma
    + int_r^inf sigma/r'."""
    from scipy.integrate import cumulative_trapezoid as ctz
    Qin = ctz(sigma, r, initial=0.0)                    # int_0^r sigma
    tail_full = ctz(sigma/r, r, initial=0.0)
    tail = tail_full[-1] - tail_full                    # int_r^rmax sigma/r'
    return Qin/r + tail

def scf(Z, occ, Rw, qnet, rmax=25.0, npts=3500, beta=0.3, tol=1e-6, maxit=300, watson=True):
    r = np.linspace(rmax/npts, rmax, npts); h = r[1]-r[0]
    lmax = max(l for (n,l) in occ)
    nmax_per_l = {l: (max((n for (n,ll) in occ if ll==l), default=l) - l) for l in range(lmax+1)}
    Vws = np.where(r <= Rw, -abs(qnet)/Rw, -abs(qnet)/r) if (watson and qnet < 0) else np.zeros_like(r)
    V = -Z/r + Vws                                     # start bare + Watson
    nold = None
    for it in range(maxit):
        orb = solve_ks(r, h, V, lmax, nmax_per_l)
        sigma = np.zeros_like(r)
        for (n,l),o in occ.items():
            if (n,l) in orb: sigma += o*orb[(n,l)][1]**2
        n_e = sigma/(4*math.pi*r**2)
        VH = hartree(r, sigma)
        Vnew = -Z/r + VH + vxc_lda(n_e) + Vws
        if nold is not None:
            dn = np.trapezoid(np.abs(sigma-nold), r)
            if dn < tol: V = (1-beta)*V + beta*Vnew; break
        V = (1-beta)*V + beta*Vnew
        nold = sigma
    return r, V, orb, occ

def solve_inhom(r, h, rhs, lp, V, eps):
    """Solve [-1/2 d2/dr2 + lp(lp+1)/(2r^2) + V - eps] w = rhs on the uniform
    grid (u=0 outside), general tridiagonal."""
    from scipy.linalg import solve_banded
    n = len(r)
    d = 1.0/h**2 + lp*(lp+1)/(2.0*r**2) + V - eps
    e = np.full(n-1, -0.5/h**2)
    ab = np.zeros((3, n)); ab[0,1:] = e; ab[1,:] = d; ab[2,:-1] = e
    return solve_banded((1,1), ab, rhs)

def alpha_from(r, V, orb, occ):
    h = r[1]-r[0]
    occ_by_l = {}                                        # occupied orbitals per l (Pauli block)
    for (n,l),o in occ.items():
        if o > 0 and (n,l) in orb:
            occ_by_l.setdefault(l, []).append(orb[(n,l)][1])
    a = 0.0
    for (n,l),o in occ.items():
        if o <= 0 or (n,l) not in orb: continue
        eps, u = orb[(n,l)]
        for lp, A in [(l+1,(l+1)/(2.0*l+1.0)), (l-1, l/(2.0*l+1.0))]:
            if lp < 0: continue
            rhs = -r*u
            for uj in occ_by_l.get(lp, []):              # project RHS off occupied lp (Pauli)
                rhs = rhs - np.trapezoid(uj*rhs, r)*uj
            w = solve_inhom(r, h, rhs, lp, V, eps)
            for uj in occ_by_l.get(lp, []):              # keep response in the virtual space
                w = w - np.trapezoid(uj*w, r)*uj
            a += (2.0/3.0)*o*A*(-np.trapezoid(u*r*w, r))
    return a

def compute(sym, q, wfcdir, rmax=25.0, npts=3500, watson=True):
    Z = G.ELEM.index(sym); N = Z - q
    if N < 1 or N > len(G._CONF): return None, f"unsupported N={N}", None
    from gen_ion_alpha_watson import rmoment_neutral
    Rw = rmoment_neutral(sym, wfcdir, p=3)
    occ = parse_config(G._CONF[N-1])
    r, V, orb, occ = scf(Z, occ, Rw, qnet=q, rmax=rmax, npts=npts, watson=watson)
    return alpha_from(r, V, orb, occ), G._CONF[N-1], Rw

def bare_hydrogen(npts=6000, rmax=60.0):
    r = np.linspace(rmax/npts, rmax, npts); h = r[1]-r[0]
    V = -1.0/r
    orb = solve_ks(r, h, V, 0, {0:1})
    eps, u = orb[(1,0)]
    return alpha_from(r, V, orb, {(1,0):1.0}), eps

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("symbol", nargs="?"); ap.add_argument("charge", nargs="?", type=int)
    ap.add_argument("--wfcdir", default=os.path.expanduser("~/critic2_dev/dat/wfc"))
    ap.add_argument("--rmax", type=float, default=25.0); ap.add_argument("--npts", type=int, default=3500)
    ap.add_argument("--no-watson", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--gridscan", action="store_true")
    a = ap.parse_args()
    if a.check:
        al, eps = bare_hydrogen()
        print(f"bare-H: eps={eps:.4f} Ha (exact -0.5), alpha={al:.4f} a0^3 (exact 4.5)")
        return
    if a.gridscan:
        for rmax in (18,22,26,30):
            al,_,Rw = compute(a.symbol, a.charge, a.wfcdir, rmax=rmax, npts=int(rmax*140),
                              watson=not a.no_watson)
            print(f"  rmax={rmax}  R_W={Rw:.3f}  alpha={al:.4f} a0^3")
        return
    al, conf, Rw = compute(a.symbol, a.charge, a.wfcdir, a.rmax, a.npts, watson=not a.no_watson)
    print(f"{a.symbol} q{a.charge:+d}  R_W={Rw:.3f}  alpha(0)={al:.4f} a0^3  [{conf}]")

if __name__ == "__main__":
    main()
