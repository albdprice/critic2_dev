#!/usr/bin/env python3
"""Build an AP-format XDM JSON (base_energy, coords[bohr], c6,c8,c10,rc) from a
critic2 xdm stdout + the species fchk + base_energy. Coords come from the fchk
(native bohr), aligned to critic2's atom order (= input/fchk order)."""
import sys, json, re, numpy as np

def parse_critic2_coeffs(text):
    """Return (n, c6,c8,c10,rc) NxN symmetric from the 'coefficients (a.u.)' block."""
    lines = text.splitlines()
    try:
        start = next(k for k, l in enumerate(lines) if l.strip().startswith("coefficients (a.u.)"))
    except StopIteration:
        return 0, None, None, None, None
    pairs = []
    nmax = 0
    for l in lines[start+2:]:
        s = l.strip()
        if s.startswith("#") or not s:
            break
        t = s.split()
        if len(t) < 6:
            continue
        i, j = int(t[0]), int(t[1])
        c6, c8, c10, rc = float(t[2]), float(t[3]), float(t[4]), float(t[5])
        pairs.append((i, j, c6, c8, c10, rc))
        nmax = max(nmax, i, j)
    n = nmax
    C6 = np.zeros((n, n)); C8 = np.zeros((n, n)); C10 = np.zeros((n, n)); RC = np.zeros((n, n))
    for i, j, c6, c8, c10, rc in pairs:
        a, b = i-1, j-1
        if a == b:
            continue                      # self term, unused by BJ sum (k=1)
        C6[a, b] = C6[b, a] = c6
        C8[a, b] = C8[b, a] = c8
        C10[a, b] = C10[b, a] = c10
        RC[a, b] = RC[b, a] = rc
    return n, C6, C8, C10, RC

def fchk_coords_bohr(fchk):
    """Read 'Current cartesian coordinates' (bohr) from a Gaussian fchk."""
    vals = []; grab = False; need = None
    for l in open(fchk):
        if l.startswith("Current cartesian coordinates"):
            need = int(l.split()[-1]); grab = True; continue
        if grab:
            for tok in l.split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    grab = False; break
            if len(vals) >= need:
                grab = False
    a = np.array(vals[:need]).reshape(-1, 3)
    return a

if __name__ == "__main__":
    out_txt, fchk, edft_file, out_json = sys.argv[1:5]
    text = open(out_txt).read()
    n, C6, C8, C10, RC = parse_critic2_coeffs(text)
    coords = fchk_coords_bohr(fchk)
    base = float(open(edft_file).read().strip())
    if n == 0:                            # single atom: no dispersion
        d = {"base_energy": base, "coords": coords.tolist(),
             "c6": [], "c8": [], "c10": [], "rc": []}
    else:
        assert coords.shape[0] == n, f"{out_json}: {coords.shape[0]} coords vs {n} coeff atoms"
        d = {"base_energy": base, "coords": coords.tolist(),
             "c6": C6.tolist(), "c8": C8.tolist(), "c10": C10.tolist(), "rc": RC.tolist()}
    json.dump(d, open(out_json, "w"))
    print(f"wrote {out_json}  n={n}")
