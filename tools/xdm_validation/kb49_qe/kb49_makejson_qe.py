#!/usr/bin/env python3
"""Build AP-format XDM JSON from a critic2 GRID-XDM stdout (+ xyz for coords[bohr],
+ base_energy[Ha]). Parses the '+ Dispersion coefficients' block (i j C6 C8 C10 Rc).
usage: kb49_makejson_qe.py <critic2_out> <xyz> <base_energy_Ha> <out.json>"""
import sys, json, numpy as np
BOHR = 0.52917720859

def parse_grid_coeffs(text):
    lines = text.splitlines()
    try:
        s = next(k for k, l in enumerate(lines) if "Dispersion coefficients" in l)
    except StopIteration:
        return 0, None, None, None, None
    # find the column header '# i   j ... C6 ...'
    h = next(k for k in range(s, len(lines)) if lines[k].lstrip().startswith("# i") and "C6" in lines[k])
    pairs = []; nmax = 0
    for l in lines[h+1:]:
        st = l.strip()
        if not st or st.startswith("#") or st.startswith("+"):
            break
        t = st.split()
        if len(t) < 6:
            break
        try:
            i, j = int(t[0]), int(t[1])
        except ValueError:
            break
        c6, c8, c10, rc = (float(t[2]), float(t[3]), float(t[4]), float(t[5]))
        pairs.append((i, j, c6, c8, c10, rc)); nmax = max(nmax, i, j)
    n = nmax
    C6 = np.zeros((n, n)); C8 = np.zeros((n, n)); C10 = np.zeros((n, n)); RC = np.zeros((n, n))
    for i, j, c6, c8, c10, rc in pairs:
        a, b = i-1, j-1
        if a == b:
            continue                      # self term unused by BJ sum (k=1)
        C6[a,b]=C6[b,a]=c6; C8[a,b]=C8[b,a]=c8; C10[a,b]=C10[b,a]=c10; RC[a,b]=RC[b,a]=rc
    return n, C6, C8, C10, RC

def xyz_coords_bohr(xyz):
    L = open(xyz).read().split("\n"); n = int(L[0]); c = []
    for ln in L[2:2+n]:
        t = ln.split()
        if len(t) >= 4:
            c.append([float(t[1])/BOHR, float(t[2])/BOHR, float(t[3])/BOHR])
    return np.array(c)

if __name__ == "__main__":
    out_txt, xyz, base, out_json = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
    n, C6, C8, C10, RC = parse_grid_coeffs(open(out_txt).read())
    coords = xyz_coords_bohr(xyz)
    if n == 0:                            # single atom
        d = {"base_energy": base, "coords": coords.tolist(), "c6": [], "c8": [], "c10": [], "rc": []}
    else:
        assert coords.shape[0] == n, f"{out_json}: {coords.shape[0]} coords vs {n} coeff atoms"
        d = {"base_energy": base, "coords": coords.tolist(),
             "c6": C6.tolist(), "c8": C8.tolist(), "c10": C10.tolist(), "rc": RC.tolist()}
    json.dump(d, open(out_json, "w"))
    print(f"wrote {out_json} n={n}")
