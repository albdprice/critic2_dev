#!/usr/bin/env python3
"""Build a provenance SUMMARY index of the Hirshfeld-I anion/ion reference .rho set
by parsing each file's header. Writes markdown."""
import os, re, glob, sys

BASE = sys.argv[1] if len(sys.argv) > 1 else \
    "/tank/research/xdm_chargeaware/data/reference_densities/hirshfeld_proatoms"
OUT = sys.argv[2] if len(sys.argv) > 2 else \
    "/tank/research/xdm_chargeaware/data/reference_densities/SUMMARY.md"

def field(pat, s, cast=str, default="-"):
    m = re.search(pat, s)
    return cast(m.group(1)) if m else default

rows = []
for sub in sorted(os.listdir(BASE)):
    d = os.path.join(BASE, sub)
    if not os.path.isdir(d): continue
    for f in glob.glob(os.path.join(d, "*.rho")):
        h = []
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"): h.append(line.strip())
                else: break
        hdr = " ".join(h)
        el   = field(r"element\s+(\S+)", hdr)
        Z    = field(r"\bZ\s+(\d+)", hdr, int, 0)
        q    = field(r"\bq\s+(-?\d+)", hdr, int, 0)
        nel  = field(r"nelec\s+(\S+)", hdr)
        meth = field(r"method\s+(\S+)", hdr)
        rmax = field(r"rmax\s+([0-9.]+)", hdr)
        stepped = "stepped-down" in hdr
        integ = field(r"integrated_electrons\s+([0-9.]+)", hdr)
        conf = field(r"config\s+'([^']*)'", hdr, str, "-")
        rows.append(dict(route=sub, el=el, Z=Z, q=q, nel=nel, conf=conf,
                         meth=meth, rmax=rmax, stepped=stepped, integ=integ))

rows.sort(key=lambda r: (r["route"], r["Z"], -r["q"]))

lines = []
lines.append("# Hirshfeld-I reference densities — provenance SUMMARY\n")
lines.append(f"Auto-generated index of `reference_densities/hirshfeld_proatoms/**/*.rho` "
             f"({len(rows)} files). Each row = one (element, charge) confined reference; "
             "`box` flags whether the standard 3.6·R99(neutral) box bound the ion or was "
             "stepped down (Route 2), and `∫e⁻` is the integrated electron count before "
             "rescaling to the exact N (a self-consistency check on the confinement).\n")
# per-route counts + stepped-down tally
from collections import Counter
c = Counter(r["route"] for r in rows); sd = Counter(r["route"] for r in rows if r["stepped"])
lines.append("| route | files | stepped-down | charges |")
lines.append("|---|---|---|---|")
for route in sorted(c):
    chs = sorted({r["q"] for r in rows if r["route"]==route})
    lines.append(f"| `{route}` | {c[route]} | {sd.get(route,0)} | {', '.join(map(str,chs))} |")
lines.append("")
lines.append("## Full table")
lines.append("| route | ion | Z | q | N e⁻ | config | method | rmax (a₀) | box | ∫e⁻ (pre-rescale) |")
lines.append("|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    ion = f"{r['el'].capitalize()}{'' if r['q']==0 else ('%+d'%r['q'])}"
    box = "**stepped-down**" if r["stepped"] else "standard"
    lines.append(f"| `{r['route']}` | {ion} | {r['Z']} | {r['q']} | {r['nel']} | "
                 f"`{r['conf']}` | {r['meth']} | {r['rmax']} | {box} | {r['integ']} |")

with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"wrote {OUT}  ({len(rows)} rows; stepped-down: {sum(sd.values())})")
