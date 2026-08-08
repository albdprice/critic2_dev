# Charge-aware XDM reference generator — anion/ion densities & polarizabilities

Self-contained kit to (re)generate the **Hirshfeld-I ion reference densities** and the
**charge-aware polarizability tables** used by critic2's charge-aware XDM (`hirshfeld_i alpharef …`),
and to reproduce the ionic-solid comparisons. Intended so a collaborator can rebuild everything from
scratch. All pre-generated outputs are already committed in this repo (see below), so you can also
just *use* them without regenerating.

Repo: https://github.com/albdprice/critic2_dev — branch `research/xdm-hirshfeld-i`.

---

## 1. Why this exists (the anion problem)
XDM scales a free-atom reference polarizability by the atom-in-molecule volume. For **ions** the
free *neutral* reference is wrong (a cation is not a compressed neutral atom). We replace it with a
**charge-matched** reference tied to the iterative-Hirshfeld (Hirshfeld-I) charge Q. That needs
reference **densities** for fractional/integer ions — including anions.

**Free anions are electronically unbound at semilocal DFT** (O²⁻/S²⁻ etc. don't hold the extra
electron; a plain SCF diverges / "convergence not achieved"). The reference density is obtained by
**confinement** — the standard stand-in for the crystal Madelung field that binds the anion in a solid.

## 2. The anion/ion reference densities (`.rho`)
Two independent routes, both writing critic2's radial `.rho` format:

**Route 2 — confined numerical ld1.x (primary; consistent with critic2's shipped neutral/cation set)**
```
gen_anion_rho_ld1.py <Sym> <q>            # e.g.  gen_anion_rho_ld1.py O -2
gen_anion_rho_ld1.py O -2 --rmax 5        # override the confinement box
./batch_anions_ld1.sh <outdir> [alpha]    # whole periodic table, q=-1 (+ p-block q=-2)
```
- All-electron scalar-relativistic PBE via **QE ld1.x** (must be in `PATH`).
- Box standardized to `rmax = 3.6 · R99(neutral)` (R99 from the neutral `dat/wfc/*.wfc`); if that
  doesn't bind, rmax steps down ×0.8 until ld1.x converges — flagged in the `.rho` header and in the
  batch `SUMMARY.txt`.
- Ion (Z, q): N = Z−q electrons, nuclear charge zed=Z, config of the isoelectronic neutral.

**Route 1 — Gaussian-basis confinement (independent cross-check)**
```
source <psi4 env>                          # Psi4 required
gen_anion_rho.py O -1 --lot PBE  --outdir rho_pbe
gen_anion_rho.py O -1 --lot PBE0 --outdir rho_pbe0
```
- Psi4 PBE/PBE0 def2-TZVP; the finite Gaussian basis is itself confining. Density evaluated on a
  Lebedev×log grid, spherically averaged → `.rho`.

**Output format** (`.rho`): 3 comment/header lines (element, Z, q, nelec, config, method, rmax, box
status, integrated electrons) + npts + `r  rho(r)` rows; normalized to exactly N electrons.

## 3. Pre-generated set (already in this repo)
```
dat/hirshfeld_proatoms/ld1_pbe/*.rho     # Route 2: 142 files, charges -1..-4
dat/hirshfeld_proatoms/pbe/*.rho         # Route 1 PBE:  38 files, -1/-2
dat/hirshfeld_proatoms/pbe0/*.rho        # Route 1 PBE0: 38 files, -1/-2
dat/wfc/*.wfc                            # 118 neutral references (source of R99)
```
Provenance index (which ions bound at the standard box vs stepped-down): regenerate with
`make_rho_summary.py <proatoms_dir> <out.md>` (a copy of the index ships in the tank archive as
`reference_densities/SUMMARY.md`).

## 4. Using them in critic2 (running the comparisons)
Periodic (grid) charge-aware XDM, one option at a time:
```
xdm grid rho <rho.cube> elf <elf.cube> xa1 <a1> xa2 <a2> \
    hirshfeld_i alpharef <route> wfcdir dat/hirshfeld_proatoms/ld1_pbe
```
`<route>` ∈ `gould | scale | compute | stern | sternws` (omit `hirshfeld_i alpharef …` for the
`neutral` baseline). Molecular path: `xdm … hirshfeld_i alpharef <route> wfcdir <dir>` after
`load x.fchk`. The a1/a2 (BJ damping) are the same across routes (trained on neutral dimers).

A worked 6-option per-solid comparison (inputs + outputs, incl. `c_sternws.out`) is archived at
`tank:/tank/research/xdm_chargeaware/data/headnode_cmp/{MgO,MgS,CaO}/` — the `c_<route>.cri` files
are directly runnable templates.

## 5. Companion: charge-aware polarizability tables (built into critic2)
The α references consumed by `alpharef stern`/`sternws` are compiled into `src/param.F90`
(`rstern_p1/p2/m1`, `aws_m{1,2,3}`). To regenerate:
```
gen_ion_alpha_sternheimer.py <Sym> <q>          # confined-ld1 Sternheimer alpha (stern)
batch_sternheimer.py 54                          # -> rstern_p1/p2/m1
gen_ion_alpha_watson_scf.py O -2 [--gridscan]    # Watson-sphere self-consistent alpha (sternws)
batch_watson.py 54                               # -> aws_m{1,2,3}, vws_m{1,2,3}
```
(`gen_ion_alpha_watson_scf.py --check` runs the bare-H α=4.5 validation.)

## 6. Dependencies
- **Route 2 + α (stern/sternws):** Quantum ESPRESSO **ld1.x** in `PATH`; Python 3 + numpy + scipy.
- **Route 1:** Psi4 (def2-TZVP).
- Running the solid comparisons: critic2 (this repo) + the `.cube` density/ELF inputs (from QE pp.x).

## 7. More detail
- `doc/research/ION_REFERENCE_GENERATION.md` — cation(→+4)/anion(−2/−3) coverage & recipe.
- `doc/HIRSHFELD_I.md` — the confinement convention + Hirshfeld-I background.
- `doc/research/PAPER_METHODS.md` — the methods write-up (all six α routes).
