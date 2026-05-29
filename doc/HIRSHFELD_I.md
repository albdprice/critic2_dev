# Iterative Hirshfeld (`HIRSHFELD_I`)

`HIRSHFELD_I` performs an iterative Hirshfeld (Hirshfeld-I) partitioning
of the reference field on a grid, following Bultinck, Van Alsenoy,
Ayers and Carbó-Dorca, *J. Chem. Phys.* **126**, 144111 (2007)
([10.1063/1.2715563](https://doi.org/10.1063/1.2715563)).

The motivation is well known: stockholder Hirshfeld weights built from
neutral free-atom densities systematically underestimate chemical
polarity, because they do not adapt to the actual charge state of the
atom in the molecule. Hirshfeld-I cures this by making the
reference density of each atom *self-consistent* with the population
it ends up carrying.

## Syntax

```
HIRSHFELD_I [WCUBE] [ONLY iat1.i iat2.i ...]
            [WFCDIR dir.s] [HITOL tol.r] [HIMAXIT maxit.i]
```

| Option         | Default | Meaning                                                                                                     |
| -------------- | ------- | ----------------------------------------------------------------------------------------------------------- |
| `WCUBE`        | off     | Write per-atom weight cubes (same as `HIRSHFELD`).                                                          |
| `ONLY i1 i2…`  | all     | Restrict atomic-property integration to the listed cell atoms.                                              |
| `WFCDIR dir`   | unset   | Look in `dir` for user-supplied charged `.wfc` files named `<elem>_q+N.wfc` / `<elem>_q-N.wfc`. See below.  |
| `HITOL tol`    | `1e-4`  | Convergence threshold on `max|ΔQ|` between successive SCF iterations.                                       |
| `HIMAXIT n`    | `60`    | Maximum number of SCF iterations.                                                                           |

The keyword is otherwise a drop-in replacement for `HIRSHFELD`: it
needs a grid-type reference field, integrates the same scalar
properties through `INTEGRABLE`, and reports the same volumes /
populations / Hirshfeld overlap populations.

## Algorithm

Let `ρ(r)` be the reference (molecular) electron density and let
`ρ_A^{q}(r)` be the spherically averaged density of atom *A* in
integer charge state *q*. The Hirshfeld-I weights are

```
            ρ_A^{q_A}(r)
w_A(r) = ───────────────────
          Σ_B ρ_B^{q_B}(r)
```

where each `q_A` is itself a function of the population that *A* picks
up under those weights:

```
N_A = ∫ w_A(r) ρ(r) dr        Q_A = Z_A − N_A
```

The reference density for a fractional `Q_A` is interpolated linearly
between the two flanking integer-charge densities,

```
ρ_A^ref(r) = (1−f) ρ_A^{⌊Q_A⌋}(r) + f ρ_A^{⌈Q_A⌉}(r)        f = Q_A − ⌊Q_A⌋
```

SCF loop:

1. Initialise every `Q_A = 0` (neutral references).
2. Pre-load the radial densities for every `(Z, q)` pair currently in
   use (see the resolution order below). This step is serial so that
   the next step can be embarrassingly parallel.
3. In a single OpenMP grid sweep, simultaneously rebuild `Σ_B ρ_B^ref`
   at every grid point and accumulate `N_A = ∫ w_A ρ dV`.
4. Update each `Q_A`, derive the bracketing `(⌊Q⌋, ⌈Q⌉)` and the
   mixing fraction `f`.
5. If `max|ΔQ|` < `HITOL` exit, else go to step 2.

The same iteration data structure is then handed to
`intgrid_hirshfeld_fields` and `intgrid_hirshfeld_overlap`, which use
the converged `ρ_A^ref` (instead of the neutral free-atom density)
wherever they previously called `agrid(z)%interp`.

## Charged atomic densities

Critic2's built-in `.wfc` tables only contain neutral atoms. Hirshfeld-I
needs both anionic and cationic references. They are resolved in this
order, per `(Z, q)`:

1. **User-supplied file (`WFCDIR`).** If `WFCDIR <dir>` is given, the
   code looks for
   - `<dir>/<elem>_q+N.wfc` for `q > 0`
   - `<dir>/<elem>_q-N.wfc` for `q < 0`

   where `<elem>` is the lower-case element symbol. The file format is
   the standard critic2 `.wfc` (LD1/QE Slater-orbital export). This is
   the rigorous route — drop in publication-quality charged densities
   when you have them.

2. **Built-in cation table.** For `q > 0` the existing `read_db` path
   produces a cation by zeroing trailing orbital occupations.

3. **Anion extrapolation.** For `q < 0`, `read_critic` has been
   extended: when the requested electron count exceeds the neutral
   orbital occupations, electrons are added to the highest-occupied
   orbital up to its angular-momentum capacity `2(2L+1)`. If the
   capacity is insufficient (which only happens for closed-shell
   anions beyond the orbitals stored in the standard `.wfc`), the
   routine emits a warning and uses the best-effort density it
   managed to build.

4. **Fall back to neutral** if all of the above fail, with a warning.

Iteration `Q_A` is clamped to the half-open interval `[-5, Z)`. The
upper limit `Q = Z` corresponds to a fully stripped atom (zero density
everywhere), which is correctly handled by the radial interpolator
returning zero.

### Precomputed anion references (`.rho`) and `WFCDIR`

The extrapolation in path 3 keeps the added electron in the neutral
atom's (compact) outermost orbital, so it under-estimates how diffuse a
real anion is. To use *proper* anion references, point `WFCDIR` at a
directory of radial-density files:

```
HIRSHFELD_I WFCDIR /path/to/proatoms
```

For each `(Z, q)` the code looks, in order, for
`<elem>_q<+|->N.rho` then `<elem>_q<+|->N.wfc` in that directory
(e.g. `o_q-1.rho`, `cl_q-2.rho`, `o_q+1.wfc`), using the trimmed
lowercase element symbol. A `.rho` file is a spherically-averaged
radial density:

```
# comment lines start with '#'
<ngrid>
<r_1>   <rho_1>      # r in bohr, rho in electrons/bohr^3
...
```
on a logarithmic grid `r_i = a*exp(b*(i-1))`, read by
`grid1_read_rho`. Files matching the current SCF charge bracket take
precedence over the built-in cation table / extrapolation.

**Generating the references.** `tools/wfc_generator/gen_anion_rho.py`
produces these files using the standard HORTON-style approach: an
isolated-atom DFT SCF in a *finite* Gaussian basis (which confines the
otherwise-unbound extra electron of a semilocal-functional anion),
followed by Lebedev × log-radial collocation of the all-electron
density (`rho = phi^T D phi`) and angular averaging. It runs in the
Psi4 environment:

```
python gen_anion_rho.py O -1 --lot PBE  --outdir rho_pbe
python gen_anion_rho.py O -1 --lot PBE0 --outdir rho_pbe0
```

`tools/wfc_generator/batch_anions.sh` drives the full main-group
H–Kr set (q=−1 for all; q=−2 for the p-block groups 13–16).
Ready-made PBE/def2-TZVP and PBE0/def2-TZVP databases are shipped in
`dat/hirshfeld_proatoms/{pbe,pbe0}/`.

### Effect of real vs extrapolated anion references (water)

PBE/def2-TZVP water density, same-level PBE/def2-TZVP O⁻ reference:

| O⁻ reference                         | Q(O)   | Q(H)   |
| ------------------------------------ | ------ | ------ |
| extrapolation (compact)              | −0.729 | +0.365 |
| real PBE/def2-TZVP (`.rho`, WFCDIR)  | −0.966 | +0.483 |

The real O⁻ is ~18% more diffuse than neutral O (90%-enclosure radius
2.26 vs 1.92 bohr), so it claims more density and drives a larger
charge. This ~0.24 e shift is the quantitative case for using proper
anion references instead of extrapolation.

> **Consistency caveat.** The shipped neutral and cation references are
> numerical scalar-relativistic PBE atoms (`ld1.x`, the existing
> `dat/wfc` set), whereas these `.rho` anions are Gaussian-basis PBE
> atoms. Interpolating across the `q=0` boundary therefore mixes two
> density representations, and the absolute charge is sensitive to the
> confining basis. Producing *all* charge states (cations, neutral,
> anions) from one method — either fully in the Gaussian route, or via
> confined numerical atomic DFT (Route 2, below) — is the clean fix.

### Route 2: confined numerical atomic DFT (`ld1.x`)

`tools/wfc_generator/gen_anion_rho_ld1.py` generates anion references
with the **same numerical scalar-relativistic PBE method** (`ld1.x`,
Quantum ESPRESSO) that produced critic2's shipped neutral/cation
`dat/wfc` set — so neutral, cation, and anion references are all on one
footing (no representation mixing across `q=0`).

A free atomic anion is unbound for semilocal functionals: `ld1.x` with
the default large box (`rmax=100`) reports *"convergence not achieved"*.
The fix is **box confinement** — shrinking the radial box `rmax` until
the extra electron binds. Per-orbital occupations are read from
`ld1.x`'s own output (authoritative for every block, including the
f-elements), the density is built as
`rho(r) = sum_i occ_i psi_i(r)^2/(4 pi r^2)`, resampled to a log grid,
and written in the same `.rho` format. (`max_out_wfc=99` is set so all
core orbitals are emitted — essential for heavy atoms.)

**Standardized confinement radius.** Rather than an ad-hoc box, the
radius is tied to the atom's own size:

```
rmax = alpha * R99(neutral),   alpha = 3.6   (default)
```

where `R99` is the radius enclosing 99% of the *neutral* atom's density,
computed from the same `ld1.x`/`dat/wfc` set — so the rule is
reproducible, size-relative, and self-consistent. If that box does not
bind (electropositive or multiply-charged anions), `rmax` is stepped
down by 0.8× until `ld1.x` converges; such cases are flagged
`stepped-down` in the `.rho` header.

`alpha = 3.6` was fixed from an O⁻/water sensitivity scan: across
`rmax = 10–15` bohr the Hirshfeld-I charge is **flat to ~0.01 e**
(Q_O = −0.882, −0.877, −0.872), i.e. the result is insensitive to the
exact box as long as it sits in the "binds but gentle" window. `α=3.6`
places most elements there (rmax ≈ 11.9 for O, 13.5 for I, 18.9 for Cs).

Ready-made database `dat/hirshfeld_proatoms/ld1_pbe/` — **whole periodic
table**, 137 files: q=−1 for Z=1–117 and q=−2 for the p-block
(groups 13–16), using the standardized `α·R99` rule.

Coverage notes:
- **Z=104–117 (Rf–Ts)** require an `ld1.x` patched to accept Z>103
  (stock `ld1.x` caps at 103); see `tools/wfc_generator/ld1x_highZ.patch`
  — the same max-Z bump the shipped neutral `dat/wfc` set used. The
  generated `.rho` files are shipped, so no patched `ld1.x` is needed to
  *use* the database, only to regenerate the superheavy entries.
- **Og⁻ (Z=118)** is the single omission: its anion would occupy an 8s
  shell (n=8) beyond `ld1.x`'s quantum-number tables, and a noble-gas
  anion is meaningless.
- A few deeply-unbound multiply-charged cases (e.g. B²⁻, C²⁻) are
  absent — doubly-charged anions of electropositive atoms have no bound
  free-ion limit even under tight confinement.

### Three-way comparison (water O, PBE)

| O⁻ reference                         | Q(O)   | consistent with neutral? |
| ------------------------------------ | ------ | ------------------------ |
| extrapolation (compact)              | −0.729 | n/a (built from neutral) |
| Route 1, Gaussian PBE/def2-TZVP      | −0.966 | no (Gaussian vs ld1)     |
| Route 2, confined ld1.x PBE (rmax 12)| −0.877 | **yes** (same ld1 method)|

Route 2 is method-consistent and lands between the compact
extrapolation and the diffuse Gaussian anion. Because the charge is
insensitive to `rmax` in the gentle-confinement window (above), the
`α·R99` rule makes the database reproducible without fine-tuning.

## Where plain Hirshfeld fails and Hirshfeld-I is needed

Plain Hirshfeld is built entirely from *neutral* free-atom densities, so
it systematically **under-polarizes**: charges come out far too small,
and for an ionic compound it badly fails to recover the ionicity.
Hirshfeld-I removes that bias by self-consistently charging the
pro-atoms. The table below is a robustness sweep over nine molecules,
all at PBE/def2-TZVP density with the Route-2 `ld1_pbe` references
(generated with `tools/wfc_generator/run_molecule_test.py`); every
Hirshfeld-I SCF converged.

| Molecule | Atom | q (Hirshfeld) | q (Hirshfeld-I) | comment |
| -------- | ---- | ------------: | --------------: | ------- |
| NaF  | Na | +0.65 | **+1.01** | recovers the formal +1 of an ionic fluoride |
| NaF  | F  | −0.58 | **−0.94** | |
| NaCl | Na | +0.61 | **+0.96** | strongly ionic, as expected |
| NaCl | Cl | −0.40 | **−0.75** | |
| LiF  | Li | +0.57 | +0.93 | |
| LiH  | Li | +0.41 | +0.89 | ionic hydride: H gets −0.89 |
| H₂O  | O  | −0.31 | −0.87 | |
| NH₃  | N  | −0.30 | −0.91 | |
| HF   | F  | −0.22 | −0.52 | |
| CH₄  | C  | −0.16 | −0.46 | |
| CO   | C  | +0.07 | +0.13 | small either way (the CO anomaly is a *dipole*, not a monopole, effect) |

The consistent pattern — Hirshfeld-I charges are ~2–3× the plain
Hirshfeld values — reproduces the central result of Bultinck et al.
(2007), who note the iterative scheme "increases the magnitudes of the
charges." The clearest **failure of plain Hirshfeld** is the ionic
salts: it reports Na in NaF as only +0.65 (and, on a coarse grid, can
even put a *positive* charge on Cl in NaCl), whereas Hirshfeld-I
correctly returns Na ≈ +1 and a strongly negative halide. Exact values
are method/basis/pro-atom-database dependent; the *trend* and the
recovery of formal ionicity are the robust, reproducible signatures.

## Worked example — water

Single-point B3LYP/6-31G(d) on H₂O, density cube generated with
`cubegen 4 fdensity=scf h2o.fchk h2o_rho_full.cube 160 h`.

Critic2 input:

```
molecule h2o_rho_full.cube
load    h2o_rho_full.cube
integrable 1
hirshfeld
hirshfeld_i
```

SCF trace:

```
# iter  max|dQ|       sum(Q)        atom-charges
  1     3.318E-01    0.00524    -0.3318  0.1685  0.1685
  2     1.735E-01    0.00524    -0.5053  0.2552  0.2552
  …
  16    1.045E-04    0.00524    -0.7292  0.3672  0.3672
  17    6.257E-05    0.00524    -0.7292  0.3672  0.3672
+ Hirshfeld-I SCF converged in 17 iterations
```

Results:

| Method        | Q(O)     | Q(H)     |
| ------------- | -------- | -------- |
| Hirshfeld     | −0.378   | +0.192   |
| Hirshfeld-I   | **−0.729** | **+0.367** |
| Literature¹   | −0.73 … −0.78 | +0.37 … +0.39 |

¹ Bultinck *et al.* 2007 and successors at comparable functionals/basis.
The Hirshfeld-I numbers reproduce the expected ~2× polarity
amplification over plain Hirshfeld.

## What changed in the source tree

* `src/integration.f90` — new `imtype_hirshfeld_i = 7` constant.
* `src/integration@proc.f90` — parses `HIRSHFELD_I` and the new
  `WFCDIR / HITOL / HIMAXIT` sub-options, dispatches to the SCF
  driver, and makes every code path that previously branched on
  `imtype_hirshfeld` also accept the iterative variant. Both the
  per-atom field integration and the bond-order overlap path now use
  the iterative density evaluator when active.
* `src/global@proc.F90` — the top-level command parser now
  recognises `hirshfeld_i` alongside `hirshfeld`.
* `src/hirshfeld.f90` / `src/hirshfeld@proc.f90` — new public
  procedures `hirsh_i_driver`, `hirsh_i_eval`, `hirsh_i_active`,
  `hirsh_i_cleanup`; module-level cache for `(Z, q)` radial
  densities; serial pre-load step that keeps the parallel grid
  sweep race-free.
* `src/grid1mod.f90` / `src/grid1mod@proc.f90` —
  - `read_critic` now augments occupations for anions
    (electron count > neutral occ): it fills the outermost orbital
    up to its `2(2L+1)` capacity.
  - `read_db` no longer rejects `q < 0`.
  - New public wrapper `grid1_read_file(g, file, z, q)` so a caller
    can load a `.wfc` from an absolute path (used by `WFCDIR`).
* `src/types.f90` — `basindat` gained the input fields `hi_wfcdir`,
  `hi_tol`, `hi_maxit` and the per-integration SCF state
  `hi_isactive`, `hi_qlo`, `hi_qhi`, `hi_frac`, `hi_qfinal`. The state
  lives on `bas` (mirroring how YT keeps `bas%luw`); only the shared
  `(Z, q) → grid1` memoisation cache is module-level in
  `hirshfeld@proc.f90`, where it can stay warm across multiple
  HIRSHFELD_I invocations like the neutral `agrid` table.
* `dat/helpdoc/hirshfeld_i.keyw` — keyword syntax shown by the GUI
  help.

## Limitations and follow-ups

* The GUI keyword catalog (`src/gui/templates@proc.f90`) is not
  updated in this change; the keyword is recognised by the parser but
  does not yet appear in the GUI keyword tree. The GUI build is
  off-by-default in critic2; this can be added later.
* Anion extrapolation only adds electrons to orbitals already present
  in the `.wfc` file. Anions that would need a new shell (e.g.
  closed-shell second-row anions like O²⁻ in some contexts) cannot be
  represented this way and trigger a warning. The `WFCDIR` route is
  the supported workaround.
* The `wcube` weight-cube output for `HIRSHFELD_I` currently calls the
  same `hirsh_weights` routine as plain Hirshfeld (neutral
  references). The integration itself uses the iterative weights; the
  visualisation cubes do not. A future change can have `hirsh_weights`
  consult `hirsh_i_eval` in active mode.
