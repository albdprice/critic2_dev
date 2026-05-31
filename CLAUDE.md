# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

critic2 — a Fortran program for analysis of quantum-mechanical calculations in
molecules and solids (structures, scalar fields, Bader AIM, NCI/ELF, Hirshfeld
& iterative-Hirshfeld partitioning, XDM dispersion). This checkout is the
**`albdprice/critic2_dev` fork**, not upstream.

## Editing/build topology (READ FIRST — non-obvious)

- **`/root/critic2_development` (this cwd) is an editing MIRROR, not a git
  repo** (`git` commands here fail). Edit files here, then `scp` the changed
  file(s) to dev-srv and build/test there.
- **The git repo + build live on dev-srv:** `ssh albd@10.10.49.104`, repo
  `~/critic2_dev`, out-of-source build `~/critic2_dev/build`. Commit and push
  on dev-srv, never locally. Sudo pw on the VMs: `800127828`.
- Typical loop: edit `src/foo@proc.f90` locally →
  `scp src/foo@proc.f90 albd@10.10.49.104:/home/albd/critic2_dev/src/` →
  `ssh albd@... 'cd ~/critic2_dev/build && make -j$(nproc) critic2'`.
- Heredoc/`scp` gotcha: `git commit -m` / heredocs with `(` or `$` over `ssh`
  get mangled — write the commit message to a file and use `git commit -F`.

## Repo conventions (hard rules for this fork)

- **Push only to the fork** (`origin` = `albdprice/critic2_dev`). NEVER push to
  upstream `aoterodelaroza/critic2`. Active work is on branch
  `research/xdm-hirshfeld-i`.
- **Commit attribution: the user only.** Do NOT add "implemented with Claude"
  / Co-Authored-By lines to commits on this fork.

## Build & test

```bash
# build (on dev-srv, out-of-source)
cd ~/critic2_dev/build && make -j$(nproc) critic2     # or: cmake .. first time
# (deps: gfortran, LAPACK/BLAS, libxc, OpenMP; official docs:
#  https://aoterodelaroza.github.io/critic2/installation/)

# run critic2 on an input script
~/critic2_dev/build/src/critic2 input.cri output.cro

# tests (ctest, from the build dir)
cd ~/critic2_dev/build && ctest                 # all
ctest -R hirshfeld                              # subset by regex
ctest -R "022_hirshfeld_i.cro" --rerun-failed --output-on-failure   # single, show diff
```

- Tests live in `tests/NNN_category/<name>.cri` (input) with an expected
  `<name>.cro`. Each `.cri` declares its check inline:
  `## check: <name>.cro -a<abs-tol>` and `## labels: <labels>`. A `-RUN-` test
  generates the output; a `CHECK:` test diffs it against the reference.
- **Large reference data is gated:** `tests/zz_source.tar.xz` (geometries,
  wavefunctions for many tests) is downloaded separately and is often ABSENT
  on the dev clone — so e.g. molecular-XDM tests that need `tests/zz_source/…`
  cannot run there. Don't assume they're runnable.
- After an intentional algorithm change, a `.cro` may legitimately differ
  (e.g. iteration count): re-baseline by copying
  `build/tests/<cat>/<name>.cro` → `tests/<cat>/ref/<name>.cro` only once
  you've confirmed the physically-meaningful numbers are unchanged.

## Code architecture

Fortran with a strict **submodule pattern**: each `src/<mod>.f90` declares the
module (types + `interface` blocks for its `module subroutine`/`function`s),
and `src/<mod>@proc.f90` (or `.F90`) holds the implementations. ~54 modules,
~59 `@proc` submodules. **To change behavior you almost always edit the
`@proc` file; touch the bare `.f90` only to change a signature/type** (and then
update both, or the build errors).

Entry point: `src/critic2.F90` (command dispatcher). The big-picture data flow:

- **`systemmod`** — the central `system` type; the module global `sy` is "the
  current system" used pervasively (`sy%c` = crystal, `sy%f(:)` = loaded
  fields, `sy%iref` = reference field index). Most analysis routines read `sy`.
- **`crystalmod`** (`@complex`, `@env`, `@proc`, …) — crystal/molecule
  structure, neighbor lists, promolecular & core densities
  (`promolecular_atom`, with `zpsp` → core reconstruction for pseudopotential
  densities).
- **`fieldmod` / `grid3mod` / `grid1mod`** — scalar fields. `grid1mod` holds the
  radial atomic densities: `agrid(Z)` (neutral all-electron) and `cgrid(Z,q)`
  (charged/core), loaded from `dat/`.
- **`param.F90`** — physical constants AND tabulated reference data compiled in
  as `parameter` arrays (e.g. `alpha_free` free-atom polarizabilities in Å³
  ÷-converted to bohr³; `frevol*` free volumes; and for the charge-aware XDM
  work, `alpha_gb_{m1,0,p1}`, `pprime_gb`, `rstern_{p1,m1}`). Add new tabulated
  data here.
- **`hirshfeld` / `hirshfeld@proc`** — Hirshfeld and iterative-Hirshfeld
  (Hirshfeld-I) partitioning; the HI SCF (`hirsh_i_driver` grid path,
  `hirsh_i_prepare`/`hirsh_i_refrho`/`hirsh_i_eval` evaluators,
  `hirsh_i_qfloor` element-aware anion clamp).
- **`xdm` / `xdm@proc`** — Becke–Johnson XDM dispersion. Two paths:
  `xdm_wfn` (molecular, Becke/Franchini **mesh** — cusp-safe) and `xdm_grid`
  (periodic, uniform/FFT **grid**). `calc_coefs` builds C6/C8/C10 from volumes
  & exchange-hole moments. Both consume `sy`.

## Active research: charge-aware (Hirshfeld-I) XDM

The current project adds charge-matched references to XDM dispersion. **The
source of truth is `doc/research/xdm_hirshfeld_i_notebook.md` — read its §0
"CURRENT STATE & HANDOFF" first.** Background partitioning doc: `doc/HIRSHFELD_I.md`.

- Keyword surface (gated; default XDM unchanged):
  `xdm <a1> <a2> <chf> hirshfeld_i [volonly] [volref] [alpharef gould|scale|compute|stern] [wfcdir <dir>]`
  (molecular); the periodic form is `xdm grid rho … elf … hirshfeld_i alpharef … wfcdir …`.
- Anion/charged reference densities: `dat/hirshfeld_proatoms/{pbe,pbe0,ld1_pbe}`;
  generators + the Sternheimer/validation tooling in `tools/wfc_generator/`.

## Gotchas (verified)

- Read wavefunctions as a field via **fchk** (`load x.fchk`); the
  wfx-as-structure reader is buggy.
- Inside arithmetic/chemical functions, reference a field by **name**
  (`gkin($mol)`), not `$1`, once several fields are loaded.
- For periodic XDM on a pseudopotential density, set `ZPSP` on the loaded field
  (`load … zpsp Na 9 Cl 7`) so critic2 reconstructs the all-electron density;
  the all-electron density on a uniform grid mis-integrates the core cusp.
- A source-built QE 7.2 (pw.x/pp.x/ld1.x) lives at `/tmp/qe-build` on dev-srv
  (the apt `pw.x` is `_FORTIFY_SOURCE`-broken). Atomic refs use the patched
  `ld1.x`.
