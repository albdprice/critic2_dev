#!/usr/bin/env python3
# Psi4 PBE/def2-TZVP on one GMTKN55-style xyz -> fchk + cached energy sidecar.
# xyz format: line1=natoms, line2="charge mult", then atoms (angstrom).
# usage: python psi4_species.py <xyz> <out_fchk> [mem_gb] [nthreads]
import sys, os, psi4
xyz, outfchk = sys.argv[1], sys.argv[2]
mem = sys.argv[3] if len(sys.argv) > 3 else '7'
nthr = int(sys.argv[4]) if len(sys.argv) > 4 else 4
# DF scratch -> /data (34 TB) so it never fills root (/ was 100% full -> write errors)
_scratch = os.environ.get('PSI_SCRATCH', '/data/XDM_Psi4/psi_scratch_kb49')
os.makedirs(_scratch, exist_ok=True)
psi4.core.IOManager.shared_object().set_default_path(_scratch)
L = open(xyz).read().splitlines()
nat = int(L[0]); tok = L[1].split(); chg, mult = int(tok[0]), int(tok[1])
atoms = "\n".join(L[2:2+nat])
psi4.set_memory(f'{mem} GB'); psi4.core.set_num_threads(nthr)
psi4.core.set_output_file(outfchk.replace('.fchk', '.psi4.out'), False)
mol = psi4.geometry(f"{chg} {mult}\n{atoms}\nunits angstrom\nno_reorient\nno_com\nsymmetry c1")
# match AP's KB49 worker convention exactly (Gaussian ultrafine grid + sad guess);
# two-stage convergence: clean defaults first, SOSCF retry only on failure
psi4.set_options({'basis': 'def2-tzvp', 'scf_type': 'df', 'guess': 'sad',
                  'reference': 'rks' if mult == 1 else 'uks', 'freeze_core': 'true',
                  'dft_spherical_points': 590, 'dft_radial_points': 99,
                  'e_convergence': 1e-6, 'd_convergence': 1e-6, 'maxiter': 80})
try:
    e, wfn = psi4.energy('pbe', return_wfn=True)
except Exception:
    psi4.set_options({'soscf': True})
    e, wfn = psi4.energy('pbe', return_wfn=True)
psi4.fchk(wfn, outfchk)
with open(outfchk + '.edft', 'w') as f:
    f.write(f"{e:.10f}\n")
print(f"EDFT {xyz} {e:.10f} Ha  (chg={chg} mult={mult} nat={nat})")
