# Papers used in the charge-aware (Hirshfeld-I) XDM development

`references.bib` has machine-readable entries. PDFs that are open-access (arXiv)
are in `pdf/`. The rest are paywalled journal articles — click the DOI to pull them
through your institution (they are NOT on arXiv). ✓ = PDF in `pdf/`.

## XDM foundations (Becke–Johnson; Otero-de-la-Roza & Johnson)
- Johnson & Becke, *J. Chem. Phys.* **123**, 024101 (2005) — post-HF exchange-hole
  dipole model of intermolecular interactions (the seed of XDM).
  https://doi.org/10.1063/1.1949201
- Becke & Johnson, *J. Chem. Phys.* **127**, 124108 (2007) — unified DFT treatment
  of dynamical/nondynamical/dispersion correlation (XDM as used in critic2).
  https://doi.org/10.1063/1.2768530
- Otero-de-la-Roza & Johnson, *J. Chem. Phys.* **136**, 174109 (2012) — **XDM in
  solids / cohesive energies** (the periodic-XDM benchmark we target).
  https://doi.org/10.1063/1.4705760
- Otero-de-la-Roza & Johnson, *J. Chem. Phys.* **138**, 204109 (2013) — XDM with
  hybrid / range-separated functionals. https://doi.org/10.1063/1.4807330

## Charge-aware atoms-in-molecules dispersion (the literature we emulate)
- Tkatchenko & Scheffler, *Phys. Rev. Lett.* **102**, 073005 (2009) — TS dispersion
  (neutral free-atom reference; the baseline our work generalizes).
  https://doi.org/10.1103/PhysRevLett.102.073005
- Bučko, Lebègue, Hafner, Ángyán, *J. Chem. Theory Comput.* **9**, 4293 (2013) —
  TS + **iterative Hirshfeld** (charge-aware references in periodic DFT).
  https://doi.org/10.1021/ct400694h
- Bučko, Lebègue, Ángyán, Hafner, *J. Chem. Phys.* **141**, 034114 (2014) —
  extending TS via iterative-Hirshfeld partitioning (the free-ion confinement Q).
  https://doi.org/10.1063/1.4890003
- **Gould, Lebègue, Ángyán, Bučko, *J. Chem. Theory Comput.* **12**, 5920 (2016) —
  fractionally-ionic (FI) MBD — THE method we emulate.** ✓ `pdf/Gould2016_FI_MBD.pdf`
  · arXiv:1703.08786 · https://doi.org/10.1021/acs.jctc.6b00925
- Bučko, Lebègue, Gould, Ángyán, *J. Mol. Model.* (2017) — periodic many-body
  dispersion, reciprocal-space implementation. https://doi.org/10.1007/s00894-017-3514-6

## Ion-reference data (what our charge-aware routes ingest)
- **Gould & Bučko, *J. Chem. Theory Comput.* **12**, 3603 (2016) — C₆ & dipole
  polarizabilities for all atoms + many ions (rows 1–6); the ion-α DB behind the
  "gould" route.** ✓ `pdf/GouldBucko2016_ion_polarizability_DB.pdf`
  · arXiv:1604.02751 · https://doi.org/10.1021/acs.jctc.6b00361
- **Gould, *J. Chem. Phys.* **145**, 084308 (2016) — how α and C₆ vary with atomic
  volume (p′≈p−0.615); the basis of the "scale" route.** ✓ `pdf/Gould2016_volume_scaling.pdf`
  · arXiv:1608.04161 · https://doi.org/10.1063/1.4961643

## Partitioning & benchmarks
- Bultinck, Van Alsenoy, Ayers, Carbó-Dorca, *J. Chem. Phys.* **126**, 144111 (2007)
  — the Hirshfeld-I method. https://doi.org/10.1063/1.2715563
- Kannemann & Becke, *J. Chem. Theory Comput.* **6**, 1081 (2010) — **KB49** vdW
  benchmark (the XDM a1/a2 training set). https://doi.org/10.1021/ct900699r
- Goerigk, Hansen, Bauer, Ehrlich, Najibi, Grimme, *Phys. Chem. Chem. Phys.* **19**,
  32184 (2017) — **GMTKN55** (the IL16/AHB21/CHB6 ionic subsets).
  https://doi.org/10.1039/c7cp04913g
- Manz, *RSC Adv.* **9** (2019) — MCLF charge-aware polarizability/C₆ scaling laws
  (alternative no-ion-reference approach). https://doi.org/10.1039/c9ra03003d

---
**Status:** 3/16 open-access PDFs staged in `pdf/` (the core charge-aware
methodology papers: FI-MBD, the ion-α database, and the volume-scaling law). The
other 13 are paywalled and not on arXiv — use the DOI links above for institutional
access. (Auto arXiv title-search was blocked by API rate-limiting; verified IDs only.)
