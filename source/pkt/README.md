# qed algorithm

## Overview

The qed algorithm of OSIRIS includes both first-principles and reduced descriptions of multiple QED effects

* Inverse Compton scattering (Classical and QED)
* Nonlinear Breit-Wheeler pair production in intense EM fields
* Linear Compton scattering
* Bremsstrahlung emission of gamma-rays in atomic Coulomb field
* Bethe-Heitler pair production in atomic Coulomb field
* Trident pair production in atomic Coulomb field

They can be used in 1/2/3D cartesian, as well as in quasi-3D (except Linear Compton scattering that is not yet implemented in quasi-3D).

## Physics

(TBD)

## Input parameters

* The QED algorithm is currently implemented as a simulation mode. To activate it the following parameter must be added to the ``simulation`` section
  
  * `algorithm = "qed",`

* Species interacting through QED processes are part of the same QED group. The number of QED groups is controlled in the ``particles`` section, using the ``num_qed`` option. A QED group consists in 3 species : an electron species, a positron species and a photon species (in this exact order).

```
  particles 
  { 
    num_species = 0, 
    num_qed = 1,
    dlb_fld2_thresh = 2500.,
    dlb_part_weight = 2.0,
  }
```

In this section, the relevant parameters are:
* ``num_qed`` - The number of qed groups present
* ``dlb_fld2_thresh`` - A real number (default is ``1.0e6``) used for dynamic load balancing.  For each particle, the load is determined by first computing ``|E|^2 + |B|^2`` at the particle. If this quantity is greater than ``dlb_fld2_thresh``, the particle is assigned a load weight of ``dlb_part_weight``.  Otherwise, the load weight is given by ``1 + (|E|^2 + |B|^2) * norm_schw2``, where ``norm_schw2`` is the inverse Schwinger field squared.
* ``dlb_part_weight`` - Weighting to give particles with fields present above the Schwinger limit (see ``dlb_fld2_thresh`` above).  Real number, defaults to ``2.0``.

For each QED group, a ``qed_group`` section can be included to control the parameters for photon emission and pair production within this group:
  
```
  qed_group
  {
    pairprod = "qed",
    if_damp_classical = .true.,
    if_damp_qed = .true.,
    if_pairprod = .true.,
    
    zeta_p = 1.0,
    zeta_g = 1.0,
    qed_g_cutoff = 10.0,
    p_emit_cutoff = 2.0,
  }
```


  In this section, the relevant parameters are:
  * ``pairprod`` - Controls the pair production model used in the simulation. The models allowed are ``"qed"`` (__ab initio__, default) and ``"gthr"`` (energy based model)
  * ``if_damp_classical`` and ``if_damp_qed`` - Control whether electrons and positrons lose energy due to classical radiation reaction / QED photon emission
  * ``if_pairprod`` - Controls whether photons emit electron-positron pairs
  * ``qed_g_cutoff`` - Sets the threshold energy between classical radiation reaction and QED photon emission
  * ``p_emit_cutoff`` - Sets a minimum energy for which photons emitted due to QED processes will be tracked
  * ``zeta_p`` and ``zeta_g`` - Sets the factors by which the quantum parameters of leptons and photons should be multiplied, respectively (for fully __ab initio__ calculations, use ``zeta_p = 1`` and ``zeta_g = 1``)
  * ``pairprod_gthr`` and ``pairprod_gpair`` - Control the parameters of the energy based pair production model (``pairprod = "gthr"``); leptons emit pairs with combined, equally split energy ``pairprod_gpair`` whenever they reach an energy ``pairprod_gthr``

Following each ``qed_group`` section, a section ``qed_group_diag`` can be included to control the parameters for diagnostics relative to the QED group:

```
qed_group_diag
{
  ndump_fac_rad = 1,
  ndump_fac_pairs = 1,
  ndump_fac_pairs_from_elec = 1,

  ndump_fac_radspect = 10,
  radspect_emin = 1.957d-5,
  radspect_emax = 1.957d5,
  radspect_bins = 400,

  ndump_fac_chi_emit = 10,
  chi_emit_min = 1.d-4,
  chi_emit_max = 1.d2,
  chi_emit_nbins = 240,

  ndump_fac_raddetector = 1,
  raddetector_nbins = 100,
  raddetector_emin = 0,
  raddetector_emax = 1.d15,
  raddetector_phimin = -0.05,
  raddetector_phimax = 0.05,
  raddetector_thetamin = -0.05, 
  raddetector_thetamax = 0.05,
  detector_direction = "x1_positive",
  
  ndump_fac_radsphere = 1,
  radsphere_nbins = 100, 
  radsphere_emin = 0.0,
  radsphere_emax = 1.d15,
  if_add_classical_radsphere = .false.,

  ndump_fac_n_emit = 1,
  ndump_fac_chi = 1,
}
```

In this section, the relevant parameters are:

  * ``ndump_fac_rad`` - Controls the dump frequency of the global radiated energy diagnostic, which includes independent columns for energy radiated by charged particles in the classical and QED regimes (for tracked and non-tracked photons) and also for Bremsstrahlung emission; this diagnostic produces one file for each electron and positron species in the QED group (see ``HIST/par01_rad`` and ``HIST/par02_rad``)
  * ``ndump_fac_pairs`` - Controls the dump frequency of a global diagnostic that reports the total number, energy and weight of pairs emitted by the photon species in the QED group via nonlinear Breit-Wheeler and Bethe Heitler (see HIST/pho01_ene and HIST/pho01_pairs). The diagnostic is stored in ``HIST/pho01_ene/`` and ``HIST/pho01_pairs/``.
  * ``ndump_fac_pairs_from_elec`` - Controls the dump frequency of a global diagnostic that reports the total number, energy and weight of pairs emitted by the electron/positron species in the QED group by Coulomb Trident (see HIST/par01_pairs). The diagnostic is stored in ``HIST/par01_pairs/``.
  * ``ndump_fac_radspect``, ``radspect_emin``, ``radspect_emax`` and ``radspect_bins`` - Respectively control the dump frequency, minimum and maximum energies, and number of bins of a runtime diagnostic of the energy spectrum of the photons emitted by each electron and positron species; the binning is only done logarithmically (see QED/RADSPECT/). This diagnostic is normalised in J/MeV for a 3D simulation, and in J/MeV/m for a 2D simulation. The min/max bounds and number of bins are optional and have the default values indicated in the example. The diagnostic is stored in ``QED/radspect/``.
  * ``ndump_fac_chi_emit``, ``chi_emit_min``, ``chi_emit_max`` and ``chi_emit_nbins`` - Respectively control the dump frequency, minimum and maximum chi parameters, and number of bins of a runtime diagnostic of the quantum parameter spectrum of both the electrons/positrons emitting and the photons emitted; the binning is only done logarithmically (see QED/chi_emission/). The min/max bounds and number of bins are optional and have the default values indicated in the example. The diagnostic is stored in ``QED/chi_emission/``.
  * ``ndump_fac_raddetector`` - Controls the frequency of a "detector-like" 2D diagnostic for photons emitted via NLIC, stored in ``QED/rad-detector/``. All optional parameters are given default values showed in the example. These optional parameters are the number of bins ``raddetector_nbins``, the min/max photon energies selected ``raddetector_emin``, ``raddetector_emax``, the min/max azimuthal angles selected ``raddetector_phimin``, ``raddetector_phimax``, the min/max polar angles selected ``raddetector_thetamin``, ``raddetector_thetamax``, and the direction of the detector ``detector_direction``.
  * ``ndump_fac_radsphere`` - Controls the frequency of a photon energy diagnostic on a spherical surface, stored in ``QED/radsphere/``. It is done via six 2D diagnostics. The six grids correspond to the six directions +/-x1, +/-x2 and +/-x3. They enable to reconstruct a spherical 4 pi steraian detector . Optional parameter are given default values whichcan be tuned. ``radsphere_nbins`` is the number of bins, ``radsphere_emin``, ``radsphere_emax`` are the min and max energies of photons selected and ``if_add_classical_radsphere`` is a boolean to add the contribution of photons radiated classically. 
  * ``ndump_fac_n_emit`` - Controls the dump frequency of a grid diagnostic with the number of photon emission/pair production events per species in each PIC cell.  The diagnostic is stored in ``QED/n_emit/``.
  * ``ndump_fac_chi`` - Controls the dump frequency of a grid diagnostic with the weighted-average quantum parameter per species in each PIC cell. To get this diagnostic for photons, one needs to turn on the pair production flag. The diagnostic is stored in ``QED/chi_grid/``.

After these sections, a species section (and corresponding mandatory and/or optional udist, profile, spe_bound and diag_species sections) should follow for the electron, positron and photon species in the QED group (in this order). For the moment, we recommend using the standard (Boris) particle pusher is allowed for charged particles, so those species' sections should include ``push_type="standard"``.

After the blocs with the initialization of qed groups, one can set up QED processes in a Coulomb field. To do so, the following blocs should be added in this order:
* blocs to initialize ion species
* blocs to set up Coulomb processes

```
qed_bremsstrahlung
{
  if_bremsstrahlung = .true.,
  Z_ion = 29,
  proba_mult = 1.e5,
  energy_damp = .true.,
  i_ion(1) = 1,
}
```
In this section, the relevant parameters are:

  * ``if_bremsstrahlung`` - Activates / deactivate the Bremsstrahlung process
  * ``Z_ion`` - Atomic number of the ion species. There are only two atomic numbers Z=13 and Z=29 (tables are large and unconvenient).
  * ``proba_mult`` - Multiplicative factor for the Bremsstrahlung probability. To compensate for this statistical increase, the weight of macro-photons created is divided by the same factor.
  * ``energy_damp`` - Activate / deactive the energy damping of electrons experiencing Bremsstrahlung (for test purpose only).
  * ``i_ion(1:num_qed)`` - Array that links each qed group with an ion species. In this case, qed group 1 is linked and experiences Bremsstrahlung with ion species 1. 

  ```
qed_betheheitler
{
  if_betheheitler = .true.,
  Z_ion = 13,
  proba_mult = 1.e6,
  i_ion(1) = 1,
}
```
In this section, the relevant parameters are:

  * ``if_betheheitler`` - Activates / deactivate the Bethe Heitler process
  * ``Z_ion`` - Atomic number of the ion species. There are only two atomic numbers Z=13 and Z=29 (tables are large and unconvenient).
  * ``proba_mult`` - Multiplicative factor for the Bethe Heitler probability. To compensate for this statistical increase, the weight of macro-pairs created is divided by the same factor.
  * ``i_ion(1:num_qed)`` - Array that links each qed group with an ion species. In this case, qed group 1 is linked and experiences Bethe Heitler with ion species 1. 

```
qed_trident_coul
{
  if_trident_coul = .true.,
  Z_ion = 13,
  proba_mult = 1.e6,
  i_ion(1) = 1,
}
```
In this section, the relevant parameters are:

  * ``if_trident_coul`` - Activates / deactivate the Coulomb Trident process
  * ``Z_ion`` - Atomic number of the ion species. There are only two atomic numbers Z=13 and Z=29 (tables are large and unconvenient).
  * ``proba_mult`` - Multiplicative factor for the Coulomb Trident probability. To compensate for this statistical increase, the weight of macro-pairs created is divided by the same factor.
  * ``i_ion(1:num_qed)`` - Array that links each qed group with an ion species. In this case, qed group 1 is linked and experiences Coulomb Trident with ion species 1. 
