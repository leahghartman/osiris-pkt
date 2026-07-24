# Quasi-3D Algorithm (cylindrical modes)

## Overview

3D effects lead to important qualitative and quantitative physics differences, but simulations are computationally expensive. Fortunately, there exists a class of problems with sufficient azimuthal symmetry such that the dominant effects can be captured in just a few azimuthal modes. The quasi-3d algorithm takes advantage of this by using a truncated azimuthal mode expansion: electromagnetic fields, charge densities and current densities, written in cylindrical coordinates (*z*,*r*,*ϕ*), are expanded into a Fourier series in *ϕ*. In the current implementation, up to 15 of these modes can be kept. In this way, the quasi-3D algorithm obtains 3D physics at orders of magnitude speedup over a comparable, honest 3D simulation.

For more background information, find the full paper on the quasi-3d algorithm, "Implementation of a hybrid particle code with a PIC description in *r*–*z* and a gridless description in *ϕ* into OSIRIS" by Davidson et al. here: https://www.sciencedirect.com/science/article/pii/S0021999114007529 (DOI 10.1016/j.jcp.2014.10.064).

Note that the expansion in the paper is in exp(*imϕ*), where *i* is the imaginary unit, *m* is the summation index, and *ϕ* is the azimuthal coordinate. However, the code assumes the expansion is in exp(–*imϕ*). Both conventions are valid but give differnt signs in front of the sin modes. See below for an example of expressing the cartesian unit vectors in terms of azimuthal field components in the convention used in the code.

## Interpreting the data and coordinate system

The 2D coordinate system is defined such that *x*<sub>1</sub> and *x*<sub>2</sub> correspond to *z* and *r*, respectively.  We define *ϕ* such that *ϕ* = 0 and *ϕ* = *π*/2 correspond to the *x* and *y* planes, respectively.  To determine the value of a given quantity, say charge density ρ, at a certain value of *ϕ*, we use the real (ρ<sup>*m*</sup><sub>re</sub>) and imaginary (ρ<sup>*m*</sup><sub>im</sub>) parts of mode *m*:

ρ(*z*,*r*,*ϕ*) = ρ<sup>0</sup><sub>re</sub>(*z*,*r*) + ρ<sup>1</sup><sub>re</sub>(*z*,*r*)\*cos(*ϕ*) + ρ<sup>1</sup><sub>im</sub>(*z*,*r*)\*sin(*ϕ*) + ρ<sup>2</sup><sub>re</sub>(*z*,*r*)\*cos(2*ϕ*) + ρ<sup>2</sup><sub>im</sub>(*z*,*r*)\*sin(2*ϕ*) + ...

However, care must be taken when working with vector quantities.  For example, the electric field vectors in OSIRIS are defined such that *E*<sub>1</sub>, *E*<sub>2</sub> and *E*<sub>3</sub> correspond to *E*<sub>*z*</sub>, *E*<sub>*r*</sub> and *E*<sub>*ϕ*</sub>, respectively.  Thus at *ϕ* = 0, *E*<sub>2</sub> and *E*<sub>3</sub> vectors correspond to the cartesian vectors *E*<sub>*x*</sub> and *E*<sub>*y*</sub>, respectively, but respectively correspond to *E*<sub>*y*</sub> and –*E*<sub>*x*</sub> at *ϕ* = *π*/2.  So if one wanted to look at a laser field linearly polarized in the *x*-direction, the principal electric fields would be the *E*<sub>2</sub><sup>1</sup><sub>re</sub> and *E*<sub>3</sub><sup>1</sup><sub>im</sub> fields (along with the divergence-corrected *E*<sub>1</sub><sup>1</sup><sub>re</sub> field).

## Input parameters

* The quasi-3D algorithm is currently implemented as a simulation mode. To activate it the following parameter must be added to the `simulation` section:
  
  * `algorithm = "quasi-3D",`

* The following parameter has been added to the ``grid`` section

  * ``n_cyl_modes``, integer, default = 0

    This integer specifies how many azimuthal modes to use, ranging between 0 and 15.  Using 0 modes will assume complete azimuthal symmetry and run with only the real part of mode 0 (no imaginary part).

* The following parameters have been added to the ``el_mag_fld`` section
  
  * ``init_b0_re, init_e0_re, ext_b0_re, ext_e0_re``(0:``n_cyl_modes``, 1:``p_f_dim``), real, default =  0.0
  * ``init_b0_im, init_e0_im, ext_b0_im, ext_e0_im``(1:``n_cyl_modes``, 1:``p_f_dim``), real, default =  0.0
    
    Arrays for real and imaginary components of electric or magnetic fields. The first dimension dimension specifies the mode in the azimuthal expansion (e.g., mode 0 corresponds to index 0). The second dimension specifies the field component. Note that the order of the field components is (*z*,*r*,*ϕ*).

    Below is example of an external *B* field in quasi-3D in the *x*, *y* and *z* cartesian directions:
    ```
    el_mag_fld
    {
      ext_fld = 'static',
      type_ext_b(1:3) = 'uniform', 'uniform', 'uniform',

      ! b ext pointing in x-direction
      ext_b0_re(1,1:3) = 0.0, 1.0,  0.0,
      ext_b0_im(1,1:3) = 0.0, 0.0, -1.0,

      ! b ext pointing in y-direction
      ! ext_b0_re(1,1:3) = 0.0, 0.0, 1.0,
      ! ext_b0_im(1,1:3) = 0.0, 1.0, 0.0,

      ! b ext pointing in z-direction
      ! ext_b0_re(0,1:3) = 1.0, 0.0, 0.0,
    }
    ```

  * ``init_b_re_mfunc, init_e_re_mfunc, ext_b_re_mfunc, ext_e_re_mfunc``(0:``n_cyl_modes``, 1:``p_f_dim``), character(\*), default = "NO_FUNCTION_SUPPLIED!"
  * ``init_b_im_mfunc, init_e_im_mfunc, ext_b_im_mfunc, ext_e_im_mfunc``(1:``n_cyl_modes``, 1:``p_f_dim``), character(\*), default = "NO_FUNCTION_SUPPLIED!"

    Arrays for the characters specifying the math functions to use for the components of electric or magnetic fields.

  See the input deck decks/cyl_modes/os-stdin-emf-ext for examples on using constant or function external fields.

* The quasi-3D simluation mode has a unique set of reports available in the ``diag_emf`` section. These are:

  * ``e1_cyl_m``, ``e2_cyl_m``, ``e3_cyl_m``
  * ``b1_cyl_m``, ``b2_cyl_m``, ``b3_cyl_m``
  * ``ext_e1_cyl_m ``, ``ext_e2_cyl_m ``, ``ext_e3_cyl_m ``
  * ``ext_b1_cyl_m ``, ``ext_b2_cyl_m ``, ``ext_b3_cyl_m ``
  * ``part_e1_cyl_m``, ``part_e2_cyl_m``, ``part_e3_cyl_m``
  * ``part_b1_cyl_m``, ``part_b2_cyl_m``, ``part_b3_cyl_m``
  * ``ene_e1_cyl_m``, ``ene_e2_cyl_m``, ``ene_e3_cyl_m``
  * ``ene_b1_cyl_m``, ``ene_b2_cyl_m``, ``ene_b3_cyl_m``
  * ``ene_e_cyl_m``, ``ene_b_cyl_m``, ``ene_emf_cyl_m``
  * ``div_e_cyl_m``, ``div_b_cyl_m``, ``psi_cyl_m``
  * ``s1_int_cyl_m``, ``s2_int_cyl_m``, ``s3_int_cyl_m``

  Reports with suffix "_cyl_m" are written to disk by azimuthal mode.

  The reports ``s1_int_cyl_m``, ``s2_int_cyl_m``, ``s3_int_cyl_m`` are defined as poynting flux components averaged over the azimuthal angle *ϕ*, i.e. ``si_int_cyl_m`` = integral_0^{2*π*} si(*z*,*r*,*ϕ*) d*ϕ* /(2*π*). This quantity is of course defined in 2 spatial dimensions, *z* and *r*, and so it is dumped into a single 2 dimensional grid.

* The following parameters have been added to the ``species`` section

  * ``num_par_theta``, integer, default = –1

    This parameter specifies how many particles to space in rings in *ϕ*.  We recommend this number be at least 8\*``n_cyl_modes``.

  * ``rand_theta``, logical, default = .false.

    This flag specifies if a random offset should be added in theta to particle rings in each cell (one random number generated per cell if .true.).  This can help sample modes more efficiently, but may also add noise to the simulation.

  * ``shift_cell_spokes``, logical, default = .true.

    This flag specifies if the individual rings of particles within a given cell for ``num_par`` > 1 will be offset from one another (.true.) or all aligned at the same values of theta (.false.).  The default is .true. to sample more of the theta space.

* The following parameters have been added to the ``profile`` section

  * ``aspect_ratio``, real, default = 1.0

    When using a Gaussian particle distribution in *r*, this parameter specifies the ratio of sigma_x/sigma_y.  When ``aspect_ratio`` is not equal to 1, the gauss_sigma(2) value always specifies the sigma_x value [sigma_y = gauss_sigma(2) * ``aspect_ratio``].

* The following parameters have been added to the ``diag_species`` section

  New reports: ``charge_cyl_m`` (reports all modes for charge density).

  * ``ps_xcmin``(1:3), real, default = simulation minimum
  * ``ps_xcmax``(1:3), real, default = simulation maximum
  * ``ps_nxc``(1:3), integer, default = 64
  * ``ps_nxc_3D``(1:3), integer, default = ps_nxc

  New phasespaces: ``xc2`` and ``xc3``

  The above parameters control the ``xc2`` and ``xc3`` phasespaces, which are x-cartesian data.  In other words, since the particles exist with cartesian and cylindrical coordinates, these phasespaces collect particle data in those coordinates (*x* and *y* directions).  The first entry in each of the ``ps_xc``\* phasespace parameters is meaningless, and the ``xc1`` phasespace does not exist (since it would be the same as ``x1``).  See the input deck in decks/cyl_modes/os-stdin for usage examples.

* The following reports have been added to the ``diag_current`` section

  * ``j1_cyl_m``, ``j2_cyl_m``, ``j3_cyl_m``

