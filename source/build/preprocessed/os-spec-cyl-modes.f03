# 1 "cyl_modes/os-spec-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-spec-cyl-modes.f03" 2
! m_species_cyl_modes module
!
! Handles the species object for the quasi-3D algorithm

# 1 "./os-config.h" 1
! Configuration file for osiris

! ----------------------------------------------------------------------------------------
! Algorithm options
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! System options
! ----------------------------------------------------------------------------------------

! MPI supports MPI_IN_PLACE in global operations
!#define 1

! Use OpenMP
!#define 1

! Use SIMD optimized code
!#define SIMD
!#define SIMD_SSE
!#define SIMD_AVX
!#define SIMD_BGQ
!#define SIMD_MIC

! Use SION for checkpointing (not available on all systems)
!#define __RST_IO__ = __RST_SION__

! Use log files
!#define __USE_LOG__

! Use MPE for logging and profiling
!#define __USE_MPE__

! Compiler does not fully support the sizeof intrinsic which is a Fortran 2003 feature
!#define __NO_SIZEOF__

! Use the PAPI library for profiling
!#define __USE_PAPI__


! ----------------------------------------------------------------------------------------
! Distribution Options
! ----------------------------------------------------------------------------------------
!
! Optional modules to be removed for distribution
! These are all turned on by default unless the __DISTRO__ preprocessor macro is defined
!
! ----------------------------------------------------------------------------------------



! Include Ionization module


! Include binary collisions module


! Include particle tracking module


! Include perfectly matched layers boundary conditions for EMF


! Include spin advance module which is disabled by default
! #define __HAS_SPIN__

! Include a debug flag for cylindrical modes simulations
! #define __CYL_MODES_DEBUG__

! Compile functions for reporting B fields in the RaDiO module
! #define __HAS_RAD_BFLD__
# 6 "cyl_modes/os-spec-cyl-modes.f03" 2
# 1 "./os-preprocess.fpp" 1
!
! File: os-preprocess.fpp
!
! A set of preprocessing macros for osiris, and the fpp/cpp preprocessor
!
!



! Macros for the IBM CPP/GNU preprocessor




! Assertion macro



! Debug functions
# 33 "./os-preprocess.fpp"
! SCR hostname functions







! LOG functions




! ERROR functions






! WARNING functions




! functions for restart io
# 7 "cyl_modes/os-spec-cyl-modes.f03" 2

module m_species_cyl_modes

# 1 "./memory/memory.h" 1
! Include file for the memory module
! Files must include this file using #include "memory/memory.h" rather than just using the module
! the #include must be placed where the module use statement would usually be i.e.
!
! module module2
!
! use module1
! #include "memory/memory.h"
!
# 19 "./memory/memory.h"
use memory
# 11 "cyl_modes/os-spec-cyl-modes.f03" 2
use m_system
use m_parameters
use stringutil, only: replace_blanks
use m_fparser, only: p_max_expr_len, setup
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_species_define, only: p_max_phasespace_namelen, p_max_n_ene_bins, &
                                      t_track_set, t_species, t_spe_bound, p_thermal, &
                                      p_waterbag, p_random_dir
use m_species_diagnostics
use m_species_phasespace, only: setup
use m_space, only: t_space, nx_move
use m_grid_define, only: t_grid
use m_node_conf, only: t_node_conf
use m_restart, only: t_restart_handle, restart_io_write, restart_io_read
use m_time_step, only: t_time_step
use m_input_file, only: t_input_file, get_namelist
use m_cyl_modes, only: t_cyl_modes, update_boundary, setup, p_real, p_imag, &
                                      cleanup, copy_component
use m_species_cyl_modes_define, only: t_species_cyl_modes
use m_species_current_cyl_modes,only: get_coef
use m_species_current, only: tmp_xold, tmp_xnew, tmp_p, tmp_rg, tmp_q, tmp_dxi, &
                                      tmp_ix
use m_current_define, only: t_current
use m_vdf_define, only: t_vdf, t_vdf_report, p_max_reports_len, p_n_report_type
use m_vdf_define, only: p_max_reports, p_full, p_savg, p_senv, p_line, p_slice
use m_vdf_report, only: new
use m_emf_interpolate, only: get_emf
use m_emf_interpolate_cm, only: get_emf_cyl_modes
use m_species_phasespace_cyl_modes, only: t_phasespace_diag_cyl_modes
use m_diagnostic_utilities, only: p_diag_prec
use m_input_file, only: disp_out
use m_species_charge_cyl_modes, only: norm_charge_cyl_m
use m_vdf_comm, only: t_vdf_msg, update_boundary
use m_emf_psi, only: get_psi
use m_vdf_interpolate, only: interpolate
use m_random, only: rng
use m_boosted_diag_cyl_modes, only: t_boosted_diag_cyl_modes, p_boosted_4vector
use m_logprof, only: create_event



  use m_species_tracks



implicit none

private

interface deposit_quant_cyl_modes
  module procedure deposit_quant_cyl_modes
end interface

interface add_track_set_data_cyl_modes
  module procedure add_track_set_data_cyl_modes
end interface

interface interpolate_cyl_modes
  module procedure interpolate_cyl_modes
end interface

interface thermal_bath_cyl_modes
  module procedure thermal_bath_cyl_modes
end interface


type, extends( t_diag_species ) :: t_diag_species_cyl_modes

  integer :: n_cyl_modes

  contains

  procedure :: allocate_objs => allocate_objs_diag_spec_cyl_modes
  procedure :: avail_report_quants => avail_report_quants_species_cyl_modes
  procedure :: init_report_quants => init_report_quants_species_cyl_modes
  procedure :: init => init_diag_spec_cyl_modes
  procedure :: read_input => read_input_diag_species_cyl_modes

end type t_diag_species_cyl_modes

type, extends( t_spe_bound ) :: t_spe_bound_cyl_modes

  ! No class-specific member data

contains

  procedure, nopass :: init_tmp_buf_current => init_tmp_buf_current_spe_bnd_cyl

end type t_spe_bound_cyl_modes

! diagnostic quantities specific to the quasi-3D geometry
character(len=12), dimension(9), parameter, public :: p_report_quants_spec_cyl_modes = &
   (/ 'charge_cyl_m', &
      'mass_cyl_m  ', &
      'ene_cyl_m   ', &
      'q1_cyl_m    ', &
      'q2_cyl_m    ', &
      'q3_cyl_m    ', &
      'j1_cyl_m    ', &
      'j2_cyl_m    ', &
      'j3_cyl_m    ' /)

character(len=18), dimension(4), parameter, public :: &
  p_boosted_report_quants_spec_cyl_modes = (/ &
                                'charge_boost_cyl_m', &
                                'j1_boost_cyl_m    ', &
                                'j2_boost_cyl_m    ', &
                                'j3_boost_cyl_m    ' /)

! cyl_modes specific diagnostics
integer, parameter, public :: p_charge_cyl_m = p_charge
integer, parameter, public :: p_mass_cyl_m = p_mass
integer, parameter, public :: p_ene_cyl_m = p_ene
integer, parameter, public :: p_q1_cyl_m = p_q1
integer, parameter, public :: p_q2_cyl_m = p_q2
integer, parameter, public :: p_q3_cyl_m = p_q3
integer, parameter, public :: p_j1_cyl_m = p_j1
integer, parameter, public :: p_j2_cyl_m = p_j2
integer, parameter, public :: p_j3_cyl_m = p_j3

public :: t_diag_species_cyl_modes, t_spe_bound_cyl_modes
public :: deposit_quant_cyl_modes, get_current_4vector_cyl_modes
public :: interpolate_cyl_modes, add_track_set_data_cyl_modes, thermal_bath_cyl_modes
! -----------------------------------------------------------------------------
contains

!-----------------------------------------------------------------------------------------
! Polymorphic diagnostic functions
!-----------------------------------------------------------------------------------------

subroutine allocate_objs_diag_spec_cyl_modes( this )

  class( t_diag_species_cyl_modes ), intent(inout) :: this

  if ( .not. associated( this%phasespaces ) ) then
    allocate( t_phasespace_diag_cyl_modes :: this%phasespaces )
  endif

  if ( .not. associated( this%boosted_diag ) ) then
    allocate( t_boosted_diag_cyl_modes :: this%boosted_diag )
  endif

  call this % t_diag_species % allocate_objs()

end subroutine

function avail_report_quants_species_cyl_modes( this )

  class( t_diag_species_cyl_modes ), intent(in) :: this
  integer :: avail_report_quants_species_cyl_modes

  avail_report_quants_species_cyl_modes = size( p_report_quants_spec_cyl_modes )

end function avail_report_quants_species_cyl_modes


subroutine init_report_quants_species_cyl_modes( this )

  class( t_diag_species_cyl_modes ), intent(inout) :: this

  if ( .not. associated( this % report_quants ) ) then
    allocate( this % report_quants( this % avail_report_quants( ) ) )
  endif

  ! These reports, and only these reports
  this % report_quants( 1 : size(p_report_quants_spec_cyl_modes) ) = p_report_quants_spec_cyl_modes


  ! Boosted diagnostic reports
  if ( .not. associated( this % boosted_report_quants ) ) then
    allocate( this % boosted_report_quants( size(p_boosted_report_quants_spec_cyl_modes) ) )
  endif

  this % boosted_report_quants( 1 : size( p_boosted_report_quants_spec_cyl_modes ) ) = p_boosted_report_quants_spec_cyl_modes

end subroutine init_report_quants_species_cyl_modes


subroutine init_diag_spec_cyl_modes( this, spec_name, n_x_dim, ndump_fac, interpolation, &
                                     restart, restart_handle )

  implicit none

  class ( t_diag_species_cyl_modes ), intent( inout ) :: this
  character( len=* ), intent(in) :: spec_name
  integer, intent(in) :: n_x_dim
  integer, intent(in) :: ndump_fac
  integer, intent(in) :: interpolation
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  integer :: ierr, izero = ichar('0')

!-----------------------------------------------------------------------------------------
! copy/paste of init_diag_species
!-----------------------------------------------------------------------------------------

  if ( ndump_fac > 0 ) then

    ! setup density diagnostics
    call setup_report_list( this%reports, spec_name, 'DENSITY', p_rep_type_dens, &
                            interpolation )

    ! setup cell average diagnostics
    call setup_report_list( this%rep_cell_avg, spec_name, 'CELL_AVG', p_rep_type_dens, &
                            interpolation, label_prefix = 'Cell average')

    ! setup u dist diagnostics
    call setup_report_list( this%rep_udist, spec_name, 'UDIST', p_rep_type_udist, &
                            interpolation )

    ! setup phasespaces
    call setup( this%phasespaces )

    ! set variable list for raw_math_expr
    if (this%raw_math_expr /= '') then

      !! This is the only difference in the quasi-3D code...
      call setup(this%raw_func, trim(this%raw_math_expr), &
                (/'x1', 'x2', 'x3', 'x4' ,'p1', 'p2', 'p3', 'g ', 't '/), ierr)

      ! check if function compiled ok
      if (ierr /= 0) then
        if (mpi_node() == 0) then
          write(0,*) "(*error*) Unable to compile 'raw_math_expr' function:"
          write(0,*) "(*error*) '", trim(this%raw_math_expr), "'"
          write(0,*) "(*error*) aborting..."
        endif
        call abort_program(-1)
      endif
    endif



    ! modify nfields now that we know the number of modes
    ! we need one for the total field, and an additional entry for each mode
    ! (including m=0)
    this%tracks%nfields = count(this%tracks%ifdmp_tracks_efl) + &
                          count(this%tracks%ifdmp_tracks_bfl)
    this%tracks%nfields = this%tracks%nfields * ( this%n_cyl_modes + 2 )

    if ( this%tracks%ifdmp_tracks_psi ) this%tracks%nfields = this%tracks%nfields + 1

    ! setup tracking diagnostic
    call setup( this%tracks, spec_name, n_x_dim, ndump_fac*this%ndump_fac_tracks, restart, restart_handle )



  endif


  contains

  subroutine setup_report_list( list, spec_name, path, rep_type, interpolation, label_prefix )

    type( t_vdf_report ), pointer :: list
    character(len=*), intent(in) :: spec_name, path
    integer, intent(in) :: rep_type
    integer, intent(in) :: interpolation
    character(len=*), optional, intent(in) :: label_prefix

    character(len=64) :: prefix
    type( t_vdf_report ), pointer :: report
    real(p_double) :: interp_offset

    if ( present(label_prefix )) then
      prefix = trim(label_prefix) // " "
    else
      prefix = ''
    endif

    select case( interpolation )
    case( p_linear, p_cubic )
      interp_offset = 0.0_p_double
    case( p_quadratic, p_quartic )
      interp_offset = 0.5_p_double
    case default
      interp_offset = 0.0_p_double
      write(err_buf__,*) 'Interpolation value not supported';call err__("cyl_modes/os-spec-cyl-modes.f03",289)
      call abort_program(p_err_invalid)
    end select

    report => list
    do
      if ( .not. associated( report ) ) exit

      report%fileLabel = replace_blanks(trim(spec_name))
      report%basePath = trim(path_mass) // trim(path) // p_dir_sep // &
                          trim(report%fileLabel) // p_dir_sep

      report%xname = (/'x1', 'x2', 'x3'/)
      report%xlabel = (/'x_1', 'x_2', 'x_3'/)
      report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

      ! these are just dummy values for now
      report%time_units = '1 / \omega_p'
      report%dt = 1.0

      ! set default units
      report%units = 'n_0'

      ! set units and labels
      select case ( rep_type )
      case ( p_rep_type_dens )

        select case ( report%quant )
        case(p_charge_cyl_m)
          report%units = 'e \omega_p^'//char(izero+p_x_dim)// &
                          '/ c^'//char(izero+p_x_dim)
          report%label = '\rho'
          report%offset_t = 0.0_p_double
        case(p_mass_cyl_m)
          report%label = 'm'
          report%offset_t = 0.0_p_double
        case(p_ene_cyl_m)
          report%label = 'Kinetic Energy'
          report%offset_t = -0.5_p_double
        case(p_q1_cyl_m)
          report%label = 'q_1'
          report%offset_t = -0.5_p_double
        case(p_q2_cyl_m)
          report%label = 'q_2'
          report%offset_t = -0.5_p_double
        case(p_q3_cyl_m)
          report%label = 'q_3'
          report%offset_t = -0.5_p_double
        case(p_j1_cyl_m)
          report%label = 'j_1'
          report%offset_t = -0.5_p_double
        case(p_j2_cyl_m)
          report%label = 'j_2'
          report%offset_t = -0.5_p_double
        case(p_j3_cyl_m)
          report%label = 'j_3'
          report%offset_t = -0.5_p_double
        end select

      case ( p_rep_type_udist )

        select case ( report%quant )
        case(p_ufl1)
          report%label = 'u_{fl1}'
          report%units = 'c'
        case(p_ufl2)
          report%label = 'u_{fl2}'
          report%units = 'c'
        case(p_ufl3)
          report%label = 'u_{fl3}'
          report%units = 'c'
        case(p_uth1)
          report%label = 'u_{th1}'
          report%units = 'c'
        case(p_uth2)
          report%label = 'u_{th2}'
          report%units = 'c'
        case(p_uth3)
          report%label = 'u_{th3}'
          report%units = 'c'
        case(p_vfl1)
          report%label = 'v_{fl1}'
          report%units = 'c'
        case(p_vfl2)
          report%label = 'v_{fl2}'
          report%units = 'c'
        case(p_vfl3)
          report%label = 'v_{fl3}'
          report%units = 'c'
        case(p_nufl1)
          report%label = 'nu_{fl1}'
          report%units = 'n_0c'
        case(p_nufl2)
          report%label = 'nu_{fl2}'
          report%units = 'n_0c'
        case(p_nufl3)
          report%label = 'nu_{fl3}'
          report%units = 'n_0c'
        case(p_T11)
          report%label = 'T_{11}'
          report%units = 'm c^2'
        case(p_T22)
          report%label = 'T_{22}'
          report%units = 'm c^2'
        case(p_T33)
          report%label = 'T_{33}'
          report%units = 'm c^2'
        case(p_T12)
          report%label = 'T_{12}'
          report%units = 'm c^2'
        case(p_T13)
          report%label = 'T_{13}'
          report%units = 'm c^2'
        case(p_T23)
          report%label = 'T_{23}'
          report%units = 'm c^2'
        case(p_P11)
          report%label = 'P_{11}'
          report%units = 'm c^2'
        case(p_P22)
          report%label = 'P_{22}'
          report%units = 'm c^2'
        case(p_P33)
          report%label = 'P_{33}'
          report%units = 'm c^2'
        case(p_P12)
          report%label = 'P_{12}'
          report%units = 'm c^2'
        case(p_P13)
          report%label = 'P_{13}'
          report%units = 'm c^2'
        case(p_P23)
          report%label = 'P_{23}'
          report%units = 'm c^2'
        end select

        report%offset_t = -0.5_p_double

      end select

      report % label = trim(prefix) // trim(report % label)
      report%offset_x = (/ 0.0_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)

      report => report%next
    enddo

end subroutine setup_report_list

!-----------------------------------------------------------------------------------------
! end copy/paste of init_diag_species
!-----------------------------------------------------------------------------------------

end subroutine


!-----------------------------------------------------------------------------------------
subroutine read_input_diag_species_cyl_modes( this, input_file, gamma )
!-----------------------------------------------------------------------------------------
! read necessary information from input file
!-----------------------------------------------------------------------------------------

  implicit none

  class ( t_diag_species_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  real(p_double), intent(in) :: gamma

  integer :: ndump_fac_pha, ndump_fac_pha_tavg, &
         ndump_fac_ene, ndump_fac_heatflux, ndump_fac_raw, &
         ndump_fac, ndump_fac_ave, &
         ndump_fac_lineout, ndump_fac_temp

  integer, dimension(3) :: n_ave

  ! physical range for phasespace data dumps
  real(p_diag_prec), dimension(p_x_dim) :: ps_xmin, ps_xmax
  real(p_diag_prec), dimension(p_p_dim) :: ps_pmin, ps_pmax
  real(p_diag_prec), dimension(p_p_dim) :: ps_lmin, ps_lmax
  real(p_diag_prec), dimension(p_max_dim) :: ps_xcmin, ps_xcmax

  ! switch for autorange of ps_p
  logical, dimension(p_p_dim) :: if_ps_p_auto

  ! switch for autorange of l for phasespaces
  logical, dimension(p_p_dim) :: if_ps_l_auto

  ! physical range for 1D momenta phasespace data dumps
  real(p_diag_prec) :: ps_gammamin, ps_gammamax
  real(p_diag_prec) :: ps_kemin, ps_kemax

  ! switch for autorange of 1D momenta phasespace data dumps
  logical :: if_ps_gamma_auto
  logical :: if_ps_ke_auto

  ! resolutions for phasespace data dumps
  integer, dimension(p_x_dim) :: ps_nx, ps_nx_3D
  integer, dimension(p_p_dim) :: ps_np, ps_np_3D
  integer, dimension(p_p_dim) :: ps_nl, ps_nl_3D
  integer, dimension(p_max_dim) :: ps_nxc, ps_nxc_3D

  ! resolution for 1D momenta phasespace data dump
  integer :: ps_ngamma
  integer :: ps_nke

  ! parameters for raw data dump
  real(p_diag_prec) :: raw_gamma_limit
  real(p_diag_prec) :: raw_fraction
  character(len = p_max_expr_len) :: raw_math_expr
  logical, dimension(p_x_dim) :: raw_if_pos_ref_box

  ! parameters for energy binned diagnostics
  integer :: n_ene_bins
  real(p_diag_prec), dimension(p_max_n_ene_bins) :: ene_bins

  integer :: ndump_fac_tracks, niter_tracks, n_start_tracks
  character(len = p_max_filename_len) :: file_tags

  logical, dimension(p_f_dim) :: ifdmp_tracks_efl
  logical, dimension(p_f_dim) :: ifdmp_tracks_bfl
  logical :: ifdmp_tracks_psi

  ! phasespace input
  integer, parameter :: p_max_phasespaces = 64
  character(len = p_max_phasespace_namelen), &
      dimension(p_max_phasespaces) :: phasespaces
  character(len = p_max_phasespace_namelen), &
      dimension(p_max_phasespaces) :: pha_ene_bin
  character(len = p_max_phasespace_namelen), &
      dimension(p_max_phasespaces) :: pha_cell_avg
  character(len = p_max_phasespace_namelen), &
      dimension(p_max_phasespaces) :: pha_time_avg

  ! number of time steps to average over for time averaged dumps
  integer :: n_tavg

  ! New reports
  character( len = p_max_reports_len ), dimension( p_max_reports ) :: &
     reports, rep_cell_avg, rep_udist

  integer :: prec

  namelist /nl_diag_species/ &
       ndump_fac_pha, ndump_fac_pha_tavg, ndump_fac_lineout, ndump_fac_ene, ndump_fac_heatflux, &
       ndump_fac_temp, ndump_fac_raw, ndump_fac, ndump_fac_ave, n_ave, prec, &
       ps_xmin, ps_xmax, ps_pmin, ps_pmax, ps_lmin, ps_lmax, ps_xcmin, ps_xcmax, &
       if_ps_p_auto, if_ps_l_auto, &
       ps_gammamin, ps_gammamax, if_ps_gamma_auto, &
       ps_kemin, ps_kemax, if_ps_ke_auto, &
       ps_nx, ps_nx_3D, ps_np, ps_np_3D, ps_nl, ps_nl_3D, ps_nxc, ps_nxc_3D, ps_ngamma, &
       ps_nke, raw_gamma_limit, raw_fraction, raw_math_expr, raw_if_pos_ref_box, &
       n_ene_bins, ene_bins, &
       ndump_fac_tracks, n_start_tracks, niter_tracks, file_tags, ifdmp_tracks_efl, ifdmp_tracks_bfl, ifdmp_tracks_psi, &
       phasespaces, pha_ene_bin, pha_cell_avg, pha_time_avg, n_tavg, &
       reports, rep_cell_avg, rep_udist

  integer, dimension( p_n_report_type ) :: ndump_fac_all
  integer :: ierr,i

  reports = "-"
  rep_cell_avg = "-"
  rep_udist = "-"

  ndump_fac = 0
  ndump_fac_ave = 0
  ndump_fac_lineout = 0
  n_ave = 0
  n_tavg = 0
  prec = p_diag_prec

  ndump_fac_ene = 0
  ndump_fac_heatflux = 0
  ndump_fac_temp = 0
  ndump_fac_raw = 0

  ndump_fac_pha = 0
  ndump_fac_pha_tavg = 0

  ps_xmin = 0
  ps_xmax = 0
  ps_pmin = 0
  ps_pmax = 0
  ps_lmin = 0
  ps_lmax = 0
  ps_xcmin = 0
  ps_xcmax = 0

  if_ps_p_auto = .false.
  if_ps_l_auto = .false.

  ps_gammamin = 1.0_p_diag_prec
  ps_gammamax = 0
  ps_kemin = 0.000001_p_diag_prec
  ps_kemax = 0
  if_ps_gamma_auto = .false.
  if_ps_ke_auto = .false.
  ps_nx = 64
  ps_nx_3D = -1
  ps_np = 64
  ps_np_3D = -1

  ps_nl = 64
  ps_nl_3D = -1

  ps_nxc = 64
  ps_nxc_3D = -1

  ps_ngamma = 64
  ps_nke = 64

  raw_gamma_limit = 0
  raw_fraction = 1.0_p_diag_prec

  ! default values for math function
  raw_math_expr = ''

  raw_if_pos_ref_box = .false.

  n_ene_bins = 0
  ene_bins = 0

  ndump_fac_tracks = 0
  n_start_tracks = -1
  niter_tracks = 1
  file_tags = ''

  ifdmp_tracks_efl = .false.
  ifdmp_tracks_bfl = .false.
  ifdmp_tracks_psi = .false.

  phasespaces = "-"
  pha_ene_bin = "-"
  pha_cell_avg = "-"
  pha_time_avg = "-"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_diag_species", ierr )

  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_diag_species, iostat = ierr)
    if (ierr /= 0) then
      print *, ""
      print *, "   Error reading diag_species parameters "
      print *, "   aborting..."
      stop
    endif
  else
    if (disp_out(input_file)) then
      if (mpi_node()==0) print *,"   - no diagnostics specified"
    endif
  endif

  this%ndump_fac_ene = ndump_fac_ene
  this%ndump_fac_temp = ndump_fac_temp
  this%ndump_fac_heatflux = ndump_fac_heatflux
  this%ndump_fac_raw = ndump_fac_raw

  ! Phasespace diagnostics data
  this%phasespaces%ndump_fac = ndump_fac_pha
  this%phasespaces%ndump_fac_tavg = ndump_fac_pha_tavg

  this%phasespaces%xmin = ps_xmin
  this%phasespaces%xmax = ps_xmax

  do i = 1, p_x_dim
    if ( ps_nx(i) < 2 ) then
      print "(A,I0,A)", " (*warning*) invalid ps_nx(",i,"), setting to default value (64)"
      ps_nx(i) = 64
    endif
  enddo
  this%phasespaces%nx = ps_nx

  this%phasespaces%pmin = ps_pmin
  this%phasespaces%pmax = ps_pmax
  this%phasespaces%if_p_auto = if_ps_p_auto

  do i = 1, p_p_dim
    if ( ps_np(i) < 2 ) then
      print "(A,I0,A)", " (*warning*) invalid ps_np(",i,"), setting to default value (64)"
      ps_np(i) = 64
    endif
  enddo
  this%phasespaces%np = ps_np

  this%phasespaces%lmin = ps_lmin
  this%phasespaces%lmax = ps_lmax
  this%phasespaces%if_l_auto = if_ps_l_auto

  do i = 1, p_p_dim
    if ( ps_nl(i) < 2 ) then
      print "(A,I0,A)", " (*warning*) invalid ps_nl(",i,"), setting to default value (64)"
      ps_nl(i) = 64
    endif
  enddo
  this%phasespaces%nl = ps_nl

  this%phasespaces%if_gamma_auto = if_ps_gamma_auto
  if ( ps_gammamin < 1.0_p_diag_prec ) then
    print "(A)", " (*warning*) invalid ps_gammamin, setting to 1.0"
    ps_gammamin = 1.0_p_diag_prec
  endif
  this%phasespaces%gammamin = ps_gammamin
  this%phasespaces%gammamax = ps_gammamax


  if ( ps_ngamma < 2 ) then
    print "(A)", " (*warning*) invalid ps_ngamma, setting to default value (64)"
    ps_ngamma = 64
  endif
  this%phasespaces%ngamma = ps_ngamma

  if (n_ene_bins > p_max_n_ene_bins) then
    print *, "(*error*) n_ene_bins must be <= ", p_max_n_ene_bins
    stop
  endif

  this%phasespaces%if_ke_auto = if_ps_ke_auto
  if ( ps_kemin <= 0.0_p_diag_prec ) then
    print "(A)", " (*warning*) invalid ps_kemin, setting to 1e-6"
    ps_kemin = 0.000001_p_diag_prec
  endif
  this%phasespaces%kemin = ps_kemin
  this%phasespaces%kemax = ps_kemax


  if ( ps_nke < 2 ) then
    print "(A)", " (*warning*) invalid ps_nke, setting to default value (64)"
    ps_nke = 64
  endif
  this%phasespaces%nke = ps_nke

  this%phasespaces%n_ene_bins = n_ene_bins
  do i=1, n_ene_bins
    this%phasespaces%ene_bins(i) = ene_bins(i)
  enddo

  if ( ps_nx_3D(1) == -1 ) then
    this%phasespaces%nx_3D = ps_nx
  else
    this%phasespaces%nx_3D = ps_nx_3D
  endif
  if ( ps_np_3D(1) == -1 ) then
    this%phasespaces%np_3D = ps_np
  else
    this%phasespaces%np_3D = ps_np_3D
  endif
  if ( ps_nl_3D(1) == -1 ) then
    this%phasespaces%nl_3D = ps_nl
  else
    this%phasespaces%nl_3D = ps_nl_3D
  endif

  this%phasespaces%n_tavg = n_tavg


  !---------------------------------------------------------
  ! New parameters for the t_phasespace_diag_cyl_modes class
  !---------------------------------------------------------
  select type( phasespaces_cyl => this%phasespaces )
  class is( t_phasespace_diag_cyl_modes )

    phasespaces_cyl%xcmin = ps_xcmin
    phasespaces_cyl%xcmax = ps_xcmax

    do i = 1, p_x_dim
      if ( ps_nxc(i) < 2 ) then
        print "(A,I0,A)", " (*warning*) invalid ps_nxc(",i,"), setting to default value (64)"
        ps_nxc(i) = 64
      endif
    enddo
    phasespaces_cyl%nxc = ps_nxc

    if ( ps_nxc_3D(1) == -1 ) then
      phasespaces_cyl%nxc_3D = ps_nxc
    else
      phasespaces_cyl%nxc_3D = ps_nxc_3D
    endif

  end select
  !---------------------------------------------------------


  ! Raw diagnostics data
  this%raw_gamma_limit = raw_gamma_limit
  this%raw_fraction = raw_fraction

  ! set math func variables
  this%raw_math_expr = trim(raw_math_expr)

  ! raw_if_pos_ref_box is a boolean. If .true., position in raw diagnostics will be saved
  ! as the position in reference to the box edge.
  ! Turn this on in the moving window direction to prevent roundoff errors.
  this%raw_if_pos_ref_box = raw_if_pos_ref_box


  ! store particle tracking data


  this%ndump_fac_tracks = ndump_fac_tracks
  this%n_start_tracks = n_start_tracks
  if ( this%ndump_fac_tracks > 0 ) then
    this%tracks%niter = niter_tracks
    this%tracks%file_tags = file_tags
    if ( .not. file_exists(file_tags)) then
      print *, "(*error*) ndump_fac_tracks is set but tags file '"//trim(file_tags)// &
                "' cannot be found."
      stop
    endif

    this%tracks%ifdmp_tracks_efl = ifdmp_tracks_efl
    this%tracks%ifdmp_tracks_bfl = ifdmp_tracks_bfl
    this%tracks%nfields = count(ifdmp_tracks_efl) + count(ifdmp_tracks_bfl)

    this%tracks%ifdmp_tracks_psi = ifdmp_tracks_psi
    if ( ifdmp_tracks_psi ) this%tracks%nfields = this%tracks%nfields + 1

  endif
# 815 "cyl_modes/os-spec-cyl-modes.f03"
  ! process phasespaces
  call this%phasespaces%init_list( this%phasespaces%phasespace_list, phasespaces, "phasespace" )

  ! process energy binned phasespaces
  call this%phasespaces%init_list( this%phasespaces%pha_ene_bin_list, pha_ene_bin, "pha_ene_bin" )

  ! process cell average phasespace
  call this%phasespaces%init_list( this%phasespaces%pha_cell_avg_list, pha_cell_avg, "pha_cell_avg" )

  ! process time average phasespace
  call this%phasespaces%init_list( this%phasespaces%pha_time_avg_list, pha_time_avg, "pha_time_avg", &
                                    time_average = .true. )

  ! process new density reports
  ndump_fac_all(p_full) = ndump_fac
  ndump_fac_all(p_savg) = ndump_fac_ave
  ndump_fac_all(p_senv) = ndump_fac_ave
  ndump_fac_all(p_line) = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  call this % init_report_quants( )

  call new( this%reports, reports, this % report_quants, &
            ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
    write(0,*) "(*error*) diag_species section:"
    write(0,*) "(*error*) Invalid report, aborting..."
    stop
  endif

  call new( this%rep_cell_avg, rep_cell_avg, this % report_quants, &
            ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
    write(0,*) "(*error*) diag_species section:"
    write(0,*) "(*error*) Invalid cell average report, aborting..."
    stop
  endif

  call new( this%rep_udist, rep_udist, p_rep_udist, &
            ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
    write(0,*) "(*error*) diag_species section:"
    write(0,*) "(*error*) Invalid cell average report, aborting..."
    stop
  endif

  select type( boosted_diag => this%boosted_diag ); class is ( t_boosted_diag_cyl_modes )
    boosted_diag%n_cyl_modes = this%n_cyl_modes
  end select
  call this % boosted_diag % read_input( input_file, this%boosted_report_quants, gamma, &
                                         p_boosted_4vector, this%if_use_boosted_diag, &
                                         this%if_use_boosted_raw, p_x_dim )

end subroutine read_input_diag_species_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_tmp_buf_current_spe_bnd_cyl()
!-----------------------------------------------------------------------------------------
! initialize buffers for current deposition with extra space
!-----------------------------------------------------------------------------------------

  implicit none

  if ( .not. associated( tmp_xold ) ) then
    call alloc(tmp_xold, (/ p_x_dim+2, p_cache_size /),"cyl_modes/os-spec-cyl-modes.f03",883)
    call alloc(tmp_xnew, (/ p_x_dim+2, p_cache_size /),"cyl_modes/os-spec-cyl-modes.f03",884)
    call alloc(tmp_p, (/ p_p_dim, p_cache_size /),"cyl_modes/os-spec-cyl-modes.f03",885)
    call alloc(tmp_rg, (/ p_cache_size /),"cyl_modes/os-spec-cyl-modes.f03",886)
    call alloc(tmp_q, (/ p_cache_size /),"cyl_modes/os-spec-cyl-modes.f03",887)
    call alloc(tmp_dxi, (/ p_x_dim, p_cache_size /),"cyl_modes/os-spec-cyl-modes.f03",888)
    call alloc(tmp_ix, (/ p_x_dim, p_cache_size /),"cyl_modes/os-spec-cyl-modes.f03",889)
  endif

end subroutine init_tmp_buf_current_spe_bnd_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Deposit species charge on cylindrical mode grids - mostly copied from deposit_quant in
! os-sec-diagnostics.f90.
!-----------------------------------------------------------------------------------------
subroutine deposit_quant_cyl_modes( spec, grid, no_co, charge_cyl_m, charge, quant, send_msg, recv_msg )

  implicit none

  class( t_species_cyl_modes ), intent(in) :: spec
  class( t_grid ),intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_cyl_modes ), intent(inout) :: charge_cyl_m
  type( t_vdf ), intent(inout) :: charge
  integer, intent(in) :: quant
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! local variables
  integer, parameter :: p_part_block = 4096
  real(p_k_part), dimension(p_part_block) :: q, coeff_re, coeff_im
  integer :: i, i1, i2

  integer, dimension (2, p_x_dim) :: gc_num
  integer, dimension(p_x_dim) :: move_num
  integer :: lquant, n_modes, mode
  type(t_vdf), pointer :: vdf

  n_modes = spec%n_cyl_modes

  ! The number of guard cells required depends on the interpolation level
  gc_num(p_lower,:) = spec%interpolation
  gc_num(p_upper,:) = spec%interpolation + 1

  ! note that this automatically sets the vdf value to 0.0
  call charge % new( p_x_dim, 1, grid%my_nx(3,:), gc_num, spec%dx, .true. )

  ! create modes
  call setup( charge_cyl_m, n_modes, charge )

  ! deposit given quantity
  if ( quant == p_norm ) then
    lquant = p_charge
  else
    lquant = quant
  endif

  do i1 = 1, spec%num_par, p_part_block
    i2 = i1 + p_part_block - 1
    if ( i2 > spec%num_par ) i2 = spec%num_par
    call spec % get_quant( i1, i2, lquant, q )

    call spec % deposit_density( charge, i1, i2, q )
    do mode = 1, n_modes
      ! weigh the charges according to their modes
      call get_coef( coeff_re, coeff_im, spec%x(3:4,i1:i2), (i2-i1+1), mode )
      call spec % deposit_density( charge_cyl_m%pf_re(mode), i1, i2, q(:(i2-i1+1))*coeff_re(:(i2-i1+1)))
      call spec % deposit_density( charge_cyl_m%pf_im(mode), i1, i2, q(:(i2-i1+1))*coeff_im(:(i2-i1+1)))
    enddo

  enddo

  move_num = 0
  call update_boundary( charge, p_vdf_add, no_co, send_msg, recv_msg, move_num )
  call update_boundary( charge_cyl_m, p_vdf_add, no_co, send_msg, recv_msg )

  ! normalize charge for cylindrical coordinates
  call norm_charge_cyl_m( charge, grid%my_nx( 1, p_r_dim ), real( spec%dx( p_r_dim ), p_k_fld ), 0 )
  do mode = 1, n_modes
    call norm_charge_cyl_m( charge_cyl_m%pf_re(mode), grid%my_nx( 1, p_r_dim ), real( spec%dx( p_r_dim ), p_k_fld ), mode )
    call norm_charge_cyl_m( charge_cyl_m%pf_im(mode), grid%my_nx( 1, p_r_dim ), real( spec%dx( p_r_dim ), p_k_fld ), mode )
  enddo

  ! get p_norm quantity requires special treatment
  if ( quant == p_norm ) then
    do mode = 0, n_modes
      do i = p_real, p_imag
        if (i==p_real) then
          vdf => charge_cyl_m%pf_re(mode)
        else
          if (mode==0) cycle
          vdf => charge_cyl_m%pf_im(mode)
        endif
        do i2 = lbound( vdf%f2, 3 ), ubound( vdf%f2, 3 )
          do i1 = lbound( vdf%f2, 2 ), ubound( vdf%f2, 2 )
            if ( vdf%f2( 1, i1, i2 ) /= 0.0 ) then
              vdf%f2( 1, i1, i2 ) = abs( 1.0 / vdf%f2( 1, i1, i2 ) )
            else
              vdf%f2( 1, i1, i2 ) = 1.0
            endif
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine deposit_quant_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Construct the current density 4 vector; only construct
! components required for boosted diags.
!-----------------------------------------------------------------------------------------
subroutine get_current_4vector_cyl_modes( spec, grid, no_co, current_cyl_m, current, send_msg, recv_msg )

  implicit none

  class( t_species_cyl_modes ), intent(in) :: spec
  class( t_grid ),intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_cyl_modes ), intent(inout) :: current_cyl_m
  type( t_vdf ), intent(inout) :: current
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer, dimension (2, p_x_dim) :: gc_num
  integer :: quant, n, n_modes
  integer, dimension(4), parameter :: quants = (/ p_charge_cyl_m, p_j1_cyl_m, p_j3_cyl_m, p_j2_cyl_m /)
  type(t_vdf) :: tmp_vdf
  type(t_cyl_modes) :: tmp_cyl_m
  integer :: component

  n_modes = spec%n_cyl_modes

  gc_num(p_lower,:) = spec%interpolation
  gc_num(p_upper,:) = spec%interpolation + 1

  call current % new( p_x_dim, 4, grid%my_nx(3,:), gc_num , spec%dx, .true.)
  call setup( current_cyl_m, n_modes, current )

  do n = 1, spec%diag%boosted_diag%num_diag_components_to_buffer

    component = spec%diag%boosted_diag%diag_components_to_buffer( n )

    quant = quants( component )

    call deposit_quant_cyl_modes( spec, grid, no_co, tmp_cyl_m, tmp_vdf, quant, send_msg, recv_msg )

    call copy_component( current_cyl_m, tmp_cyl_m, 1, component )

    call cleanup( tmp_cyl_m )
    call tmp_vdf % cleanup()

  enddo

end subroutine get_current_4vector_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Polymorphic diagnostic functions end
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function ntrim(x)
!-----------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
! This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
! routine that would break for x = nearest( 0.5, -1.0 )
!-----------------------------------------------------------------------------------------
  implicit none

  real(p_k_part), intent(in) :: x
  integer :: ntrim, a, b

  if ( x < -.5 ) then
    a = -1
  else
    a = 0
  endif

  if ( x >= .5 ) then
    b = +1
  else
    b = 0
  endif

  ntrim = a+b

end function ntrim
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! auxiliary function for interpolating field values of particles
!-----------------------------------------------------------------------------------------
subroutine interpolate_cyl_modes(emf, nptotal, interpolation, track_set)

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  class( t_emf_cyl_modes ), intent(in) :: emf
  integer, intent(in) :: nptotal
  integer, intent(in) :: interpolation

  real(p_k_part), dimension(p_cache_size) :: cos_th, sin_th
  real(p_k_part), dimension(p_cache_size) :: coeff_re, coeff_im
  real(p_k_part), dimension(0:ubound( emf%e_cyl_m%pf_re, 1 ),p_p_dim,p_cache_size) :: ep_re, ep_im, bp_re, bp_im

  integer :: n_modes, i, mode, np, ptrcur, k, stride


  ! n_modes
  n_modes = emf % n_cyl_modes
  stride = n_modes + 2 ! For total field and m=0

  do ptrcur = 1, nptotal, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > nptotal ) then
      np = nptotal - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! start with mode 0
    mode = 0
    call get_emf( emf, bp_re(mode,:,:), ep_re(mode,:,:), track_set%ipos_buffer(1:2,ptrcur:), &
                  track_set%pos_buffer(1:2,ptrcur:), np, interpolation )

    ! get cos(theta) and sin(theta) for particles being processed
    call get_coef( cos_th, sin_th, track_set%pos_buffer(3:4,ptrcur:), np, 1 )

    ! the fields need to be converted to cartesian coordinates
    do i = 1, np
      k = 1
      if(track_set%ifdmp_tracks_efl(1)) then
        ! all 6 components follow the following strategy
        ! first set the buffer value corresponding to mode 0
        track_set%field_buffer(i,k+mode+1) = ep_re(mode,1,i)
        ! next set the total field value, which will be added to in the next loop
        track_set%field_buffer(i,k) = track_set%field_buffer(i,k+mode+1)
        k = k + stride
      endif

      if(track_set%ifdmp_tracks_efl(2)) then
        track_set%field_buffer(i,k+mode+1) = ep_re(mode,2,i) * cos_th( i ) &
                                           - ep_re(mode,3,i) * sin_th( i )
        track_set%field_buffer(i,k) = track_set%field_buffer(i,k+mode+1)
        k = k + stride
      endif

      if(track_set%ifdmp_tracks_efl(3)) then
        track_set%field_buffer(i,k+mode+1) = ep_re(mode,3,i) * cos_th( i ) &
                                           + ep_re(mode,2,i) * sin_th( i )
        track_set%field_buffer(i,k) = track_set%field_buffer(i,k+mode+1)
        k = k + stride
      endif

      if(track_set%ifdmp_tracks_bfl(1)) then
        track_set%field_buffer(i,k+mode+1) = bp_re(mode,1,i)
        track_set%field_buffer(i,k) = track_set%field_buffer(i,k+mode+1)
        k = k + stride
      endif

      if(track_set%ifdmp_tracks_bfl(2)) then
        track_set%field_buffer(i,k+mode+1) = bp_re(mode,2,i) * cos_th( i ) &
                                           - bp_re(mode,3,i) * sin_th( i )
        track_set%field_buffer(i,k) = track_set%field_buffer(i,k+mode+1)
        k = k + stride
      endif

      if(track_set%ifdmp_tracks_bfl(3)) then
        track_set%field_buffer(i,k+mode+1) = bp_re(mode,3,i) * cos_th( i ) &
                                           + bp_re(mode,2,i) * sin_th( i )
        track_set%field_buffer(i,k) = track_set%field_buffer(i,k+mode+1)
      endif

    enddo

    do mode = 1, n_modes
      call get_emf_cyl_modes( emf, mode, bp_re(mode,:,:), ep_re(mode,:,:), &
                              bp_im(mode,:,:), ep_im(mode,:,:), &
                              track_set%ipos_buffer(:,ptrcur:), &
                              track_set%pos_buffer(:,ptrcur:), np, &
                              interpolation )
      call get_coef( coeff_re, coeff_im, track_set%pos_buffer(3:4,ptrcur:), np, mode )

      do i = 1, np
        k = 1
        if(track_set%ifdmp_tracks_efl(1)) then
          ! all 6 components follow the following strategy
          ! first set the buffer value corresponding to the appropriate mode
          track_set%field_buffer(i,k+mode+1) = ep_re(mode,1,i) * coeff_re( i ) &
                                              + ep_im(mode,1,i) * coeff_im( i )
          ! next add this value to the total field value
          track_set%field_buffer(i,k) = track_set%field_buffer(i,k) &
                                      + track_set%field_buffer(i,k+mode+1)
          k = k + stride
        endif

        if(track_set%ifdmp_tracks_efl(2)) then
          track_set%field_buffer(i,k+mode+1) = ep_re(mode,2,i) * coeff_re( i ) * cos_th( i ) &
                                              - ep_re(mode,3,i) * coeff_re( i ) * sin_th( i ) &
                                              + ep_im(mode,2,i) * coeff_im( i ) * cos_th( i ) &
                                              - ep_im(mode,3,i) * coeff_im( i ) * sin_th( i )
          track_set%field_buffer(i,k) = track_set%field_buffer(i,k) &
                                      + track_set%field_buffer(i,k+mode+1)
          k = k + stride
        endif

        if(track_set%ifdmp_tracks_efl(3)) then
          track_set%field_buffer(i,k+mode+1) = ep_re(mode,3,i) * coeff_re( i ) * cos_th( i ) &
                                              + ep_re(mode,2,i) * coeff_re( i ) * sin_th( i ) &
                                              + ep_im(mode,3,i) * coeff_im( i ) * cos_th( i ) &
                                              + ep_im(mode,2,i) * coeff_im( i ) * sin_th( i )
          track_set%field_buffer(i,k) = track_set%field_buffer(i,k) &
                                      + track_set%field_buffer(i,k+mode+1)
          k = k + stride
        endif

        if(track_set%ifdmp_tracks_bfl(1)) then
          track_set%field_buffer(i,k+mode+1) = bp_re(mode,1,i) * coeff_re( i ) &
                                              + bp_im(mode,1,i) * coeff_im( i )
          track_set%field_buffer(i,k) = track_set%field_buffer(i,k) &
                                      + track_set%field_buffer(i,k+mode+1)
          k = k + stride
        endif

        if(track_set%ifdmp_tracks_bfl(2)) then
          track_set%field_buffer(i,k+mode+1) = bp_re(mode,2,i) * coeff_re( i ) * cos_th( i ) &
                                              - bp_re(mode,3,i) * coeff_re( i ) * sin_th( i ) &
                                              + bp_im(mode,2,i) * coeff_im( i ) * cos_th( i ) &
                                              - bp_im(mode,3,i) * coeff_im( i ) * sin_th( i )
          track_set%field_buffer(i,k) = track_set%field_buffer(i,k) &
                                      + track_set%field_buffer(i,k+mode+1)
          k = k + stride
        endif

        if(track_set%ifdmp_tracks_bfl(3)) then
          track_set%field_buffer(i,k+mode+1) = bp_re(mode,3,i) * coeff_re( i ) * cos_th( i ) &
                                              + bp_re(mode,2,i) * coeff_re( i ) * sin_th( i ) &
                                              + bp_im(mode,3,i) * coeff_im( i ) * cos_th( i ) &
                                              + bp_im(mode,2,i) * coeff_im( i ) * sin_th( i )
          track_set%field_buffer(i,k) = track_set%field_buffer(i,k) &
                                      + track_set%field_buffer(i,k+mode+1)
        endif
      enddo
    enddo
  enddo

end subroutine interpolate_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Add current values to track data for cyl_modes
!-----------------------------------------------------------------------------------------
subroutine add_track_set_data_cyl_modes( track_set, spec, no_co, emf, n, t, send_msg, recv_msg )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  class( t_species_cyl_modes ), intent(inout) :: spec
  class( t_node_conf ), intent(in) :: no_co
  class( t_emf_cyl_modes ), intent(inout) :: emf
  integer, intent(in) :: n
  real( p_double ), intent(in) :: t
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, idx, n_x_dim, k

  type( t_vdf ), pointer :: psi

  psi => null()

  n_x_dim = spec % get_n_x_dims()

  if ( track_set%ifdmp_tracks_psi ) then
    call get_psi( emf, n, no_co, psi, send_msg, recv_msg )
  endif

  if ( track_set%npresent > 0 ) then

    if ( track_set%nfields > 0 ) then

      if ( track_set%npresent > track_set%size_buffer ) then
        if ( track_set%size_buffer > 0 ) then
          call freemem(track_set%field_buffer,"cyl_modes/os-spec-cyl-modes.f03",1274)
          call freemem(track_set%pos_buffer,"cyl_modes/os-spec-cyl-modes.f03",1275)
          call freemem(track_set%ipos_buffer,"cyl_modes/os-spec-cyl-modes.f03",1276)
        endif

        ! grow buffer in sizes of 1024 tracks
        track_set%size_buffer = ceiling( track_set%npresent / 1024.0 ) * 1024

        call alloc(track_set%pos_buffer, (/ n_x_dim, track_set%size_buffer /),"cyl_modes/os-spec-cyl-modes.f03",1282)
        call alloc(track_set%ipos_buffer, (/ p_x_dim, track_set%size_buffer /),"cyl_modes/os-spec-cyl-modes.f03",1283)

        ! If interpolating any fields also grow field_buffer
        if ( track_set%nfields > 0 ) &
          call alloc(track_set%field_buffer, (/ track_set%size_buffer, track_set%nfields /),"cyl_modes/os-spec-cyl-modes.f03",1287)

      endif

      ! Get positions of particles being tracked
      do i = 1, track_set%npresent

        idx = track_set%tracks( track_set%present(i) )%part_idx
        track_set% ipos_buffer(1:p_x_dim, i) = spec%ix( 1:p_x_dim, idx )
        track_set% pos_buffer(1:n_x_dim, i) = spec%x( 1:n_x_dim, idx )
      enddo

      call interpolate_cyl_modes( emf, track_set%npresent, spec%interpolation, track_set )

      if ( track_set%ifdmp_tracks_psi ) then
        k = track_set%nfields
        call interpolate( psi, 1, 0, track_set%ipos_buffer, track_set%pos_buffer, track_set%npresent, &
                          spec%interpolation, track_set%field_buffer( :, k ) )
      endif

      ! Add track data including interpolated E and B fields, and Psi diagnostic
      do i = 1, track_set%npresent
        call add_single_track( track_set%tracks( track_set%present(i) ), &
                                spec, n, t, track_set%field_buffer( i, : ) )
      enddo

    else

      ! Add track data
      do i = 1, track_set%npresent
        call add_single_track( track_set%tracks( track_set%present(i) ), &
                                spec, n, t )
      enddo

    endif

  endif

end subroutine add_track_set_data_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Set the momenta of particles coming from thermal bath boundary.
!-----------------------------------------------------------------------------------------
subroutine thermal_bath_cyl_modes( this, i_dim, bnd_idx, u, cos_sin, npar, tmp_ix )

  implicit none

  class ( t_spe_bound ),intent(in) :: this
  integer, intent(in) :: i_dim, bnd_idx
  real(p_k_part), intent(inout), dimension(:,:) :: u
  real(p_k_part), intent(in), dimension(:,:) :: cos_sin
  integer, intent(in) :: npar
  integer, pointer, dimension(:,:) :: tmp_ix

  integer :: k
  real(p_k_part) :: usign, u1, u2, u3, u_tmp, upar, uperp, theta
  real(p_k_part) :: a, b, gth

  real(p_double), parameter :: pi = 3.14159265358979323846264_p_double

  ! Make sure the component normal to the boundary points inward
  select case( bnd_idx )
  case( p_lower )
    usign = 1
  case( p_upper )
    usign = -1
  case default
    usign = 0
  end select

  u1 = this%uth_bnd( 1, bnd_idx, i_dim )
  u2 = this%uth_bnd( 2, bnd_idx, i_dim )
  u3 = this%uth_bnd( 3, bnd_idx, i_dim )

  select case ( this%thermal_type(bnd_idx, i_dim) )

  case( p_thermal )

    u_tmp = usign * this%uth_bnd(i_dim,bnd_idx,i_dim)

    do k = 1, npar
      u(1,k) = u1 * real( rng % genrand_gaussian( tmp_ix(:, k)), p_k_part )
      u(2,k) = u2 * real( rng % genrand_gaussian( tmp_ix(:, k)), p_k_part )
      u(3,k) = u3 * real( rng % genrand_gaussian( tmp_ix(:, k)), p_k_part )

      ! This is required to maintain a null flux across the boundary
      ! Higer momentum particles are more likely to exit the boundary, so re-injecting
      ! from a normal gaussian distribution will generally yield lower energy particles
      ! and cool down the plasma
      u(i_dim, k) = u_tmp * real(rng % genrand_half_max(tmp_ix(:, k)), p_k_part)
    enddo

  case ( p_waterbag )
    do k = 1, npar
      u(1,k) = u1 * real( rng % genrand_real3( tmp_ix(:, k)) - 0.5, p_k_part )
      u(2,k) = u2 * real( rng % genrand_real3( tmp_ix(:, k)) - 0.5, p_k_part )
      u(3,k) = u3 * real( rng % genrand_real3( tmp_ix(:, k)) - 0.5, p_k_part )

      u(i_dim, k) = usign * abs( u(i_dim, k) )
    enddo

  case ( p_random_dir )

    u_tmp = sqrt( u1**2 + u2**2 + u3**2 )

    do k = 1, npar
      ! get p3
      upar = real(2 * rng % genrand_real3( tmp_ix(:, k)) - 1.0_p_double, p_k_part)

      ! get p2 and p1
      uperp = sqrt(1 - upar)*u_tmp
      theta = real(2 * pi * rng % genrand_real3( tmp_ix(:, k)), p_k_part)
      u(1,k) = cos(theta) * uperp
      u(2,k) = sin(theta) * uperp
      u(3,k) = upar * u_tmp

      u(i_dim, k) = usign * abs( u(i_dim, k) )
    enddo

  end select

  if ( i_dim == 2 ) then
    ! u2 and u3 correspond to u_r and u_theta, respectively. Need to rotate to cartesian
    do k = 1, npar
      u2 = u(2,k)
      u3 = u(3,k)
      u(2,k) = u2*cos_sin(1,k) - u3*cos_sin(2,k)
      u(3,k) = u2*cos_sin(2,k) + u3*cos_sin(1,k)
    enddo
  endif

  ! Add fluid momenta if required
  u1 = this%ufl_bnd( 1, bnd_idx, i_dim )
  u2 = this%ufl_bnd( 2, bnd_idx, i_dim )
  u3 = this%ufl_bnd( 3, bnd_idx, i_dim )

  if ( u1 /= 0.0 .and. u2 /= 0.0 .and. u3 /= 0.0 ) then
    a = 1.0_p_k_part / (1.0_p_k_part + sqrt( 1 + u1**2 + u2**2 + u3**2 ) )

    do k = 1, npar
      gth = sqrt( 1 + u(1,k)**2 + u(2,k)**2 + u(3,k)**2 )
      b = gth + a * (u1*u(1,k) + u2*u(2,k) + u3*u(3,k))
      u(1,k) = u(1,k) + b * u1
      u(2,k) = u(2,k) + b * u2
      u(3,k) = u(3,k) + b * u3
    enddo
  endif

end subroutine thermal_bath_cyl_modes
!-----------------------------------------------------------------------------------------


end module m_species_cyl_modes


!-----------------------------------------------------------------------------------------
! Populate bp, ep with the cartesian components of the magnetic and electric fields
! at the specified particle positions
!-----------------------------------------------------------------------------------------
subroutine get_emf_spec_cyl_modes( this, emf, bp, ep, np, ptrcur, t )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_emf_define, only : t_emf
  use m_emf_cyl_modes, only : t_emf_cyl_modes
  use m_emf_interpolate, only : get_emf
  use m_parameters
  use m_species_current_cyl_modes, only : get_coef
  use m_emf_interpolate_cm, only : get_emf_cyl_modes

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  class( t_emf ), intent(in), target :: emf
  real(p_k_part), dimension(:,:), intent(out) :: bp, ep
  integer, intent(in) :: np, ptrcur
  real(p_double), intent(in), optional :: t

  integer :: n_modes, mode, i
  real(p_k_part), dimension(0:this%n_cyl_modes,p_p_dim,p_cache_size) :: ep_im, bp_re, ep_re, bp_im
  real(p_k_part), dimension(p_cache_size) :: cos_th, sin_th, coeff_re, coeff_im

  n_modes = this%n_cyl_modes

  ! Interpolate fields adding up contribution from all modes
  ! get mode 0 (only real part meaningfull)
  call get_emf( emf, bp_re(0,:,:), ep_re(0,:,:), this%ix(1:2,ptrcur:), this%x(1:2,ptrcur:), np, &
                this%interpolation )

  ! get cos(theta) and sin(theta) for particles being processed
  call get_coef( cos_th, sin_th, this%x(3:4,ptrcur:), np, 1 )

  ! the fields need to be converted to cartesian coordinates
  do i = 1, np
    ep(1,i) = ep_re(0,1,i)

    ep(2,i) = ep_re(0,2,i) * cos_th( i ) &
            - ep_re(0,3,i) * sin_th( i )

    ep(3,i) = ep_re(0,3,i) * cos_th( i ) &
            + ep_re(0,2,i) * sin_th( i )

    bp(1,i) = bp_re(0,1,i)

    bp(2,i) = bp_re(0,2,i) * cos_th( i ) &
            - bp_re(0,3,i) * sin_th( i )

    bp(3,i) = bp_re(0,3,i) * cos_th( i ) &
            + bp_re(0,2,i) * sin_th( i )
  enddo

  ! Add contribution of all remaining modes
  select type( emf )
  class is( t_emf_cyl_modes )
    do mode = 1, n_modes
      call get_emf_cyl_modes( emf, mode, bp_re(mode,:,:), ep_re(mode,:,:), &
                              bp_im(mode,:,:), ep_im(mode,:,:), &
                              this%ix(:,ptrcur:), this%x(:,ptrcur:), np, &
                              this%interpolation )
      call get_coef( coeff_re, coeff_im, this%x(3:4,ptrcur:), np, mode )

      do i = 1, np
        ep(1,i) = ep(1,i) + ep_re(mode,1,i) * coeff_re( i ) &
                          + ep_im(mode,1,i) * coeff_im( i )
        ep(2,i) = ep(2,i) + ep_re(mode,2,i) * coeff_re( i ) * cos_th( i ) &
                          - ep_re(mode,3,i) * coeff_re( i ) * sin_th( i ) &
                          + ep_im(mode,2,i) * coeff_im( i ) * cos_th( i ) &
                          - ep_im(mode,3,i) * coeff_im( i ) * sin_th( i )
        ep(3,i) = ep(3,i) + ep_re(mode,3,i) * coeff_re( i ) * cos_th( i ) &
                          + ep_re(mode,2,i) * coeff_re( i ) * sin_th( i ) &
                          + ep_im(mode,3,i) * coeff_im( i ) * cos_th( i ) &
                          + ep_im(mode,2,i) * coeff_im( i ) * sin_th( i )

        bp(1,i) = bp(1,i) + bp_re(mode,1,i) * coeff_re( i ) &
                          + bp_im(mode,1,i) * coeff_im( i )
        bp(2,i) = bp(2,i) + bp_re(mode,2,i) * coeff_re( i ) * cos_th( i ) &
                          - bp_re(mode,3,i) * coeff_re( i ) * sin_th( i ) &
                          + bp_im(mode,2,i) * coeff_im( i ) * cos_th( i ) &
                          - bp_im(mode,3,i) * coeff_im( i ) * sin_th( i )
        bp(3,i) = bp(3,i) + bp_re(mode,3,i) * coeff_re( i ) * cos_th( i ) &
                          + bp_re(mode,2,i) * coeff_re( i ) * sin_th( i ) &
                          + bp_im(mode,3,i) * coeff_im( i ) * cos_th( i ) &
                          + bp_im(mode,2,i) * coeff_im( i ) * sin_th( i )
      enddo
    enddo
  end select

end subroutine get_emf_spec_cyl_modes

!-----------------------------------------------------------------------------------------
! Return the number of spatial dimensions on the x buffer
! - For the quasi-3D calculations this will be 4 (z,r,x,y)
!-----------------------------------------------------------------------------------------
pure function get_n_x_dims_cyl_modes( this )

  use m_species_cyl_modes_define, only : t_species_cyl_modes

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer :: get_n_x_dims_cyl_modes

  get_n_x_dims_cyl_modes = 4

end function get_n_x_dims_cyl_modes


subroutine allocate_objs_spec_cyl_modes( this )
!-----------------------------------------------------------------------------------------
! Allocate any objects contained within the species (t_species_cyl_modes) object
!-----------------------------------------------------------------------------------------
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_species_cyl_modes, only : t_diag_species_cyl_modes, t_spe_bound_cyl_modes

  implicit none

  class( t_species_cyl_modes ), intent(inout) :: this

  ! Allocate t_diag_species_cyl_modes object class
  if ( .not. associated( this%diag ) ) then
    allocate( t_diag_species_cyl_modes :: this%diag )
  endif

  ! Allocate t_spe_bound_cyl_modes object class
  if ( .not. associated( this%bnd_con ) ) then
    allocate( t_spe_bound_cyl_modes :: this%bnd_con )
  endif

  ! call superclass allocation
  call this % t_species % allocate_objs( )

end subroutine


subroutine init_particle_buffer_cyl_modes( this, num_par_req )
!-----------------------------------------------------------------------------------------
! This is where we allocate the particle position buffers in the species object.
! In the quasi-3D geometry there are 4 coordinates in species%x
! (z(cell coord),r(cell coord),x(box coord),y(box coord))
!-----------------------------------------------------------------------------------------

  use m_system
  use m_parameters
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_species_memory

  implicit none

  class(t_species_cyl_modes), intent(inout) :: this
  integer, intent(in) :: num_par_req

  ! The buffer size must always be a multiple of vector width in size because of SIMD code
  this%num_par_max = (( num_par_req + p_vecwidth - 1 ) / p_vecwidth) * p_vecwidth

  ! Under some compilers / configurations (e.g. gfortran 4.9.1, OS X, single precision)
  ! if the total q buffer size is below 1Kb it may not be allocated to a 32bit boundary
  ! which is required by AVX code
  if ( this%num_par_max < 256 ) this%num_par_max = 256

  ! setup position buffer
  ! This is the main difference between the species object in quasi-3D and in
  ! 2D cylindrical : the number of dimensions allocated in species%x
  call freemem(this%x,"cyl_modes/os-spec-cyl-modes.f03",1609)
  call alloc(this%x, (/ p_x_dim + 2, this%num_par_max /),"cyl_modes/os-spec-cyl-modes.f03",1610)

  ! initialize particle cell information
  call freemem(this%ix,"cyl_modes/os-spec-cyl-modes.f03",1613)
  call alloc(this%ix, (/ p_x_dim, this%num_par_max /),"cyl_modes/os-spec-cyl-modes.f03",1614)

  ! setup momenta buffer
  call freemem(this%p,"cyl_modes/os-spec-cyl-modes.f03",1617)
  call alloc(this%p, (/ p_p_dim, this%num_par_max /),"cyl_modes/os-spec-cyl-modes.f03",1618)
# 1628 "cyl_modes/os-spec-cyl-modes.f03"
  ! setup particle charge buffer
  call freemem(this%q,"cyl_modes/os-spec-cyl-modes.f03",1629)
  call alloc(this%q, (/ this%num_par_max /),"cyl_modes/os-spec-cyl-modes.f03",1630)

  ! initialize tracking data if necessary
  if ( this%add_tag ) then
    call freemem(this%tag,"cyl_modes/os-spec-cyl-modes.f03",1634)
    call alloc(this%tag, (/ 2, this%num_par_max /),"cyl_modes/os-spec-cyl-modes.f03",1635)
  endif

end subroutine


!-----------------------------------------------------------------------------------------
! Read information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_species_cyl_modes( this, input_file, def_name, periodic, if_move, grid, &
                             dt, read_prof, sim_options )

  use m_system
  use m_parameters
  use m_species_define, only: p_max_spname_len
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_input_file, only : t_input_file, disp_out, get_namelist
  use m_piston, only : read_nml_piston, alloc, check_piston
  use m_species_boundary, only : read_nml, type
  use m_species_udist, only : read_nml



  use m_grid_define, only : t_grid
  use m_profile_cyl_modes, only : new_source_cyl_modes
  use m_grid_cyl_modes, only : t_grid_cyl_modes
  use m_species_cyl_modes, only : t_diag_species_cyl_modes

  implicit none

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  character(len = *), intent(in) :: def_name
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid

  type(t_options), intent(in) :: sim_options

  real(p_double), intent(in) :: dt
  logical, intent(in) :: read_prof

  character(len = p_max_spname_len) :: name
  character(len = 10 ) :: push_type

  logical :: if_use_delayed_injection
  real(p_double) :: delayed_injection_time

  integer :: num_par_max
  integer :: n_sort
  real(p_k_part) :: rqm

  real(p_k_part) :: q_real



  integer, dimension(p_x_dim) :: num_par_x
  integer(p_int64), dimension(p_x_dim) :: tot_par_x
  integer :: num_par_theta ! main difference in quasi-3D
  logical :: rand_theta, shift_cell_spokes

  integer :: num_pistons

  logical :: add_tag
  logical :: free_stream

  real(p_double) :: push_start_time

  logical :: if_collide
  logical :: if_like_collide

  logical :: init_fields

  character(len = 16) :: init_type

  real(p_double) :: iter_tol
  logical :: rad_react

  ! gravity
  real(p_k_part), dimension(p_f_dim) :: gravity ! gravitational acceleration
# 1728 "cyl_modes/os-spec-cyl-modes.f03"
  namelist /nl_species/ name, num_par_max, n_sort, rqm, q_real, &
                        num_par_x, tot_par_x, num_par_theta, rand_theta, push_type, &
                        shift_cell_spokes, push_start_time, num_pistons, &
                        add_tag, free_stream, init_fields, &
                        if_use_delayed_injection, delayed_injection_time, &
                        if_collide, if_like_collide, init_type, iter_tol, rad_react, &
                        gravity


  integer :: i, ierr, piston_id

  ! Copy in relevant simulation parameters
  this%dt = dt
  this%g_nx( 1:p_x_dim ) = grid%g_nx( 1:p_x_dim )
  this%coordinates = grid%coordinates

  push_type = "standard"

  ! Initialization type
  init_type = "standard"

  num_par_max = 0
  rqm = 0
  q_real = 0
  num_par_x = -1
  tot_par_x = -1
  num_par_theta = -1
  rand_theta = .false.
  shift_cell_spokes = .true.
  n_sort = 25

  push_start_time = -1.0_p_double
  name = trim(adjustl(def_name))

  num_pistons = 0

  add_tag = .false.

  free_stream = .false.




  if_collide = .false.
  if_like_collide = .false.
  init_fields = .false.

  if_use_delayed_injection = .false.
  delayed_injection_time = 0

  ! convergence threshold of root finding routines in exact pusher
  iter_tol = 1.0d-3

  ! if include radiation reaction (only valid for exact pusher)
  rad_react = .false.

  ! gravity
  gravity = 0.0_p_k_fld

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_species", ierr )
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading species parameters"
      else
        write(0,*) "Error: species parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif


  read (input_file%nml_text, nml = nl_species, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error reading species parameters"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  if (disp_out(input_file)) then
    if (mpi_node()==0) print *,"   Species name : ", trim(name)
  endif

  ! -------------- Allocate any sub-objects needed --------------------
  call this%allocate_objs()

  ! Fill in n_cyl_modes values
  select type( grid_cyl => grid )
  class is( t_grid_cyl_modes )
    this%n_cyl_modes = grid_cyl%n_cyl_modes
  end select

  select type( diag => this%diag )
  class is ( t_diag_species_cyl_modes )
    diag%n_cyl_modes = this%n_cyl_modes
  end select

  this%name = name
  this%free_stream = free_stream

  this%num_par_max = num_par_max
  this%rqm = rqm

  ! Validate num_par_x and tot_par_x
  do i = 1, p_x_dim
    if ( (num_par_x(i) <= 0) .and. (tot_par_x(i) <= 0) ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "   Error reading species parameters"
        write(0,*) "   Number of particles per cell must be >= 1 in all directions:"
        write(0,"(A,I0,A,I0,A,I0,A,I0)") "    -> num_par_x(",i,") = ", &
                 num_par_x(i), ", tot_par_x(",i,") = ", tot_par_x(i)
        write(0,*) "   aborting..."
      endif
      stop
    endif

    ! num_par_x overrides tot_par_x
    if ( num_par_x(i) <= 0 ) then
      this%num_par_x(i) = tot_par_x(i) / this % g_nx(i)
      this%tot_par_x(i) = tot_par_x(i)
    else
      this%num_par_x(i) = num_par_x(i)
      this%tot_par_x(i) = num_par_x(i) * this % g_nx(i)
    endif
  enddo

  if ( num_par_theta <= 0 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error reading species parameters"
      write(0,*) "   Number of particles per in theta must be >= 1:"
      write(0,"(A,I0)") "    -> num_par_theta = ", num_par_theta
      write(0,*) "   aborting..."
    endif
    stop
  endif

  if ( num_par_theta < 8*this%n_cyl_modes ) then
    if ( mpi_node() == 0 ) then
      write(0,"(A)") " (*WARNING*) num_par_theta is recommended to be"
      write(0,"(A,I0)") " (*WARNING*) at least 8*n_cyl_modes = ", 8*this%n_cyl_modes
      write(0,"(A,I0)") " (*WARNING*) current value = ", num_par_theta
    endif
  endif

  this%num_par_theta = num_par_theta
  this%rand_theta = rand_theta
  this%shift_cell_spokes= shift_cell_spokes
  this%n_sort = n_sort

  ! Pusher type
  select case ( trim( push_type ) )
    case('standard')
       this%push_type = p_std
    case('vay')
      this%push_type = p_vay
    case('fullrot')
      this%push_type = p_fullrot
    case('euler')
      this%push_type = p_euler
    case('cond_vay')
      this%push_type = p_cond_vay
    case('cary')
      this%push_type = p_cary
    case('exact', 'analytic')
      this%push_type = p_exact
      this%iter_tol = iter_tol
    case('exact-rr', 'analytic-rr')
      this%push_type = p_exact_rr
      this%iter_tol = iter_tol
    case('gravity')
      this%push_type = p_gravity
      this%gravity = gravity
      if (disp_out(input_file)) then
        if (mpi_node()==0) print *,"  ::::::::::::::::::::::::::"
        if (mpi_node()==0) print *,"  Using gravity"
        if (mpi_node()==0) print *,"  gravity = ", this%gravity
        if (mpi_node()==0) print *,"  ::::::::::::::::::::::::::"
      endif

    case default
      if ( mpi_node() == 0) then
        write(0,*) "   Invalid push type:'", trim(push_type), "'"
        write(0,*) "   Allowed values are 'standard', 'vay', 'fullrot', 'euler', 'cond_vay',"
        write(0,*) "   'cary',' 'gravity', 'exact'/'analytic' and 'exact-rr'/'analytic-rr'"
        write(0,*) "   Error reading species parameters"
        write(0,*) "   aborting..."
      endif
      stop
  end select

  ! parameters for RR
  this%rad_react = rad_react
  if ( this%rad_react ) then
    if ( sim_options%omega_p0 <= 0. ) then
      if ( mpi_node() == 0) then
        write(0,*) "   Error: omega_p0 missing or invalid"
        write(0,*) "   When choosing pusher with radiation reaction you must"
        write(0,*) "   also set omega_p0 in the simulation section of the input file"
      endif
      stop
    else
      ! the leading coefficient is a characteristic time 2*e^2/(3*m_e*c^3)
      ! (e = elementary charge, m_e = electron static mass) in cgs unit.
      ! see Eq. (16.3) in Chapter 16 of Jackson
      this%k_rr = 6.266424752633625d-24 * sim_options%omega_p0
    endif

    if ( this%push_type /= p_std .and. this%push_type /= p_exact .and. &
         this%push_type /= p_gravity ) then
      if ( mpi_node() == 0) then
        write(0,*) "   Error: radiation reaction correction is only implemented for"
        write(0,*) "   'standard' and 'exact'/'analytic' push types"
      endif
      stop
    endif
  endif

  ! Free streaming is not implemented in simd code
  if ( free_stream .and. this%push_type == p_simd ) then
    this%push_type = p_std
  endif

  this%push_start_time = push_start_time

  this%add_tag = add_tag

  this%num_pistons = num_pistons

  ! this is also used by some people even without collisions
  ! q_real and rqm having different signs is inconsistent, give rqm precedence
  this%q_real = sign(q_real, rqm)

  if ( (this%push_type == p_exact .and. this%rad_react) .or. &
    this%push_type == p_exact_rr ) then
    if ( this%q_real == 0. ) then
      if ( mpi_node() == 0) then
        write(0,*) "   Error: q_real missing or invalid"
        write(0,*) "   When choosing push_type = 'exact' with radiation reaction"
        write(0,*) "   or push_type = 'exact-rr', you must also set q_real"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  endif



  this%if_collide = if_collide
  this%if_like_collide = if_like_collide

  if (( if_collide .or. if_like_collide ) .and. q_real == 0.0 ) then
     if ( mpi_node() == 0) then
       write(0,*) "   Error reading species parameters"
       write(0,*) "   When using collisions the user must also set q_real"
       write(0,*) "   aborting..."
     endif
     stop
  endif
# 2004 "cyl_modes/os-spec-cyl-modes.f03"
  ! read momentum distribution
  call read_nml( this%udist, input_file )
# 2015 "cyl_modes/os-spec-cyl-modes.f03"
  ! Initialize fields from initial density / momentum profile
  this%init_fields = init_fields

  this%if_use_delayed_injection = if_use_delayed_injection
  this%delayed_injection_time = delayed_injection_time

  ! read profile definition unless read_prof was set to false
  if ( read_prof ) then
    this % source => new_source_cyl_modes( init_type, input_file, this%coordinates )
    if ( .not. associated( this%source )) then
      if ( mpi_node() == 0) then
        write(0,*) "   Error reading species parameters"
        write(0,*) "   Invalid 'init_type' value (", trim(init_type), ")"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
    this % source => null()
  endif

  ! read boundary condtions
  ! If the particle is a photon, characterized by this%rqm=0, we add it to the list of arguments
  if ( this % rqm /= 0.0_p_k_part ) then
    call read_nml( this%bnd_con, input_file, periodic, if_move, grid%coordinates )
  else
    call read_nml( this%bnd_con, input_file, periodic, if_move, grid%coordinates, this%rqm )
  endif

  ! don't allow different thermal temperatures in x/y at r_max
  if ( type( this%bnd_con, p_upper, p_r_dim ) == p_bc_thermal ) then
    if ( this%bnd_con%uth_bnd( 2, p_upper, p_r_dim ) /= &
         this%bnd_con%uth_bnd( 3, p_upper, p_r_dim ) ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "   Error reading species parameters"
        write(0,*) "   Thermal temperatures in x and y (dims 2 and 3)"
        write(0,*) "   must be identical for thermal boundary at r_max"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  endif

  ! read pistons
  if ( num_pistons > 0 ) then
    call alloc(this%pistons, (/ num_pistons /),"cyl_modes/os-spec-cyl-modes.f03",2060)
    do piston_id=1, num_pistons
      call read_nml_piston( this%pistons(piston_id), input_file )
    enddo

    call check_piston( this%pistons, dt )

  endif

  ! read diagnostics
  call this%diag%read_input( input_file, sim_options%gamma )

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Process physical boundaries in cylindrical modes
!-----------------------------------------------------------------------------------------
subroutine phys_boundary_spec_cyl_modes( this, current, dt, i_dim, bnd_idx, par_idx, npar )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_current_define, only : t_current
  use m_current_cyl_modes, only : t_current_cyl_modes
  use m_parameters
  use m_vdf_define, only : t_vdf
  use m_species_current, only : tmp_xold, tmp_ix, tmp_q, tmp_p, tmp_rg, tmp_xnew, tmp_dxi
  use m_random, only : rng
  use m_species_cyl_modes, only : thermal_bath_cyl_modes
  use m_species_current, only : deposit_current_1d, deposit_current_2d, deposit_current_3d
  use m_cyl_modes, only : t_cyl_modes
  use m_species_current_cyl_modes
  use m_species_boundary, only : thermal_bath
  use m_species_define, only : p_cell_near

  implicit none

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_current ), intent(inout) :: current
  real(p_double), intent(in) :: dt

  integer, intent(in) :: i_dim, bnd_idx
  integer, dimension(:), intent(in) :: par_idx
  integer, intent(in) :: npar

  ! local variables
  real(p_k_part), dimension(p_x_dim) :: dt_dx
  real(p_k_part), dimension(4) :: pos
  real(p_k_part) :: shrink
  integer :: k, i, l, nbuf, bc_type

  integer :: shift, ipos
  real(p_k_part) :: delta, vperp, pperp, xpos, x, y, r2, px, py, cos_2th, sin_2th
  real(p_k_part) :: r_new, r_old
  real(p_double) :: r_shift, dr

  ! dt / dx
  dt_dx(1:p_x_dim) = real( dt/this%dx(1:p_x_dim), p_k_part )

  ! Boundary condition type
  bc_type = this%bnd_con%type( bnd_idx, i_dim )

  ! On a moving window run, the upper boundary behaves as specular
  if ( bnd_idx == p_upper .and. bc_type == p_bc_move_c ) bc_type = p_bc_specular

  select case ( bc_type )

  case( p_bc_periodic ) ! Single node periodic boundaries

    ! check if the species is a photon : this%rqm = 0
    if ( (this%rqm==0.0_p_k_part) .and. (i_dim==2) .and. (bnd_idx == p_upper) ) then

      if ( this%pos_type == p_cell_near ) then
        r_shift = 0.5_p_double
      else
        r_shift = 0.0_p_double
      endif
      dr = this%dx(p_r_dim)

      ! theta becomes pi - theta
      do k = 1, npar
        r_old = ( this%x(2,par_idx(k)) + ( this%ix(2,par_idx(k)) + this%my_nx_p(p_lower, 2) - 2 - r_shift) ) * dr
        this%ix( 2, par_idx(k) ) = this%ix( 2, par_idx(k) ) - 1
        this%x( 2, par_idx(k) ) = - this%x( 2, par_idx(k) )
        r_new = ( this%x(2,par_idx(k)) + ( this%ix(2,par_idx(k)) + this%my_nx_p(p_lower, 2) - 2 - r_shift) ) * dr
        this%x( 3, par_idx(k) ) = - r_new * this%x( 3, par_idx(k) ) / r_old
        this%x( 4, par_idx(k) ) = - r_new * this%x( 4, par_idx(k) ) / r_old
      enddo

    else

      ! this should never actually happen
      if ( bnd_idx == p_lower ) then
        shift = + this%my_nx_p( 3, i_dim )
      else
        shift = - this%my_nx_p( 3, i_dim )
      endif

      do k = 1, npar
        this%ix( i_dim , par_idx(k) ) = this%ix( i_dim , par_idx(k) ) + shift
      enddo

    endif


  case ( p_bc_specular ) ! specular reflection

    if ( bnd_idx == p_lower ) then
      ipos = 1
    else
      ipos = this%my_nx_p( 3, i_dim )
    endif

    if ( i_dim == 1 ) then

      ! process boundary like normal
      do k = 1, npar
        ! flip momentum
        this%p(i_dim,par_idx(k) ) = -this%p(i_dim,par_idx(k) )

        ! position after reflection
        this%ix(i_dim,par_idx(k) ) = ipos
        this%x(i_dim,par_idx(k) ) = - this%x(i_dim,par_idx(k))
      enddo

    else ! i_dim == 2

      ! specularly reflect in r
      do k = 1, npar
        ! flip momentum in r only by setting
        ! px_new = -[ px*cos(2theta) + py*sin(2theta) ]
        ! py_new = py*cos(2theta) - px*sin(2theta)
        x = this%x(3,par_idx(k)); y = this%x(4,par_idx(k))
        px = this%p(2,par_idx(k)); py = this%p(3,par_idx(k))
        r2 = x**2 + y**2
        cos_2th = ( x**2 - y**2 ) / r2
        sin_2th = 2 * x * y / r2
        this%p(2,par_idx(k)) = - ( px * cos_2th + py * sin_2th )
        this%p(3,par_idx(k)) = py * cos_2th - px * sin_2th

        ! position after reflection
        this%ix(i_dim,par_idx(k)) = ipos
        this%x(i_dim,par_idx(k)) = - this%x(i_dim,par_idx(k))

        ! shrink x and y positions by the amount that r shrinks
        call this%get_position(par_idx(k),pos) ! r_new = pos(2)
        shrink = pos(2) / sqrt(r2)
        this%x(3,par_idx(k)) = this%x(3,par_idx(k)) * shrink
        this%x(4,par_idx(k)) = this%x(4,par_idx(k)) * shrink

      enddo

    endif

    ! There is no current being deposited from the reflected particles here.
    ! If the user is using some form of conducting boundaries for the fields then
    ! the electric current object will reflect the current deposited out of the
    ! simulation box back into box

  case ( p_bc_thermal ) ! thermal bath

    ! Get position of boundary
    if ( bnd_idx == p_lower ) then
      ipos = 1
      xpos = -0.5_p_k_part
    else
      ipos = this%my_nx_p( 3, i_dim )
      ! to inject inside the current grid cell, xpos should be samller than 0.5_p_k_part
      xpos = nearest(+0.5_p_k_part, -1.0_p_k_part)
    endif

    do i = 1, npar, p_cache_size

      if( i + p_cache_size > npar ) then
        nbuf = npar - i + 1
      else
        nbuf = p_cache_size
      endif

      if ( i_dim == 2 ) then
        ! Store r, cosine and sine values to properly transform momentum
        ! in thermal_bath_cyl_modes
        call this % get_position( 2, par_idx(i:i+nbuf-1), nbuf, tmp_xold(2,:) )
        call this % get_position( 3, par_idx(i:i+nbuf-1), nbuf, tmp_xold(3,:) )
        call this % get_position( 4, par_idx(i:i+nbuf-1), nbuf, tmp_xold(4,:) )
        do k = 1, nbuf
          ! tmp_xold(3,:) = cos_th
          ! tmp_xold(4,:) = sin_th
          tmp_xold(3,k) = tmp_xold(3,k) / tmp_xold(2,k)
          tmp_xold(4,k) = tmp_xold(4,k) / tmp_xold(2,k)
        enddo
      endif

      do k = 1, nbuf
        this%ix(i_dim,par_idx(i+k-1) ) = ipos
        tmp_ix(:,k) = this%ix(:,par_idx(i+k-1) )
      end do

      if ( i_dim == 1 ) then

        ! Perform standard thermal boundary handling
        call thermal_bath( this%bnd_con, i_dim, bnd_idx, tmp_p, nbuf, tmp_ix )

        do k = 1, nbuf

          ! Store new momenta
          this%p (1,par_idx(i+k-1) ) = tmp_p(1,k)
          this%p (2,par_idx(i+k-1) ) = tmp_p(2,k)
          this%p (3,par_idx(i+k-1) ) = tmp_p(3,k)

          tmp_rg(k) = 1.0 / sqrt( 1.0 + tmp_p(1,k)**2+ tmp_p(2,k)**2+ tmp_p(3,k)**2 )

          call rng % harvest_real3( delta, this%ix(:,par_idx(i+k-1)))

          vperp = tmp_rg(k) * tmp_p(i_dim,k) * dt_dx(i_dim)

          this%x (i_dim,par_idx(i+k-1) ) = xpos + vperp * delta

          ! store values for current deposition, change the sign of charge and
          ! reverse the in-plane momentum
          tmp_q( k ) = -this%q( par_idx(i+k-1) )
          ! we don't really need to change the sign of tmp_p because
          ! the in-plane momentum is not used in deposit_current_{1,2,3}d directly.
          ! The particle motion is reversed in rev_current_cyl_modes below
          ! do l = p_p_dim-p_x_dim, p_p_dim
          ! tmp_p(l,k) = -tmp_p(l,k)
          ! enddo

          do l = 1, p_x_dim
            ! the following line is commented out since we do this
            ! operation already above (right before calling 'thermal_bath')
            ! since the IX data is sometimes needed for random number generation.
            !tmp_ix( l, k ) = this%ix(l, par_idx(i+k-1) )
            tmp_xold( l, k ) = this%x(l, par_idx(i+k-1) )
            tmp_xold( l+2, k ) = this%x(l+2, par_idx(i+k-1) )
          enddo

        enddo

        ! deposit current as if particles were coming from another node
        select type( current )
        class is( t_current_cyl_modes )
          call rev_current_cyl_modes( this, current%jay_cyl_m, tmp_xold, tmp_ix, tmp_q, &
                                      tmp_p, tmp_rg, nbuf, dt_dx, dt )
        end select

      else ! i_dim == 2

        call thermal_bath_cyl_modes( this%bnd_con, i_dim, bnd_idx, tmp_p, &
                                     tmp_xold(3:4,:), nbuf, tmp_ix )

        do k = 1, nbuf

          ! Store new momenta
          this%p (1,par_idx(i+k-1) ) = tmp_p(1,k)
          this%p (2,par_idx(i+k-1) ) = tmp_p(2,k)
          this%p (3,par_idx(i+k-1) ) = tmp_p(3,k)

          tmp_rg(k) = 1.0 / sqrt( 1.0 + tmp_p(1,k)**2+ tmp_p(2,k)**2+ tmp_p(3,k)**2 )

          call rng % harvest_real3( delta, this%ix(:,par_idx(i+k-1)))

          pperp = tmp_p(2,k)*tmp_xold(3,k) + tmp_p(3,k)*tmp_xold(4,k)
          vperp = tmp_rg(k) * pperp * dt_dx(i_dim)

          this%x (i_dim,par_idx(i+k-1) ) = xpos + vperp * delta

        enddo

        ! tmp_xold(2,:) = r_new
        call this % get_position( 2, par_idx(i:i+nbuf-1), nbuf, tmp_xold(2,:) )

        do k = 1, nbuf

          ! x = cos_th * r_new
          ! y = sin_th * r_new
          this%x( 3, par_idx(i+k-1) ) = tmp_xold(3,k) * tmp_xold(2,k)
          this%x( 4, par_idx(i+k-1) ) = tmp_xold(4,k) * tmp_xold(2,k)

          ! store values for current deposition, change the sign of charge and
          ! reverse the in-plane momentum
          tmp_q( k ) = -this%q( par_idx(i+k-1) )
          ! we don't really need to change the sign of tmp_p because
          ! the in-plane momentum is not used in deposit_current_{1,2,3}d directly.
          ! The particle motion is reversed in rev_current_cyl_modes below
          ! do l = 1, p_p_dim
          ! tmp_p(l,k) = -tmp_p(l,k)
          ! enddo

          do l = 1, p_x_dim
            ! the following line is commented out since we do this
            ! operation already above (right before calling 'thermal_bath')
            ! since the IX data is sometimes needed for random number generation.
            !tmp_ix( l, k ) = this%ix(l, par_idx(i+k-1) )
            tmp_xold( l, k ) = this%x(l, par_idx(i+k-1) )
            tmp_xold( l+2, k ) = this%x(l+2, par_idx(i+k-1) )
          enddo

        enddo

        ! deposit current as if particles were coming from another node
        select type( current )
        class is( t_current_cyl_modes )
          call rev_current_cyl_modes( this, current%jay_cyl_m, tmp_xold, tmp_ix, tmp_q, &
                                      tmp_p, tmp_rg, nbuf, dt_dx, dt )
        end select

      endif

    enddo

  case default

    write(err_buf__,*) "invalid particle bc";call err__("cyl_modes/os-spec-cyl-modes.f03",2372)
    print *, 'type(', bnd_idx,', ', i_dim, ' ) = ', bc_type

  end select


  contains

  subroutine rev_current_cyl_modes( this, jay_cyl_m, xold, ix, q, p, rg, nbuf, dt_dx, dt )

    implicit none

    class( t_species_cyl_modes ), intent(inout) :: this
    type( t_cyl_modes ), intent(inout) :: jay_cyl_m
    real(p_k_part), dimension(:,:), intent(inout) :: xold
    integer, dimension(:,:), intent(inout) :: ix
    real(p_k_part), dimension( : ), intent(inout) :: q
    real(p_k_part), dimension(:,:), intent(inout) :: p
    real(p_k_part), dimension( : ), intent(inout) :: rg
    integer, intent(in) :: nbuf
    real(p_k_part), dimension(:), intent(in) :: dt_dx
    real(p_double), intent(in) :: dt

    integer :: k, shift_ix2, mode, gix2
    real(p_k_part) :: p_x, p_y, r_new, cos_th, sin_th
    real(p_double) :: rdr, r_shift



    shift_ix2 = this%my_nx_p(p_lower, 2) - 2
    rdr = 1.0_p_double / this%dx(p_r_dim)
    if ( this%pos_type == p_cell_near ) then
      r_shift = 0.5_p_double
    else
      r_shift = 0.0_p_double
    endif

    ! get previous position
    do k = 1, nbuf

      ! process first dimension like normal
      tmp_xnew(1,k) = xold(1,k) - p(1,k)*rg(k)*dt_dx(1)
      tmp_dxi(1,k) = ntrim( tmp_xnew(1,k) )

      ! process x and y together
      tmp_xnew(3,k) = xold(3,k) - p(2,k)*rg(k) * real(dt,p_k_part)
      tmp_xnew(4,k) = xold(4,k) - p(3,k)*rg(k) * real(dt,p_k_part)

      ! calculate r from x and y
      gix2 = ix(2,k) + shift_ix2
      r_new = sqrt( tmp_xnew(3,k)**2 + tmp_xnew(4,k)**2 )
      tmp_xnew(2,k) = real( r_new * rdr - ( gix2 - r_shift ) , p_k_part )
      tmp_dxi(2,k) = ntrim( tmp_xnew(2,k) )

      ! convert momentum into cylindrical coordinates for current deposition
      cos_th = tmp_xnew(3,k) / r_new
      sin_th = tmp_xnew(4,k) / r_new
      p_x = p(2,k); p_y = p(3,k)
      p(2,k) = p_x*cos_th + p_y*sin_th
      p(3,k) = -p_x*sin_th + p_y*cos_th

    enddo
    ! Only handle 2-dimensional case
    call deposit_current_2d( this, jay_cyl_m%pf_re(0), tmp_dxi, tmp_xnew, ix, xold, q, &
                              rg, p, nbuf, dt)

    ! handle high order modes - should use special cylindrical mode deposit to remain charge conserving
    do mode = 1, this%n_cyl_modes

      ! still need to port cylindrical mode current deposition functions
      select case (this%interpolation)
      case( p_linear )

        call getjr_cyl_m_s1( jay_cyl_m%pf_re(mode), jay_cyl_m%pf_im(mode), tmp_dxi, &
                              tmp_xnew, ix, xold, q, rg, p, &
                              nbuf, dt, shift_ix2, mode )

      case( p_quadratic )

        call getjr_cyl_m_s2( jay_cyl_m%pf_re(mode), jay_cyl_m%pf_im(mode), tmp_dxi, &
                              tmp_xnew, ix, xold, q, rg, p, &
                              nbuf, dt, shift_ix2, mode )

      case( p_cubic )

        call getjr_cyl_m_s3( jay_cyl_m%pf_re(mode), jay_cyl_m%pf_im(mode), tmp_dxi, &
                              tmp_xnew, ix, xold, q, rg, p, &
                              nbuf, dt, shift_ix2, mode )

      case default
        write(err_buf__,*) 'Not implemented yet';call err__("cyl_modes/os-spec-cyl-modes.f03",2462)
        call abort_program( p_err_notimplemented )
      end select

    enddo ! mode
# 2477 "cyl_modes/os-spec-cyl-modes.f03"
  end subroutine rev_current_cyl_modes

  function ntrim(x)

    implicit none

    real(p_k_part), intent(in) :: x
    integer :: ntrim, a, b

    if ( x < -.5 ) then
    a = -1
    else
      a = 0
    endif

    if ( x >= .5 ) then
    b = +1
    else
      b = 0
    endif

    ntrim = a+b

  end function ntrim

end subroutine phys_boundary_spec_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Report species diagnostic
!-----------------------------------------------------------------------------------------
subroutine report_species_cyl_modes( this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )

  use m_parameters
  use m_vdf_math, only : mult
  use m_vdf_define, only : t_vdf, t_vdf_report
  use m_vdf_comm, only : t_vdf_msg
  use m_species_udist, only : spatial_ufl, spatial_boost, spatial_uth, spatial_uth_tensor
  use m_species_udist, only : spatial_momentum_flux_tensor, spatial_nufl
  use m_emf_define, only : t_emf
  use m_space, only : t_space
  use m_grid_define, only : t_grid
  use m_node_conf, only : t_node_conf
  use m_time_step, only : t_time_step, test_if_report, n, ndump
  use m_species_diagnostics, only : p_ufl1, p_ufl2, p_ufl3, p_uth1, p_uth2, p_uth3, p_norm
  use m_species_diagnostics, only : p_T11, p_T22, p_T33, p_T12, p_T13, p_T23
  use m_species_diagnostics, only : p_vfl1, p_vfl2, p_vfl3, p_P11, p_P22, p_P33
  use m_species_diagnostics, only : p_P12, p_P13, p_P23, p_nufl1, p_nufl2, p_nufl3
  use m_species_diagnostics, only : deposit_quant, report_temperature, report_heat_flux
  use m_species_diagnostics, only : diag_part_boost_ev, diag_part_boost_extract_ev
  use m_vdf_report, only : if_report, report_vdf
  use m_vdf_comm, only : t_vdf_msg
  use m_species_tracks, only : write_tracks, create_file
  use m_species_phasespace, only : report
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_species_rawdiag, only : write_raw, extract_boosted_raw, write_boosted_raw
  use m_cyl_modes, only : t_cyl_modes, report_cyl_modes, cleanup, mult
  use m_species_cyl_modes, only : add_track_set_data_cyl_modes, deposit_quant_cyl_modes
  use m_species_cyl_modes, only : get_current_4vector_cyl_modes
  use m_emf_cyl_modes, only : t_emf_cyl_modes
  use m_logprof, only: begin_event, end_event
  use m_boosted_diag_cyl_modes, only : t_boosted_diag_cyl_modes

  implicit none

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_emf ), intent(inout) :: emf
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t, tmin
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  type( t_vdf ) :: charge, norm, ufl, uth, vboost, Tdiag, Tcross, current_4vector
  type( t_vdf ) :: Pdiag, Pcross, nufl
  type( t_cyl_modes ) :: charge_cyl_m, norm_cyl_m, current_4vector_cyl_m
  type( t_vdf_report ), pointer :: rep

  ! Density diagnostics
  rep => this%diag%reports
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then

      ! deposit quantity (this also creates the required charge vdf)
      call deposit_quant_cyl_modes( this, grid, no_co, charge_cyl_m, charge, rep%quant, send_msg, recv_msg )

      ! do all reports for current quantity
      call report_cyl_modes( rep, charge_cyl_m , 1, g_space, grid, no_co, tstep, t )

      ! cleanup temp. grid
      call cleanup( charge_cyl_m )
      call charge % cleanup()

    endif

    rep => rep%next
  enddo

  ! boosted diagnostics
  if( this%diag%if_use_boosted_diag ) then

    call begin_event( diag_part_boost_extract_ev )

    call get_current_4vector_cyl_modes( this, grid, no_co, current_4vector_cyl_m, current_4vector, send_msg, recv_msg )

    select type( boosted_diag => this%diag%boosted_diag ); class is ( t_boosted_diag_cyl_modes )
      call boosted_diag%extract_cyl_modes( tstep, t, grid, g_space, current_4vector_cyl_m )
    end select

    call cleanup( current_4vector_cyl_m )
    call current_4vector % cleanup()

    call end_event( diag_part_boost_extract_ev )

    call begin_event( diag_part_boost_ev )

    call this%diag%boosted_diag%report( tstep, grid )

    call end_event( diag_part_boost_ev )

  endif

  if( this%diag%if_use_boosted_raw ) then

    call begin_event( diag_part_boost_extract_ev )

    select type( boosted_diag => this%diag%boosted_diag ); class is ( t_boosted_diag_cyl_modes )
      call extract_boosted_raw( this, boosted_diag%pf_re(0), tstep, t, grid, g_space )
    end select

    call end_event( diag_part_boost_extract_ev )

    if ( test_if_report( tstep, this%diag%boosted_diag%ndump_fac_raw ) ) then
      call begin_event( diag_part_boost_ev )

      select type( boosted_diag => this%diag%boosted_diag ); class is ( t_boosted_diag_cyl_modes )
        call write_boosted_raw( this, boosted_diag%pf_re(0), no_co, grid, n(tstep) )
      end select

      call end_event( diag_part_boost_ev )
    endif

  endif


  ! Momentum distribution diagnostics
  call norm % cleanup()
  call ufl % cleanup()
  call vboost % cleanup()
  call uth % cleanup()
  call Tdiag % cleanup()
  call Tcross % cleanup()
  call Pdiag % cleanup()
  call Pcross % cleanup()
  call nufl % cleanup()

  rep => this%diag%rep_udist
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then
      select case ( rep%quant )

      case ( p_ufl1, p_ufl2, p_ufl3 )

        ! Initialize fluid velocity
        if ( ufl%x_dim_ < 0 ) then
          ! deposit normalization value
          if ( norm%x_dim_ < 0 ) then
            call deposit_quant( this, grid, no_co, norm, p_norm, send_msg, recv_msg )
          endif

          ! Get fluid velocity
          call ufl % new( norm, f_dim = 3 )
          call spatial_ufl( this, no_co, norm, ufl, send_msg, recv_msg )
        endif

        ! do all reports for current quantity
        call report_vdf( rep, ufl, rep%quant-p_ufl1+1, g_space, grid, no_co, tstep, t )

      case ( p_uth1, p_uth2, p_uth3 )

        ! Initialize thermal velocity
        if ( uth%x_dim_ < 0 ) then
          ! deposit normalization value
          if ( norm%x_dim_ < 0 ) then
            call deposit_quant( this, grid, no_co, norm, p_norm, send_msg, recv_msg )
          endif

          ! Get fluid velocity
          if ( ufl%x_dim_ < 0 ) then
            call ufl % new( norm, f_dim = 3 )
            call spatial_ufl( this, no_co, norm, ufl, send_msg, recv_msg )
          endif

          ! Get thermal velocity
          call uth % new( ufl )
          call spatial_uth( this, no_co, norm, ufl, uth, send_msg, recv_msg )
        endif

        ! do all reports for current quantity
        call report_vdf( rep, uth, rep%quant-p_uth1+1, g_space, grid, no_co, tstep, t )

      case ( p_T11, p_T22, p_T33, p_T12, p_T13, p_T23)

        ! Initialize thermal velocity tensor
        if ( Tdiag%x_dim_ < 0 ) then
          ! deposit normalization value
          if ( norm%x_dim_ < 0 ) then
            call deposit_quant( this, grid, no_co, norm, p_norm, send_msg, recv_msg )
          endif

          ! Get boost velocity
          call vboost % new( norm, f_dim = 3 )
          call spatial_boost( this, no_co, norm, vboost, send_msg, recv_msg )

          ! Get thermal velocity tensor
          call Tdiag % new( vboost )
          call Tcross % new( vboost )
          call spatial_uth_tensor( this, no_co, norm, vboost, Tdiag, Tcross, send_msg, recv_msg )
        endif

        ! do all reports for current quantity
        select case ( rep%quant )
        case ( p_T11, p_T22, p_T33 )
          !Initialize the diagonal terms of the temperature tensor
          call report_vdf( rep, Tdiag, rep%quant-p_T11+1, g_space, grid, no_co, tstep, t )
        case ( p_T12, p_T13, p_T23 )
          !Initialize the off-diagonal terms of the temperature tensor
          call report_vdf( rep, Tcross, rep%quant-p_T12+1, g_space, grid, no_co, tstep, t )
        end select

      case (p_vfl1, p_vfl2, p_vfl3, p_P11, p_P22, p_P33, p_P12, p_P13, p_P23, p_nufl1, p_nufl2, p_nufl3)

        ! Initialize thermal velocity tensor
        if ( Pdiag%x_dim_ < 0 ) then
          ! deposit normalization value
          if ( norm%x_dim_ < 0 ) then
            call deposit_quant( this, grid, no_co, norm, p_norm, send_msg, recv_msg )
          endif

          ! Get boost velocity
          call vboost % new( norm, f_dim = 3 )
          call spatial_boost( this, no_co, norm, vboost, send_msg, recv_msg )

          ! Get thermal velocity tensor
          call Pdiag % new( vboost )
          call Pcross % new( vboost )
          call spatial_momentum_flux_tensor( this, no_co, Pdiag, Pcross, send_msg, recv_msg )

          call nufl % new( vboost )
          call spatial_nufl( this, no_co, nufl, send_msg, recv_msg )
        endif

        ! do all reports for current quantity
        select case ( rep%quant )
          case ( p_vfl1, p_vfl2, p_vfl3 )
            call report_vdf( rep, vboost, rep%quant-p_vfl1+1, g_space, grid, no_co, tstep, t )
          case ( p_P11, p_P22, p_P33 )
            !Initialize the diagonal terms of the temperature tensor
            call report_vdf( rep, Pdiag, rep%quant-p_P11+1, g_space, grid, no_co, tstep, t )
          case ( p_P12, p_P13, p_P23 )
            !Initialize the off-diagonal terms of the temperature tensor
            call report_vdf( rep, Pcross, rep%quant-p_P12+1, g_space, grid, no_co, tstep, t )
          case ( p_nufl1, p_nufl2, p_nufl3 )
            call report_vdf( rep, nufl, rep%quant-p_nufl1+1, g_space, grid, no_co, tstep, t )
        end select

      end select
    endif

    rep => rep%next
  enddo

  ! Clear ufl and uth vdfs, but keep norm vdf if available
  call ufl % cleanup()
  call uth % cleanup()
  call vboost % cleanup()
  call Tdiag % cleanup()
  call Tcross % cleanup()
  call Pdiag % cleanup()
  call Pcross % cleanup()
  call nufl % cleanup()

  ! Cell average diagnostics
  rep => this%diag%rep_cell_avg
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then

      ! deposit normalization value
      if ( norm % x_dim_ < 1 ) then
        call deposit_quant_cyl_modes( this, grid, no_co, norm_cyl_m, norm, p_norm, send_msg, recv_msg )
      endif

      ! deposit quantity and normalize it
      call deposit_quant_cyl_modes( this, grid, no_co, charge_cyl_m, charge, rep%quant, send_msg, recv_msg )
      call mult( charge_cyl_m, norm_cyl_m )

      ! do all reports for current quantity
      call report_cyl_modes( rep, charge_cyl_m , 1, g_space, grid, no_co, tstep, t )

      ! cleanup temp. grid
      call cleanup( charge_cyl_m )
      call charge % cleanup()

    endif

    rep => rep%next
  enddo
  call cleanup( norm_cyl_m )
  call norm % cleanup()

  ! Phasespaces
  call report( this%diag%phasespaces, this, no_co, g_space, tstep, t )

  ! Raw particle data diagnostics
  if (test_if_report( tstep, this%diag%ndump_fac_raw )) then
     ! new hdf5 routine
     call write_raw( this, no_co, n(tstep), t, n(tstep) / ndump(tstep) )
  endif

  ! Temperature and heat flux
  call report_temperature( this, no_co, tstep, t, tmin )
  call report_heat_flux( this, no_co, tstep, t, tmin )

  ! Particle tracking diagnostics



  if (this%diag%ndump_fac_tracks > 0) then
    if ( mod( n(tstep), this%diag%tracks%niter ) == 0 ) then
      if ( n( tstep ) > this%diag%n_start_tracks ) then

        select type( emf )
        class is (t_emf_cyl_modes)
          call add_track_set_data_cyl_modes( this%diag%tracks, this, no_co, emf, n(tstep), t, send_msg, recv_msg )
        end select
      endif
    endif
  endif

  if ( test_if_report( tstep, this%diag%ndump_fac_tracks )) then
     if ( n(tstep) == 0 ) then
        ! create the diagnostics file
        call create_file( this%diag%tracks, this, ndump(tstep) )
     else
        ! write the data
        call write_tracks( this%diag%tracks, no_co, this )
     endif
  endif



end subroutine report_species_cyl_modes

!-----------------------------------------------------------------------------------------
! Write metadata for track files
!-----------------------------------------------------------------------------------------
subroutine enumerate_quants_spec_cyl_modes( this, diagFile, track_set )

  use m_parameters
  use m_diagnostic_utilities, only : t_diag_file
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_species_define, only : t_track_set

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  class( t_diag_file ),intent(inout) :: diagFile
  type( t_track_set ), intent(inout) :: track_set
  integer :: i, j, k

  ! enumerate quants
  j = 0

  j = j+1
  diagFile % tracks % quants(j) = 'n'
  diagFile % tracks % qlabels(j) = 'n'
  diagFile % tracks % qunits(j) = ''
  diagFile % tracks % offset_t(j)= 0.0_p_double

  j = j+1
  diagFile % tracks % quants(j) = 't'
  diagFile % tracks % qlabels(j) = 't'
  diagFile % tracks % qunits(j) = '1/\omega_p'
  diagFile % tracks % offset_t(j)= 0.0_p_double

  j = j+1
  diagFile % tracks % quants(j) = 'q'
  diagFile % tracks % qlabels(j) = 'q'
  diagFile % tracks % qunits(j) = 'e'
  diagFile % tracks % offset_t(j)= 0.0_p_double

  j = j+1
  diagFile % tracks % quants(j) = 'ene'
  diagFile % tracks % qlabels(j) = 'Ene'
  diagFile % tracks % qunits(j) = 'm c^2'
  diagFile % tracks % offset_t(j)= -0.5_p_double

  do i = 1, this % get_n_x_dims()
    j = j+1
    diagFile % tracks % quants(j) = 'x'//(char(iachar('0')+i))
    diagFile % tracks % qlabels(j) = 'x_'//(char(iachar('0')+i))
    diagFile % tracks % qunits(j) = 'c/\omega_p'
    diagFile % tracks % offset_t(j)= 0.0_p_double
  enddo

  do i = 1, p_p_dim
    j = j+1
    diagFile % tracks % quants(j) = 'p'//(char(iachar('0')+i))
    diagFile % tracks % qlabels(j) = 'p_'//(char(iachar('0')+i))
    diagFile % tracks % qunits(j) = 'm c'
    diagFile % tracks % offset_t(j)= -0.5_p_double
  enddo
# 2911 "cyl_modes/os-spec-cyl-modes.f03"
  do i = 1, p_f_dim
    if (track_set%ifdmp_tracks_efl(i)) then
      j = j+1
      diagFile % tracks % quants(j) = 'E'//(char(iachar('0')+i))
      diagFile % tracks % qlabels(j) = 'E_'//(char(iachar('0')+i))
      diagFile % tracks % qunits(j) = 'm_e c \omega_p/e'
      diagFile % tracks % offset_t(j)= 0.0_p_double

      do k = 0, this % n_cyl_modes
        j = j+1
        diagFile % tracks % quants(j) = 'E'//(char(iachar('0')+i))// &
                                        ' MODE '//(char(iachar('0')+k))
        diagFile % tracks % qlabels(j) = 'E_'//(char(iachar('0')+i))// &
                                        ' MODE '//(char(iachar('0')+k))
        diagFile % tracks % qunits(j) = 'm_e c \omega_p/e'
        diagFile % tracks % offset_t(j)= 0.0_p_double
      enddo
    endif
  enddo

  do i = 1, p_f_dim
    if (track_set%ifdmp_tracks_bfl(i)) then
      j = j+1
      diagFile % tracks % quants(j) = 'B'//(char(iachar('0')+i))
      diagFile % tracks % qlabels(j) = 'B_'//(char(iachar('0')+i))
      diagFile % tracks % qunits(j) = 'm_e c \omega_p/e'
      diagFile % tracks % offset_t(j)= 0.0_p_double

      do k = 0, this % n_cyl_modes
        j = j+1
        diagFile % tracks % quants(j) = 'B'//(char(iachar('0')+i))// &
                                        ' MODE '//(char(iachar('0')+k))
        diagFile % tracks % qlabels(j) = 'B_'//(char(iachar('0')+i))// &
                                        ' MODE '//(char(iachar('0')+k))
        diagFile % tracks % qunits(j) = 'm_e c \omega_p/e'
        diagFile % tracks % offset_t(j)= 0.0_p_double
      enddo
    endif
  enddo

  if (track_set%ifdmp_tracks_psi) then
    j = j+1
    diagFile % tracks % quants(j) = 'psi'
    diagFile % tracks % qlabels(j) = 'Psi'
    diagFile % tracks % qunits(j) = ''
    diagFile % tracks % offset_t(j)= 0.0_p_double
  endif

  diagFile % tracks % nquants = j

end subroutine enumerate_quants_spec_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Position functions begin
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_idx_comp_4_cyl_modes( this, comp, idx, np, pos )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_parameters, only : p_single, p_x_dim, p_lower

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_single ), dimension(:), intent(out) :: pos

  integer :: j, ixmin
  real( p_single ) :: xmin, dx

  if (comp <= p_x_dim) then
    xmin = real( this%g_box( p_lower, comp ), p_single )
    ixmin = this%my_nx_p( 1, comp ) - 2
    dx = real( this%dx( comp ), p_single )
  endif

  if (comp > p_x_dim) then ! extra coordinates do not have a corresponding ix
    do j = 1, np
      pos(j) = real(this%x(comp, idx(j)), p_single)
    enddo
  else
    do j = 1, np
      pos(j) = xmin + dx * real((this%ix( comp, idx(j) ) + ixmin) + &
                                 this%x ( comp, idx(j) ) , p_single )
    enddo
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_idx_comp_8_cyl_modes( this, comp, idx, np, pos )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_parameters, only : p_double, p_x_dim, p_lower

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_double ), dimension(:), intent(out) :: pos

  integer :: j, ixmin
  real( p_double ) :: xmin, dx

  if (comp <= p_x_dim) then
    xmin = this%g_box( p_lower, comp )
    ixmin = this%my_nx_p( 1, comp ) - 2
    dx = this%dx( comp )
  endif

  if (comp > p_x_dim) then ! extra coordinates do not have a corresponding ix
    do j = 1, np
      pos(j) = real(this%x(comp, idx(j)), p_double)
    enddo
  else
    do j = 1, np
      pos(j) = xmin + dx * (real(( this%ix( comp, idx(j) ) + ixmin), p_double ) + &
                                    this%x( comp, idx(j) ) )
    enddo
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_ref_box_idx_comp_4_cyl_modes( this, comp, idx, np, pos, if_pos_ref_box )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_parameters, only : p_single, p_x_dim, p_lower

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_single ), dimension(:), intent(out) :: pos
  logical, intent(in) :: if_pos_ref_box

  integer :: j, ixmin
  real( p_single ) :: xmin, dx

  if (comp <= p_x_dim) then
    if (if_pos_ref_box) then
      xmin = 0.0_p_single
    else
      xmin = real( this%g_box( p_lower, comp ), p_single )
    endif
    ixmin = this%my_nx_p( 1, comp ) - 2
    dx = real( this%dx( comp ), p_single )
  endif

  if (comp > p_x_dim) then ! extra coordinates do not have a corresponding ix
    do j = 1, np
      pos(j) = real(this%x(comp, idx(j)), p_single)
    enddo
  else
    do j = 1, np
      pos(j) = xmin + dx * real((this%ix( comp, idx(j) ) + ixmin) + &
                                 this%x ( comp, idx(j) ) , p_single )
    enddo
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_ref_box_idx_comp_8_cyl_modes( this, comp, idx, np, pos, if_pos_ref_box )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_parameters, only : p_double, p_x_dim, p_lower

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_double ), dimension(:), intent(out) :: pos
  logical, intent(in) :: if_pos_ref_box

  integer :: j, ixmin
  real( p_double ) :: xmin, dx

  if (comp <= p_x_dim) then
    if (if_pos_ref_box) then
      xmin = 0.0_p_double
    else
      xmin = this%g_box( p_lower, comp )
    endif
    ixmin = this%my_nx_p( 1, comp ) - 2
    dx = this%dx( comp )
  endif

  if (comp > p_x_dim) then ! extra coordinates do not have a corresponding ix
    do j = 1, np
      pos(j) = real(this%x(comp, idx(j)), p_double)
    enddo
  else
    do j = 1, np
      pos(j) = xmin + dx * (real(( this%ix( comp, idx(j) ) + ixmin), p_double ) + &
                                    this%x( comp, idx(j) ) )
    enddo
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_single_4_cyl_modes( this, idx, pos )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_parameters, only : p_single, p_x_dim, p_lower

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: idx
  real( p_single ), dimension(:), intent(out) :: pos

  integer :: j, n_x_dim

  n_x_dim = size( this%x, 1 )

  do j = 1, p_x_dim
    pos(j) = real( ( (this%ix( j, idx ) + this%my_nx_p(p_lower, j) - 2) + &
                      this%x( j, idx ) ) * this%dx( j ) + &
                      this%g_box( p_lower, j ), p_single )
  enddo

  if (n_x_dim > p_x_dim) then ! extra coordinates do not have a corresponding ix
    do j=p_x_dim+1, n_x_dim
      pos(j) = real(this%x(j,idx))
    enddo
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_single_8_cyl_modes( this, idx, pos )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_parameters, only : p_double, p_x_dim, p_lower

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: idx
  real( p_double ), dimension(:), intent(out) :: pos

  integer :: j, n_x_dim

  n_x_dim = size( this%x, 1 )

  do j = 1, p_x_dim
    pos(j) = this%g_box( p_lower, j ) + this%dx( j ) * &
            ( real((this%ix( j, idx ) + this%my_nx_p(p_lower, j) - 2), p_double) + &
                    this%x( j, idx ) )
  enddo

  if (n_x_dim > p_x_dim) then ! extra coordinates do not have a corresponding ix
    do j=p_x_dim+1, n_x_dim
      pos(j) = real(this%x(j,idx))
    enddo
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_range_comp_4_cyl_modes( this, comp, idx0, idx1, pos )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_parameters, only : p_single, p_x_dim, p_lower

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp, idx0, idx1
  real( p_single ), dimension(:), intent(out) :: pos

  integer :: j, ixmin
  real( p_single ) :: xmin, dx

  if (comp <= p_x_dim) then

    xmin = real( this%g_box( p_lower, comp ), p_single )
    ixmin = this%my_nx_p( 1, comp ) - 2
    dx = real( this%dx( comp ), p_single )

    do j = idx0, idx1
      pos(j-idx0+1) = xmin + dx * (real(( this%ix( comp, j ) + ixmin ) + &
                                          this%x( comp, j ), p_single ) )
    enddo

  else

    do j = idx0, idx1
      pos(j-idx0+1) = real( this%x( comp, j ), p_single )
    enddo

  endif

end subroutine position_range_comp_4_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_range_comp_8_cyl_modes( this, comp, idx0, idx1, pos )

  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_parameters, only : p_double, p_x_dim, p_lower

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp, idx0, idx1
  real( p_double ), dimension(:), intent(out) :: pos

  integer :: j, ixmin
  real( p_double ) :: xmin, dx

  if (comp <= p_x_dim) then

    xmin = this%g_box( p_lower, comp )
    ixmin = this%my_nx_p( 1, comp ) - 2
    dx = this%dx( comp )

    do j = idx0, idx1
      pos(j-idx0+1) = xmin + dx * (real(( this%ix( comp, j ) + ixmin ), p_double ) + &
                                          this%x( comp, j ) )
    enddo

  else

    do j = idx0, idx1
      pos(j-idx0+1) = real( this%x( comp, j ), p_double )
    enddo

  endif

end subroutine position_range_comp_8_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine create_particle_single_cell_cyl( this, ix, x, q )
!-----------------------------------------------------------------------------------------
! creates a single particle in this this
! does not initialize particle momentum (to be done with set_momentum)
! ix and x are expected to match the this pos_type:
! - cell_lower => 0 <= x < 1
! - cell_near => -0.5 <= x < 0.5
!-----------------------------------------------------------------------------------------

  use m_parameters
  use m_species_define, only : p_spec_buf_block
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_species_tag, only : set_tags

  implicit none

  class(t_species_cyl_modes), intent(inout) :: this
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x
  real(p_k_part), intent(in) :: q

  integer :: bseg_size

  if ( this%num_par + 1 > this%num_par_max ) then
    ! buffer is full, resize it
    bseg_size = max(this%num_par_max / 16, p_spec_buf_block )
    call this%grow_buffer( this%num_par_max + bseg_size )
  endif

  this%num_par = this%num_par + 1

  ! if required set tags of particles
  if (this%add_tag) then
    call set_tags( this, this%num_par )
  endif

  this%x(:,this%num_par) = x
  this%ix(1:p_x_dim,this%num_par) = ix(1:p_x_dim)

  this%q( this%num_par) = q
  this%num_created = this%num_created + 1


 end subroutine create_particle_single_cell_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine create_particle_single_cell_p_cyl( this, ix, x, p, q )
!-----------------------------------------------------------------------------------------
! creates a single particle in this this and initialized momentum to the given value
! ix and x are expected to match the this pos_type:
! - cell_lower => 0 <= x < 1
! - cell_near => -0.5 <= x < 0.5
!-----------------------------------------------------------------------------------------

  use m_parameters
  use m_species_define, only : p_spec_buf_block
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_species_tag, only : set_tags

  implicit none

  class(t_species_cyl_modes), intent(inout) :: this
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x
  real(p_k_part), dimension(:), intent(in) :: p
  real(p_k_part), intent(in) :: q

  integer :: bseg_size
# 3354 "cyl_modes/os-spec-cyl-modes.f03"
  if ( this%num_par + 1 > this%num_par_max ) then
    ! buffer is full, resize it
    bseg_size = max(this%num_par_max / 16, p_spec_buf_block )
    call this%grow_buffer( this%num_par_max + bseg_size )
  endif

  this%num_par = this%num_par + 1

  ! if required set tags of particles
  if (this%add_tag) then
    call set_tags( this, this%num_par )
  endif

  this%x(:,this%num_par) = x
  this%ix(1:p_x_dim,this%num_par) = ix(1:p_x_dim)
  this%q(this%num_par) = q
  this%p(1,this%num_par) = p(1)
  this%p(2,this%num_par) = p(2)
  this%p(3,this%num_par) = p(3)

  this%num_created = this%num_created + 1


 end subroutine create_particle_single_cell_p_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return species kinetic energy on local node
!-----------------------------------------------------------------------------------------
function get_energy_spec_cyl_modes( this )

  use m_parameters
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_species_define, only : p_ene_recalc
  use m_math, only : pi

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  real( p_double ) :: get_energy_spec_cyl_modes

  integer :: i
  real( p_double ) :: energy, u2, gamma, kin, cell_volume

  ! Check if energy was calculated during the push (time-centered), otherwise
  ! calculate it here
  if ( this % energy(1) == p_ene_recalc ) then
    ! Energy not available, recalculate it
    energy = 0

    !$omp parallel do private(u2, gamma, kin) reduction(+ : energy)
    do i = 1, this%num_par
      u2 = this%p(1,i)**2 + this%p(2,i)**2 + this%p(3,i)**2
      gamma = sqrt( u2 + 1 )
      kin = u2 / (gamma + 1)
      energy = energy + this%q(i) * kin
    enddo
    !$omp end parallel do

    ! Store value so it can be reused
    this%energy(1) = energy

    do i = 2, ubound(this % energy, 1)
      this % energy(i) = 0
    enddo


  else

    ! Get energy calculated previously
    energy = this % energy(1)


    ! If doing a multi-threaded run, add contribution from all threads
    do i = 2, ubound(this % energy, 1)
      energy = energy + this % energy(i)
    enddo


  endif

  ! Normalize to cell size and charge over mass ratio
  ! Include factor of 2*pi due to cylindrical geometry
  cell_volume = this % dx(1) * this % dx(2) * 2.0_p_double * pi

  get_energy_spec_cyl_modes = energy * this%rqm * cell_volume


end function get_energy_spec_cyl_modes
!-----------------------------------------------------------------------------------------
