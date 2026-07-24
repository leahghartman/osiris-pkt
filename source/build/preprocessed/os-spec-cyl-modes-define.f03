# 1 "cyl_modes/os-spec-cyl-modes-define.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-spec-cyl-modes-define.f03" 2

!-----------------------------------------------------------------------------------------
! Species quasi-3D definition module
!
! This file contains the class definition for the following classes:
!
! t_species_cyl_modes
!-----------------------------------------------------------------------------------------


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
# 12 "cyl_modes/os-spec-cyl-modes-define.f03" 2
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
# 13 "cyl_modes/os-spec-cyl-modes-define.f03" 2

module m_species_cyl_modes_define

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
# 17 "cyl_modes/os-spec-cyl-modes-define.f03" 2


use m_system
use m_parameters
use m_time_step, only : t_time_step
use m_space, only : t_space
use m_grid_define, only : t_grid
use m_node_conf, only : t_node_conf
use m_restart, only : t_restart_handle
use m_input_file, only : t_input_file
use m_emf_define, only : t_emf
use m_emf_cyl_modes, only : t_emf_cyl_modes
use m_vdf_define, only : t_vdf, t_vdf_report
use m_vdf_comm, only : t_vdf_msg
use m_current_define, only : t_current
use m_species_define, only : t_species, t_part_idx, t_spec_msg, t_track_set
use m_diagnostic_utilities, only : p_diag_prec, t_diag_file
use m_cyl_modes, only : t_cyl_modes

implicit none

private

type, extends( t_species ) :: t_species_cyl_modes

  integer :: num_par_theta
  logical :: rand_theta = .false.
  logical :: shift_cell_spokes = .true. ! spokes in a single cell will be offset for
                                        ! num_par > 1
  integer :: n_cyl_modes = -1

  ! Parameters for charge correction near the axis
  real(p_k_part), dimension(4) :: alpha_const ! Standard coefficients, 4 for up to quartic
  real(p_k_part) :: dx2_thresh ! Threshold in x2 below which we modify particle charge
  ! Parameters for 1 ppc in x2
  integer :: n_alpha ! number of charge coefficients used
  real(p_k_part), dimension(6) :: alpha_single ! Coefficients for 1 ppcell in x2

contains

  procedure :: read_input => read_input_species_cyl_modes
  procedure :: allocate_objs => allocate_objs_spec_cyl_modes
  procedure :: init_buffer => init_particle_buffer_cyl_modes
  procedure :: init => init_species_cyl_modes
  procedure :: get_phasespace_axis => get_phasespace_axis_cyl_modes
  procedure :: phys_boundary => phys_boundary_spec_cyl_modes

  ! Can't override push since we need jay_cyl_m instead of just the vdf
  procedure :: push_cyl => push_species_cyl_modes
  procedure :: get_emf => get_emf_spec_cyl_modes
  procedure :: create_particle_single_cell_p_cyl => create_particle_single_cell_p_cyl
  procedure :: create_particle_single_cell_cyl => create_particle_single_cell_cyl
  procedure :: position_single_4 => position_single_4_cyl_modes
  procedure :: position_single_8 => position_single_8_cyl_modes
  procedure :: position_range_comp_4 => position_range_comp_4_cyl_modes
  procedure :: position_range_comp_8 => position_range_comp_8_cyl_modes
  procedure :: position_idx_comp_4 => position_idx_comp_4_cyl_modes
  procedure :: position_idx_comp_8 => position_idx_comp_8_cyl_modes
  procedure :: position_ref_box_idx_comp_4 => position_ref_box_idx_comp_4_cyl_modes
  procedure :: position_ref_box_idx_comp_8 => position_ref_box_idx_comp_8_cyl_modes

  ! diagnostic
  procedure :: report => report_species_cyl_modes
  procedure :: enumerate_quants => enumerate_quants_spec_cyl_modes

  procedure :: get_n_x_dims => get_n_x_dims_cyl_modes
  procedure :: get_energy => get_energy_spec_cyl_modes

  ! charge normalization near the axis
  procedure :: norm_charge_axis => norm_charge_axis_spec_cyl_modes

end type t_species_cyl_modes

! this extra type necessary to make an array of pointers in fortran
type :: t_spec_cyl_arr

  class( t_species_cyl_modes ), pointer :: s => null()

end type t_spec_cyl_arr

public :: t_species_cyl_modes, t_spec_cyl_arr


!-----------------------------------------------------------------------------------------
! Species quasi-3D interfaces
!
! methods defined in os-spec-cyl-modes.f03
!-----------------------------------------------------------------------------------------

interface
  subroutine get_emf_spec_cyl_modes( this, emf, bp, ep, np, ptrcur, t )

    import t_species_cyl_modes, t_emf, p_k_part, p_double

    class( t_species_cyl_modes ), intent(in) :: this
    class( t_emf ), intent(in), target :: emf
    real(p_k_part), dimension(:,:), intent(out) :: bp, ep
    integer, intent(in) :: np, ptrcur
    real(p_double), intent(in), optional :: t

  end subroutine
end interface

interface
  subroutine read_input_species_cyl_modes( this, input_file, def_name, periodic, if_move, grid, &
                             dt, read_prof, sim_options )

  import t_species_cyl_modes, t_input_file, t_options, p_double, t_grid

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  character(len = *), intent(in) :: def_name
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  type(t_options), intent(in) :: sim_options
  real(p_double), intent(in) :: dt
  logical, intent(in) :: read_prof


  end subroutine
end interface

interface
  subroutine allocate_objs_spec_cyl_modes( this )

  import t_species_cyl_modes

  class( t_species_cyl_modes ), intent(inout) :: this

  end subroutine
end interface

interface
  subroutine init_particle_buffer_cyl_modes( this, num_par_req )

  import t_species_cyl_modes

  class(t_species_cyl_modes), intent(inout) :: this
  integer, intent(in) :: num_par_req

  end subroutine
end interface

! Add interface to superclass init function to avoid type casting
interface
  subroutine init_species( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
              no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
              recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )

  import t_species, t_emf, t_current, t_space, t_grid, t_node_conf, t_restart_handle, &
         t_options, t_vdf_msg, t_part_idx, t_spec_msg, p_double, t_time_step

  class( t_species ), intent(inout) :: this

  integer, intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in) :: grid
  type( t_space ), intent(in) :: g_space
  class( t_emf ), intent(inout) :: emf
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_vdf, recv_vdf
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_spec, recv_spec
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type(t_options), intent(in) :: sim_options
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax

  end subroutine
end interface

interface
  subroutine get_phasespace_axis_cyl_modes( spec, xp, l, lp, x_or_p, xp_dim )

  import t_species_cyl_modes, p_diag_prec

  class(t_species_cyl_modes), intent(in) :: spec
  real(p_diag_prec), dimension(:), intent(out) :: xp
  integer, intent(in) :: l, lp
  integer, intent(in) :: x_or_p
  integer, intent(in) :: xp_dim

  end subroutine
end interface

interface
  subroutine phys_boundary_spec_cyl_modes( this, current, dt, i_dim, bnd_idx, par_idx, npar )

    import t_species_cyl_modes, t_current, p_double

    class( t_species_cyl_modes ), intent(inout) :: this
    class( t_current ), intent(inout) :: current
    real(p_double), intent(in) :: dt
    integer, intent(in) :: i_dim, bnd_idx
    integer, dimension(:), intent(in) :: par_idx
    integer, intent(in) :: npar

  end subroutine
end interface

interface
  subroutine push_species_cyl_modes( this, emf, jay_cyl_m, t, tstep, tid, n_threads, options )

  import t_species_cyl_modes, t_emf_cyl_modes, t_cyl_modes, t_time_step, t_options, p_double

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_emf_cyl_modes ), intent( inout ) :: emf
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  type( t_options ), intent(in) :: options
  integer, intent(in) :: tid
  integer, intent(in) :: n_threads

  end subroutine
end interface

interface
  subroutine create_particle_single_cell_cyl( this, ix, x, q )

  import t_species_cyl_modes , p_k_part

  class(t_species_cyl_modes), intent(inout) :: this
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x
  real(p_k_part), intent(in) :: q

  end subroutine
end interface

interface
  subroutine create_particle_single_cell_p_cyl( this, ix, x, p, q )

  import t_species_cyl_modes, p_k_part

  class(t_species_cyl_modes), intent(inout) :: this
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x
  real(p_k_part), dimension(:), intent(in) :: p
  real(p_k_part), intent(in) :: q

  end subroutine
end interface

interface
  subroutine position_single_4_cyl_modes( this, idx, pos )

  import t_species_cyl_modes, p_single

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: idx
  real( p_single ), dimension(:), intent(out) :: pos

  end subroutine
end interface

interface
  subroutine position_single_8_cyl_modes( this, idx, pos )

  import t_species_cyl_modes, p_double

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: idx
  real( p_double ), dimension(:), intent(out) :: pos

  end subroutine
end interface

interface
  subroutine position_range_comp_4_cyl_modes( this, comp, idx0, idx1, pos )

  import t_species_cyl_modes, p_single

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp, idx0, idx1
  real( p_single ), dimension(:), intent(out) :: pos

  end subroutine
end interface

interface
  subroutine position_range_comp_8_cyl_modes( this, comp, idx0, idx1, pos )

  import t_species_cyl_modes, p_double

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp, idx0, idx1
  real( p_double ), dimension(:), intent(out) :: pos

  end subroutine
end interface

interface
  subroutine position_idx_comp_4_cyl_modes( this, comp, idx, np, pos )

  import t_species_cyl_modes, p_single

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_single ), dimension(:), intent(out) :: pos

  end subroutine
end interface

interface
  subroutine position_idx_comp_8_cyl_modes( this, comp, idx, np, pos )

  import t_species_cyl_modes, p_double

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_double ), dimension(:), intent(out) :: pos

  end subroutine
end interface

interface
  subroutine position_ref_box_idx_comp_4_cyl_modes( this, comp, idx, np, pos, if_pos_ref_box )

  import t_species_cyl_modes, p_single

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_single ), dimension(:), intent(out) :: pos
  logical, intent(in) :: if_pos_ref_box

  end subroutine
end interface

interface
  subroutine position_ref_box_idx_comp_8_cyl_modes( this, comp, idx, np, pos, if_pos_ref_box )

  import t_species_cyl_modes, p_double

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_double ), dimension(:), intent(out) :: pos
  logical, intent(in) :: if_pos_ref_box

  end subroutine
end interface

interface
  subroutine report_species_cyl_modes( this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )

  import t_species_cyl_modes, t_emf, t_space, t_grid, t_node_conf, t_time_step, p_double, t_vdf_msg

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_emf ), intent(inout) :: emf
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t, tmin
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine enumerate_quants_spec_cyl_modes( this, diagFile, track_set )

  import t_species_cyl_modes, t_diag_file, t_track_set

  class( t_species_cyl_modes ), intent(in) :: this
  class( t_diag_file ),intent(inout) :: diagFile
  type( t_track_set ), intent(inout) :: track_set

  end subroutine
end interface

interface
  INTEGER pure function get_n_x_dims_cyl_modes( this )

  import t_species_cyl_modes

  class( t_species_cyl_modes ), intent(in) :: this

  end function
end interface

interface
  function get_energy_spec_cyl_modes( this )

  import t_species_cyl_modes, p_double

  class( t_species_cyl_modes ), intent(in) :: this
  real( p_double ) :: get_energy_spec_cyl_modes

  end function
end interface


contains

! The below functions are defined here to (1) include the interface statement to the
! superclass init_species function and (2) to put the charge normalization routines
! near one another

!-----------------------------------------------------------------------------------------
subroutine init_species_cyl_modes( this, sp_id, interpolation, grid_center, grid, &
                          g_space, emf, jay, no_co, send_vdf, recv_vdf, bnd_cross, &
                          node_cross, send_spec, recv_spec, ndump_fac, restart, &
                          restart_handle, sim_options, tstep, tmin, tmax )
!-----------------------------------------------------------------------------------------
! sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species_cyl_modes ), intent(inout) :: this

  integer, intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in) :: grid
  type( t_space ), intent(in) :: g_space
  class( t_emf ), intent(inout) :: emf
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_vdf, recv_vdf
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_spec, recv_spec
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type(t_options), intent(in) :: sim_options
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax

  ! number of particles per cell
  integer :: ppc2
  real(p_k_part) :: ppc2_r

  ! Initialize parameters for charge correction near the axis
  ! The number of coefficients for one particle varies based on interpolation, e.g.,
  ! q ~ r * ( alpha(1) + alpha(2) * r + alpha(3) * r**2 ) for cubic interpolation.
  ! Cubic interpolation is the maximum for which the coefficients are derived,
  ! but we have 4 alpha_const parameters to at least run (without correction) for quartic.
  ! Also, coefficients are different with only 1 ppc in the r direction. We technically
  ! need coefficients for each cell, but we only find the first few (max 6) because
  ! they approach 1 quickly, and because we only have to modify the node nearest the axis.

  ! These coefficients were derived by Dan Gordon and Kyle Miller, assuming that the
  ! charge at each cell is equal.

  ppc2 = this % num_par_x(2)
  ppc2_r = real( ppc2, p_k_part )

  select case( interpolation )
  case(p_linear)
    ! if ppc2 > 1, no singularity, so we always calculate it
    if ( mod( ppc2, 2 ) == 0 ) then
      this%alpha_const(1) = 2.0 / 3.0 - 2.0 / (3.0*ppc2_r**2)
    else
      this%alpha_const(1) = 2.0 / 3.0
    endif

    this % dx2_thresh = emf%e%dx(p_r_dim) / 2

    ! Linear interpolation with 1 ppc is exactly correct already
    this % n_alpha = 0

  case(p_quadratic)
    ! avoid singularity
    if ( ppc2 > 1 ) then
      this%alpha_const(1) = ( 7.0 - 22.0*ppc2_r**2 ) / ( 24.0*ppc2_r**2 * ( ppc2_r**2 - 1.0 ) )
      this%alpha_const(2) = 5.0 * ( 4.0*ppc2_r**4 - 1.0 ) / ( 4.0 * ( 4.0*ppc2_r**4 - 5.0*ppc2_r**2 + 1.0 ) )
    endif

    this % dx2_thresh = emf%e%dx(p_r_dim)

    ! if ppc2 == 1
    this % n_alpha = 4
    this % alpha_single(1) = 985.0/1393.0
    this % alpha_single(2) = 607.0/597.0
    this % alpha_single(3) = 6953.0/6965.0
    this % alpha_single(4) = 9753.0/9751.0

  case(p_cubic)
    ! avoid singularity
    if ( ppc2 > 1 ) then
      if ( mod( ppc2, 2 ) == 0 ) then
        this%alpha_const(1) = ( 8680.0 - 294460.0*ppc2_r**2 - 6.0*ppc2_r**4 * (-107514.0 + &
                          134135.0*ppc2_r**2 + 367533.0*ppc2_r**4 - 895320.0*ppc2_r**6 + &
                          450036.0*ppc2_r**8))/(27.0*ppc2_r**2 * (-2.0 + 3.0*ppc2_r) * (2.0 + &
                          3.0*ppc2_r) * (3232.0 - 2236.0*ppc2_r**2 + 15225.0*ppc2_r**4 - &
                          26199.0*ppc2_r**6 + 9558.0*ppc2_r**8))
        this%alpha_const(2) = 2.0 * (-3496.0 + 22114.0*ppc2_r**2 - 25518.0*ppc2_r**4 + &
                          605937.0*ppc2_r**6 - 1151361.0*ppc2_r**8 + 547074.0*ppc2_r**10)/ &
                          (3.0 * (-12928.0 + 38032.0*ppc2_r**2 - 81024.0*ppc2_r**4 + &
                          241821.0*ppc2_r**6 - 274023.0*ppc2_r**8 + 86022.0*ppc2_r**10))
        this%alpha_const(3) = -(56.0 * (-196.0 + 2926.0*ppc2_r**2 + 5997.0*ppc2_r**4 + &
                          47673.0*ppc2_r**6 - 137376.0*ppc2_r**8 + 80676.0*ppc2_r**10))/ &
                          (27.0 * (-12928.0 + 38032.0*ppc2_r**2 - 81024.0*ppc2_r**4 + &
                          241821.0*ppc2_r**6 - 274023.0*ppc2_r**8 + 86022.0*ppc2_r**10))
      else
        this%alpha_const(1) = (2996.0 - 46310.0*ppc2_r**2 + 138738.0*ppc2_r**4 + 263520.0*ppc2_r**6 &
                          - 300024.0*ppc2_r**8)/(27.0 * (341.0 - 6457.0*ppc2_r**2 + &
                          34601.0*ppc2_r**4 - 38043.0*ppc2_r**6 + 9558.0*ppc2_r**8))
        this%alpha_const(2) = (2.0 * (-339.0 - 11847.0*ppc2_r**2 + 169111.0*ppc2_r**4 - &
                          346023.0*ppc2_r**6 + 182358.0*ppc2_r**8))/(9.0 * (1.0 - &
                          10.0*ppc2_r**2 + 9.0*ppc2_r**4) * (341.0 - 3047.0*ppc2_r**2 + &
                          1062.0*ppc2_r**4))
        this%alpha_const(3) = -(56.0*ppc2_r**2 * (301.0 + 4601.0*ppc2_r**2 - 13536.0*ppc2_r**4 + &
                          8964.0*ppc2_r**6))/(27.0 * (1.0 - 10.0*ppc2_r**2 + 9.0*ppc2_r**4) * &
                          (341.0 - 3047.0*ppc2_r**2 + 1062.0*ppc2_r**4))
      endif
    endif

    this % dx2_thresh = 3 * emf%e%dx(p_r_dim) / 2

    ! if ppc2 == 1
    this % n_alpha = 6
    this % alpha_single(1) = 10121.0/11087.0
    this % alpha_single(2) = 11592.0/11087.0
    this % alpha_single(3) = 32249.0/33261.0
    this % alpha_single(4) = 11340.0/11087.0
    this % alpha_single(5) = 10885.0/11087.0
    this % alpha_single(6) = 11248.0/11087.0

  ! If desired, the below code block could be enabled to at least allow the code
  ! to run with quartic interpolation. The charge near the axis will just be
  ! slightly too large.
  ! case(p_quartic)
  ! ! if ppc2 > 1
  ! this%alpha_const(1) = 1.0_p_k_part
  ! this%alpha_const(2) = 0.0_p_k_part
  ! this%alpha_const(3) = 0.0_p_k_part
  ! this%alpha_const(4) = 0.0_p_k_part

  ! this % dx2_thresh = 2 * emf%e%dx(p_r_dim)

  ! ! if ppc2 == 1
  ! this % n_alpha = 4
  ! this % alpha_single(1) = 1.0_p_k_part
  ! this % alpha_single(2) = 1.0_p_k_part
  ! this % alpha_single(3) = 1.0_p_k_part
  ! this % alpha_single(4) = 1.0_p_k_part

  case default
    write(err_buf__,*) 'Not a valid interpolation order';call err__("cyl_modes/os-spec-cyl-modes-define.f03",572)
    call abort_program(p_err_notimplemented)
  end select

  ! Call superclass init
  call init_species( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
                      no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
                      recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, &
                      tmin, tmax )

end subroutine init_species_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function norm_charge_axis_spec_cyl_modes( this, ix2, r, ppc2_in )
!-----------------------------------------------------------------------------------------
! returns a coefficient that correctly normalizes injected charge near the axis
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species_cyl_modes ), intent(in) :: this
  integer, intent(in) :: ix2 ! integer position of particle
  real(p_k_part), intent(in) :: r ! absolute radial position (not within cell)
  integer, intent(in), optional :: ppc2_in ! force calculation based on specified ppc2
  real(p_k_part) :: norm_charge_axis_spec_cyl_modes

  ! number of particles per cell
  integer :: ppc2, j
  real(p_k_part) :: alpha
  real(p_double) :: dr

  if (present(ppc2_in)) then
    ppc2 = ppc2_in
    if ( ppc2_in /= 1 ) then
      write(err_buf__,*) 'Specifying unique ppc2 /= 1 is not yet supported.';call err__("cyl_modes/os-spec-cyl-modes-define.f03",607)
      write(err_buf__,*) 'alpha_const coefficients would need to be recalculated.';call err__("cyl_modes/os-spec-cyl-modes-define.f03",608)
      write(err_buf__,*) 'aborting...';call err__("cyl_modes/os-spec-cyl-modes-define.f03",609)
      call abort_program(p_err_notimplemented)
    endif
  else
    ppc2 = this % num_par_x(p_r_dim)
  endif

  dr = this % dx(p_r_dim)

  alpha = 1.0_p_k_part
  if ( ppc2 == 1 ) then
    if ( (r <= this % n_alpha * dr) .and. (r > 0.0_p_k_part) ) then
      alpha = this%alpha_single(ix2-1)
    endif
  else
    if ( (r <= this % dx2_thresh ) .and. (r > 0.0_p_k_part) ) then
      alpha = 0.0_p_k_part
      do j = 1, this % interpolation
        alpha = alpha + this%alpha_const(j) * ( r / dr ) ** (j-1)
      enddo
    endif
  endif

  norm_charge_axis_spec_cyl_modes = alpha

end function norm_charge_axis_spec_cyl_modes
!-----------------------------------------------------------------------------------------


end module m_species_cyl_modes_define
