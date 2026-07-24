! module defining the PGC module for advancing the envelope
#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_pgc_define

#include "memory/memory.h"

use m_system
use m_parameters

use m_node_conf,   only : t_node_conf
use m_input_file,  only : t_input_file
use m_restart,     only : t_restart_handle
use m_time_step,   only : t_time_step

use m_grid_define, only : t_grid
use m_space,       only : t_space

use m_vdf_define, only : t_vdf
use m_cmplx_vdf,  only : t_cmplx_vdf
use m_vdf_smooth, only : t_smooth

use m_emf_define,     only : t_emf
use m_current_define, only : t_current

use m_tridiag, only : t_tridiag
use m_vdf_comm, only : t_vdf_msg

implicit none

private

! IBM XL compilers
! This has to be explicitly made public so that you can access the superclass
! of t_emf_pgc
public :: t_emf

! string to id restart data
character(len=*), parameter :: p_emf_pgc_rst_id = "emf_pgc rst data - 0x0001"

!-------------------------------------------------------------------------------
! parameters for pulse profiles
!-------------------------------------------------------------------------------
! longitudinal profiles
integer, parameter :: p_pgc_lon_gaussian       = 1   ! gaussian
integer, parameter :: p_pgc_lon_polynomial     = 2   ! polynomial
integer, parameter :: p_pgc_lon_sin2           = 3   ! sin^2
! perpendicular profiles
integer, parameter :: p_pgc_per_plane          = 4   ! plane wave
integer, parameter :: p_pgc_per_gaussian       = 5   ! gaussian
integer, parameter :: p_pgc_per_gaussian_asym  = 6   ! gaussian assymetric

!-------------------------------------------------------------------------------
! boundary mask for longitudinal envelope
! - p_env_bnd_overlap - condition for envelope being intialized below this value
! - p_mask_length - length of the applied boundary mask
! - p_mask_values - mask values array with length of `p_mask_length`
!-------------------------------------------------------------------------------
real, parameter :: p_env_bnd_overlap = 1e-3
integer, parameter :: p_mask_length = 11
real(p_double), parameter, dimension(p_mask_length) :: p_mask_values = &
  (/ &
    0.00000000000000d+00, 2.44717418524232d-02, 9.54915028125263d-02, &
    2.06107373853763d-01, 3.45491502812526d-01, 5.00000000000000d-01, &
    6.54508497187474d-01, 7.93892626146236d-01, 9.04508497187474d-01, &
    9.75528258147577d-01, 1.00000000000000d+00                        &
  /)

!-------------------------------------------------------------------------------
! t_emf_pgc class definition
!-------------------------------------------------------------------------------
! field components constants - representing different time steps
integer, parameter :: p_env_nm = 1    ! t = n-1/2
integer, parameter :: p_env_n  = 2    ! t = n
integer, parameter :: p_env_np = 3    ! t = n+1/2
!
type, extends( t_emf ) :: t_emf_pgc

  !> pulse evolution variables <
  ! laser envelope
  type( t_cmplx_vdf ) :: env
  ! laser envelope squared
  type( t_vdf ) :: sqrt_nm
  type( t_vdf ) :: sqrt_n
  type( t_vdf ) :: sqrt_np
  type( t_vdf ) :: env_mod
  ! ponderomotive force parameters
  type( t_vdf ) :: f_np   ! fp: n+1/2
  type( t_vdf ) :: f_n    ! fp: n
  ! coupling parameter
  type( t_vdf ), dimension(:), pointer :: chi => null()
  ! boundary variables
  logical :: lb_apply, ub_apply
  integer :: bnd_lb_low, bnd_lb_upp, bnd_ub_low, bnd_ub_upp
  ! boundary mask values
  real(p_double), pointer, dimension(:) :: bd_val_lb
  real(p_double), pointer, dimension(:) :: bd_val_ub

  !> TDMA solver variables <
  ! discritization matrices
  type( t_tridiag ), allocatable, dimension(:) :: disc_mat
  ! solver type for matrix
  integer :: solver_type
  ! error cleaning for tridiagonal solver
  logical :: error_cleaning
  logical :: pre_calc_smoothing, post_calc_smoothing
  ! source term - right-hand-side
  type( t_cmplx_vdf ) :: source

  !> parameters for different laser types <
  ! general laser envelope variables
  real(p_k_fld) :: a0               ! laser amplitude
  real(p_k_fld) :: omega            ! laser frequency
  logical       :: free_stream      ! free steaming laser option
  integer       :: per_type         ! perpendicular laser profile type
  integer       :: lon_type         ! longitudinal laser profile type

  ! longitudinal laser parameters
  real(p_k_fld) :: lon_center       ! central position of pulse (gauss/sin2)
  real(p_k_fld) :: lon_duration     ! FWHM for gaussian envelope
  real(p_k_fld) :: lon_range        ! full length of gaussian pulse
  real(p_k_fld) :: lon_rise         ! length of pulse rise region (poly/sin2)
  real(p_k_fld) :: lon_fall         ! length of pulse fall region (poly/sin2)
  real(p_k_fld) :: lon_flat         ! length of flat region (poly/sin2)
  real(p_k_fld) :: lon_start        ! begining of laser pulse (poly)

  ! perpendicular laser parameters
  real(p_k_fld), dimension(2) :: per_center       ! position of optical axis
  real(p_k_fld)               :: w0               ! spot size at focal plane
  real(p_k_fld)               :: per_focus        ! position of focal plane
  real(p_k_fld), dimension(2) :: w0_asym          ! spot size of focal planes
  real(p_k_fld), dimension(2) :: per_focus_asym   ! positions of focal planes

  !> smoothing parameters for PGC <
  ! coupling parameter 'chi' smoothing
  logical          :: if_smooth_chi
  type( t_smooth ) :: chi_smooth
  ! ponderomotive force 'fp' smoothing
  logical          :: if_smooth_fp
  integer          :: fp_smooth_niter, fp_smooth_nmax
  type( t_smooth ) :: fp_smooth
  ! envelope smoothing
  logical          :: if_smooth_env
  type( t_smooth ) :: env_smooth

  ! Parallel configuration
  class( t_node_conf ), pointer :: no_co => null()

contains

    procedure :: read_input           => read_input_pgc
    procedure :: init                 => init_pgc
    procedure :: advance              => advance_pgc
    procedure :: report               => report_diag_pgc
    procedure :: cleanup              => cleanup_pgc
    procedure :: update_particle_fld  => update_particle_fld_pgc
    procedure :: write_checkpoint     => write_checkpoint_pgc
    procedure :: allocate_objs        => allocate_objs_pgc

    ! local methods
    procedure :: reset_chi => reset_chi_pgc

end type t_emf_pgc

!-------------------------------------------------------------------------------
! interfaces for t_emf_pgc
!-------------------------------------------------------------------------------

interface
  subroutine read_input_pgc( this, input_file, periodic, if_move, grid, dx, dt, gamma )
    import :: t_emf_pgc, t_input_file, t_grid, p_double
    class( t_emf_pgc ), intent(inout)   :: this
    class( t_input_file ), intent(inout) :: input_file
    logical, dimension(:), intent(in)   :: periodic, if_move
    class( t_grid ), intent(in)         :: grid
    real(p_double), dimension(:), intent(in) :: dx
    real(p_double), intent(in) :: dt
    real(p_double), intent(in) :: gamma 
  end subroutine
end interface

interface
  subroutine init_pgc( this, part_grid_center, part_interpolation, &
                       g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                       no_co, send_msg, recv_msg, restart, restart_handle, sim_options )
    import :: t_emf_pgc, t_space, t_grid, p_double, t_node_conf, t_vdf_msg
    import :: t_restart_handle, t_options, t_time_step
    class( t_emf_pgc ), intent( inout ), target :: this
    logical, intent(in)                         :: part_grid_center
    integer, intent(in)                         :: part_interpolation
    type( t_space ), intent(in)                 :: g_space
    class( t_grid ), intent(in)                 :: grid
    integer, dimension(:,:), intent(in)         :: gc_min
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: tmin, tmax
    real( p_double ), dimension(:), intent(in)  :: dx
    class( t_node_conf ), intent(in), target    :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    logical, intent(in)                         :: restart
    type( t_restart_handle ), intent(in)        :: restart_handle
    type( t_options ), intent(in)               :: sim_options
  end subroutine
end interface

interface
  subroutine advance_pgc( this, jay, dt, send_msg, recv_msg )
    import :: t_emf_pgc, t_current, p_double, t_vdf_msg
    class( t_emf_pgc ), intent( inout ) ::  this
    class( t_current ), intent(inout)   :: jay
    real( p_double ), intent(in)        :: dt
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  end subroutine
end interface

interface
  subroutine report_diag_pgc( this, g_space, grid, no_co, tstep, t, send_msg, recv_msg )
    import :: t_emf_pgc, t_space, t_grid, t_node_conf, t_time_step, p_double, t_vdf_msg
    class( t_emf_pgc ), intent(inout) :: this
    type( t_space ), intent(in)       :: g_space
    class( t_grid ), intent(in)       :: grid
    class( t_node_conf ), intent(in)   :: no_co
    type( t_time_step ), intent(in)   :: tstep
    real( p_double ), intent(in)      :: t
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  end subroutine
end interface

interface
  subroutine cleanup_pgc( this )
    import :: t_emf_pgc
    class( t_emf_pgc ), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine update_particle_fld_pgc( this, n, t, dt, space, nx_p_min )
    import :: t_emf_pgc, p_double, t_space
    class( t_emf_pgc ), intent( inout ) :: this
    integer, intent(in)                 :: n
    real( p_double ), intent(in)        :: t, dt
    type( t_space ), intent(in)         :: space
    integer, dimension(:), intent(in)   :: nx_p_min
  end subroutine
end interface

interface
  subroutine write_checkpoint_pgc( this, restart_handle )
    import :: t_emf_pgc, t_restart_handle
    class( t_emf_pgc ), intent(in)          :: this
    type( t_restart_handle ), intent(inout) :: restart_handle
  end subroutine
end interface

interface
  subroutine allocate_objs_pgc( this )
    import t_emf_pgc
    class( t_emf_pgc ), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine reset_chi_pgc( this )
    import :: t_emf_pgc
    class( t_emf_pgc ), intent(inout) :: this
  end subroutine
end interface


!-------------------------------------------------------------------------------
! exported symbols
!-------------------------------------------------------------------------------
public :: p_emf_pgc_rst_id
public :: t_emf_pgc, p_env_nm, p_env_n, p_env_np
public :: p_env_bnd_overlap, p_mask_length, p_mask_values
public :: p_pgc_lon_gaussian, p_pgc_lon_sin2, p_pgc_lon_polynomial
public :: p_pgc_per_plane, p_pgc_per_gaussian, p_pgc_per_gaussian_asym

end module m_emf_pgc_define
