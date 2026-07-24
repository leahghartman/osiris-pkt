# 1 "gr/os-emf-define-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-emf-define-gr.f03" 2
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
# 2 "gr/os-emf-define-gr.f03" 2
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
# 3 "gr/os-emf-define-gr.f03" 2

module m_emf_define_gr

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
# 7 "gr/os-emf-define-gr.f03" 2

use m_system
use m_parameters

use m_node_conf, only : t_node_conf
use m_input_file, only : t_input_file

use m_grid_define, only : t_grid
use m_space, only : t_space

use m_restart, only : t_restart_handle

use m_vdf_define, only : t_vdf
use m_vdf_comm, only : t_vdf_msg
use m_emf_define, only : t_emf
use m_current_define, only : t_current
use m_geometry_gr, only : t_geometry_gr
use m_time_step, only : t_time_step


implicit none

private

! string to id restart data
character(len=*), parameter, public :: p_emf_gr_rst_id = "emf_gr rst data - 0x0006"

! IBM XL compilers
! This has to be explicitly made public so that you can access the superclass
! of t_emf_gr
public :: t_emf


!-------------------------------------------------------------------------------
! t_emf_gr class definition
!-------------------------------------------------------------------------------

type, extends( t_emf ) :: t_emf_gr

  class( t_geometry_gr ), pointer :: geometry

contains

  procedure :: read_input => read_input_gr
  procedure :: allocate_objs => allocate_objs_gr
  procedure :: init => init_gr
  procedure :: advance => advance_gr
  ! procedure :: cleanup => cleanup_pulsar
  procedure :: write_checkpoint => write_checkpoint_gr
  procedure :: restart_read => restart_read_gr

  procedure :: min_gc => min_gc_emf_gr

end type t_emf_gr

!-------------------------------------------------------------------------------
! interfaces for t_emf_gr
!-------------------------------------------------------------------------------

interface
  subroutine read_input_gr( this, input_file, periodic, if_move, &
                                     grid, dx, dt, gamma )
    import :: t_emf_gr, t_input_file, t_grid, p_double
    class( t_emf_gr ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ), intent(in) :: grid
    real(p_double), dimension(:), intent(in) :: dx
    real(p_double), intent(in) :: dt
    real(p_double), intent(in) :: gamma
  end subroutine
end interface

interface
  subroutine allocate_objs_gr( this )
    import t_emf_gr
    class( t_emf_gr ), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine init_gr( this, part_grid_center, part_interpolation, &
                       g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                       no_co, send_msg, recv_msg, restart, restart_handle, sim_options )
    import :: t_emf_gr, t_space, t_grid, p_double, t_node_conf, t_vdf_msg
    import :: t_restart_handle, t_options, t_time_step
    class( t_emf_gr ), intent( inout ), target :: this
    logical, intent(in) :: part_grid_center
    integer, intent(in) :: part_interpolation
    type( t_space ), intent(in) :: g_space
    class( t_grid ), intent(in) :: grid
    integer, dimension(:,:), intent(in) :: gc_min
    real( p_double ), dimension(:), intent(in) :: dx
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: tmin, tmax
    class( t_node_conf ), intent(in), target :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    logical, intent(in) :: restart
    type( t_restart_handle ), intent(in) :: restart_handle
    type( t_options ), intent(in) :: sim_options
  end subroutine
end interface

interface
  subroutine advance_gr( this, jay, dt, send_msg, recv_msg )
    import :: t_emf_gr, t_current, p_double, t_vdf_msg
    class( t_emf_gr ), intent( inout ) :: this
    class( t_current ), intent(inout) :: jay
    real(p_double), intent(in) :: dt
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  end subroutine
end interface

interface
  subroutine write_checkpoint_gr( this, restart_handle )
    import :: t_emf_gr, t_restart_handle
    class( t_emf_gr ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle
  end subroutine
end interface

interface
  subroutine restart_read_gr( this, restart_handle )
  import t_emf_gr, t_restart_handle
  class( t_emf_gr ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle
  end subroutine
end interface

! exported symbols
public :: t_emf_gr

contains

function min_gc_emf_gr( this )

  use m_vdf_smooth, only: smooth_order

  implicit none

  class( t_emf_gr ), intent(in) :: this
  integer, dimension(2,p_x_dim) :: min_gc_emf_gr
  integer, dimension(2,p_x_dim) :: gc_num_temp
  integer :: i

  ! number of guard cells in radial direction
  min_gc_emf_gr(p_lower,1) = 2
  min_gc_emf_gr(p_upper,1) = 2

  ! number of guard cells in polar direction
  min_gc_emf_gr(p_lower,2) = 1
  min_gc_emf_gr(p_upper,2) = 2

  ! Get number of guard cells required by smoothing and choose the maximum
  gc_num_temp(p_lower,:) = smooth_order(this%smooth)
  gc_num_temp(p_upper,:) = gc_num_temp(p_lower,:)+1

  do i = 1, p_x_dim
    min_gc_emf_gr(p_lower,i)=max(min_gc_emf_gr(p_lower,i),gc_num_temp(p_lower,i))
    min_gc_emf_gr(p_upper,i)=max(min_gc_emf_gr(p_upper,i),gc_num_temp(p_upper,i))
  enddo

end function min_gc_emf_gr

end module m_emf_define_gr
