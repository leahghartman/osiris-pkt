# 1 "zpulse/os-zpulse.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "zpulse/os-zpulse.f03" 2
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
# 2 "zpulse/os-zpulse.f03" 2
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
# 3 "zpulse/os-zpulse.f03" 2

module m_zpulse

use m_system

use m_zpulse_std
use m_zpulse_wall
use m_zpulse_mov_wall
use m_zpulse_point
use m_zpulse_flyfoc
use m_zpulse_speckle
use m_zpulse_planewaves

use m_input_file

use m_space
use m_emf_define, only: t_emf_bound, t_emf
use m_grid_define
use m_node_conf


implicit none

private

!-----------------------------------------------------------------------------------------
type :: t_zpulse_list

  class( t_zpulse ), pointer :: list => null()

contains

  procedure :: read_input => read_input_zpulse_list
  procedure :: init => init_zpulse_list
  procedure :: cleanup => cleanup_zpulse_list

  procedure :: launch => launch_zpulse_list

end type t_zpulse_list
!-----------------------------------------------------------------------------------------


public :: t_zpulse_list

contains


!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_list( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  implicit none


  class( t_zpulse_list ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  class( t_zpulse ), pointer :: tail, item

  tail => null()

  do
    select case (trim(input_file % get_section_name()))
    case('nl_zpulse')
        allocate( t_zpulse :: item )
    case('nl_zpulse_wall')
        allocate( t_zpulse_wall :: item )
    case('nl_zpulse_mov_wall')
        allocate( t_zpulse_mov_wall :: item )
    case('nl_zpulse_point')
        allocate( t_zpulse_point :: item )
    case('nl_zpulse_flyfoc')
        allocate( t_zpulse_flyfoc :: item )
    case('nl_zpulse_speckle')
        allocate( t_zpulse_speckle :: item )
    case('nl_zpulse_planewaves')
        allocate( t_zpulse_planewaves :: item )
    case default
        ! Not a zpulse section, exit
        exit
    end select

    ! Read input for zpulse
    call item % read_input( input_file, g_space, bnd_con, periodic, grid, sim_options )
    call item % check_dimensionality()
    item % next => null()

    if (associated(tail)) then
      tail % next => item
    else
      this % list => item
    endif

    tail => item

  enddo

end subroutine read_input_zpulse_list
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine init_zpulse_list( this, restart, t, dt, emf, g_space, nx_p_min, no_co, &
                             interpolation )

  implicit none

  class( t_zpulse_list ), intent(inout) :: this
  logical, intent(in) :: restart
  real( p_double ), intent(in) :: t, dt
  class( t_emf ), intent(in) :: emf
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: interpolation

  class( t_zpulse ), pointer :: zpulse

  zpulse => this%list

  do
    if ( .not. associated(zpulse)) exit
    call zpulse % init( restart, t, dt, emf%b, g_space, nx_p_min, no_co, interpolation )

    zpulse => zpulse % next
  enddo

end subroutine init_zpulse_list
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine cleanup_zpulse_list( this )

  implicit none

  class( t_zpulse_list ), intent(inout) :: this
  class( t_zpulse ), pointer :: zpulse, next

  zpulse => this % list

  do
    if ( .not. associated(zpulse)) exit
    next => zpulse % next
    call zpulse % cleanup()
    deallocate( zpulse )
    zpulse => next
  enddo

end subroutine cleanup_zpulse_list
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_list( this, emf, g_space, &
                                nx_p_min, g_nx, t, dt, no_co )

  class( t_zpulse_list ), intent(inout) :: this
  class( t_emf ) , intent(inout) :: emf

  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  class( t_zpulse ), pointer :: zpulse


  zpulse => this % list

  do
    if (.not. associated(zpulse)) exit
    call zpulse % launch( emf, g_space, nx_p_min, g_nx, t, dt, no_co )
    zpulse => zpulse % next
  enddo

end subroutine launch_zpulse_list
!-----------------------------------------------------------------------------------------


end module m_zpulse
