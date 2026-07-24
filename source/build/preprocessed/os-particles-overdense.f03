# 1 "overdense/os-particles-overdense.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "overdense/os-particles-overdense.f03" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Particles module for particle damper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
# 6 "overdense/os-particles-overdense.f03" 2
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
# 7 "overdense/os-particles-overdense.f03" 2

module m_particles_overdense

use m_system
use m_parameters
use m_particles_define, only : t_particles, pushev
use m_spec_overdense
use m_spec_splitting, only : if_split, t_recombine
use m_species_define, only : t_species
use m_species_memory
use m_emf_define, only : t_emf
use m_time_step, only : t_time_step, n
use m_current_define, only : t_current
use m_logprof, only : begin_event, end_event, create_event
use m_node_conf, only : t_node_conf
use m_grid_define, only : t_grid
use m_restart, only : t_restart_handle
use m_space, only : t_space
use m_bnd, only : t_bnd
use m_zpulse, only : t_zpulse_list

implicit none

private

type, extends( t_particles ) :: t_particles_overdense

contains

  procedure :: allocate_objs => allocate_objs_part_overdense
  procedure :: init => init_part_overdense
  procedure :: advance_deposit => advance_deposit_part_overdense

end type t_particles_overdense

public :: t_particles_overdense

integer, public :: split_ev = 0, recombine_ev, prep_damper_ev

contains

!-----------------------------------------------------------------------------------------
! Allocate particle species, cathodes, etc.
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_overdense( this )

  implicit none

  class( t_particles_overdense ), intent(inout) :: this

  integer :: i
  class(t_species), pointer :: spec, tail


  if ( this%num_species > 0 ) then

    allocate( t_species_overdense :: spec )
    this % species => spec
    do i = 2, this%num_species
      tail => spec
      allocate( t_species_overdense :: spec )
      tail % next => spec
    enddo

    if ( this%num_cathode > 0 ) then
      allocate( this%cathode( this%num_cathode ) )
    endif


    if ( this%num_neutral > 0 ) then
      allocate( this%neutral( this%num_neutral ) )
    endif


  endif

end subroutine allocate_objs_part_overdense
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Sets up this data structure from the given information
!-----------------------------------------------------------------------------------------
subroutine init_part_overdense( this, g_space, jay, emf, grid, no_co, bnd, zpulse_list, &
                                ndump_fac, restart, restart_handle, t, tstep, tmin, tmax, &
                                sim_options )

  implicit none

  class( t_particles_overdense ), intent(inout) :: this

  type( t_space ), intent(in) :: g_space
  class( t_emf) ,intent(inout) :: emf
  class( t_current ) ,intent(inout) :: jay
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  class( t_bnd ), intent(inout) :: bnd
  class( t_zpulse_list ), intent(in) :: zpulse_list
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  real(p_double), intent(in) :: t
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax
  type( t_options ), intent(in) :: sim_options

  call this % t_particles % init( g_space, jay, emf, grid, no_co, bnd, zpulse_list, &
                                  ndump_fac, restart, restart_handle, t, tstep, tmin, &
                                  tmax, sim_options )

  ! setup events in this file
  if (split_ev==0) then
    split_ev = create_event('particle splitting')
    recombine_ev = create_event('particle recombination')
    prep_damper_ev = create_event('particle damping preparation')
  endif

end subroutine init_part_overdense

!-----------------------------------------------------------------------------------------
! Prepare fluid moment arrays if damping
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_part_overdense( this, emf, jay, tstep, t, no_co, options )

  implicit none

  class( t_particles_overdense ), intent(inout) :: this
  class( t_emf ), intent( inout ) :: emf
  class( t_current ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  class(t_species), pointer :: species
  type(t_recombine), pointer :: recombine

  ! prepare the damper if required
  call begin_event( pushev )
  species => this % species
  do
    if (.not. associated(species)) exit
    select type(species)
    class is (t_species_overdense)
      if ( t >= species%push_start_time ) then

        call begin_event( prep_damper_ev )
        call species % prep_damper( n(tstep), t )
        call end_event( prep_damper_ev )

      endif
    end select
    species => species % next
  enddo
  call end_event( pushev )

  call this % t_particles % advance_deposit( emf, jay, tstep, t, no_co, options )

  ! split and recombine particles if required
  call begin_event( pushev )
  species => this % species
  do
    if (.not. associated(species)) exit
    select type(species)
    class is (t_species_overdense)
      if ( t >= species%push_start_time ) then

        call begin_event( split_ev )
        call species % split( n(tstep) )
        call end_event( split_ev )

        call begin_event( recombine_ev )
        recombine => species%splitter%recombine
        do
          if (.not. associated(recombine)) exit
          call species % recombine( recombine, n(tstep), t )
          recombine => recombine%next
        enddo
        call end_event( recombine_ev )

      endif
    end select
    species => species % next
  enddo
  call end_event( pushev )

end subroutine advance_deposit_part_overdense
!-----------------------------------------------------------------------------------------

end module m_particles_overdense
