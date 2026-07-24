# 1 "pgc/os-particles-pgc.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "pgc/os-particles-pgc.f03" 2
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
# 2 "pgc/os-particles-pgc.f03" 2
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
# 3 "pgc/os-particles-pgc.f03" 2

module m_particles_pgc

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
# 7 "pgc/os-particles-pgc.f03" 2

use m_system
use m_parameters
use m_space
use m_grid_define
use m_node_conf
use m_restart
use m_time_step
use m_input_file

use m_particles_define
use m_species_define
use m_current_define

implicit none

private

type, extends( t_particles ) :: t_particles_pgc

  ! no additional class members
  ! all species must be of type species_pgc
  ! all neutrals must be of type neutral_pgc

contains

  procedure :: advance_deposit => advance_deposit_pgc
  procedure :: allocate_objs => allocate_objs_part_pgc

end type t_particles_pgc

public :: t_particles_pgc

contains

!-----------------------------------------------------------------------------------------
! Overrides allocate_objs method and allocates t_species_pgc and t_neutral_pgc objects
! instead
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_pgc( this )

  use m_species_define
  use m_species_pgc
  use m_neutral_pgc
  use m_cathode

  implicit none

  class( t_particles_pgc ), intent(inout) :: this

  integer :: i
  class( t_species_pgc ), pointer :: spec, tail

  if ( this % num_species > 0 ) then

    allocate(spec)
    this % species => spec
    do i = 2, this % num_species
      tail => spec
      allocate(spec)
      tail % next => spec
    enddo

    if ( this % num_cathode > 0 ) then
      allocate( this%cathode( this%num_cathode ) )
    endif


   if ( this%num_neutral > 0 ) then
      allocate( t_neutral_pgc :: this%neutral( this%num_neutral ) )
    endif


  endif

end subroutine allocate_objs_part_pgc
!-----------------------------------------------------------------------------------------

subroutine advance_deposit_pgc( this, emf, jay, tstep, t, no_co, options )

  use m_emf_define, only : t_emf
  use m_emf_pgc_define, only : t_emf_pgc
  use m_vdf_define
  use m_vdf_math
  use m_time_step
  use m_node_conf
  use m_current_define

  implicit none

  class( t_particles_pgc ), intent(inout) :: this
  class( t_emf ), intent( inout ) :: emf
  class( t_current ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep

  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  !print *, '[pgc] In advance_deposit_pgc'

  ! before calling the pusher we must reset chi

  select type( emf )

  class is ( t_emf_pgc )

 ! reset_chi() only works with t_emf_pgc objects, and since emf was declared as
 ! t_emf we must use a select type construct to check that the object is indeed a
 ! t_emf_pgc instance

    call emf % reset_chi()

  class default

    ! This must never happen, advance_deposit_pgc must always be called with a
    ! t_emf_pgc object
    call abort_program( p_err_invalid )

  end select

  ! call superclass pusher
  call this % t_particles % advance_deposit( emf, jay, tstep, t, no_co, options )

  ! reduce chi
  select type( emf )
  class is ( t_emf_pgc )

    call reduce( emf%chi, 1 )

  class default

    ! This must never happen, advance_deposit_pgc must always be called with a
    ! t_emf_pgc object
    call abort_program( p_err_invalid )

  end select

end subroutine advance_deposit_pgc

end module m_particles_pgc
