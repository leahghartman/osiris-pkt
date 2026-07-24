# 1 "os-bnd.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "os-bnd.f03" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary class
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
# 6 "os-bnd.f03" 2
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
# 7 "os-bnd.f03" 2

module m_bnd

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
# 11 "os-bnd.f03" 2

use m_system
use m_parameters
use m_species_define, only : t_part_idx, t_spec_msg
use m_vdf_comm
use m_species_comm

implicit none

private

! In each t_sim_tiles
type :: t_bnd

  ! indexes of particles crossing lower/upper boundary
  type( t_part_idx ), dimension(2) :: bnd_cross
  ! indexes of particles leaving the node
  type( t_part_idx ) :: node_cross
  ! Replaces module variables for vdf comm
  type( t_vdf_msg ), dimension(2) :: send_vdf, recv_vdf
  ! Replaces module variables for part comm
  type( t_spec_msg ), dimension(2) :: send_spec, recv_spec

  contains

  procedure :: init => init_bnd
  procedure :: cleanup => cleanup_bnd

end type t_bnd

public :: t_bnd

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_bnd( this, buffer_size )

  implicit none

  class( t_bnd ), intent(inout) :: this
  integer, intent(in), optional :: buffer_size

  if (present(buffer_size)) then
    call init_spec_comm( this%send_spec, this%recv_spec, buffer_size )
  else
    call init_spec_comm( this%send_spec, this%recv_spec )
  endif

end subroutine init_bnd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cleanup_bnd( this )

  implicit none

  class( t_bnd ), intent(inout) :: this

  integer :: i

  call freemem(this%node_cross % idx,"os-bnd.f03",73)
  do i = p_lower, p_upper
    call freemem(this%bnd_cross( i ) % idx,"os-bnd.f03",75)
  enddo

  call cleanup( this%send_spec )
  call cleanup( this%recv_spec )
  call cleanup( this%send_vdf )
  call cleanup( this%recv_vdf )


end subroutine cleanup_bnd
!-----------------------------------------------------------------------------------------

end module m_bnd
