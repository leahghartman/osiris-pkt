# 1 "spec/os-spec-tag.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/os-spec-tag.f03" 2
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
# 2 "spec/os-spec-tag.f03" 2
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
# 3 "spec/os-spec-tag.f03" 2

module m_species_tag

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
# 7 "spec/os-spec-tag.f03" 2

use m_species_define


use m_species_tracks


implicit none

private

interface set_tags
  module procedure set_tags_range
  module procedure set_tags_single
end interface

public :: set_tags

contains

!-------------------------------------------------------------------------------
subroutine set_tags_range( this, idx0, idx1 )
!-------------------------------------------------------------------------------
! Sets tags of particle range
!-------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: this
  integer, intent(in) :: idx0, idx1

  integer :: i,tag

  ! version 1 of tags

  tag = this%num_created

  do i = idx0, idx1
    tag = tag + 1
    this%tag(1,i) = this%ngp_id
    this%tag(2,i) = tag
  enddo



  ! if the tracking diagnostic is on and tags are in the missing list
  ! move them to the present list and mark the particle as being tracked
  if ( this%diag%ndump_fac_tracks > 0) then
     call new_particles( this%diag%tracks, this, idx0, idx1 )
  endif



end subroutine set_tags_range
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine set_tags_single( this, idx )
!-------------------------------------------------------------------------------
! Sets tag of single particle
!-------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: this
  integer, intent(in) :: idx

  integer :: tag

  ! version 1 of tags

  tag = this%num_created + 1
  this%tag(1,idx) = this%ngp_id
  this%tag(2,idx) = tag



  ! if the tracking diagnostic is on and tags are in the missing list
  ! move them to the present list and mark the particle as being tracked
  if ( this%diag%ndump_fac_tracks > 0) then
     call new_particles( this%diag%tracks, this%tag(:,idx), idx )
  endif



end subroutine set_tags_single
!-------------------------------------------------------------------------------

end module m_species_tag
