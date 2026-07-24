# 1 "qed/os-bremsstrahlung.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "qed/os-bremsstrahlung.f03" 2
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
# 2 "qed/os-bremsstrahlung.f03" 2
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
# 3 "qed/os-bremsstrahlung.f03" 2

module m_bremsstrahlung

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
# 7 "qed/os-bremsstrahlung.f03" 2
 !use m_random_class

  use m_parameters
  use m_system
  use m_node_conf
  use m_random
  use m_vdf_define
  use m_vdf
  use m_species_define


implicit none

private

! -----
! CLASS DEFINITIONS
! -----

type :: t_bremsstrahlung_info

  logical :: if_bremsstrahlung
  integer :: Z_ion
  type( t_vdf ) :: rho_ion

  ! multiplicative factor applied on the probability
  real( p_double ) :: proba_mult

  ! Type of cross section
  integer :: cs_model
  logical :: energy_damp
  real( p_double ) :: p_emit_cutoff

  ! pointer to the ion species for bremsstrahlung
  class( t_species ), pointer :: ion

  integer, dimension(:), allocatable :: i_ion

  contains

  procedure :: read_input => read_input_bremsstrahlung

end type

public :: t_bremsstrahlung_info

contains

  !-----------------------------------------------------------------------------------------
  ! Reads info for bremsstrahlung
  !-----------------------------------------------------------------------------------------

  subroutine read_input_bremsstrahlung( this, input_file, num_qed)

  use m_input_file, only : t_input_file, get_namelist

  implicit none

  class( t_bremsstrahlung_info ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: num_qed

  logical :: if_bremsstrahlung, energy_damp
  integer :: ierr, Z_ion
  integer :: cs_model
  real(p_double) :: proba_mult
  integer, dimension(:), allocatable :: i_ion

  namelist /nl_qed_bremsstrahlung/ if_bremsstrahlung, Z_ion, proba_mult, cs_model, energy_damp, i_ion

  if_bremsstrahlung = .false.
  Z_ion = -1
  proba_mult = 1.0
  cs_model = -1
  energy_damp = .true.
  allocate(i_ion(num_qed))
  i_ion = 0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_qed_bremsstrahlung", ierr )

  if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_qed_bremsstrahlung, iostat = ierr)
      if (ierr /= 0) then
          print *, ""
          print *, "   Error reading qed_bremsstrahlung parameters "
          print *, "   aborting..."
          stop
      endif
  else
      if (mpi_node()==0) print *,"   - no particle bremsstrahlung specified"
  endif

  this%if_bremsstrahlung = if_bremsstrahlung
  this%Z_ion = Z_ion
  this%proba_mult = proba_mult
  this%cs_model = cs_model
  this%energy_damp = energy_damp
  this%i_ion = i_ion

  end subroutine read_input_bremsstrahlung

end module m_bremsstrahlung
