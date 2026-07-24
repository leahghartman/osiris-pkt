# 1 "spec/os-spec-initfields.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/os-spec-initfields.f90" 2
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
# 2 "spec/os-spec-initfields.f90" 2
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
# 3 "spec/os-spec-initfields.f90" 2

module m_species_initfields

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
# 7 "spec/os-spec-initfields.f90" 2

use m_system
use m_parameters

use m_species_define, only: t_species
use m_emf_define, only: t_emf
use m_grid_define, only: t_grid
use m_node_conf, only: t_node_conf
use m_vdf_define, only: t_vdf
use m_vdf_comm, only: update_boundary, t_vdf_msg
use m_vdf_math, only: add
use m_species_charge, only: deposit_rho
use m_emf_es_solver, only: es_solver, es_solver_beam

use m_vpml, only: update_beam

implicit none

private

interface init_fields
  module procedure init_fields_species
end interface

public :: init_fields

contains

!-----------------------------------------------------------------------------------------
! Initialize fields from initial charge / momentum distribution
!-----------------------------------------------------------------------------------------
subroutine init_fields_species( this, emf, no_co, grid, send_msg, recv_msg )

  implicit none

  class( t_species ), intent(in) :: this
  class( t_emf ), intent(inout) :: emf
  class( t_node_conf ), intent(in) :: no_co
  class( t_grid ), intent(in) :: grid
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  type( t_vdf ) :: rho, ef, bf

  if (mpi_node()==0) print *,' - Initializing EM fields for species "', trim(this%name),'"'

  ! Get charge density from species

  ! Create new charge grid with same dimensions as E field grid and 1 field component
  call rho%new( emf%e, f_dim = 1 )

  ! zero the array
  call rho%zero( )

  ! Deposit the charge density

  call deposit_rho( this, rho )

  ! Get contributions from neighboring parallel nodes
  call update_boundary( rho, p_vdf_add, no_co, send_msg, recv_msg )

  ! Find E field from initial charge distribution

  ! Create temporary E field grid
  call ef%new( emf%e )

  if ( this%udist%ufl(1) == 0.0 ) then

     ! Call ES solver to find E field from charge density
     call es_solver( rho, ef, grid, no_co, send_msg, recv_msg )

     ! Add ES field to global simulation field
     call add( emf%e, ef )

  else

     call bf % new( emf%b )

     ! Call ES solver to find E field from beam charge density
     call es_solver_beam( rho, ef, bf, this%udist%ufl(1), grid, no_co, send_msg, recv_msg )

     ! Add static fields to global simulation field
     call add( emf%e, ef )
     call add( emf%b, bf )

     ! Update PML boundaries, if any
     call update_beam( emf%bnd_con%vpml_all, ef, bf, no_co, send_msg, recv_msg )

     ! cleanup
     call bf % cleanup()

  endif

  ! cleanup
  call ef % cleanup()
  call rho % cleanup()

end subroutine

end module m_species_initfields
