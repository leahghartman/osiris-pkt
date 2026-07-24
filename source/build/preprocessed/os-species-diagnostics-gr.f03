# 1 "gr/os-species-diagnostics-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-species-diagnostics-gr.f03" 2
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
# 2 "gr/os-species-diagnostics-gr.f03" 2
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
# 3 "gr/os-species-diagnostics-gr.f03" 2

module m_species_diagnostics_gr

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
# 7 "gr/os-species-diagnostics-gr.f03" 2

use m_system
use m_parameters

use m_species_define_gr, only : t_species_gr

implicit none

private

interface report_nemit_gr
    module procedure report_nemit_gr
end interface report_nemit_gr

public :: report_nemit_gr

contains

subroutine report_nemit_gr( this, g_space, grid, no_co, tstep, t )

  use stringutil
  use m_vdf_define, only : t_vdf_report
  use m_space, only : t_space
  use m_grid_define, only : t_grid
  use m_node_conf, only : t_node_conf
  use m_time_step

  implicit none

  class( t_species_gr ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real( p_double ) :: t

  type( t_vdf_report ) :: info

  ! Prepare the metadata
  info%name = 'nemit'
  info%ndump = ndump(tstep)

  info%xname = (/'x1', 'x2', 'x3'/)

  info%xlabel = (/'x_1', 'x_2', 'x_3'/)
  info%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

  info%time_units = '1 / \omega_p'
  info%dt = dt(tstep)

  info%label = '# emissions'
  info%units = ''

  info%n = n(tstep)
  info%t = t

  info%fileLabel = ''
  info%path = trim(path_mass) // &
              'GROUP' // p_dir_sep // &
              trim(this%name) // p_dir_sep // &
              'nemit' // p_dir_sep

  info%filename = 'nemit-' // trim(this%name) // '-' // idx_string( n(tstep)/ndump(tstep), 6 )

  ! write the file
  call this % nemit % write( info, 1, g_space, grid, no_co )

  ! reset vdf values
  call this % nemit % zero()

end subroutine report_nemit_gr

end module m_species_diagnostics_gr
