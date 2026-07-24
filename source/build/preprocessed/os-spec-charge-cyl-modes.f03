# 1 "cyl_modes/os-spec-charge-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-spec-charge-cyl-modes.f03" 2
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
# 2 "cyl_modes/os-spec-charge-cyl-modes.f03" 2
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
# 3 "cyl_modes/os-spec-charge-cyl-modes.f03" 2

module m_species_charge_cyl_modes

use m_species_define
use m_system
use m_parameters
use m_input_file
use m_vdf_define

private

interface norm_charge_cyl_m
  module procedure norm_charge_cyl_m
end interface

public :: norm_charge_cyl_m

! -----------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
! Normalize charge in cylindrical coordinates.
!-------------------------------------------------------------------------------
subroutine norm_charge_cyl_m( rho, gir_pos, dr, cyl_m )

  implicit none

  type( t_vdf ), intent(inout) :: rho
  integer, intent(in) :: gir_pos ! local grid position on global grid
  real( p_k_fld ), intent(in) :: dr
  integer, intent(in) :: cyl_m

  integer :: ir, iz
  real( p_k_fld ) :: r, vol_fac_cor

  if ( cyl_m == 0 ) then
    vol_fac_cor = 1.0_p_k_fld
  else
    ! Modes > 0 need an extra factor of 2, akin to that in norm_ring_grid_cyl_modes
    vol_fac_cor = 2.0_p_k_fld
  endif

  ! normalize for 'ring' particles
  ! Note that the lower radial spatial boundary of a cylindrical geometry simulation is
  ! always -dr/2, where dr is the radial cell size. Turns out that charge below the r = 0
  ! boundary should always be added back in (verified with charge conservation test).
  ! The final values in the guard cells are really meaningless, so we just set them to the
  ! same value as for |r|.

  do ir = lbound( rho%f2, 3 ), ubound( rho%f2, 3 )
    ! r = ( (ir+ gir_pos - 2) - 0.5_p_k_fld )* dr
    r = ( ABS((ir+ gir_pos - 2) - 0.5_p_k_fld ))* dr
    do iz = lbound( rho%f2, 2 ), ubound( rho%f2, 2 )
      rho%f2(1,iz,ir) = rho%f2(1,iz,ir) / r * vol_fac_cor
    enddo
  enddo

  ! Fold axial guard cells back into simulation space
  if ( gir_pos == 1 ) then
    do ir = 0, (1 - lbound( rho%f2, 3))
      do iz = lbound( rho%f2, 2 ), ubound( rho%f2, 2 )
        rho%f2(1,iz,ir+2) = rho%f2(1,iz,ir+2) + rho%f2(1,iz,1-ir)
        rho%f2(1,iz,1-ir) = rho%f2(1,iz,ir+2) ! won't really affect anything
      enddo
    enddo
  endif


end subroutine norm_charge_cyl_m
!-------------------------------------------------------------------------------


end module m_species_charge_cyl_modes
