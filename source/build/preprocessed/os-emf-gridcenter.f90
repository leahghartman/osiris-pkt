# 1 "emf/os-emf-gridcenter.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "emf/os-emf-gridcenter.f90" 2
!-----------------------------------------------------------------------------------------
! Centers field values on the corner of the grid
!
!-----------------------------------------------------------------------------------------

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
# 7 "emf/os-emf-gridcenter.f90" 2
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
# 8 "emf/os-emf-gridcenter.f90" 2

module m_emf_gridcenter

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
# 12 "emf/os-emf-gridcenter.f90" 2

use m_parameters
use m_emf_define

implicit none

private

interface grid_center
  module procedure grid_center
end interface

public :: grid_center

contains

!-----------------------------------------------------------------------------------------
! Converts from a staggered Yee mesh to a grid with field values centered on the corner
! of the cell
!-----------------------------------------------------------------------------------------
subroutine grid_center( this )
!-----------------------------------------------------------------------------------------

  use m_vdf_define
  use m_vdf

  implicit none

  class( t_emf ), intent(inout) :: this
  integer :: i1, i2, i3

  ! temporary vdf buffer for centering calculations
  type( t_vdf ) :: buf

  call buf % new( this%e )

  select case ( p_x_dim )

    case (1)
      do i1 = lbound( buf%f1, 2 ) + 1, ubound( buf%f1, 2 )
        buf%f1( 1, i1 ) = 0.5 * ( this%e_part%f1( 1, i1 ) + &
                                  this%e_part%f1( 1, i1-1) )
        buf%f1( 2, i1 ) = this%e_part%f1( 2, i1 )
        buf%f1( 3, i1 ) = this%e_part%f1( 3, i1 )
      enddo

      this%e_part = buf

      do i1 = lbound( buf%f1, 2 ) + 1, ubound( buf%f1, 2 )
        buf%f1( 1, i1 ) = this%b_part%f1( 1, i1 )
        buf%f1( 2, i1 ) = 0.5 * ( this%b_part%f1( 2, i1 ) + &
                                   this%b_part%f1( 2, i1-1) )
        buf%f1( 3, i1 ) = 0.5 * ( this%b_part%f1( 3, i1 ) + &
                                   this%b_part%f1( 3, i1-1) )
      enddo

      this%b_part = buf

    case (2)

      do i2 = lbound( buf%f2, 3 ) + 1, ubound( buf%f2, 3 )
        do i1 = lbound( buf%f2, 2 ) + 1, ubound( buf%f2, 2 )
          buf%f2( 1, i1, i2 ) = 0.5 * ( this%e_part%f2( 1, i1 , i2 ) + &
                                        this%e_part%f2( 1, i1-1, i2 ) )
          buf%f2( 2, i1, i2 ) = 0.5 * ( this%e_part%f2( 2, i1 , i2 ) + &
                                        this%e_part%f2( 2, i1 , i2-1) )
          buf%f2( 3, i1, i2 ) = this%e_part%f2( 3, i1 , i2 )
        enddo
      enddo

      this%e_part = buf

      do i2 = lbound( buf%f2, 3 ) + 1, ubound( buf%f2, 3 )
        do i1 = lbound( buf%f2, 2 ) + 1, ubound( buf%f2, 2 )
          buf%f2( 1, i1, i2 ) = 0.5 * ( this%b_part%f2( 1, i1 , i2 ) + &
                                         this%b_part%f2( 1, i1 , i2-1) )

          buf%f2( 2, i1, i2) = 0.5 * ( this%b_part%f2( 2, i1 , i2 ) + &
                                         this%b_part%f2( 2, i1-1, i2 ) )

          buf%f2( 3, i1, i2 ) = 0.25 * ( this%b_part%f2( 3, i1 , i2 ) + &
                                         this%b_part%f2( 3, i1-1, i2 ) + &

                                         this%b_part%f2( 3, i1 , i2-1) + &
                                         this%b_part%f2( 3, i1-1, i2-1) )
        enddo
      enddo

      this%b_part = buf

    case (3)

      do i3 = lbound( buf%f3, 4 ) + 1, ubound( buf%f3, 4 )
        do i2 = lbound( buf%f3, 3 ) + 1, ubound( buf%f3, 3 )
          do i1 = lbound( buf%f3, 2 ) + 1, ubound( buf%f3, 2 )
            buf%f3( 1, i1, i2, i3 ) = 0.5 * ( this%e_part%f3( 1, i1 , i2 , i3 ) + &
                                              this%e_part%f3( 1, i1-1, i2 , i3 ) )
            buf%f3( 2, i1, i2, i3 ) = 0.5 * ( this%e_part%f3( 2, i1 , i2 , i3 ) + &
                                              this%e_part%f3( 2, i1 , i2-1, i3 ) )
            buf%f3( 3, i1, i2, i3 ) = 0.5 * ( this%e_part%f3( 3, i1 , i2 , i3 ) + &
                                              this%e_part%f3( 3, i1 , i2 , i3-1 ) )
          enddo
        enddo
      enddo

      this%e_part = buf

      do i3 = lbound( buf%f3, 4 ) + 1, ubound( buf%f3, 4 )
        do i2 = lbound( buf%f3, 3 ) + 1, ubound( buf%f3, 3 )
          do i1 = lbound( buf%f3, 2 ) + 1, ubound( buf%f3, 2 )
            buf%f3( 1, i1, i2, i3 ) = 0.25 * ( this%b_part%f3( 1, i1 , i2 , i3 ) + &
                                               this%b_part%f3( 1, i1 , i2-1, i3 ) + &

                                               this%b_part%f3( 1, i1 , i2 , i3-1 ) + &
                                               this%b_part%f3( 1, i1 , i2-1, i3-1 ) )

            buf%f3( 2, i1, i2, i3 ) = 0.25 * ( this%b_part%f3( 2, i1 , i2 , i3 ) + &
                                               this%b_part%f3( 2, i1-1, i2 , i3 ) + &

                                               this%b_part%f3( 2, i1 , i2 , i3-1 ) + &
                                               this%b_part%f3( 2, i1-1, i2 , i3-1 ) )

            buf%f3( 3, i1, i2, i3 ) = 0.25 * ( this%b_part%f3( 3, i1 , i2 , i3 ) + &
                                               this%b_part%f3( 3, i1-1, i2 , i3 ) + &

                                               this%b_part%f3( 3, i1 , i2-1, i3 ) + &
                                               this%b_part%f3( 3, i1-1, i2-1, i3 ) )
          enddo
        enddo
      enddo

      this%b_part = buf

  end select

  call buf % cleanup()

end subroutine grid_center
!-----------------------------------------------------------------------------------------

end module m_emf_gridcenter
