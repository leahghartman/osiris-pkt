# 1 "current/os-current-boundary.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "current/os-current-boundary.f03" 2
module m_current_boundary

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
# 4 "current/os-current-boundary.f03" 2
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
# 5 "current/os-current-boundary.f03" 2
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
# 6 "current/os-current-boundary.f03" 2

use m_parameters
use m_vdf_define

implicit none

private

interface spec_bc_1d
  module procedure spec_bc_1d
end interface

interface spec_bc_2d
  module procedure spec_bc_2d
end interface

interface spec_bc_3d
  module procedure spec_bc_3d
end interface

public :: spec_bc_1d, spec_bc_2d, spec_bc_3d

contains

!-----------------------------------------------------------------------------------------
! Specular boundary for currents in 1D
!-----------------------------------------------------------------------------------------
subroutine spec_bc_1d( jay, bnd, interpolation )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ), intent(inout) :: jay
  integer, intent(in) :: bnd, interpolation

  integer :: ipos

  select case ( interpolation )
  case ( p_linear )
    ! Reflection at the cell edge
    select case ( bnd )
      case ( p_lower )
        jay%f1( 1, 1 ) = jay%f1( 1, 1 ) - jay%f1( 1, 0 )
        jay%f1( 2, 1 ) = jay%f1( 2, 1 ) * 2
        jay%f1( 3, 1 ) = jay%f1( 3, 1 ) * 2

        jay%f1( 2, 2 ) = jay%f1( 2, 2 ) + &
                    jay%f1( 2, 0 )
        jay%f1( 3, 2 ) = jay%f1( 3, 2 ) + &
                    jay%f1( 3, 0 )

      case ( p_upper )
         ipos = jay % nx_( 1 )

         jay%f1( 1, ipos ) = jay%f1( 1, ipos ) - &
                        jay%f1( 1, ipos + 1 )

         jay%f1( 2, ipos ) = jay%f1( 2, ipos ) + &
                         jay%f1( 2, ipos+2 )

         jay%f1( 3, ipos ) = jay%f1( 3, ipos ) + &
                         jay%f1( 3, ipos+2 )

         ! This point is outside of the simulation volume
         ! jay%f1( 1, ipos+1 ) = ...

         ! These are on the edge of the simulation volume
         jay%f1( 2, ipos+1 ) = jay%f1( 2, ipos+1 ) * 2
         jay%f1( 3, ipos+1 ) = jay%f1( 3, ipos+1 ) * 2

    end select

  case ( p_quadratic )
    ! Reflection at the cell center
    select case ( bnd )
      case ( p_lower )
    jay%f1( 1, 0 ) = 0

    jay%f1( 1, 1 ) = jay%f1( 1, 1 ) - jay%f1( 1, -1 )
    jay%f1( 2, 1 ) = jay%f1( 2, 1 ) + jay%f1( 2, 0 )
    jay%f1( 3, 1 ) = jay%f1( 3, 1 ) + jay%f1( 3, 0 )

    jay%f1( 2, 2 ) = jay%f1( 2, 2 ) + jay%f1( 2, -1 )
    jay%f1( 3, 2 ) = jay%f1( 3, 2 ) + jay%f1( 3, -1 )

      case ( p_upper )
     ipos = jay%nx_( 1 )

     jay%f1( 1, ipos-1 ) = jay%f1( 1, ipos-1 ) - &
                                   jay%f1( 1, ipos+1 )
     jay%f1( 2, ipos-1 ) = jay%f1( 2, ipos-1 ) + &
                     jay%f1( 2, ipos+2 )
     jay%f1( 3, ipos-1 ) = jay%f1( 3, ipos-1 ) + &
                     jay%f1( 3, ipos+2 )

     jay%f1( 1, ipos ) = 0
     jay%f1( 2, ipos ) = jay%f1( 2, ipos ) + &
                   jay%f1( 2, ipos + 1 )
     jay%f1( 3, ipos ) = jay%f1( 3, ipos ) + &
                   jay%f1( 3, ipos + 1 )

    end select

  case ( p_cubic )

    select case ( bnd )
      case ( p_lower )
    jay%f1( 1, 1 ) = jay%f1( 1, 1 ) - jay%f1( 1, 0 )
    jay%f1( 2, 1 ) = jay%f1( 2, 1 ) * 2
    jay%f1( 3, 1 ) = jay%f1( 3, 1 ) * 2

    jay%f1( 1, 2 ) = jay%f1( 1, 2 ) - &
                jay%f1( 1, -1 )
    jay%f1( 2, 2 ) = jay%f1( 2, 2 ) + &
                jay%f1( 2, 0 )
    jay%f1( 3, 2 ) = jay%f1( 3, 2 ) + &
                jay%f1( 3, 0 )

    jay%f1( 2, 3 ) = jay%f1( 2, 3 ) + &
                jay%f1( 2, -1 )
    jay%f1( 3, 3 ) = jay%f1( 3, 3 ) + &
                jay%f1( 3, -1 )

      case ( p_upper )
     ipos = jay%nx_( 1 )
     jay%f1( 1, ipos-1 ) = jay%f1( 1, ipos - 1 ) - &
                    jay%f1( 1, ipos + 2 )
     jay%f1( 2, ipos-1 ) = jay%f1( 2, ipos - 1 ) + &
                    jay%f1( 2, ipos + 3 )
     jay%f1( 3, ipos-1 ) = jay%f1( 3, ipos - 1 ) + &
                    jay%f1( 3, ipos + 3 )

     jay%f1( 1, ipos ) = jay%f1( 1, ipos ) - &
                    jay%f1( 1, ipos + 1 )
     jay%f1( 2, ipos ) = jay%f1( 2, ipos ) + &
                    jay%f1( 2, ipos + 2 )
     jay%f1( 3, ipos ) = jay%f1( 3, ipos ) + &
                    jay%f1( 3, ipos + 2 )

     ! This point is outside of the simulation volume
     ! jay%f1( 1, ipos+1 ) = ...

     ! These are on the edge of the simulation volume
     jay%f1( 2, ipos+1 ) = jay%f1( 2, ipos+1 ) * 2
     jay%f1( 3, ipos+1 ) = jay%f1( 3, ipos+1 ) * 2

    end select

  case ( p_quartic )
    ! Reflection at the cell center
    select case ( bnd )
      case ( p_lower )
    jay%f1( 1, 0 ) = 0

    jay%f1( 1, 1 ) = jay%f1( 1, 1 ) - jay%f1( 1, -1 )
    jay%f1( 2, 1 ) = jay%f1( 2, 1 ) + jay%f1( 2, 0 )
    jay%f1( 3, 1 ) = jay%f1( 3, 1 ) + jay%f1( 3, 0 )

    jay%f1( 1, 2 ) = jay%f1( 1, 2 ) - jay%f1( 1, -2 )
    jay%f1( 2, 2 ) = jay%f1( 2, 2 ) + jay%f1( 2, -1 )
    jay%f1( 3, 2 ) = jay%f1( 3, 2 ) + jay%f1( 3, -1 )

    jay%f1( 2, 3 ) = jay%f1( 2, 3 ) + jay%f1( 2, -2 )
    jay%f1( 3, 3 ) = jay%f1( 3, 3 ) + jay%f1( 3, -2 )

      case ( p_upper )
     ipos = jay%nx_( 1 )

     jay%f1( 1, ipos-2 ) = jay%f1( 1, ipos-2 ) - &
                                   jay%f1( 1, ipos+2 )
     jay%f1( 2, ipos-2 ) = jay%f1( 2, ipos-2 ) + &
                     jay%f1( 2, ipos+3 )
     jay%f1( 3, ipos-2 ) = jay%f1( 3, ipos-2 ) + &
                     jay%f1( 3, ipos+3 )

     jay%f1( 1, ipos-1 ) = jay%f1( 1, ipos-1 ) - &
                                   jay%f1( 1, ipos+1 )
     jay%f1( 2, ipos-1 ) = jay%f1( 2, ipos-1 ) + &
                     jay%f1( 2, ipos+2 )
     jay%f1( 3, ipos-1 ) = jay%f1( 3, ipos-1 ) + &
                     jay%f1( 3, ipos+2 )

     jay%f1( 1, ipos ) = 0
     jay%f1( 2, ipos ) = jay%f1( 2, ipos ) + &
                   jay%f1( 2, ipos + 1 )
     jay%f1( 3, ipos ) = jay%f1( 3, ipos ) + &
                   jay%f1( 3, ipos + 1 )

    end select

  end select

end subroutine spec_bc_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Specular boundary for currents in 2D
!-----------------------------------------------------------------------------------------
subroutine spec_bc_2d( jay, dim, bnd, interpolation )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ), intent(inout) :: jay
  integer, intent(in) :: dim, bnd, interpolation

  integer :: i1, i2, ipos

  select case ( interpolation )

  case (p_linear)
    ! linear interpolation
  select case (dim)
    case(1)
    select case ( bnd )
      case ( p_lower ) ! x1 lower boundary
      do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )
         jay%f2( 1, 1, i2 ) = jay%f2( 1, 1, i2 ) - &
                       jay%f2( 1, 0, i2 )
         jay%f2( 2, 1, i2 ) = jay%f2( 2, 1, i2 ) * 2
         jay%f2( 3, 1, i2 ) = jay%f2( 3, 1, i2 ) * 2

         jay%f2( 2, 2, i2 ) = jay%f2( 2, 2, i2 ) + &
                       jay%f2( 2, 0, i2 )
         jay%f2( 3, 2, i2 ) = jay%f2( 3, 2, i2 ) + &
                       jay%f2( 3, 0, i2 )

      enddo

      case ( p_upper ) ! x1 upper boundary
      ipos = jay%nx_( 1 )
      do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )

         jay%f2( 1, ipos, i2 ) = jay%f2( 1, ipos , i2 ) - &
                        jay%f2( 1, ipos+1, i2 )
         jay%f2( 2, ipos, i2 ) = jay%f2( 2, ipos , i2 ) + &
                        jay%f2( 2, ipos+2, i2 )
         jay%f2( 3, ipos, i2 ) = jay%f2( 3, ipos , i2 ) + &
                        jay%f2( 3, ipos+2, i2 )

         ! This point is outside of the simulation volume
         !jay%f2( 1, ipos+1, i2 ) = (...)

         ! These points are on the edge of the simulation volume
         jay%f2( 2, ipos+1, i2 ) = jay%f2( 2, ipos+1, i2 ) * 2
         jay%f2( 3, ipos+1, i2 ) = jay%f2( 3, ipos+1, i2 ) * 2
      enddo

    end select
    case(2)
    select case ( bnd )
      case ( p_lower ) ! x2 lower boundary
      do i1 = lbound( jay%f2, 2 ), ubound( jay%f2, 2 )

         jay%f2( 1, i1, 1 ) = jay%f2( 1, i1, 1 ) * 2 ! l
         jay%f2( 2, i1, 1 ) = jay%f2( 2, i1, 1 ) - & ! p
                       jay%f2( 2, i1, 0 )
         jay%f2( 3, i1, 1 ) = jay%f2( 3, i1, 1 ) * 2 ! l

         jay%f2( 1, i1, 2 ) = jay%f2( 1, i1, 2 ) + & ! l
                       jay%f2( 1, i1, 0 )
         jay%f2( 3, i1, 2 ) = jay%f2( 3, i1, 2 ) + & ! l
                       jay%f2( 3, i1, 0 )

      enddo

      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 2 )

      do i1 = lbound( jay%f2, 2 ), ubound( jay%f2, 2 )
         jay%f2( 1, i1, ipos ) = jay%f2( 1, i1, ipos ) + & ! l
                        jay%f2( 1, i1, ipos + 2 )

         jay%f2( 2, i1, ipos ) = jay%f2( 2, i1, ipos ) - & ! p
                        jay%f2( 2, i1, ipos + 1 )

         jay%f2( 3, i1, ipos ) = jay%f2( 3, i1, ipos ) + & ! l
                        jay%f2( 3, i1, ipos + 2 )

         ! This point in on the edge of the simulation box
         jay%f2( 1, i1, ipos+1 ) = jay%f2( 1, i1, ipos+1 ) * 2 ! l
         ! This point is outside the simulation box
         ! jay%f2( 2, i1, ipos+1 ) = (...)
         ! This point in on the edge of the simulation box
         jay%f2( 3, i1, ipos+1 ) = jay%f2( 3, i1, ipos+1 ) * 2 ! l
      enddo

    end select
  end select

  case (p_quadratic)
    ! quadratic interpolation, reflection at cell center
  select case (dim)
    case(1)
    select case ( bnd )
      case ( p_lower ) ! x1 lower boundary
      do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )

         jay%f2( 1, 0, i2 ) = 0

         jay%f2( 1, 1, i2 ) = jay%f2( 1, 1, i2 ) - &
                       jay%f2( 1, -1, i2 )
         jay%f2( 2, 1, i2 ) = jay%f2( 2, 1, i2 ) + &
                                     jay%f2( 2, 0, i2 )
         jay%f2( 3, 1, i2 ) = jay%f2( 3, 1, i2 ) + &
                                     jay%f2( 3, 0, i2 )

         jay%f2( 2, 2, i2 ) = jay%f2( 2, 2, i2 ) + &
                       jay%f2( 2, -1, i2 )
         jay%f2( 3, 2, i2 ) = jay%f2( 3, 2, i2 ) + &
                       jay%f2( 3, -1, i2 )

      enddo

      case ( p_upper ) ! x1 upper boundary
      ipos = jay%nx_( 1 )
      do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )

         jay%f2( 1, ipos-1, i2 ) = jay%f2( 1, ipos-1, i2 ) - &
                          jay%f2( 1, ipos+1, i2 )
         jay%f2( 2, ipos-1, i2 ) = jay%f2( 2, ipos-1, i2 ) + &
                          jay%f2( 2, ipos+2, i2 )
         jay%f2( 3, ipos-1, i2 ) = jay%f2( 3, ipos-1, i2 ) + &
                          jay%f2( 3, ipos+2, i2 )

         jay%f2( 1, ipos, i2 ) = 0
         jay%f2( 2, ipos, i2 ) = jay%f2( 2, ipos , i2 ) + &
                        jay%f2( 2, ipos+1, i2 )
         jay%f2( 3, ipos, i2 ) = jay%f2( 3, ipos , i2 ) + &
                        jay%f2( 3, ipos+1, i2 )

      enddo

    end select
    case(2)
    select case ( bnd )
      case ( p_lower ) ! x2 lower boundary
      do i1 = lbound( jay%f2, 2 ), ubound( jay%f2, 2 )
               jay%f2( 2, i1, 0 ) = 0

         jay%f2( 1, i1, 1 ) = jay%f2( 1, i1, 1 ) + &
                                     jay%f2( 1, i1, 0 )
         jay%f2( 2, i1, 1 ) = jay%f2( 2, i1, 1 ) - &
                       jay%f2( 2, i1, -1 )
         jay%f2( 3, i1, 1 ) = jay%f2( 3, i1, 1 ) + &
                                     jay%f2( 3, i1, 0 )

         jay%f2( 1, i1, 2 ) = jay%f2( 1, i1, 2 ) + &
                       jay%f2( 1, i1, -1 )
         jay%f2( 3, i1, 2 ) = jay%f2( 3, i1, 2 ) + &
                       jay%f2( 3, i1, -1 )

      enddo

      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 2 )

      do i1 = lbound( jay%f2, 2 ), ubound( jay%f2, 2 )
         jay%f2( 1, i1, ipos-1) = jay%f2( 1, i1, ipos-1 ) + &
                         jay%f2( 1, i1, ipos+2 )

         jay%f2( 2, i1, ipos-1 ) = jay%f2( 2, i1, ipos-1 ) - &
                        jay%f2( 2, i1, ipos+1 )

         jay%f2( 3, i1, ipos-1 ) = jay%f2( 3, i1, ipos-1 ) + &
                        jay%f2( 3, i1, ipos+2 )

         jay%f2( 1, i1, ipos ) = jay%f2( 1, i1, ipos ) + &
                        jay%f2( 1, i1, ipos+1 )

         jay%f2( 2, i1, ipos ) = 0

         jay%f2( 3, i1, ipos ) = jay%f2( 3, i1, ipos ) + &
                        jay%f2( 3, i1, ipos+1 )

      enddo

    end select
  end select

  case (p_cubic)
    ! cubic interpolation
  select case (dim)
    case(1)
    select case ( bnd )
      case ( p_lower ) ! x1 lower boundary
      do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )
         jay%f2( 1, 1, i2 ) = jay%f2( 1, 1, i2 ) - &
                       jay%f2( 1, 0, i2 )
         jay%f2( 2, 1, i2 ) = jay%f2( 2, 1, i2 ) * 2
         jay%f2( 3, 1, i2 ) = jay%f2( 3, 1, i2 ) * 2


         jay%f2( 1, 2, i2 ) = jay%f2( 1, 2, i2 ) - &
                       jay%f2( 1, -1, i2 )
         jay%f2( 2, 2, i2 ) = jay%f2( 2, 2, i2 ) + &
                       jay%f2( 2, 0, i2 )
         jay%f2( 3, 2, i2 ) = jay%f2( 3, 2, i2 ) + &
                       jay%f2( 3, 0, i2 )

         jay%f2( 2, 3, i2 ) = jay%f2( 2, 3, i2 ) + &
                       jay%f2( 2, -1, i2 )
         jay%f2( 3, 3, i2 ) = jay%f2( 3, 3, i2 ) + &
                       jay%f2( 3, -1, i2 )
      enddo

      case ( p_upper ) ! x1 upper boundary
      ipos = jay%nx_( 1 )
      do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )
         jay%f2( 1, ipos-1, i2 ) = jay%f2( 1, ipos-1, i2 ) - &
                        jay%f2( 1, ipos+2, i2 )
         jay%f2( 2, ipos-1, i2 ) = jay%f2( 2, ipos-1, i2 ) + &
                        jay%f2( 2, ipos+3, i2 )
         jay%f2( 3, ipos-1, i2 ) = jay%f2( 3, ipos-1, i2 ) + &
                        jay%f2( 3, ipos+3, i2 )

         jay%f2( 1, ipos, i2 ) = jay%f2( 1, ipos , i2 ) - &
                        jay%f2( 1, ipos+1, i2 )
         jay%f2( 2, ipos, i2 ) = jay%f2( 2, ipos , i2 ) + &
                        jay%f2( 2, ipos+2, i2 )
         jay%f2( 3, ipos, i2 ) = jay%f2( 3, ipos , i2 ) + &
                        jay%f2( 3, ipos+2, i2 )

         ! This point is outside of the simulation volume
         !jay%f2( 1, ipos+1, i2 ) = (...)

         ! These points are on the edge of the simulation volume
         jay%f2( 2, ipos+1, i2 ) = jay%f2( 2, ipos+1, i2 ) * 2
         jay%f2( 3, ipos+1, i2 ) = jay%f2( 3, ipos+1, i2 ) * 2
      enddo

    end select

    case(2)
    select case ( bnd )
      case ( p_lower ) ! x2 lower boundary
      do i1 = lbound( jay%f2, 2 ), ubound( jay%f2, 2 )

         jay%f2( 1, i1, 1 ) = jay%f2( 1, i1, 1 ) * 2 ! l
         jay%f2( 2, i1, 1 ) = jay%f2( 2, i1, 1 ) - & ! p
                       jay%f2( 2, i1, 0 )
         jay%f2( 3, i1, 1 ) = jay%f2( 3, i1, 1 ) * 2 ! l

         jay%f2( 1, i1, 2 ) = jay%f2( 1, i1, 2 ) + & ! l
                       jay%f2( 1, i1, 0 )
         jay%f2( 2, i1, 2 ) = jay%f2( 2, i1, 2 ) - & ! p
                       jay%f2( 2, i1, -1 )
         jay%f2( 3, i1, 2 ) = jay%f2( 3, i1, 2 ) + & ! l
                       jay%f2( 3, i1, 0 )

         jay%f2( 1, i1, 3 ) = jay%f2( 1, i1, 3 ) + & ! l
                       jay%f2( 1, i1,-1 )
         jay%f2( 3, i1, 3 ) = jay%f2( 3, i1, 2 ) + & ! l
                       jay%f2( 3, i1,-1 )
      enddo

      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 2 )

      do i1 = lbound( jay%f2, 2 ), ubound( jay%f2, 2 )
         jay%f2( 1, i1, ipos-1 ) = jay%f2( 1, i1, ipos - 1 ) + & ! l
                        jay%f2( 1, i1, ipos + 3 )

         jay%f2( 2, i1, ipos-1 ) = jay%f2( 2, i1, ipos - 1 ) - & ! p
                        jay%f2( 2, i1, ipos + 2 )

         jay%f2( 3, i1, ipos-1 ) = jay%f2( 3, i1, ipos - 1 ) + & ! l
                        jay%f2( 3, i1, ipos + 3 )

         jay%f2( 1, i1, ipos ) = jay%f2( 1, i1, ipos ) + & ! l
                        jay%f2( 1, i1, ipos + 2 )

         jay%f2( 2, i1, ipos ) = jay%f2( 2, i1, ipos ) - & ! p
                        jay%f2( 2, i1, ipos + 1 )

         jay%f2( 3, i1, ipos ) = jay%f2( 3, i1, ipos ) + & ! l
                        jay%f2( 3, i1, ipos + 2 )

         ! This point in on the edge of the simulation box
         jay%f2( 1, i1, ipos+1 ) = jay%f2( 1, i1, ipos+1 ) * 2 ! l
         ! This point is outside the simulation box
         ! jay%f2( 2, i1, ipos+1 ) = (...)
         ! This point in on the edge of the simulation box
         jay%f2( 3, i1, ipos+1 ) = jay%f2( 3, i1, ipos+1 ) * 2 ! l
      enddo

    end select

  end select

  case (p_quartic)
    ! 4th order interpolation, reflection at the cell center
  select case (dim)
    case(1)
    select case ( bnd )
      case ( p_lower ) ! x1 lower boundary
      do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )

         jay%f2( 1, 0, i2 ) = 0

         jay%f2( 1, 1, i2 ) = jay%f2( 1, 1, i2 ) - &
                       jay%f2( 1, -1, i2 )
         jay%f2( 2, 1, i2 ) = jay%f2( 2, 1, i2 ) + &
                                     jay%f2( 2, 0, i2 )
         jay%f2( 3, 1, i2 ) = jay%f2( 3, 1, i2 ) + &
                                     jay%f2( 3, 0, i2 )

         jay%f2( 1, 2, i2 ) = jay%f2( 1, 2, i2 ) - &
                       jay%f2( 1, -2, i2 )
         jay%f2( 2, 2, i2 ) = jay%f2( 2, 2, i2 ) + &
                       jay%f2( 2, -1, i2 )
         jay%f2( 3, 2, i2 ) = jay%f2( 3, 2, i2 ) + &
                       jay%f2( 3, -1, i2 )

         jay%f2( 2, 3, i2 ) = jay%f2( 2, 3, i2 ) + &
                       jay%f2( 2, -2, i2 )
         jay%f2( 3, 3, i2 ) = jay%f2( 3, 3, i2 ) + &
                       jay%f2( 3, -2, i2 )

      enddo

      case ( p_upper ) ! x1 upper boundary
      ipos = jay%nx_( 1 )
      do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )

         jay%f2( 1, ipos-2, i2 ) = jay%f2( 1, ipos-2, i2 ) - &
                          jay%f2( 1, ipos+2, i2 )
         jay%f2( 2, ipos-2, i2 ) = jay%f2( 2, ipos-2, i2 ) + &
                          jay%f2( 2, ipos+3, i2 )
         jay%f2( 3, ipos-2, i2 ) = jay%f2( 3, ipos-2, i2 ) + &
                          jay%f2( 3, ipos+3, i2 )

         jay%f2( 1, ipos-1, i2 ) = jay%f2( 1, ipos-1, i2 ) - &
                          jay%f2( 1, ipos+1, i2 )
         jay%f2( 2, ipos-1, i2 ) = jay%f2( 2, ipos-1, i2 ) + &
                          jay%f2( 2, ipos+2, i2 )
         jay%f2( 3, ipos-1, i2 ) = jay%f2( 3, ipos-1, i2 ) + &
                          jay%f2( 3, ipos+2, i2 )

         jay%f2( 1, ipos, i2 ) = 0
         jay%f2( 2, ipos, i2 ) = jay%f2( 2, ipos , i2 ) + &
                        jay%f2( 2, ipos+1, i2 )
         jay%f2( 3, ipos, i2 ) = jay%f2( 3, ipos , i2 ) + &
                        jay%f2( 3, ipos+1, i2 )

      enddo

    end select
    case(2)
    select case ( bnd )
      case ( p_lower ) ! x2 lower boundary
      do i1 = lbound( jay%f2, 2 ), ubound( jay%f2, 2 )
               jay%f2( 2, i1, 0 ) = 0

         jay%f2( 1, i1, 1 ) = jay%f2( 1, i1, 1 ) + &
                                     jay%f2( 1, i1, 0 )
         jay%f2( 2, i1, 1 ) = jay%f2( 2, i1, 1 ) - &
                       jay%f2( 2, i1, -1 )
         jay%f2( 3, i1, 1 ) = jay%f2( 3, i1, 1 ) + &
                                     jay%f2( 3, i1, 0 )

         jay%f2( 1, i1, 2 ) = jay%f2( 1, i1, 2 ) + &
                       jay%f2( 1, i1, -1 )
         jay%f2( 2, i1, 2 ) = jay%f2( 2, i1, 2 ) - &
                       jay%f2( 2, i1, -2 )
         jay%f2( 3, i1, 2 ) = jay%f2( 3, i1, 2 ) + &
                       jay%f2( 3, i1, -1 )

         jay%f2( 1, i1, 3 ) = jay%f2( 1, i1, 3 ) + &
                       jay%f2( 1, i1, -2 )
         jay%f2( 3, i1, 3 ) = jay%f2( 3, i1, 3 ) + &
                       jay%f2( 3, i1, -2 )

      enddo

      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 2 )

      do i1 = lbound( jay%f2, 2 ), ubound( jay%f2, 2 )
         jay%f2( 1, i1, ipos-2) = jay%f2( 1, i1, ipos-2 ) + &
                         jay%f2( 1, i1, ipos+3 )

         jay%f2( 2, i1, ipos-2 ) = jay%f2( 2, i1, ipos-2 ) - &
                        jay%f2( 2, i1, ipos+2 )

         jay%f2( 3, i1, ipos-2 ) = jay%f2( 3, i1, ipos-2 ) + &
                        jay%f2( 3, i1, ipos+3 )

         jay%f2( 1, i1, ipos-1) = jay%f2( 1, i1, ipos-1 ) + &
                         jay%f2( 1, i1, ipos+2 )

         jay%f2( 2, i1, ipos-1 ) = jay%f2( 2, i1, ipos-1 ) - &
                        jay%f2( 2, i1, ipos+1 )

         jay%f2( 3, i1, ipos-1 ) = jay%f2( 3, i1, ipos-1 ) + &
                        jay%f2( 3, i1, ipos+2 )

         jay%f2( 1, i1, ipos ) = jay%f2( 1, i1, ipos ) + &
                        jay%f2( 1, i1, ipos+1 )

         jay%f2( 2, i1, ipos ) = 0

         jay%f2( 3, i1, ipos ) = jay%f2( 3, i1, ipos ) + &
                        jay%f2( 3, i1, ipos+1 )

      enddo

    end select
  end select

  end select

end subroutine spec_bc_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Specular boundary for currents in 2D
!-----------------------------------------------------------------------------------------
subroutine spec_bc_3d( jay, dim, bnd, interpolation )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ), intent(inout) :: jay
  integer, intent(in) :: dim, bnd, interpolation

  integer :: i1, i2, i3, ipos

  select case ( interpolation )

  case (p_linear)

  ! linear interpolation
  select case (dim)
    case(1)
    select case ( bnd )
      case ( p_lower ) ! x1 lower boundary
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
         jay%f3( 1, 1, i2, i3 ) = jay%f3( 1, 1, i2, i3 ) - & ! p
                         jay%f3( 1, 0, i2, i3 )
         jay%f3( 2, 1, i2, i3 ) = jay%f3( 2, 1, i2, i3 ) * 2 ! l
         jay%f3( 3, 1, i2, i3 ) = jay%f3( 3, 1, i2, i3 ) * 2 ! l


         jay%f3( 2, 2, i2, i3 ) = jay%f3( 2, 2, i2, i3 ) + & ! l
                         jay%f3( 2, 0, i2, i3 )
         jay%f3( 3, 2, i2, i3 ) = jay%f3( 3, 2, i2, i3 ) + & ! l
                         jay%f3( 3, 0, i2, i3 )
        enddo
      enddo

      case ( p_upper ) ! x1 upper boundary
      ipos = jay%nx_( 1 )

      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
         jay%f3( 1, ipos, i2, i3 ) = jay%f3( 1, ipos , i2, i3 ) - & ! p
                            jay%f3( 1, ipos+1, i2, i3 )
         jay%f3( 2, ipos, i2, i3 ) = jay%f3( 2, ipos , i2, i3 ) + & ! l
                            jay%f3( 2, ipos+2, i2, i3 )
         jay%f3( 3, ipos, i2, i3 ) = jay%f3( 3, ipos , i2, i3 ) + & ! l
                            jay%f3( 3, ipos+2, i2, i3 )

         ! This point is outside of the simulation volume
         !jay%f3( 1, ipos+1, i2, i3 ) = (...) ! p

         ! These points are on the edge of the simulation volume
         jay%f3( 2, ipos+1, i2, i3 ) = jay%f3( 2, ipos+1, i2, i3 ) * 2 ! l
         jay%f3( 3, ipos+1, i2, i3 ) = jay%f3( 3, ipos+1, i2, i3 ) * 2 ! l
        enddo
      enddo
    end select
    case(2)
    select case ( bnd )
      case ( p_lower ) ! x2 lower boundary
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )

         jay%f3( 1, i1, 1, i3 ) = jay%f3( 1, i1, 1, i3 ) * 2 ! l
         jay%f3( 2, i1, 1, i3 ) = jay%f3( 2, i1, 1, i3 ) - & ! p
                         jay%f3( 2, i1, 0, i3 )
         jay%f3( 3, i1, 1, i3 ) = jay%f3( 3, i1, 1, i3 ) * 2 ! l

         jay%f3( 1, i1, 2, i3 ) = jay%f3( 1, i1, 2, i3 ) + & ! l
                         jay%f3( 1, i1, 0, i3 )
         jay%f3( 3, i1, 2, i3 ) = jay%f3( 3, i1, 2, i3 ) + & ! l
                         jay%f3( 3, i1, 0, i3 )
        enddo
      enddo
      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 2 )

      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 1, i1, ipos , i3 ) = jay%f3( 1, i1, ipos , i3 ) + & ! l
                            jay%f3( 1, i1, ipos + 2, i3 )
         jay%f3( 2, i1, ipos , i3 ) = jay%f3( 2, i1, ipos , i3 ) - & ! p
                            jay%f3( 2, i1, ipos + 1, i3 )
         jay%f3( 3, i1, ipos , i3 ) = jay%f3( 3, i1, ipos , i3 ) + & ! l
                            jay%f3( 3, i1, ipos + 2, i3 )

         ! This point in on the edge of the simulation box
         jay%f3( 1, i1, ipos+1, i3 ) = jay%f3( 1, i1, ipos+1, i3 ) * 2 ! l
         ! This point is outside the simulation box
         ! jay%f3( 2, i1, ipos+1, i3 ) = (...)
         ! This point in on the edge of the simulation box
         jay%f3( 3, i1, ipos+1, i3 ) = jay%f3( 3, i1, ipos+1, i3 ) * 2 ! l
        enddo
      enddo

    end select
    case(3)
    select case ( bnd )
      case ( p_lower ) ! x3 lower boundary
      do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )

         jay%f3( 1, i1, i2, 1 ) = jay%f3( 1, i1, i2, 1 ) * 2 ! l
         jay%f3( 2, i1, i2, 1 ) = jay%f3( 2, i1, i2, 1 ) * 2 ! l
         jay%f3( 3, i1, i2, 1 ) = jay%f3( 3, i1, i2, 1 ) - & ! p
                         jay%f3( 3, i1, i2, 0 )

         jay%f3( 1, i1, i2, 2 ) = jay%f3( 1, i1, i2, 2 ) + & ! l
                         jay%f3( 1, i1, i2, 0 )
         jay%f3( 2, i1, i2, 2 ) = jay%f3( 2, i1, i2, 2 ) + & ! l
                         jay%f3( 2, i1, i2, 0 )
        enddo
      enddo
      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 3 )

      do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 1, i1, i2, ipos ) = jay%f3( 1, i1, i2, ipos ) + & ! l
                            jay%f3( 1, i1, i2, ipos + 2 )
         jay%f3( 2, i1, i2, ipos ) = jay%f3( 2, i1, i2, ipos ) + & ! l
                            jay%f3( 2, i1, i2, ipos + 2 )
         jay%f3( 3, i1, i2, ipos ) = jay%f3( 3, i1, i2, ipos ) - & ! p
                            jay%f3( 3, i1, i2, ipos + 1 )

         ! These points are on the edge of the simulation box
         jay%f3( 1, i1, i2, ipos+1 ) = jay%f3( 1, i1, i2, ipos+1 ) * 2 ! l
         jay%f3( 2, i1, i2, ipos+1 ) = jay%f3( 2, i1, i2, ipos+1 ) * 2 ! l
         ! This point is outside the simulation box
         ! jay%f3( 3, i1, i2, ipos+1 ) = (...)
        enddo
      enddo
    end select
  end select

  case (p_quadratic)
    ! quadratic interpolation, reflection at cell center
  select case (dim)
    case(1)
    select case ( bnd )
      case ( p_lower ) ! x1 lower boundary
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )

         jay%f3( 1, 0, i2, i3 ) = 0

         jay%f3( 1, 1, i2, i3 ) = jay%f3( 1, 1, i2, i3 ) - &
                           jay%f3( 1, -1, i2, i3 )
         jay%f3( 2, 1, i2, i3 ) = jay%f3( 2, 1, i2, i3 ) + &
                           jay%f3( 2, 0, i2, i3 )
         jay%f3( 3, 1, i2, i3 ) = jay%f3( 3, 1, i2, i3 ) + &
                           jay%f3( 3, 0, i2, i3 )

         jay%f3( 2, 2, i2, i3 ) = jay%f3( 2, 2, i2, i3 ) + &
                           jay%f3( 2, -1, i2, i3 )
         jay%f3( 3, 2, i2, i3 ) = jay%f3( 3, 2, i2, i3 ) + &
                           jay%f3( 3, -1, i2, i3 )
        enddo
      enddo

      case ( p_upper ) ! x1 upper boundary
      ipos = jay%nx_( 1 )
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )

         jay%f3( 1, ipos-1, i2, i3 ) = jay%f3( 1, ipos-1, i2, i3 ) - &
                              jay%f3( 1, ipos+1, i2, i3 )
         jay%f3( 2, ipos-1, i2, i3 ) = jay%f3( 2, ipos-1, i2, i3 ) + &
                              jay%f3( 2, ipos+2, i2, i3 )
         jay%f3( 3, ipos-1, i2, i3 ) = jay%f3( 3, ipos-1, i2, i3 ) + &
                              jay%f3( 3, ipos+2, i2, i3 )

         jay%f3( 1, ipos, i2, i3 ) = 0
         jay%f3( 2, ipos, i2, i3 ) = jay%f3( 2, ipos , i2, i3 ) + &
                            jay%f3( 2, ipos+1, i2, i3 )
         jay%f3( 3, ipos, i2, i3 ) = jay%f3( 3, ipos , i2, i3 ) + &
                            jay%f3( 3, ipos+1, i2, i3 )

        enddo
      enddo
    end select

    case(2)
    select case ( bnd )
      case ( p_lower ) ! x2 lower boundary
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 2, i1, 0, i3 ) = 0

         jay%f3( 1, i1, 1, i3 ) = jay%f3( 1, i1, 1, i3 ) + &
                         jay%f3( 1, i1, 0, i3 )
         jay%f3( 2, i1, 1, i3 ) = jay%f3( 2, i1, 1, i3 ) - &
                         jay%f3( 2, i1, -1, i3 )
         jay%f3( 3, i1, 1, i3 ) = jay%f3( 3, i1, 1, i3 ) + &
                         jay%f3( 3, i1, 0, i3 )

         jay%f3( 1, i1, 2, i3 ) = jay%f3( 1, i1, 2, i3 ) + &
                         jay%f3( 1, i1, -1, i3 )
         jay%f3( 3, i1, 2, i3 ) = jay%f3( 3, i1, 2, i3 ) + &
                         jay%f3( 3, i1, -1, i3 )
        enddo
      enddo

      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 2 )

      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 1, i1, ipos-1, i3 ) = jay%f3( 1, i1, ipos-1, i3 ) + &
                              jay%f3( 1, i1, ipos+2, i3 )

         jay%f3( 2, i1, ipos-1, i3 ) = jay%f3( 2, i1, ipos-1, i3 ) - &
                              jay%f3( 2, i1, ipos+1, i3 )

         jay%f3( 3, i1, ipos-1, i3 ) = jay%f3( 3, i1, ipos-1, i3 ) + &
                              jay%f3( 3, i1, ipos+2, i3 )

         jay%f3( 1, i1, ipos, i3 ) = jay%f3( 1, i1, ipos , i3 ) + &
                            jay%f3( 1, i1, ipos+1, i3 )

         jay%f3( 2, i1, ipos, i3 ) = 0

         jay%f3( 3, i1, ipos, i3 ) = jay%f3( 3, i1, ipos , i3 ) + &
                            jay%f3( 3, i1, ipos+1, i3 )
        enddo
      enddo

    end select

    case(3)
    select case ( bnd )
      case ( p_lower ) ! x3 lower boundary
      do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 3, i1, i2, 0 ) = 0

         jay%f3( 1, i1, i2, 1 ) = jay%f3( 1, i1, i2, 1 ) + &
                         jay%f3( 1, i1, i2, 0 )
         jay%f3( 2, i1, i2, 1 ) = jay%f3( 2, i1, i2, 1 ) + &
                         jay%f3( 2, i1, i2, 0 )
         jay%f3( 3, i1, i2, 1 ) = jay%f3( 3, i1, i2, 1 ) - &
                         jay%f3( 3, i1, i2, -1 )

         jay%f3( 1, i1, i2, 2 ) = jay%f3( 1, i1, i2, 2 ) + &
                         jay%f3( 1, i1, i2, -1 )
         jay%f3( 2, i1, i2, 2 ) = jay%f3( 2, i1, i2, 2 ) + &
                         jay%f3( 2, i1, i2, -1 )
        enddo
      enddo

      case ( p_upper ) ! x3 upper boundary
      ipos = jay%nx_( 3 )

      do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 1, i1, i2, ipos-1 ) = jay%f3( 1, i1, i2, ipos-1 ) + &
                              jay%f3( 1, i1, i2, ipos+2 )

         jay%f3( 2, i1, i2, ipos-1 ) = jay%f3( 2, i1, i2, ipos-1 ) + &
                              jay%f3( 2, i1, i2, ipos+2 )

         jay%f3( 3, i1, i2, ipos-1 ) = jay%f3( 3, i1, i2, ipos-1 ) - &
                              jay%f3( 3, i1, i2, ipos+1 )

         jay%f3( 1, i1, i2, ipos ) = jay%f3( 1, i1, i2, ipos ) + &
                              jay%f3( 1, i1, i2, ipos+1 )

         jay%f3( 2, i1, i2, ipos ) = jay%f3( 2, i1, i2, ipos ) + &
                              jay%f3( 2, i1, i2, ipos+1 )

         jay%f3( 3, i1, i2, ipos ) = 0
        enddo
      enddo

    end select
  end select


  case (p_cubic)
  ! cubic interpolation
  select case (dim)
    case(1)
    select case ( bnd )
      case ( p_lower ) ! x1 lower boundary
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
         jay%f3( 1, 1, i2, i3 ) = jay%f3( 1, 1, i2, i3 ) - & ! p
                         jay%f3( 1, 0, i2, i3 )
         jay%f3( 2, 1, i2, i3 ) = jay%f3( 2, 1, i2, i3 ) * 2 ! l
         jay%f3( 3, 1, i2, i3 ) = jay%f3( 3, 1, i2, i3 ) * 2 ! l


         jay%f3( 1, 2, i2, i3 ) = jay%f3( 1, 2, i2, i3 ) - & ! p
                         jay%f3( 1, -1, i2, i3 )
         jay%f3( 2, 2, i2, i3 ) = jay%f3( 2, 2, i2, i3 ) + & ! l
                         jay%f3( 2, 0, i2, i3 )
         jay%f3( 3, 2, i2, i3 ) = jay%f3( 3, 2, i2, i3 ) + & ! l
                         jay%f3( 3, 0, i2, i3 )

         jay%f3( 2, 3, i2, i3 ) = jay%f3( 2, 3, i2, i3 ) + & ! l
                         jay%f3( 2, -1, i2, i3 )
         jay%f3( 3, 3, i2, i3 ) = jay%f3( 3, 3, i2, i3 ) + & ! l
                         jay%f3( 3, -1, i2, i3 )

        enddo
      enddo

      case ( p_upper ) ! x1 upper boundary
      ipos = jay%nx_( 1 )

      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
         jay%f3( 1, ipos-1, i2, i3 ) = jay%f3( 1, ipos-1, i2, i3 ) - & ! p
                            jay%f3( 1, ipos+2, i2, i3 )
         jay%f3( 2, ipos-1, i2, i3 ) = jay%f3( 2, ipos-1, i2, i3 ) + & ! l
                            jay%f3( 2, ipos+3, i2, i3 )
         jay%f3( 3, ipos-1, i2, i3 ) = jay%f3( 3, ipos-1, i2, i3 ) + & ! l
                            jay%f3( 3, ipos+3, i2, i3 )

         jay%f3( 1, ipos, i2, i3 ) = jay%f3( 1, ipos , i2, i3 ) - & ! p
                            jay%f3( 1, ipos+1, i2, i3 )
         jay%f3( 2, ipos, i2, i3 ) = jay%f3( 2, ipos , i2, i3 ) + & ! l
                            jay%f3( 2, ipos+2, i2, i3 )
         jay%f3( 3, ipos, i2, i3 ) = jay%f3( 3, ipos , i2, i3 ) + & ! l
                            jay%f3( 3, ipos+2, i2, i3 )

         ! This point is outside of the simulation volume
         !jay%f3( 1, ipos+1, i2, i3 ) = (...) ! p

         ! These points are on the edge of the simulation volume
         jay%f3( 2, ipos+1, i2, i3 ) = jay%f3( 2, ipos+1, i2, i3 ) * 2 ! l
         jay%f3( 3, ipos+1, i2, i3 ) = jay%f3( 3, ipos+1, i2, i3 ) * 2 ! l
        enddo
      enddo
    end select
    case(2)
    select case ( bnd )
      case ( p_lower ) ! x2 lower boundary
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )

         jay%f3( 1, i1, 1, i3 ) = jay%f3( 1, i1, 1, i3 ) * 2 ! l
         jay%f3( 2, i1, 1, i3 ) = jay%f3( 2, i1, 1, i3 ) - & ! p
                         jay%f3( 2, i1, 0, i3 )
         jay%f3( 3, i1, 1, i3 ) = jay%f3( 3, i1, 1, i3 ) * 2 ! l

         jay%f3( 1, i1, 2, i3 ) = jay%f3( 1, i1, 2, i3 ) + & ! l
                         jay%f3( 1, i1, 0, i3 )
         jay%f3( 2, i1, 2, i3 ) = jay%f3( 2, i1, 2, i3 ) - & ! p
                         jay%f3( 2, i1, -1, i3 )
         jay%f3( 3, i1, 2, i3 ) = jay%f3( 3, i1, 2, i3 ) + & ! l
                         jay%f3( 3, i1, 0, i3 )

         jay%f3( 1, i1, 3, i3 ) = jay%f3( 1, i1, 3, i3 ) + & ! l
                         jay%f3( 1, i1, -1, i3 )
         jay%f3( 3, i1, 3, i3 ) = jay%f3( 3, i1, 3, i3 ) + & ! l
                         jay%f3( 3, i1, -1, i3 )
        enddo
      enddo
      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 2 )

      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 1, i1, ipos-1, i3 ) = jay%f3( 1, i1, ipos - 1, i3 ) + & ! l
                            jay%f3( 1, i1, ipos + 3, i3 )
         jay%f3( 2, i1, ipos-1, i3 ) = jay%f3( 2, i1, ipos - 1, i3 ) - & ! p
                            jay%f3( 2, i1, ipos + 2, i3 )
         jay%f3( 3, i1, ipos-1, i3 ) = jay%f3( 3, i1, ipos - 1, i3 ) + & ! l
                            jay%f3( 3, i1, ipos + 3, i3 )

         jay%f3( 1, i1, ipos , i3 ) = jay%f3( 1, i1, ipos , i3 ) + & ! l
                            jay%f3( 1, i1, ipos + 2, i3 )
         jay%f3( 2, i1, ipos , i3 ) = jay%f3( 2, i1, ipos , i3 ) - & ! p
                            jay%f3( 2, i1, ipos + 1, i3 )
         jay%f3( 3, i1, ipos , i3 ) = jay%f3( 3, i1, ipos , i3 ) + & ! l
                            jay%f3( 3, i1, ipos + 2, i3 )

         ! This point in on the edge of the simulation box
         jay%f3( 1, i1, ipos+1, i3 ) = jay%f3( 1, i1, ipos+1, i3 ) * 2 ! l
         ! This point is outside the simulation box
         ! jay%f3( 2, i1, ipos+1, i3 ) = (...)
         ! This point in on the edge of the simulation box
         jay%f3( 3, i1, ipos+1, i3 ) = jay%f3( 3, i1, ipos+1, i3 ) * 2 ! l
        enddo
      enddo

    end select
    case(3)
    select case ( bnd )
      case ( p_lower ) ! x3 lower boundary
      do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )

         jay%f3( 1, i1, i2, 1 ) = jay%f3( 1, i1, i2, 1 ) * 2 ! l
         jay%f3( 2, i1, i2, 1 ) = jay%f3( 2, i1, i2, 1 ) * 2 ! l
         jay%f3( 3, i1, i2, 1 ) = jay%f3( 3, i1, i2, 1 ) - & ! p
                         jay%f3( 3, i1, i2, 0 )

         jay%f3( 1, i1, i2, 2 ) = jay%f3( 1, i1, i2, 2 ) + & ! l
                         jay%f3( 1, i1, i2, 0 )
         jay%f3( 2, i1, i2, 2 ) = jay%f3( 2, i1, i2, 2 ) + & ! l
                         jay%f3( 2, i1, i2, 0 )
         jay%f3( 3, i1, i2, 2 ) = jay%f3( 3, i1, i2, 2 ) - & ! p
                         jay%f3( 3, i1, i2,-1 )

         jay%f3( 1, i1, i2, 3 ) = jay%f3( 1, i1, i2, 3 ) + & ! l
                         jay%f3( 1, i1, i2,-1 )
         jay%f3( 2, i1, i2, 3 ) = jay%f3( 2, i1, i2, 3 ) + & ! l
                         jay%f3( 2, i1, i2,-1 )
        enddo
      enddo
      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 3 )

      do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 1, i1, i2, ipos-1 ) = jay%f3( 1, i1, i2, ipos - 1 ) + & ! l
                            jay%f3( 1, i1, i2, ipos + 3 )
         jay%f3( 2, i1, i2, ipos-1 ) = jay%f3( 2, i1, i2, ipos - 1 ) + & ! l
                            jay%f3( 2, i1, i2, ipos + 3 )
         jay%f3( 3, i1, i2, ipos-1 ) = jay%f3( 3, i1, i2, ipos - 1 ) - & ! p
                            jay%f3( 3, i1, i2, ipos + 2 )

         jay%f3( 1, i1, i2, ipos ) = jay%f3( 1, i1, i2, ipos ) + & ! l
                            jay%f3( 1, i1, i2, ipos + 2 )
         jay%f3( 2, i1, i2, ipos ) = jay%f3( 2, i1, i2, ipos ) + & ! l
                            jay%f3( 2, i1, i2, ipos + 2 )
         jay%f3( 3, i1, i2, ipos ) = jay%f3( 3, i1, i2, ipos ) - & ! p
                            jay%f3( 3, i1, i2, ipos + 1 )

         ! These points are on the edge of the simulation box
         jay%f3( 1, i1, i2, ipos+1 ) = jay%f3( 1, i1, i2, ipos+1 ) * 2 ! l
         jay%f3( 2, i1, i2, ipos+1 ) = jay%f3( 2, i1, i2, ipos+1 ) * 2 ! l
         ! This point is outside the simulation box
         ! jay%f3( 3, i1, i2, ipos+1 ) = (...)
        enddo
      enddo
    end select
  end select

  case (p_quartic)
    ! quartic interpolation, reflection at cell center
  select case (dim)
    case(1)
    select case ( bnd )
      case ( p_lower ) ! x1 lower boundary
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )

         jay%f3( 1, 0, i2, i3 ) = 0

         jay%f3( 1, 1, i2, i3 ) = jay%f3( 1, 1, i2, i3 ) - &
                           jay%f3( 1, -1, i2, i3 )
         jay%f3( 2, 1, i2, i3 ) = jay%f3( 2, 1, i2, i3 ) + &
                           jay%f3( 2, 0, i2, i3 )
         jay%f3( 3, 1, i2, i3 ) = jay%f3( 3, 1, i2, i3 ) + &
                           jay%f3( 3, 0, i2, i3 )

         jay%f3( 1, 2, i2, i3 ) = jay%f3( 1, 2, i2, i3 ) - &
                           jay%f3( 1, -2, i2, i3 )
         jay%f3( 2, 2, i2, i3 ) = jay%f3( 2, 2, i2, i3 ) + &
                           jay%f3( 2, -1, i2, i3 )
         jay%f3( 3, 2, i2, i3 ) = jay%f3( 3, 2, i2, i3 ) + &
                           jay%f3( 3, -1, i2, i3 )

         jay%f3( 2, 3, i2, i3 ) = jay%f3( 2, 3, i2, i3 ) + &
                           jay%f3( 2, -2, i2, i3 )
         jay%f3( 3, 3, i2, i3 ) = jay%f3( 3, 3, i2, i3 ) + &
                           jay%f3( 3, -2, i2, i3 )

        enddo
      enddo

      case ( p_upper ) ! x1 upper boundary
      ipos = jay%nx_( 1 )
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )

         jay%f3( 1, ipos-2, i2, i3 ) = jay%f3( 1, ipos-2, i2, i3 ) - &
                              jay%f3( 1, ipos+2, i2, i3 )
         jay%f3( 2, ipos-2, i2, i3 ) = jay%f3( 2, ipos-2, i2, i3 ) + &
                              jay%f3( 2, ipos+3, i2, i3 )
         jay%f3( 3, ipos-2, i2, i3 ) = jay%f3( 3, ipos-2, i2, i3 ) + &
                              jay%f3( 3, ipos+3, i2, i3 )

         jay%f3( 1, ipos-1, i2, i3 ) = jay%f3( 1, ipos-1, i2, i3 ) - &
                              jay%f3( 1, ipos+1, i2, i3 )
         jay%f3( 2, ipos-1, i2, i3 ) = jay%f3( 2, ipos-1, i2, i3 ) + &
                              jay%f3( 2, ipos+2, i2, i3 )
         jay%f3( 3, ipos-1, i2, i3 ) = jay%f3( 3, ipos-1, i2, i3 ) + &
                              jay%f3( 3, ipos+2, i2, i3 )

         jay%f3( 1, ipos, i2, i3 ) = 0
         jay%f3( 2, ipos, i2, i3 ) = jay%f3( 2, ipos , i2, i3 ) + &
                            jay%f3( 2, ipos+1, i2, i3 )
         jay%f3( 3, ipos, i2, i3 ) = jay%f3( 3, ipos , i2, i3 ) + &
                            jay%f3( 3, ipos+1, i2, i3 )

        enddo
      enddo
    end select

    case(2)
    select case ( bnd )
      case ( p_lower ) ! x2 lower boundary
      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 2, i1, 0, i3 ) = 0

         jay%f3( 1, i1, 1, i3 ) = jay%f3( 1, i1, 1, i3 ) + &
                         jay%f3( 1, i1, 0, i3 )
         jay%f3( 2, i1, 1, i3 ) = jay%f3( 2, i1, 1, i3 ) - &
                         jay%f3( 2, i1, -1, i3 )
         jay%f3( 3, i1, 1, i3 ) = jay%f3( 3, i1, 1, i3 ) + &
                         jay%f3( 3, i1, 0, i3 )

         jay%f3( 1, i1, 2, i3 ) = jay%f3( 1, i1, 2, i3 ) + &
                         jay%f3( 1, i1, -1, i3 )
         jay%f3( 2, i1, 2, i3 ) = jay%f3( 2, i1, 2, i3 ) - &
                         jay%f3( 2, i1, -2, i3 )
         jay%f3( 3, i1, 2, i3 ) = jay%f3( 3, i1, 2, i3 ) + &
                         jay%f3( 3, i1, -1, i3 )

         jay%f3( 1, i1, 3, i3 ) = jay%f3( 1, i1, 3, i3 ) + &
                         jay%f3( 1, i1, -2, i3 )
         jay%f3( 3, i1, 3, i3 ) = jay%f3( 3, i1, 3, i3 ) + &
                         jay%f3( 3, i1, -2, i3 )

        enddo
      enddo

      case ( p_upper ) ! x2 upper boundary
      ipos = jay%nx_( 2 )

      do i3 = lbound( jay%f3, 4 ), ubound( jay%f3, 4 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 1, i1, ipos-2, i3 ) = jay%f3( 1, i1, ipos-2, i3 ) + &
                              jay%f3( 1, i1, ipos+3, i3 )

         jay%f3( 2, i1, ipos-2, i3 ) = jay%f3( 2, i1, ipos-2, i3 ) - &
                              jay%f3( 2, i1, ipos+2, i3 )

         jay%f3( 3, i1, ipos-2, i3 ) = jay%f3( 3, i1, ipos-2, i3 ) + &
                              jay%f3( 3, i1, ipos+3, i3 )

         jay%f3( 1, i1, ipos-1, i3 ) = jay%f3( 1, i1, ipos-1, i3 ) + &
                              jay%f3( 1, i1, ipos+2, i3 )

         jay%f3( 2, i1, ipos-1, i3 ) = jay%f3( 2, i1, ipos-1, i3 ) - &
                              jay%f3( 2, i1, ipos+1, i3 )

         jay%f3( 3, i1, ipos-1, i3 ) = jay%f3( 3, i1, ipos-1, i3 ) + &
                              jay%f3( 3, i1, ipos+2, i3 )

         jay%f3( 1, i1, ipos, i3 ) = jay%f3( 1, i1, ipos , i3 ) + &
                            jay%f3( 1, i1, ipos+1, i3 )

         jay%f3( 2, i1, ipos, i3 ) = 0

         jay%f3( 3, i1, ipos, i3 ) = jay%f3( 3, i1, ipos , i3 ) + &
                            jay%f3( 3, i1, ipos+1, i3 )
        enddo
      enddo

    end select

    case(3)
    select case ( bnd )
      case ( p_lower ) ! x3 lower boundary
      do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 3, i1, i2, 0 ) = 0

         jay%f3( 1, i1, i2, 1 ) = jay%f3( 1, i1, i2, 1 ) + &
                         jay%f3( 1, i1, i2, 0 )
         jay%f3( 2, i1, i2, 1 ) = jay%f3( 2, i1, i2, 1 ) + &
                         jay%f3( 2, i1, i2, 0 )
         jay%f3( 3, i1, i2, 1 ) = jay%f3( 3, i1, i2, 1 ) - &
                         jay%f3( 3, i1, i2, -1 )

         jay%f3( 1, i1, i2, 2 ) = jay%f3( 1, i1, i2, 2 ) + &
                         jay%f3( 1, i1, i2, -1 )
         jay%f3( 2, i1, i2, 2 ) = jay%f3( 2, i1, i2, 2 ) + &
                         jay%f3( 2, i1, i2, -1 )
         jay%f3( 3, i1, i2, 2 ) = jay%f3( 3, i1, i2, 2 ) - &
                         jay%f3( 3, i1, i2, -2 )

         jay%f3( 1, i1, i2, 3 ) = jay%f3( 1, i1, i2, 3 ) + &
                         jay%f3( 1, i1, i2, -2 )
         jay%f3( 2, i1, i2, 3 ) = jay%f3( 2, i1, i2, 3 ) + &
                         jay%f3( 2, i1, i2, -2 )

        enddo
      enddo

      case ( p_upper ) ! x3 upper boundary
      ipos = jay%nx_( 3 )

      do i2 = lbound( jay%f3, 3 ), ubound( jay%f3, 3 )
        do i1 = lbound( jay%f3, 2 ), ubound( jay%f3, 2 )
         jay%f3( 1, i1, i2, ipos-2 ) = jay%f3( 1, i1, i2, ipos-2 ) + &
                              jay%f3( 1, i1, i2, ipos+3 )

         jay%f3( 2, i1, i2, ipos-2 ) = jay%f3( 2, i1, i2, ipos-2 ) + &
                              jay%f3( 2, i1, i2, ipos+3 )

         jay%f3( 3, i1, i2, ipos-2 ) = jay%f3( 3, i1, i2, ipos-2 ) - &
                              jay%f3( 3, i1, i2, ipos+2 )

         jay%f3( 1, i1, i2, ipos-1 ) = jay%f3( 1, i1, i2, ipos-1 ) + &
                              jay%f3( 1, i1, i2, ipos+2 )

         jay%f3( 2, i1, i2, ipos-1 ) = jay%f3( 2, i1, i2, ipos-1 ) + &
                              jay%f3( 2, i1, i2, ipos+2 )

         jay%f3( 3, i1, i2, ipos-1 ) = jay%f3( 3, i1, i2, ipos-1 ) - &
                              jay%f3( 3, i1, i2, ipos+1 )

         jay%f3( 1, i1, i2, ipos ) = jay%f3( 1, i1, i2, ipos ) + &
                              jay%f3( 1, i1, i2, ipos+1 )

         jay%f3( 2, i1, i2, ipos ) = jay%f3( 2, i1, i2, ipos ) + &
                              jay%f3( 2, i1, i2, ipos+1 )

         jay%f3( 3, i1, i2, ipos ) = 0
        enddo
      enddo

    end select
  end select

  end select

end subroutine spec_bc_3d
!-----------------------------------------------------------------------------------------

end module m_current_boundary
