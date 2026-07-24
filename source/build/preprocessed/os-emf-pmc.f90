# 1 "emf/boundary/os-emf-pmc.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "emf/boundary/os-emf-pmc.f90" 2
!-----------------------------------------------------------------------------------------
! (P)erfect (M)agnetic (C)onductor field boundaries
!-----------------------------------------------------------------------------------------
!
! These routines implement PMC field boundaries for the FDTD solver for 2 cases:
! pmc_bc_?d_?_odd - Boundary at the cell edge (lower edge [1], right edge [nx])
! pmc_bc_?d_?_even - Boundary at the cell center (center [0], center [nx])
!-----------------------------------------------------------------------------------------

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
# 11 "emf/boundary/os-emf-pmc.f90" 2
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
# 12 "emf/boundary/os-emf-pmc.f90" 2

module m_emf_pmc

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
# 16 "emf/boundary/os-emf-pmc.f90" 2

use m_parameters
use m_vdf_define

implicit none

private

interface pmc_bc_e_1d_odd
  module procedure pmc_bc_e_1d_odd
end interface

interface pmc_bc_e_2d_odd
  module procedure pmc_bc_e_2d_odd
end interface

interface pmc_bc_e_3d_odd
  module procedure pmc_bc_e_3d_odd
end interface

interface pmc_bc_e_1d_even
  module procedure pmc_bc_e_1d_even
end interface

interface pmc_bc_e_2d_even
  module procedure pmc_bc_e_2d_even
end interface

interface pmc_bc_e_3d_even
  module procedure pmc_bc_e_3d_even
end interface

interface pmc_bc_b_1d_odd
  module procedure pmc_bc_b_1d_odd
end interface

interface pmc_bc_b_2d_odd
  module procedure pmc_bc_b_2d_odd
end interface

interface pmc_bc_b_3d_odd
  module procedure pmc_bc_b_3d_odd
end interface

interface pmc_bc_b_1d_even
  module procedure pmc_bc_b_1d_even
end interface

interface pmc_bc_b_2d_even
  module procedure pmc_bc_b_2d_even
end interface

interface pmc_bc_b_3d_even
  module procedure pmc_bc_b_3d_even
end interface

public :: pmc_bc_e_1d_odd, pmc_bc_e_2d_odd, pmc_bc_e_3d_odd
public :: pmc_bc_b_1d_odd, pmc_bc_b_2d_odd, pmc_bc_b_3d_odd
public :: pmc_bc_e_1d_even, pmc_bc_e_2d_even, pmc_bc_e_3d_even
public :: pmc_bc_b_1d_even, pmc_bc_b_2d_even, pmc_bc_b_3d_even

contains

!-----------------------------------------------------------------------------------------
! pmc_bc_?d_?_odd - Boundary at the cell edge (lower edge [1], right edge [nx])
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_e_1d_odd( e, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: bnd_idx

  integer :: i1

  select case (bnd_idx)
     case ( p_lower )

       do i1 = lbound( e%f1, 2 ), 0
         e%f1( 1, i1 ) = -e%f1( 1, 1-i1 )

         e%f1( 2, i1 ) = e%f1( 2, 2-i1 )

         e%f1( 3, i1 ) = e%f1( 3, 2-i1 )

       enddo

     case ( p_upper )

       do i1 = 1, ubound( e%f1, 2 ) - e%nx_(1)
         e%f1( 1, e%nx_(1) + i1) = -e%f1( 1, e%nx_(1) + 1 - i1 )

         e%f1( 2, e%nx_(1) + i1) = e%f1( 2, e%nx_(1) + 2 - i1 )

         e%f1( 3, e%nx_(1) + i1) = e%f1( 3, e%nx_(1) + 2 - i1 )

       enddo

  end select

end subroutine pmc_bc_e_1d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_b_1d_odd( b, bnd_idx )
!-----------------------------------------------------------------------------------------

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: bnd_idx

  integer :: i1

  select case (bnd_idx)
     case ( p_lower )

       do i1 = lbound( b%f1, 2 ), 0
         b%f1( 1, i1 ) = b%f1( 1, 2-i1 )

         b%f1( 2, i1 ) = -b%f1( 2, 1-i1 )

         b%f1( 3, i1 ) = -b%f1( 3, 1-i1 )

       enddo

     case ( p_upper )
       do i1 = 1, ubound( b%f1, 2 ) - b%nx_(1)
         b%f1( 1, b%nx_(1) + i1) = b%f1( 1, b%nx_(1) + 2 - i1 )

         b%f1( 2, b%nx_(1) + i1) = -b%f1( 2, b%nx_(1) + 1 - i1 )

         b%f1( 3, b%nx_(1) + i1) = -b%f1( 3, b%nx_(1) + 1 - i1 )

       enddo

  end select

end subroutine pmc_bc_b_1d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_e_2d_odd( e, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
              do i1 = lbound( e%f2, 2 ), 0
                e%f2( 1, i1, i2 ) = -e%f2( 1, 1-i1, i2 )

                e%f2( 2, i1, i2 ) = e%f2( 2, 2-i1, i2 )

                e%f2( 3, i1, i2 ) = e%f2( 3, 2-i1, i2 )

              enddo
            enddo

          case(2) ! x2 boundary
            do i2 = lbound( e%f2, 3 ), 0
              do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
                e%f2( 1, i1, i2 ) = e%f2( 1, i1, 2-i2 )

                e%f2( 2, i1, i2 ) = -e%f2( 2, i1, 1-i2 )

                e%f2( 3, i1, i2 ) = e%f2( 3, i1, 2-i2 )

              enddo
            enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
              do i1 = 1, ubound( e%f2, 2 ) - e%nx_(1)
                e%f2( 1, e%nx_(1) + i1, i2) = -e%f2( 1, e%nx_(1) + 1 - i1, i2 )

                e%f2( 2, e%nx_(1) + i1, i2) = e%f2( 2, e%nx_(1) + 2 - i1, i2 )

                e%f2( 3, e%nx_(1) + i1, i2) = e%f2( 3, e%nx_(1) + 2 - i1, i2 )

              enddo
            enddo

          case(2) ! x2 boundary
            do i2 = 1, ubound( e%f2, 3 ) - e%nx_(2)
              do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
                e%f2( 1, i1, e%nx_(2) + i2) = e%f2( 1, i1, e%nx_(2) + 2 - i2 )

                e%f2( 2, i1, e%nx_(2) + i2) = -e%f2( 2, i1, e%nx_(2) + 1 - i2 )

                e%f2( 3, i1, e%nx_(2) + i2) = e%f2( 3, i1, e%nx_(2) + 2 - i2 )

              enddo
            enddo

       end select
  end select

end subroutine pmc_bc_e_2d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_b_2d_odd( b, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )

              do i1 = lbound( b%f2, 2 ), 0
                b%f2( 1, i1, i2 ) = b%f2( 1, 2-i1, i2 )

                b%f2( 2, i1, i2 ) = -b%f2( 2, 1-i1, i2 )

                b%f2( 3, i1, i2 ) = -b%f2( 3, 1-i1, i2 )

              enddo
            enddo

          case(2) ! x2 boundary
            do i2 = lbound( b%f2, 3 ), 0
              do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )

                b%f2( 1, i1, i2 ) = -b%f2( 1, i1, 1-i2 )

                b%f2( 2, i1, i2 ) = b%f2( 2, i1, 2-i2 )

                b%f2( 3, i1, i2 ) = -b%f2( 3, i1, 1-i2 )

              enddo
            enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )

              do i1 = 1, ubound( b%f2, 2 ) - b%nx_(1)
                b%f2( 1, b%nx_(1) + i1, i2 ) = b%f2( 1, b%nx_(1) + 2 - i1, i2 )

                b%f2( 2, b%nx_(1) + i1, i2 ) = -b%f2( 2, b%nx_(1) + 1 - i1, i2 )

                b%f2( 3, b%nx_(1) + i1, i2 ) = -b%f2( 3, b%nx_(1) + 1 - i1, i2 )

              enddo
            enddo

          case(2) ! x2 boundary
            do i2 = 1, ubound( b%f2, 3 ) - b%nx_(2)
              do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )

                b%f2( 1, i1, b%nx_(2) + i2 ) = -b%f2( 1, i1, b%nx_(2) + 1 - i2 )

                b%f2( 2, i1, b%nx_(2) + i2 ) = b%f2( 2, i1, b%nx_(2) + 2 - i2 )

                b%f2( 3, i1, b%nx_(2) + i2 ) = -b%f2( 3, i1, b%nx_(2) + 1 - i2 )

              enddo
            enddo

       end select
  end select

end subroutine pmc_bc_b_2d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_e_3d_odd( e, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
              do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                do i1 = lbound( e%f3, 2 ), 0
                  e%f3( 1, i1, i2, i3 ) = -e%f3( 1, 1-i1, i2, i3 )

                  e%f3( 2, i1, i2, i3 ) = e%f3( 2, 2-i1, i2, i3 )

                  e%f3( 3, i1, i2, i3 ) = e%f3( 3, 2-i1, i2, i3 )

                enddo
              enddo
            enddo

          case(2) ! x2 boundary
            do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
              do i2 = lbound( e%f3, 3 ), 0
                do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                  e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, 2-i2, i3 )

                  e%f3( 2, i1, i2, i3 ) = -e%f3( 2, i1, 1-i2, i3 )

                  e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, 2-i2, i3 )

                enddo
              enddo
            enddo

          case(3) ! x3 boundary
            do i3 = lbound( e%f3, 4 ), 0
              do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                  e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, 2-i3 )

                  e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, 2-i3 )

                  e%f3( 3, i1, i2, i3 ) = -e%f3( 3, i1, i2, 1-i3 )

                enddo
              enddo
            enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
              do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                do i1 = 1, ubound( e%f3, 2 ) - e%nx_(1)
                  e%f3( 1, e%nx_(1) + i1, i2, i3 ) = -e%f3( 1, e%nx_(1) + 1 - i1, i2, i3 )

                  e%f3( 2, e%nx_(1) + i1, i2, i3 ) = e%f3( 2, e%nx_(1) + 2 - i1, i2, i3 )

                  e%f3( 3, e%nx_(1) + i1, i2, i3 ) = e%f3( 3, e%nx_(1) + 2 - i1, i2, i3 )

                enddo
              enddo
            enddo

          case(2) ! x2 boundary
            do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
              do i2 = 1, ubound( e%f3, 3 ) - e%nx_(2)
                do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                  e%f3( 1, i1, e%nx_(2) + i2, i3 ) = e%f3( 1, i1, e%nx_(2) + 2 - i2, i3 )

                  e%f3( 2, i1, e%nx_(2) + i2, i3 ) = -e%f3( 2, i1, e%nx_(2) + 1 - i2, i3 )

                  e%f3( 3, i1, e%nx_(2) + i2, i3 ) = e%f3( 3, i1, e%nx_(2) + 2 - i2, i3 )

                enddo
              enddo
            enddo

          case(3) ! x3 boundary
            do i3 = 1, ubound( e%f3, 4 ) - e%nx_(3)
              do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                  e%f3( 1, i1, i2, e%nx_(3) + i3 ) = e%f3( 1, i1, i2, e%nx_(3) + 2 - i3 )

                  e%f3( 2, i1, i2, e%nx_(3) + i3 ) = e%f3( 2, i1, i2, e%nx_(3) + 2 - i3 )

                  e%f3( 3, i1, i2, e%nx_(3) + i3 ) = -e%f3( 3, i1, i2, e%nx_(3) + 1 - i3 )

                enddo
              enddo
            enddo

       end select
  end select

end subroutine pmc_bc_e_3d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_b_3d_odd( b, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
              do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )

                do i1 = lbound( b%f3, 2 ), 0
                  b%f3( 1, i1, i2, i3 ) = b%f3( 1, 2-i1, i2, i3 )

                  b%f3( 2, i1, i2, i3 ) = -b%f3( 2, 1-i1, i2, i3 )

                  b%f3( 3, i1, i2, i3 ) = -b%f3( 3, 1-i1, i2, i3 )

                enddo
              enddo
            enddo

          case(2) ! x2 boundary
            do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
              do i2 = lbound( b%f3, 3 ), 0
                do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                  b%f3( 1, i1, i2, i3 ) = -b%f3( 1, i1, 1-i2, i3 )

                  b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, 2-i2, i3 )

                  b%f3( 3, i1, i2, i3 ) = -b%f3( 3, i1, 1-i2, i3 )

                enddo
              enddo
            enddo

          case(3) ! x3 boundary
            do i3 = lbound( b%f3, 4 ), 0
              do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                  b%f3( 1, i1, i2, i3 ) = -b%f3( 1, i1, i2, 1-i3 )

                  b%f3( 2, i1, i2, i3 ) = -b%f3( 2, i1, i2, 1-i3 )

                  b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, 2-i3 )

                enddo
              enddo
            enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
              do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )

                do i1 = 1, ubound( b%f3, 2 ) - b%nx_(1)
                  b%f3( 1, b%nx_(1) + i1, i2, i3 ) = b%f3( 1, b%nx_(1) + 2 - i1, i2, i3 )

                  b%f3( 2, b%nx_(1) + i1, i2, i3 ) = -b%f3( 2, b%nx_(1) + 1 - i1, i2, i3 )

                  b%f3( 3, b%nx_(1) + i1, i2, i3 ) = -b%f3( 3, b%nx_(1) + 1 - i1, i2, i3 )

                enddo
              enddo
            enddo

          case(2) ! x2 boundary
            do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
              do i2 = 1, ubound( b%f3, 3 ) - b%nx_(2)
                do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                  b%f3( 1, i1, b%nx_(2) + i2, i3 ) = -b%f3( 1, i1, b%nx_(2) + 1 - i2, i3 )

                  b%f3( 2, i1, b%nx_(2) + i2, i3 ) = b%f3( 2, i1, b%nx_(2) + 2 - i2, i3 )

                  b%f3( 3, i1, b%nx_(2) + i2, i3 ) = -b%f3( 3, i1, b%nx_(2) + 1 - i2, i3 )

                enddo
              enddo
            enddo

          case(3) ! x3 boundary
            do i3 = 1, ubound( b%f3, 4 ) - b%nx_(3)
              do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                  b%f3( 1, i1, i2, b%nx_(3) + i3 ) = -b%f3( 1, i1, i2, b%nx_(3) + 1 - i3 )

                  b%f3( 2, i1, i2, b%nx_(3) + i3 ) = -b%f3( 2, i1, i2, b%nx_(3) + 1 - i3 )

                  b%f3( 3, i1, i2, b%nx_(3) + i3 ) = b%f3( 3, i1, i2, b%nx_(3) + 2 - i3 )

                enddo
              enddo
            enddo

       end select
  end select

end subroutine pmc_bc_b_3d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! pmc_bc_?d_?_even - Boundary at the cell center (center [0], center [nx])
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_e_1d_even( e, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: bnd_idx

  integer :: i1

  select case (bnd_idx)
     case ( p_lower )

       do i1 = lbound( e%f1, 2 ), -1
         e%f1( 1, i1 ) = -e%f1( 1, -i1 )

         e%f1( 2, i1 ) = e%f1( 2, 1-i1 )

         e%f1( 3, i1 ) = e%f1( 3, 1-i1 )

       enddo

       e%f1( 1, 0 ) = 0
       e%f1( 2, 0 ) = e%f1( 2, 1 )
       e%f1( 3, 0 ) = e%f1( 3, 1 )

     case ( p_upper )

       e%f1( 1, e%nx_(1) ) = 0

       do i1 = 1, ubound( e%f1, 2 ) - e%nx_(1)
         e%f1( 1, e%nx_(1) + i1) = -e%f1( 1, e%nx_(1) - i1 )

         e%f1( 2, e%nx_(1) + i1) = e%f1( 2, e%nx_(1) + 1 - i1 )

         e%f1( 3, e%nx_(1) + i1) = e%f1( 3, e%nx_(1) + 1 - i1 )

       enddo

  end select

end subroutine pmc_bc_e_1d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_b_1d_even( b, bnd_idx )
!-----------------------------------------------------------------------------------------

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: bnd_idx

  integer :: i1

  select case (bnd_idx)
     case ( p_lower )

       do i1 = lbound( b%f1, 2 ), -1
         b%f1( 1, i1 ) = b%f1( 1, 1-i1 )

         b%f1( 2, i1 ) = -b%f1( 2, -i1 )

         b%f1( 3, i1 ) = -b%f1( 3, -i1 )

       enddo

       b%f1( 1, 0 ) = b%f1( 1, 1 )
       b%f1( 2, 0 ) = 0
       b%f1( 3, 0 ) = 0

     case ( p_upper )

       b%f1( 2, b%nx_(1) ) = 0

       b%f1( 3, b%nx_(1) ) = 0

       do i1 = 1, ubound( b%f1, 2 ) - b%nx_(1)
         b%f1( 1, b%nx_(1) + i1) = b%f1( 1, b%nx_(1) + 1 - i1 )

         b%f1( 2, b%nx_(1) + i1) = -b%f1( 2, b%nx_(1) - i1 )

         b%f1( 3, b%nx_(1) + i1) = -b%f1( 3, b%nx_(1) - i1 )

       enddo

  end select

end subroutine pmc_bc_b_1d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_e_2d_even( e, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
              do i1 = lbound( e%f2, 2 ), -1
                e%f2( 1, i1, i2 ) = -e%f2( 1, -i1, i2 )

                e%f2( 2, i1, i2 ) = e%f2( 2, 1-i1, i2 )

                e%f2( 3, i1, i2 ) = e%f2( 3, 1-i1, i2 )

              enddo

              e%f2( 1, 0, i2 ) = 0
              e%f2( 2, 0, i2 ) = e%f2( 2, 1, i2 )
              e%f2( 3, 0, i2 ) = e%f2( 3, 1, i2 )
            enddo

          case(2) ! x2 boundary
            do i2 = lbound( e%f2, 3 ), -1
              do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
                e%f2( 1, i1, i2 ) = e%f2( 1, i1, 1-i2 )

                e%f2( 2, i1, i2 ) = -e%f2( 2, i1, -i2 )

                e%f2( 3, i1, i2 ) = e%f2( 3, i1, 1-i2 )

              enddo
            enddo

            do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
              e%f2( 1, i1, 0 ) = e%f2( 1, i1, 1 )

              e%f2( 2, i1, 0 ) = 0

              e%f2( 3, i1, 0 ) = e%f2( 3, i1, 1 )

            enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
              e%f2( 1, e%nx_(1), i2) = 0
              do i1 = 1, ubound( e%f2, 2 ) - e%nx_(1)
                e%f2( 1, e%nx_(1) + i1, i2) = -e%f2( 1, e%nx_(1) - i1, i2 )

                e%f2( 2, e%nx_(1) + i1, i2) = e%f2( 2, e%nx_(1) + 1 - i1, i2 )

                e%f2( 3, e%nx_(1) + i1, i2) = e%f2( 3, e%nx_(1) + 1 - i1, i2 )

              enddo
            enddo

          case(2) ! x2 boundary
            do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
              e%f2( 2, i1, e%nx_(2) ) = 0

            enddo

            do i2 = 1, ubound( e%f2, 3 ) - e%nx_(2)
              do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
                e%f2( 1, i1, e%nx_(2) + i2) = e%f2( 1, i1, e%nx_(2) + 1 - i2 )

                e%f2( 2, i1, e%nx_(2) + i2) = -e%f2( 2, i1, e%nx_(2) - i2 )

                e%f2( 3, i1, e%nx_(2) + i2) = e%f2( 3, i1, e%nx_(2) + 1 - i2 )

              enddo
            enddo

       end select
  end select

end subroutine pmc_bc_e_2d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_b_2d_even( b, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )

              do i1 = lbound( b%f2, 2 ), -1
                b%f2( 1, i1, i2 ) = b%f2( 1, 1-i1, i2 )

                b%f2( 2, i1, i2 ) = -b%f2( 2, -i1, i2 )

                b%f2( 3, i1, i2 ) = -b%f2( 3, -i1, i2 )

              enddo

              b%f2( 1, 0, i2 ) = b%f2( 1, 1, i2 )
              b%f2( 2, 0, i2 ) = 0
              b%f2( 3, 0, i2 ) = 0
            enddo

          case(2) ! x2 boundary
            do i2 = lbound( b%f2, 3 ), -1
              do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )

                b%f2( 1, i1, i2 ) = -b%f2( 1, i1, -i2 )

                b%f2( 2, i1, i2 ) = b%f2( 2, i1, 1-i2 )

                b%f2( 3, i1, i2 ) = -b%f2( 3, i1, -i2 )

              enddo
            enddo

            do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )

              b%f2( 1, i1, 0 ) = 0
              b%f2( 2, i1, 0 ) = b%f2( 2, i1, 1 )
              b%f2( 3, i1, 0 ) = 0
            enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )

              b%f2( 2, b%nx_(1), i2 ) = 0

              b%f2( 3, b%nx_(1), i2 ) = 0

              do i1 = 1, ubound( b%f2, 2 ) - b%nx_(1)
                b%f2( 1, b%nx_(1) + i1, i2 ) = b%f2( 1, b%nx_(1) + 1 - i1, i2 )

                b%f2( 2, b%nx_(1) + i1, i2 ) = -b%f2( 2, b%nx_(1) - i1, i2 )

                b%f2( 3, b%nx_(1) + i1, i2 ) = -b%f2( 3, b%nx_(1) - i1, i2 )

              enddo
            enddo

          case(2) ! x2 boundary
            do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )

              b%f2( 1, i1, b%nx_(2) ) = 0
              b%f2( 3, i1, b%nx_(2) ) = 0
            enddo

            do i2 = 1, ubound( b%f2, 3 ) - b%nx_(2)
              do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )

                b%f2( 1, i1, b%nx_(2) + i2 ) = -b%f2( 1, i1, b%nx_(2) - i2 )

                b%f2( 2, i1, b%nx_(2) + i2 ) = b%f2( 2, i1, b%nx_(2) + 1 - i2 )

                b%f2( 3, i1, b%nx_(2) + i2 ) = -b%f2( 3, i1, b%nx_(2) - i2 )

              enddo
            enddo

       end select
  end select

end subroutine pmc_bc_b_2d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_e_3d_even( e, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
              do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                do i1 = lbound( e%f3, 2 ), -1
                  e%f3( 1, i1, i2, i3 ) = -e%f3( 1, -i1, i2, i3 )

                  e%f3( 2, i1, i2, i3 ) = e%f3( 2, 1-i1, i2, i3 )

                  e%f3( 3, i1, i2, i3 ) = e%f3( 3, 1-i1, i2, i3 )

                enddo
                e%f3( 1, 0, i2, i3 ) = 0
                e%f3( 2, 0, i2, i3 ) = e%f3( 2, 1, i2, i3 )
                e%f3( 3, 0, i2, i3 ) = e%f3( 3, 1, i2, i3 )
              enddo
            enddo

          case(2) ! x2 boundary
            do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
              do i2 = lbound( e%f3, 3 ), -1
                do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                  e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, 1-i2, i3 )

                  e%f3( 2, i1, i2, i3 ) = -e%f3( 2, i1, -i2, i3 )

                  e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, 1-i2, i3 )

                enddo
              enddo
              do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                 e%f3( 1, i1, 0, i3 ) = e%f3( 1, i1, 1, i3 )

                 e%f3( 2, i1, 0, i3 ) = 0

                 e%f3( 3, i1, 0, i3 ) = e%f3( 3, i1, 1, i3 )

              enddo
            enddo

          case(3) ! x3 boundary
            do i3 = lbound( e%f3, 4 ), -1
              do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                  e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, 1-i3 )

                  e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, 1-i3 )

                  e%f3( 3, i1, i2, i3 ) = -e%f3( 3, i1, i2, -i3 )

                enddo
              enddo
            enddo

            do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
              do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                e%f3( 1, i1, i2, 0 ) = e%f3( 1, i1, i2, 1 )

                e%f3( 2, i1, i2, 0 ) = e%f3( 2, i1, i2, 1 )

                e%f3( 3, i1, i2, 0 ) = 0

              enddo
            enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
              do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                e%f3( 1, e%nx_(1), i2, i3 ) = 0
                do i1 = 1, ubound( e%f3, 2 ) - e%nx_(1)
                  e%f3( 1, e%nx_(1) + i1, i2, i3 ) = -e%f3( 1, e%nx_(1) - i1, i2, i3 )

                  e%f3( 2, e%nx_(1) + i1, i2, i3 ) = e%f3( 2, e%nx_(1) + 1 - i1, i2, i3 )

                  e%f3( 3, e%nx_(1) + i1, i2, i3 ) = e%f3( 3, e%nx_(1) + 1 - i1, i2, i3 )

                enddo
              enddo
            enddo

          case(2) ! x2 boundary
            do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
              do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                e%f3( 2, i1, e%nx_(2), i3 ) = 0
              enddo

              do i2 = 1, ubound( e%f3, 3 ) - e%nx_(2)
                do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                  e%f3( 1, i1, e%nx_(2) + i2, i3 ) = e%f3( 1, i1, e%nx_(2) + 1 - i2, i3 )

                  e%f3( 2, i1, e%nx_(2) + i2, i3 ) = -e%f3( 2, i1, e%nx_(2) - i2, i3 )

                  e%f3( 3, i1, e%nx_(2) + i2, i3 ) = e%f3( 3, i1, e%nx_(2) + 1 - i2, i3 )

                enddo
              enddo
            enddo

          case(3) ! x3 boundary
            do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
              do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                e%f3( 3, i1, i2, e%nx_(3) ) = 0

              enddo
            enddo

            do i3 = 1, ubound( e%f3, 4 ) - e%nx_(3)
              do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                  e%f3( 1, i1, i2, e%nx_(3) + i3 ) = e%f3( 1, i1, i2, e%nx_(3) + 1 - i3 )

                  e%f3( 2, i1, i2, e%nx_(3) + i3 ) = e%f3( 2, i1, i2, e%nx_(3) + 1 - i3 )

                  e%f3( 3, i1, i2, e%nx_(3) + i3 ) = -e%f3( 3, i1, i2, e%nx_(3) - i3 )

                enddo
              enddo
            enddo

       end select
  end select

end subroutine pmc_bc_e_3d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pmc_bc_b_3d_even( b, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
              do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )

                do i1 = lbound( b%f3, 2 ), -1
                  b%f3( 1, i1, i2, i3 ) = b%f3( 1, 1-i1, i2, i3 )

                  b%f3( 2, i1, i2, i3 ) = -b%f3( 2, -i1, i2, i3 )

                  b%f3( 3, i1, i2, i3 ) = -b%f3( 3, -i1, i2, i3 )

                enddo

                b%f3( 1, 0, i2, i3 ) = b%f3( 1, 1, i2, i3 )
                b%f3( 2, 0, i2, i3 ) = 0
                b%f3( 3, 0, i2, i3 ) = 0
              enddo
            enddo

          case(2) ! x2 boundary
            do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
              do i2 = lbound( b%f3, 3 ), -1
                do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                  b%f3( 1, i1, i2, i3 ) = -b%f3( 1, i1, -i2, i3 )

                  b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, 1-i2, i3 )

                  b%f3( 3, i1, i2, i3 ) = -b%f3( 3, i1, -i2, i3 )

                enddo
              enddo
              do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                b%f3( 1, i1, 0, i3 ) = 0

                b%f3( 2, i1, 0, i3 ) = b%f3( 2, i1, 1, i3 )

                b%f3( 3, i1, 0, i3 ) = 0

              enddo
            enddo

          case(3) ! x3 boundary
            do i3 = lbound( b%f3, 4 ), -1
              do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                  b%f3( 1, i1, i2, i3 ) = -b%f3( 1, i1, i2, -i3 )

                  b%f3( 2, i1, i2, i3 ) = -b%f3( 2, i1, i2, -i3 )

                  b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, 1-i3 )

                enddo
              enddo
            enddo

            do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
              do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                b%f3( 1, i1, i2, 0 ) = 0

                b%f3( 2, i1, i2, 0 ) = 0

                b%f3( 3, i1, i2, 0 ) = b%f3( 3, i1, i2, 1 )

              enddo
            enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
            do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
              do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )

                b%f3( 2, b%nx_(1), i2, i3 ) = 0

                b%f3( 3, b%nx_(1), i2, i3 ) = 0

                do i1 = 1, ubound( b%f3, 2 ) - b%nx_(1)
                  b%f3( 1, b%nx_(1) + i1, i2, i3 ) = b%f3( 1, b%nx_(1) + 1 - i1, i2, i3 )

                  b%f3( 2, b%nx_(1) + i1, i2, i3 ) = -b%f3( 2, b%nx_(1) - i1, i2, i3 )

                  b%f3( 3, b%nx_(1) + i1, i2, i3 ) = -b%f3( 3, b%nx_(1) - i1, i2, i3 )

                enddo
              enddo
            enddo

          case(2) ! x2 boundary
            do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
              do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                b%f3( 1, i1, b%nx_(2), i3 ) = 0

                b%f3( 3, i1, b%nx_(2), i3 ) = 0
              enddo

              do i2 = 1, ubound( b%f3, 3 ) - b%nx_(2)
                do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                  b%f3( 1, i1, b%nx_(2) + i2, i3 ) = -b%f3( 1, i1, b%nx_(2) - i2, i3 )

                  b%f3( 2, i1, b%nx_(2) + i2, i3 ) = b%f3( 2, i1, b%nx_(2) + 1 - i2, i3 )

                  b%f3( 3, i1, b%nx_(2) + i2, i3 ) = -b%f3( 3, i1, b%nx_(2) - i2, i3 )

                enddo
              enddo
            enddo

          case(3) ! x3 boundary
            do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
              do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                b%f3( 1, i1, i2, b%nx_(3) ) = 0

                b%f3( 2, i1, i2, b%nx_(3) ) = 0
              enddo
            enddo

            do i3 = 1, ubound( b%f3, 4 ) - b%nx_(3)
              do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )

                  b%f3( 1, i1, i2, b%nx_(3) + i3 ) = -b%f3( 1, i1, i2, b%nx_(3) - i3 )

                  b%f3( 2, i1, i2, b%nx_(3) + i3 ) = -b%f3( 2, i1, i2, b%nx_(3) - i3 )

                  b%f3( 3, i1, i2, b%nx_(3) + i3 ) = b%f3( 3, i1, i2, b%nx_(3) + 1 - i3 )

                enddo
              enddo
            enddo

       end select
  end select

end subroutine pmc_bc_b_3d_even
!-----------------------------------------------------------------------------------------

end module m_emf_pmc
