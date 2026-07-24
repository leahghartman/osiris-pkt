# 1 "emf/os-emf-gridval.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "emf/os-emf-gridval.f03" 2
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
# 2 "emf/os-emf-gridval.f03" 2
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
# 3 "emf/os-emf-gridval.f03" 2


module m_emf_gridval

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
# 8 "emf/os-emf-gridval.f03" 2

use m_parameters

implicit none

private

character(len=*), parameter :: p_emf_gridval_rst_id = "emf gridval rst data - 0x0000"

interface dipole_field
  module procedure dipole_field_
end interface

public :: p_emf_gridval_rst_id, dipole_field

contains

!---------------------------------------------------------------------------------------------------
function dipole_field_( p, x0, r0, dir, x )

implicit none

! dummy variables

  real(p_k_fld) :: dipole_field_

  real(p_k_fld), dimension(p_f_dim), intent(in) :: p ! dipole moment
  real(p_k_fld), dimension(p_x_dim), intent(in) :: x0 ! dipole position
  real(p_k_fld), intent(in) :: r0 ! minimum distance to calc. field
  integer, intent(in) :: dir ! field component to return
  real(p_double), dimension(p_x_dim), intent(in) :: x ! field position

! local variables

  real(p_k_fld), dimension(p_f_dim) :: n
  real(p_k_fld) :: mod_r, p_int_n

! executable statements

  n = 0.0_p_k_fld
  n(1:p_x_dim) = real(x,p_k_fld) - x0
  mod_r = sqrt(sum(n**2))

  if (mod_r >= r0) then
    n = n/mod_r
    p_int_n = sum(n*p)

    dipole_field_ = (3 * n(dir) * p_int_n - p(dir)) / ( mod_r**3 )
  else
    dipole_field_ = 8.37758041_p_k_fld * p(dir) ! 8 * pi / 3
  endif

end function dipole_field_
!---------------------------------------------------------------------------------------------------

end module m_emf_gridval

!---------------------------------------------------------------------------------------------------
subroutine set_fld_values_emf_gridval( this, e, b, g_space, nx_p_min, t, lbin, ubin )
!---------------------------------------------------------------------------------------------------

  use m_emf_define, only: t_emf_gridval, p_emf_uniform, p_emf_math, p_emf_dipole
  use m_space, only: t_space, xmin, xmax
  use m_fparser, only: p_k_fparse, eval
  use m_parameters
  use m_vdf_define, only: t_vdf
  use m_emf_gridval
  use m_emf_define, only: p_emf_python
  use m_callpy

  implicit none

  ! dummy variables

  class(t_emf_gridval), intent( inout ) :: this
  type(t_vdf), intent(inout) :: e, b

  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real( p_double ), intent(in) :: t

  integer, dimension(:), intent(in), optional :: lbin
  integer, dimension(:), intent(in), optional :: ubin


  ! local variables
  real(p_double), dimension(p_max_dim) :: g_xmin, ldx
  real(p_k_fparse), dimension(p_max_dim+1) :: x_eval

  real(p_double), dimension(p_max_dim,p_f_dim) :: odx_b
  real(p_double), dimension(p_max_dim,p_f_dim) :: odx_e
  integer :: i, j1, j2, j3
  integer :: x_dim

  integer, dimension(p_max_dim) :: lb, ub
  integer, dimension( p_max_dim ) :: shift

  real(p_double), dimension(2,p_max_dim) :: x_bnd
  character(2) :: fld

  ! executable statements

  g_xmin(1:p_x_dim) = xmin(g_space)
  ldx(1:p_x_dim) = b%dx()

  shift(1:p_x_dim) = nx_p_min(1:p_x_dim) - 1

  select case (p_x_dim)
    case(1)
      x_dim = 1

      ! get boundaries to inject
      if (present(lbin)) then
        lb(1) = lbin(1)
      else
        lb(1) = lbound( b%f1, 2)
      endif
      if (present(ubin)) then
        ub(1) = ubin(1)
      else
        ub(1) = ubound( b%f1, 2)
      endif

      ! set the grid offsets

      odx_b(1,1) = 0.0_p_double*ldx(1)
      odx_b(1,2) = 0.5_p_double*ldx(1)
      odx_b(1,3) = 0.5_p_double*ldx(1)

      odx_e(1,1) = 0.5_p_double*ldx(1)
      odx_e(1,2) = 0.0_p_double*ldx(1)
      odx_e(1,3) = 0.0_p_double*ldx(1)

      select case( this%interpolation )
      case( p_linear, p_cubic )
        ! Nothing to change here
      case( p_quadratic, p_quartic )
        do j2 = 1, p_f_dim
          do j1 = 1, p_x_dim
            odx_b(j1,j2) = odx_b(j1,j2) + 0.5_p_double*ldx(j1)
            odx_e(j1,j2) = odx_e(j1,j2) + 0.5_p_double*ldx(j1)
          enddo
        enddo
      end select

      do i= 1, p_f_dim
        ! set the external b field

         select case (this%type_b(i))
          case(p_emf_uniform)
            ! set the external field to the supplied uniform value
            do j1 = lb(1), ub(1)
              b%f1(i,j1) = this%uniform_b0(i)
            enddo

          case(p_emf_math)

            ! set the external field to the function value
            ! also set the value in the guard cells
            x_eval(2) = t
            do j1 = lb(1), ub(1)
                ! set the proper correction to account for field positions
                x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_b(1, i)
                b%f1(i,j1) = real( eval( this%mfunc_b(i), x_eval ), p_k_fld )
            enddo

          case(p_emf_dipole)

            ! dipole field not allowed in 1D

            write(err_buf__,*) "dipole field not allowed in 1D";call err__("emf/os-emf-gridval.f03",178)
            call abort_program( p_err_invalid )

          case(p_emf_python)

            x_bnd(1,1) = g_xmin(1) + (lb(1)-1+shift(1))*ldx(1) + odx_b(1,i)

            x_bnd(2,1) = g_xmin(1) + (ub(1)-1+shift(1))*ldx(1) + odx_b(1,i)

            fld = ""
            write(fld,'(A,I1)') "b", i
            call set_state( "data", b%f1(i,lb(1):ub(1)) )
            call set_state( "fld", trim(fld) )
            call set_state( "x_bnd", x_bnd(:,1:p_x_dim) )
            call set_state( "g_xmin", g_xmin(1:p_x_dim) )
            call set_state( "g_xmax", xmax(g_space) )
            call set_state( "nx_p_min", nx_p_min(1:p_x_dim) )
            call set_state( "t", (/t/) )

            call call_function( trim(this%py_mod), trim(this%py_func) )

            call get_state( "data", b%f1(i,lb(1):ub(1)) )

        end select

        ! set the external e field

        select case (this%type_e(i))

          case(p_emf_uniform)
            ! set the external field to the supplied uniform value
            do j1 = lb(1), ub(1)
              e%f1(i,j1) = this%uniform_e0(i)
            enddo

          case(p_emf_math)

            ! set the external field to the function value
            ! also set the value in the guard cells
            x_eval(2) = t
            do j1 = lb(1), ub(1)

                ! set the proper correction to account for field positions
                x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_e(1, i)
                e%f1(i,j1) = real( eval( this%mfunc_e(i), x_eval ), p_k_fld )

            enddo

          case(p_emf_dipole)

            ! dipole field not allowed in 1D

            write(err_buf__,*) "dipole field not allowed in 1D";call err__("emf/os-emf-gridval.f03",230)
            call abort_program( p_err_invalid )

          case(p_emf_python)

            x_bnd(1,1) = g_xmin(1) + (lb(1)-1+shift(1))*ldx(1) + odx_e(1,i)

            x_bnd(2,1) = g_xmin(1) + (ub(1)-1+shift(1))*ldx(1) + odx_e(1,i)

            fld = ""
            write(fld,'(A,I1)') "e", i
            call set_state( "data", e%f1(i,lb(1):ub(1)) )
            call set_state( "fld", trim(fld) )
            call set_state( "x_bnd", x_bnd(:,1:p_x_dim) )
            call set_state( "g_xmin", g_xmin(1:p_x_dim) )
            call set_state( "g_xmax", xmax(g_space) )
            call set_state( "nx_p_min", nx_p_min(1:p_x_dim) )
            call set_state( "t", (/t/) )

            call call_function( trim(this%py_mod), trim(this%py_func) )

            call get_state( "data", e%f1(i,lb(1):ub(1)) )

        end select

      enddo ! i= 1, p_f_dim

    case(2)

      x_dim = 2

      ! get boundaries to inject
      if (present(lbin)) then
        lb(1:2) = lbin(1:2)
      else
        lb(1) = lbound( b%f2, 2)
        lb(2) = lbound( b%f2, 3)
      endif
      if (present(ubin)) then
        ub(1:2) = ubin(1:2)
      else
        ub(1) = ubound( b%f2, 2)
        ub(2) = ubound( b%f2, 3)
      endif

      ! set the grid offsets

      odx_b(1,1) = 0.0_p_double*ldx(1)
      odx_b(2,1) = 0.5_p_double*ldx(2)

      odx_b(1,2) = 0.5_p_double*ldx(1)
      odx_b(2,2) = 0.0_p_double*ldx(2)

      odx_b(1,3) = 0.5_p_double*ldx(1)
      odx_b(2,3) = 0.5_p_double*ldx(2)

      odx_e(1,1) = 0.5_p_double*ldx(1)
      odx_e(2,1) = 0.0_p_double*ldx(2)

      odx_e(1,2) = 0.0_p_double*ldx(1)
      odx_e(2,2) = 0.5_p_double*ldx(2)

      odx_e(1,3) = 0.0_p_double*ldx(1)
      odx_e(2,3) = 0.0_p_double*ldx(2)

      select case( this%interpolation )
      case( p_linear, p_cubic )
        ! Nothing to change here
      case( p_quadratic, p_quartic )
        do j2 = 1, p_f_dim
          do j1 = 1, p_x_dim
            odx_b(j1,j2) = odx_b(j1,j2) + 0.5_p_double*ldx(j1)
            odx_e(j1,j2) = odx_e(j1,j2) + 0.5_p_double*ldx(j1)
          enddo
        enddo
      end select

      do i= 1, p_f_dim
        ! set the external b field

         select case (this%type_b(i))
          case(p_emf_uniform)
            ! set the external field to the supplied uniform value
            do j2 = lb(2), ub(2)
              do j1 = lb(1), ub(1)
                 b%f2(i,j1,j2) = this%uniform_b0(i)
              enddo
            enddo

          case(p_emf_math)

            ! set the external field to the function value
            ! also set the value in the guard cells
            x_eval(3) = t
            do j2 = lb(2), ub(2)
              do j1 = lb(1), ub(1)

                ! set the proper correction to account for field positions
                x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_b(1, i)
                x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(x_dim) + odx_b(2, i)

                b%f2(i,j1,j2) = real( eval( this%mfunc_b(i), x_eval ), p_k_fld )

              enddo
            enddo

          case(p_emf_dipole)

            ! set the external field to a dipole field

            do j2 = lb(2), ub(2)
              do j1 = lb(1), ub(1)

                ! set the proper correction to account for field positions
                x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_b(1, i)
                x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_b(2, i)

                b%f2(i,j1,j2) = dipole_field( this%dipole_b_m, &
                                              this%dipole_b_x0, &
                                              this%dipole_b_r0, &
                                              i, x_eval )

              enddo
            enddo

          case(p_emf_python)

            x_bnd(1,1) = g_xmin(1) + (lb(1)-1+shift(1))*ldx(1) + odx_b(1,i)
            x_bnd(1,2) = g_xmin(2) + (lb(2)-1+shift(2))*ldx(2) + odx_b(2,i)

            x_bnd(2,1) = g_xmin(1) + (ub(1)-1+shift(1))*ldx(1) + odx_b(1,i)
            x_bnd(2,2) = g_xmin(2) + (ub(2)-1+shift(2))*ldx(2) + odx_b(2,i)

            fld = ""
            write(fld,'(A,I1)') "b", i
            call set_state( "data", b%f2(i,lb(1):ub(1),lb(2):ub(2)) )
            call set_state( "fld", trim(fld) )
            call set_state( "x_bnd", x_bnd(:,1:p_x_dim) )
            call set_state( "g_xmin", g_xmin(1:p_x_dim) )
            call set_state( "g_xmax", xmax(g_space) )
            call set_state( "nx_p_min", nx_p_min(1:p_x_dim) )
            call set_state( "t", (/t/) )

            call call_function( trim(this%py_mod), trim(this%py_func) )

            call get_state( "data", b%f2(i,lb(1):ub(1),lb(2):ub(2)) )

        end select

        ! set the external e field

        select case (this%type_e(i))

          case(p_emf_uniform)
            ! set the external field to the supplied uniform value
            do j2 = lb(2), ub(2)
              do j1 = lb(1), ub(1)
                 e%f2(i,j1,j2) = this%uniform_e0(i)
              enddo
            enddo

          case(p_emf_math)

            ! set the external field to the function value
            ! also set the value in the guard cells
            x_eval(3) = t
            do j2 = lb(2), ub(2)
              do j1 = lb(1), ub(1)

                ! set the proper correction to account for field positions
                x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_e(1, i)
                x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_e(2, i)

                e%f2(i,j1,j2) = real( eval( this%mfunc_e(i), x_eval ), p_k_fld )
              enddo
            enddo

          case(p_emf_dipole)

            ! set the external field to a dipole field
            do j2 = lb(2), ub(2)
              do j1 = lb(1), ub(1)

                ! set the proper correction to account for field positions
                x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_e(1, i)
                x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_e(2, i)

                e%f2(i,j1,j2) = 0.5_p_k_fld * dipole_field( this%dipole_e_p, &
                                                             this%dipole_e_x0, &
                                                             this%dipole_e_r0, &
                                                             i, x_eval )

              enddo
            enddo

          case(p_emf_python)

            x_bnd(1,1) = g_xmin(1) + (lb(1)-1+shift(1))*ldx(1) + odx_e(1,i)
            x_bnd(1,2) = g_xmin(2) + (lb(2)-1+shift(2))*ldx(2) + odx_e(2,i)

            x_bnd(2,1) = g_xmin(1) + (ub(1)-1+shift(1))*ldx(1) + odx_e(1,i)
            x_bnd(2,2) = g_xmin(2) + (ub(2)-1+shift(2))*ldx(2) + odx_e(2,i)

            fld = ""
            write(fld,'(A,I1)') "e", i
            call set_state( "data", e%f2(i,lb(1):ub(1),lb(2):ub(2)) )
            call set_state( "fld", trim(fld) )
            call set_state( "x_bnd", x_bnd(:,1:p_x_dim) )
            call set_state( "g_xmin", g_xmin(1:p_x_dim) )
            call set_state( "g_xmax", xmax(g_space) )
            call set_state( "nx_p_min", nx_p_min(1:p_x_dim) )
            call set_state( "t", (/t/) )

            call call_function( trim(this%py_mod), trim(this%py_func) )

            call get_state( "data", e%f2(i,lb(1):ub(1),lb(2):ub(2)) )

        end select

      enddo ! i= 1, p_f_dim

    case(3)

      x_dim = 3

      ! get boundaries to inject
      if (present(lbin)) then
        lb(1:3) = lbin(1:3)
      else
        lb(1) = lbound( b%f3, 2)
        lb(2) = lbound( b%f3, 3)
        lb(3) = lbound( b%f3, 4)
      endif
      if (present(ubin)) then
        ub(1:3) = ubin(1:3)
      else
        ub(1) = ubound( b%f3, 2)
        ub(2) = ubound( b%f3, 3)
        ub(3) = ubound( b%f3, 4)
      endif

      ! set the grid offsets

      odx_b(1,1) = 0.0_p_double*ldx(1)
      odx_b(2,1) = 0.5_p_double*ldx(2)
      odx_b(3,1) = 0.5_p_double*ldx(3)

      odx_b(1,2) = 0.5_p_double*ldx(1)
      odx_b(2,2) = 0.0_p_double*ldx(2)
      odx_b(3,2) = 0.5_p_double*ldx(3)

      odx_b(1,3) = 0.5_p_double*ldx(1)
      odx_b(2,3) = 0.5_p_double*ldx(2)
      odx_b(3,3) = 0.0_p_double*ldx(3)

      odx_e(1,1) = 0.5_p_double*ldx(1)
      odx_e(2,1) = 0.0_p_double*ldx(2)
      odx_e(3,1) = 0.0_p_double*ldx(3)

      odx_e(1,2) = 0.0_p_double*ldx(1)
      odx_e(2,2) = 0.5_p_double*ldx(2)
      odx_e(3,2) = 0.0_p_double*ldx(3)

      odx_e(1,3) = 0.0_p_double*ldx(1)
      odx_e(2,3) = 0.0_p_double*ldx(2)
      odx_e(3,3) = 0.5_p_double*ldx(3)

      select case( this%interpolation )
      case( p_linear, p_cubic )
        ! Nothing to change here
      case( p_quadratic, p_quartic )
        do j2 = 1, p_f_dim
          do j1 = 1, p_x_dim
            odx_b(j1,j2) = odx_b(j1,j2) + 0.5_p_double*ldx(j1)
            odx_e(j1,j2) = odx_e(j1,j2) + 0.5_p_double*ldx(j1)
          enddo
        enddo
      end select

      ! set the field values

      do i= 1, p_f_dim
        ! set the external b field

         select case (this%type_b(i))
          case(p_emf_uniform)
            ! set the external field to the supplied uniform value
            do j3 = lb(3), ub(3)
              do j2 = lb(2), ub(2)
                do j1 = lb(1), ub(1)
                   b%f3(i,j1,j2,j3) = this%uniform_b0(i)
                enddo
              enddo
            enddo

          case(p_emf_math)

            ! set the external field to the function value
            ! also set the value in the guard cells
            x_eval(4) = t
            do j3 = lb(3), ub(3)
              do j2 = lb(2), ub(2)
                do j1 = lb(1), ub(1)
                   ! set the proper correction to account for field positions
                   x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_b(1,i)
                   x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_b(2,i)
                   x_eval(3) = g_xmin(3) + (j3-1+shift(3))*ldx(3) + odx_b(3,i)

                   b%f3(i,j1,j2,j3) = real( eval( this%mfunc_b(i), x_eval ), p_k_fld )
                enddo
              enddo
            enddo

          case(p_emf_dipole)

            ! set the external field to a dipole field

            do j3 = lb(3), ub(3)
              do j2 = lb(2), ub(2)
                do j1 = lb(1), ub(1)

                   ! set the proper correction to account for field positions
                   x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_b(1, i)
                   x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_b(2, i)
                   x_eval(3) = g_xmin(3) + (j3-1+shift(3))*ldx(3) + odx_b(3, i)

                   b%f3(i,j1,j2,j3) = dipole_field( this%dipole_b_m, &
                                                    this%dipole_b_x0, &
                                                    this%dipole_b_r0, &
                                                    i, x_eval )
                enddo
              enddo
            enddo

          case(p_emf_python)

            x_bnd(1,1) = g_xmin(1) + (lb(1)-1+shift(1))*ldx(1) + odx_b(1,i)
            x_bnd(1,2) = g_xmin(2) + (lb(2)-1+shift(2))*ldx(2) + odx_b(2,i)
            x_bnd(1,3) = g_xmin(3) + (lb(3)-1+shift(3))*ldx(3) + odx_b(3,i)

            x_bnd(2,1) = g_xmin(1) + (ub(1)-1+shift(1))*ldx(1) + odx_b(1,i)
            x_bnd(2,2) = g_xmin(2) + (ub(2)-1+shift(2))*ldx(2) + odx_b(2,i)
            x_bnd(2,3) = g_xmin(3) + (ub(3)-1+shift(3))*ldx(3) + odx_b(3,i)

            fld = ""
            write(fld,'(A,I1)') "b", i
            call set_state( "data", b%f3(i,lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) )
            call set_state( "fld", trim(fld) )
            call set_state( "x_bnd", x_bnd(:,1:p_x_dim) )
            call set_state( "g_xmin", g_xmin(1:p_x_dim) )
            call set_state( "g_xmax", xmax(g_space) )
            call set_state( "nx_p_min", nx_p_min(1:p_x_dim) )
            call set_state( "t", (/t/) )

            call call_function( trim(this%py_mod), trim(this%py_func) )

            call get_state( "data", b%f3(i,lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) )

        end select

        ! set the external e field

        select case (this%type_e(i))

          case(p_emf_uniform)
            ! set the external field to the supplied uniform value
            do j3 = lb(3), ub(3)
              do j2 = lb(2), ub(2)
                do j1 = lb(1), ub(1)
                   e%f3(i,j1,j2,j3) = this%uniform_e0(i)
                enddo
              enddo
            enddo

          case(p_emf_math)

            ! set the external field to the function value
            ! also set the value in the guard cells
            x_eval(4) = t
            do j3 = lb(3), ub(3)
              do j2 = lb(2), ub(2)
                do j1 = lb(1), ub(1)

                   ! set the proper correction to account for field positions
                   x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_e(1, i)
                   x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_e(2, i)
                   x_eval(3) = g_xmin(3) + (j3-1+shift(3))*ldx(3) + odx_e(3, i)

                   e%f3(i,j1,j2,j3) = real( eval( this%mfunc_e(i), x_eval ), p_k_fld )

                enddo
              enddo
            enddo

          case(p_emf_dipole)

            ! set the external field to a dipole field
            do j3 = lb(3), ub(3)
              do j2 = lb(2), ub(2)
                do j1 = lb(1), ub(1)

                   ! set the proper correction to account for field positions
                   x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_e(1, i)
                   x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_e(2, i)
                   x_eval(3) = g_xmin(3) + (j3-1+shift(3))*ldx(3) + odx_e(3, i)

                   e%f3(i,j1,j2,j3) = 0.5_p_k_fld * dipole_field( this%dipole_e_p, &
                                                                  this%dipole_e_x0, &
                                                                  this%dipole_e_r0, &
                                                                  i, x_eval )

                enddo
              enddo
            enddo

          case(p_emf_python)

            x_bnd(1,1) = g_xmin(1) + (lb(1)-1+shift(1))*ldx(1) + odx_e(1,i)
            x_bnd(1,2) = g_xmin(2) + (lb(2)-1+shift(2))*ldx(2) + odx_e(2,i)
            x_bnd(1,3) = g_xmin(3) + (lb(3)-1+shift(3))*ldx(3) + odx_e(3,i)

            x_bnd(2,1) = g_xmin(1) + (ub(1)-1+shift(1))*ldx(1) + odx_e(1,i)
            x_bnd(2,2) = g_xmin(2) + (ub(2)-1+shift(2))*ldx(2) + odx_e(2,i)
            x_bnd(2,3) = g_xmin(3) + (ub(3)-1+shift(3))*ldx(3) + odx_e(3,i)

            fld = ""
            write(fld,'(A,I1)') "e", i
            call set_state( "data", e%f3(i,lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) )
            call set_state( "fld", trim(fld) )
            call set_state( "x_bnd", x_bnd(:,1:p_x_dim) )
            call set_state( "g_xmin", g_xmin(1:p_x_dim) )
            call set_state( "g_xmax", xmax(g_space) )
            call set_state( "nx_p_min", nx_p_min(1:p_x_dim) )
            call set_state( "t", (/t/) )

            call call_function( trim(this%py_mod), trim(this%py_func) )

            call get_state( "data", e%f3(i,lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) )

        end select

      enddo ! i= 1, p_f_dim
  end select

end subroutine set_fld_values_emf_gridval
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine write_checkpoint_emf_gridval( this, restart_handle )

  use m_emf_define, only: t_emf_gridval
  use m_parameters
  use m_restart
  use m_emf_gridval

  implicit none

  ! dummy variables

  class( t_emf_gridval ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  ! local variables
  character(len=*), parameter :: err_msg = 'error writing restart data for emf gridval object.'
  integer :: ierr

  ! ---

  call restart_io_write("p_emf_gridval_rst_id", p_emf_gridval_rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",699)

  call restart_io_write("this%type_b", this%type_b, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",702)

  call restart_io_write("this%type_e", this%type_e, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",705)

  call restart_io_write("this%uniform_b0", this%uniform_b0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",708)

  call restart_io_write("this%uniform_e0", this%uniform_e0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",711)

  call restart_io_write("this%mfunc_expr_b", this%mfunc_expr_b, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",714)

  call restart_io_write("this%mfunc_expr_e", this%mfunc_expr_e, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",717)

  call restart_io_write("this%dipole_b_m", this%dipole_b_m, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",720)

  call restart_io_write("this%dipole_b_x0", this%dipole_b_x0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",723)

  call restart_io_write("this%dipole_b_r0", this%dipole_b_r0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",726)

  call restart_io_write("this%dipole_e_p", this%dipole_e_p, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",729)

  call restart_io_write("this%dipole_e_x0", this%dipole_e_x0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",732)

  call restart_io_write("this%dipole_e_r0", this%dipole_e_r0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",735)


end subroutine write_checkpoint_emf_gridval
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine restart_read_emf_gridval( this, restart_handle )

  use m_emf_define, only: t_emf_gridval
  use m_parameters
  use m_restart
  use m_emf_gridval

  implicit none

  ! dummy variables

  class( t_emf_gridval ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables
  character(len=*), parameter :: err_msg = 'error reading restart data for emf gridval object.'
  character(len=len(p_emf_gridval_rst_id)) :: rst_id
  integer :: ierr

  ! ---

  call restart_io_read(rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"emf/os-emf-gridval.f03",764)

  ! check if restart file is compatible
  if ( rst_id /= p_emf_gridval_rst_id) then
    write(err_buf__,*) 'Corrupted restart file, or restart file';call err__("emf/os-emf-gridval.f03",768)
    write(err_buf__,*) 'from incompatible binary (emf_gridval)';call err__("emf/os-emf-gridval.f03",769)
    write(err_buf__,*) 'rst_id = "', rst_id,'"';call err__("emf/os-emf-gridval.f03",770)
    call abort_program(p_err_rstrd)
  endif

  call restart_io_read(this%type_b, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",775)

  call restart_io_read(this%type_e, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",778)

  call restart_io_read(this%uniform_b0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",781)

  call restart_io_read(this%uniform_e0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",784)

  call restart_io_read(this%mfunc_expr_b, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",787)

  call restart_io_read(this%mfunc_expr_e, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",790)

  call restart_io_read(this%dipole_b_m, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",793)

  call restart_io_read(this%dipole_b_x0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",796)

  call restart_io_read(this%dipole_b_r0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",799)

  call restart_io_read(this%dipole_e_p, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",802)

  call restart_io_read(this%dipole_e_x0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",805)

  call restart_io_read(this%dipole_e_r0, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"emf/os-emf-gridval.f03",808)


end subroutine restart_read_emf_gridval
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine setup_emf_gridval( this, dynamic, interpolation, restart, restart_handle )

  use m_emf_define, only: t_emf_gridval, p_emf_math
  use m_fparser, only: setup
  use m_parameters
  use m_restart

  implicit none

  class( t_emf_gridval ), intent(inout) :: this
  logical, intent(in) :: dynamic
  integer, intent(in) :: interpolation
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  integer :: i, nvars, ierr
  character(len=2), dimension(p_x_dim+1) :: vars

  ! if restarting read information from checkpoint file
  if ( restart ) then
    call this%restart_read( restart_handle )
  endif

  this%dynamic = dynamic
  this%interpolation = interpolation

  ! compile math functions if required
  do i=1, p_x_dim
    vars(i)='x'//char(ichar('0')+i)
  enddo

  if ( this%dynamic ) then
    nvars = p_x_dim + 1
    vars(nvars) = 't'
  else
    nvars = p_x_dim
  endif

  do i= 1, p_f_dim
    if (this%type_b(i) == p_emf_math) then
      call setup(this%mfunc_b(i), trim(this%mfunc_expr_b(i)), vars(1:nvars), ierr)
      if ( ierr /= 0) then
         print *, "(*error*) Error compiling supplied function :"
         print *, trim(this%mfunc_expr_b(i))
         print *, "(*error*) for external/initial B field values."
         print *, "(*error*) bailing out..."
         call abort_program()
      endif
    endif

    if (this%type_e(i) == p_emf_math) then
        call setup(this%mfunc_e(i), trim(this%mfunc_expr_e(i)), vars(1:nvars), ierr)
      if ( ierr /= 0) then
         print *, "(*error*) Error compiling supplied function :"
         print *, trim(this%mfunc_expr_e(i))
         print *, "(*error*) for external/initial E field values."
         print *, "(*error*) bailing out..."
         call abort_program()
      endif
    endif
  enddo

end subroutine setup_emf_gridval
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine cleanup_emf_gridval( this )

  use m_emf_define, only: t_emf_gridval, p_emf_math
  use m_fparser, only: cleanup
  use m_parameters

  implicit none

  class( t_emf_gridval ), intent(inout) :: this
  integer :: i

  do i= 1, p_f_dim
    if (this%type_b(i) == p_emf_math) then
        call cleanup(this%mfunc_b(i))
    endif

    if (this%type_e(i) == p_emf_math) then
        call cleanup(this%mfunc_e(i))
    endif
  enddo


end subroutine cleanup_emf_gridval
!---------------------------------------------------------------------------------------------------
