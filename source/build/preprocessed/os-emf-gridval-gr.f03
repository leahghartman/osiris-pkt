# 1 "gr/os-emf-gridval-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-emf-gridval-gr.f03" 2
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
# 2 "gr/os-emf-gridval-gr.f03" 2
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
# 3 "gr/os-emf-gridval-gr.f03" 2

module m_emf_gridval_gr

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
# 7 "gr/os-emf-gridval-gr.f03" 2

use m_system
use m_parameters

use m_emf_define, only : t_emf_gridval
use m_geometry_gr, only : t_geometry_gr

implicit none

private

! string to id restart data
character(len=*), parameter, public :: p_emf_gridval_gr_rst_id = "emf_gridval_gr rst data - 0x0006"

! IBM XL compilers
! This has to be explicitly made public so that you can access the superclass
! of t_emf_gridval_gr
public :: t_emf_gridval_gr

!-------------------------------------------------------------------------------
! t_emf_gridval_gr class definition
!-------------------------------------------------------------------------------

type, extends( t_emf_gridval ) :: t_emf_gridval_gr

  class( t_geometry_gr ), pointer :: geometry

contains

  procedure :: set_fld_values => set_fld_values_gr
  procedure :: write_checkpoint => write_checkpoint_gr
  procedure :: restart_read => restart_read_gr

end type t_emf_gridval_gr

contains

subroutine set_fld_values_gr( this, e, b, g_space, nx_p_min, t, lbin, ubin )

  use m_emf_define
  use m_vdf_define, only : t_vdf
  use m_space, only : t_space, xmin
  use m_fparser

  implicit none

  class(t_emf_gridval_gr), intent( inout ) :: this
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
  integer :: i, j1, j2
  integer :: x_dim

  integer, dimension(p_max_dim) :: lb, ub
  integer, dimension( p_max_dim ) :: shift

  ! executable statements
  g_xmin(1:p_x_dim) = xmin(g_space)
  ldx(1:p_x_dim) = b%dx()

  shift(1:p_x_dim) = nx_p_min(1:p_x_dim) - 1

  select case (p_x_dim)
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
    odx_b(1,1) = 0.0_p_double
    odx_b(2,1) = 0.5_p_double

    odx_b(1,2) = 0.5_p_double
    odx_b(2,2) = 0.0_p_double

    odx_b(1,3) = 0.5_p_double
    odx_b(2,3) = 0.5_p_double

    odx_e(1,1) = 0.5_p_double
    odx_e(2,1) = 0.0_p_double

    odx_e(1,2) = 0.0_p_double
    odx_e(2,2) = 0.5_p_double

    odx_e(1,3) = 0.0_p_double
    odx_e(2,3) = 0.0_p_double


    do i= 1, p_f_dim
      ! set the b field
      select case (this%type_b(i))
      case(p_emf_uniform)

        ! set the field to the supplied uniform value
        do j2 = lb(2), ub(2)
          do j1 = lb(1), ub(1)
              b%f2(i,j1,j2) = this%uniform_b0(i)
          enddo
        enddo

      case(p_emf_math)

        ! set the field to the function value
        ! also set the value in the guard cells
        x_eval(3) = t
        do j2 = lb(2), ub(2)
          do j1 = lb(1), ub(1)
            ! set the proper correction to account for field positions
            x_eval(1) = this%geometry%r%f1(1, j1) + odx_b(1, i) * this%geometry%dr%f1(1, j1)
            x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_b(2, i) * ldx(2)

            b%f2(i,j1,j2) = real( eval( this%mfunc_b(i), x_eval ), p_k_fld )

          enddo
        enddo
      end select

      ! set the e field
      select case (this%type_e(i))
      case(p_emf_uniform)
        ! set the field to the supplied uniform value
        do j2 = lb(2), ub(2)
          do j1 = lb(1), ub(1)
              e%f2(i,j1,j2) = this%uniform_e0(i)
          enddo
        enddo

      case(p_emf_math)

        ! set the field to the function value
        ! also set the value in the guard cells
        x_eval(3) = t
        do j2 = lb(2), ub(2)
          do j1 = lb(1), ub(1)
            ! set the proper correction to account for field positions
            x_eval(1) = this%geometry%r%f1(1, j1) + odx_e(1, i) * this%geometry%dr%f1(1, j1)
            x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_e(2, i) * ldx(2)

            e%f2(i,j1,j2) = real( eval( this%mfunc_e(i), x_eval ), p_k_fld )
          enddo
        enddo
      end select
    enddo
  end select

end subroutine

!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_gr( this, restart_handle )
!-----------------------------------------------------------------------------------------
! write object information into a restart file
!-----------------------------------------------------------------------------------------

  use m_restart
  use m_parameters

  implicit none

  class( t_emf_gridval_gr ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for emf_gridval_gr object.'
  integer :: ierr

  ! change the rst_id for the GR value
  call restart_io_write("p_emf_gridval_gr_rst_id", p_emf_gridval_gr_rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"gr/os-emf-gridval-gr.f03",198)

  ! write superclass checkpoint data first
  call this % t_emf_gridval % write_checkpoint( restart_handle )

end subroutine write_checkpoint_gr
!-----------------------------------------------------------------------------------------

subroutine restart_read_gr( this, restart_handle )

  use m_restart
  use m_parameters

  implicit none

  class( t_emf_gridval_gr ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for emf_gridval_gr object.'
  character(len=len(p_emf_gridval_gr_rst_id)) :: rst_id
  integer :: ierr

  call restart_io_read(rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"gr/os-emf-gridval-gr.f03",221)

  ! read superclass checkpoint data first
  call this % t_emf_gridval % restart_read( restart_handle )

end subroutine restart_read_gr

end module m_emf_gridval_gr
