# 1 "cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03" 2
!-----------------------------------------------------------------------------------------
! Customized (Fei) Solver
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
# 6 "cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03" 2
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
# 7 "cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03" 2

module m_emf_solver_cyl_modes_fei

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
# 11 "cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03" 2

use m_parameters
use m_emf_define
use m_emf_solver_fei
use m_vdf_define, only : t_vdf
use m_fft
use m_math, only : pi
use m_cyl_modes, only : t_cyl_modes, check_nan
use m_emf_bound_cyl_modes, only : t_emf_bound_cyl_modes

private

type, extends( t_emf_solver_fei ), public :: t_emf_solver_cyl_modes_fei

    integer :: n_cyl_modes = -1
    procedure(dedt_cyl_modes_fei), pointer :: dedt_cyl_modes => null()
    procedure(dbdt_cyl_modes_fei), pointer :: dbdt_cyl_modes => null()

    ! Current filtering
    procedure(filter_corr_current_cyl_modes_fei), pointer :: filter_corr_current_cyl_modes => null()
    real(p_k_fld), dimension(:,:,:), pointer :: f_buffer_re => null(), &
                                                f_buffer_im => null()

contains

    procedure :: init => init_cyl_modes_fei
    procedure :: advance_cyl_modes => advance_cyl_modes_fei
    procedure :: test_stability => test_stability_cyl_modes_fei
    procedure :: cleanup => cleanup_cyl_modes_fei

end type t_emf_solver_cyl_modes_fei

interface
    subroutine dedt_cyl_modes_fei( this, e_cyl_m, b_cyl_m, jay_cyl_m, dt, bnd_type )
        import t_emf_solver_cyl_modes_fei, t_cyl_modes, p_double, p_max_dim
        class( t_emf_solver_cyl_modes_fei ), intent(in) :: this
        type( t_cyl_modes ), intent(inout) :: e_cyl_m
        type( t_cyl_modes ), intent(in) :: b_cyl_m, jay_cyl_m
        real(p_double), intent(in) :: dt
        integer, intent(in), dimension(2,p_max_dim) :: bnd_type
    end subroutine
end interface

interface
    subroutine dbdt_cyl_modes_fei( this, b_cyl_m, e_cyl_m, dt, bnd_type )
        import t_emf_solver_cyl_modes_fei, t_cyl_modes, p_double, p_max_dim
        class( t_emf_solver_cyl_modes_fei ), intent(in) :: this
        type( t_cyl_modes ), intent(inout) :: b_cyl_m
        type( t_cyl_modes ), intent(in) :: e_cyl_m
        real(p_double), intent(in) :: dt
        integer, intent(in), dimension(2,p_max_dim) :: bnd_type
    end subroutine
end interface

interface
    subroutine filter_corr_current_cyl_modes_fei( this, jay_cyl_m )
        import t_emf_solver_cyl_modes_fei, t_cyl_modes
        class( t_emf_solver_cyl_modes_fei ), intent(inout) :: this
        type( t_cyl_modes ), intent(inout) :: jay_cyl_m
    end subroutine
end interface

interface
    function kernel( x, arglst )
        import p_double
        real(p_double), intent(in) :: x
        real(p_double), intent(in), dimension(:) :: arglst
        real(p_double) :: kernel
    end function
end interface

contains

!-----------------------------------------------------------------------------------------
! B-field advance
!-----------------------------------------------------------------------------------------
subroutine dbdt_2d_cyl_modes_fei( this, b_cyl_m, e_cyl_m, dt, bnd_type )

  implicit none

  class ( t_emf_solver_cyl_modes_fei ), intent(in) :: this
  type( t_cyl_modes ), intent(inout) :: b_cyl_m
  type( t_cyl_modes ), intent(in) :: e_cyl_m
  real(p_double), intent(in) :: dt
  integer, intent(in), dimension(2,p_max_dim) :: bnd_type


  real(p_k_fld) :: dtdz, dtdr, dtif, dt_r, dt_r_2
  real(p_k_fld) :: tmp, rm, rp, tmp_b, n_e, o_e, phase_factor
  integer :: i1, i2, i2_0, gshift_i2, mode, mm
  type( t_vdf ), pointer :: b_re, b_im, e_re, e_im
  integer :: j, i1_lb, i1_ub, offset, ord_hf, i1_lb_ax, i1_ub_ax


  integer :: n_modes

  n_modes = ubound(e_cyl_m%pf_re,1)

  dtdz = real(dt/b_cyl_m%pf_re(0)%dx_(1), p_k_fld)
  dtdr = real(dt/b_cyl_m%pf_re(0)%dx_(2), p_k_fld)
  dtif = real( dt, p_k_fld )
  n_e = real(9.0_p_k_fld/8.0_p_k_fld, p_k_fld)
  o_e = real(1.0_p_k_fld/8.0_p_k_fld, p_k_fld)

  i1_lb = 2 - 2*this%n_coef
  i1_ub = b_cyl_m%pf_re(0)%nx_(1) - 1 + 2*this%n_coef

  i1_lb_ax = i1_lb
  i1_ub_ax = i1_ub

  if ( this%taper_bnd ) then

    offset = this%taper_order / 2

    if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
         bnd_type(p_lower, 1) /= p_bc_other_node ) then
      i1_lb = offset
    endif

    if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
        bnd_type(p_upper, 1) /= p_bc_other_node ) then
      i1_ub = b_cyl_m%pf_re(0)%nx_(1) - offset + 1
    endif

  endif

  gshift_i2 = this%gir_pos - 2

  if ( gshift_i2 < 0 ) then
    ! This node contains the cylindrical axis so start the solver in cell 2
    i2_0 = 2
  else
    i2_0 = 0
  endif


  do mode = 0, n_modes ! look over each mode independently 0th mode is included
    b_re => b_cyl_m%pf_re(mode)
    e_re => e_cyl_m%pf_re(mode)
    if (mode > 0) then
      b_im => b_cyl_m%pf_im(mode)
      e_im => e_cyl_m%pf_im(mode)
    endif

    !$omp parallel do private(rp,rm,tmp,dt_r,dt_r_2,i1,tmp_b,j)
    do i2 = i2_0, b_re%nx_(2)+1
      rp = i2 + gshift_i2 + 0.5 ! position of the upper edge of the cell (normalized to dr) ! at the axis this is 1.5
      rm = i2 + gshift_i2 - 0.5 ! position of the lower edge of the cell (normalized to dr) ! at the axis this is 0.5
      tmp = dtdr/(i2 + gshift_i2) ! (dt/dr) / position of the middle of the cell (normalized to dr)
      dt_r = tmp
      dt_r_2 = dtdr/(i2 + gshift_i2 - 0.5)
      !dt_r = dtif/(i2 + gshift_i2) ! (dt) / position of the middle of the cell (normalized to dr) ! at the axis this is 1
      !dt_r_2 = dtif/(i2 + gshift_i2 - 0.5) ! (dt) / position of the middle of the cell (normalized to dr)

      do i1 = i1_lb, i1_ub
        ! First the Real part
        !B1
        b_re%f2( 1, i1, i2 ) = b_re%f2(1, i1, i2) &
                        - tmp * ( rp * e_re%f2( 3, i1, i2+1 ) - rm * e_re%f2( 3, i1, i2 ))
        !B2
        tmp_b = 0.0_p_k_fld
        do j = 1, this%n_coef
          tmp_b = tmp_b + this%coef_e(j) * ( e_re%f2(3, i1+j, i2) - e_re%f2(3, i1-j+1, i2) )
        enddo
        b_re%f2( 2, i1, i2 ) = b_re%f2(2, i1, i2) + dtdz * tmp_b
        !B3
        tmp_b = 0.0_p_k_fld
        do j = 1, this%n_coef
          tmp_b = tmp_b + this%coef_e(j) * ( e_re%f2( 2, i1+j, i2 ) - e_re%f2( 2, i1-j+1, i2 ) )
        enddo
        b_re%f2( 3, i1, i2 ) = b_re%f2( 3, i1, i2 ) - dtdz * tmp_b &
                                + dtdr * ( e_re%f2( 1, i1, i2+1 ) - e_re%f2( 1, i1, i2 ))

        if (mode == 0) cycle

        b_re%f2( 1, i1, i2 ) = b_re%f2( 1, i1, i2 ) + dt_r * mode * e_im%f2(2, i1, i2)
        b_re%f2( 2, i1, i2 ) = b_re%f2( 2, i1, i2 ) - dt_r_2 * mode * e_im%f2(1, i1, i2)

        ! Now the Imaginary part
        !B1
        b_im%f2( 1, i1, i2 ) = b_im%f2(1, i1, i2) &
                        - tmp*( rp * e_im%f2( 3, i1, i2+1 ) - rm * e_im%f2( 3, i1, i2 )) &
                        - dt_r * mode * e_re%f2(2, i1, i2)

        !B2
        tmp_b = 0.0_p_k_fld
        do j = 1, this%n_coef
          tmp_b = tmp_b + this%coef_e(j) * ( e_im%f2(3, i1+j, i2) - e_im%f2(3, i1-j+1, i2) )
        enddo
        b_im%f2( 2, i1, i2 ) = b_im%f2(2, i1, i2) + dtdz * tmp_b &
                                + dt_r_2 * mode * e_re%f2(1, i1, i2)

        !B3
        tmp_b = 0.0_p_k_fld
        do j = 1, this%n_coef
          tmp_b = tmp_b + this%coef_e(j) * ( e_im%f2(2, i1+j, i2) - e_im%f2(2, i1-j+1, i2) )
        enddo
        b_im%f2( 3, i1, i2 ) = b_im%f2( 3, i1, i2 ) - dtdz * tmp_b + &
                                dtdr * ( e_im%f2( 1, i1, i2+1 ) - e_im%f2( 1, i1, i2 ))

      enddo ! i1

    enddo ! i2
    !$omp end parallel do

    ! do the tapered region
    if ( this%taper_bnd ) then

      ! advance lower boundary
      if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
           bnd_type(p_lower, 1) /= p_bc_other_node ) then

        !$omp parallel do private(rp,rm,tmp,dt_r,dt_r_2,ord_hf,i1,tmp_b,j)
        do i2 = i2_0, b_re%nx_(2)+1
          rp = i2 + gshift_i2 + 0.5 ! position of the upper edge of the cell (normalized to dr) ! at the axis this is 1.5
          rm = i2 + gshift_i2 - 0.5 ! position of the lower edge of the cell (normalized to dr) ! at the axis this is 0.5
          tmp = dtdr/(i2 + gshift_i2) ! (dt/dr) / position of the middle of the cell (normalized to dr)
          dt_r = tmp
          dt_r_2 = dtdr/(i2 + gshift_i2 - 0.5)
          !dt_r = dtif/(i2 + gshift_i2) ! (dt) / position of the middle of the cell (normalized to dr) ! at the axis this is 1
          !dt_r_2 = dtif/(i2 + gshift_i2 - 0.5) ! (dt) / position of the middle of the cell (normalized to dr)

          ord_hf = 0

          do i1 = 0, i1_lb-1
            ! First the Real part

            ord_hf = ord_hf + 1
            !B1
            b_re%f2( 1, i1, i2 ) = b_re%f2(1, i1, i2) &
                        - tmp * ( rp * e_re%f2( 3, i1, i2+1 ) - rm * e_re%f2( 3, i1, i2 ))
            !B2
            tmp_b = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_re%f2(3, i1+j, i2) - e_re%f2(3, i1-j+1, i2) )
            enddo
            b_re%f2( 2, i1, i2 ) = b_re%f2(2, i1, i2) + dtdz * tmp_b
            !B3
            tmp_b = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_re%f2( 2, i1+j, i2 ) - e_re%f2( 2, i1-j+1, i2 ) )
            enddo
            b_re%f2( 3, i1, i2 ) = b_re%f2( 3, i1, i2 ) - dtdz * tmp_b &
                                + dtdr * ( e_re%f2( 1, i1, i2+1 ) - e_re%f2( 1, i1, i2 ))

            if (mode == 0) cycle

            b_re%f2( 1, i1, i2 ) = b_re%f2( 1, i1, i2 ) + dt_r * mode * e_im%f2(2, i1, i2)
            b_re%f2( 2, i1, i2 ) = b_re%f2( 2, i1, i2 ) - dt_r_2 * mode * e_im%f2(1, i1, i2)

            ! Now the Imaginary part
            !B1
            b_im%f2( 1, i1, i2 ) = b_im%f2(1, i1, i2) &
                        - tmp*( rp * e_im%f2( 3, i1, i2+1 ) - rm * e_im%f2( 3, i1, i2 )) &
                        - dt_r * mode * e_re%f2(2, i1, i2)

            !B2
            tmp_b = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_im%f2(3, i1+j, i2) - e_im%f2(3, i1-j+1, i2) )
            enddo
            b_im%f2( 2, i1, i2 ) = b_im%f2(2, i1, i2) + dtdz * tmp_b &
                                    + dt_r_2 * mode * e_re%f2(1, i1, i2)

            !B3
            tmp_b = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_im%f2(2, i1+j, i2) - e_im%f2(2, i1-j+1, i2) )
            enddo
            b_im%f2( 3, i1, i2 ) = b_im%f2( 3, i1, i2 ) - dtdz * tmp_b + &
                                   dtdr * ( e_im%f2( 1, i1, i2+1 ) - e_im%f2( 1, i1, i2 ))

          enddo ! i1

        enddo ! i2
        !$omp end parallel do

      endif

      ! advance upper boundary
      if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
           bnd_type(p_upper, 1) /= p_bc_other_node ) then

        !$omp parallel do private(rp,rm,tmp,dt_r,dt_r_2,ord_hf,i1,tmp_b,j)
        do i2 = i2_0, b_re%nx_(2)+1
          rp = i2 + gshift_i2 + 0.5 ! position of the upper edge of the cell (normalized to dr) ! at the axis this is 1.5
          rm = i2 + gshift_i2 - 0.5 ! position of the lower edge of the cell (normalized to dr) ! at the axis this is 0.5
          tmp = dtdr/(i2 + gshift_i2) ! (dt/dr) / position of the middle of the cell (normalized to dr)
          dt_r = tmp
          dt_r_2 = dtdr/(i2 + gshift_i2 - 0.5)
          !dt_r = dtif/(i2 + gshift_i2) ! (dt) / position of the middle of the cell (normalized to dr) ! at the axis this is 1
          !dt_r_2 = dtif/(i2 + gshift_i2 - 0.5) ! (dt) / position of the middle of the cell (normalized to dr)

          ord_hf = offset + 1

          do i1 = i1_ub+1, b_re%nx_(1)+1
            ! First the Real part

            ord_hf = ord_hf - 1
            !B1
            b_re%f2( 1, i1, i2 ) = b_re%f2(1, i1, i2) &
                        - tmp * ( rp * e_re%f2( 3, i1, i2+1 ) - rm * e_re%f2( 3, i1, i2 ))
            !B2
            tmp_b = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_re%f2(3, i1+j, i2) - e_re%f2(3, i1-j+1, i2) )
            enddo
            b_re%f2( 2, i1, i2 ) = b_re%f2(2, i1, i2) + dtdz * tmp_b
            !B3
            tmp_b = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_re%f2( 2, i1+j, i2 ) - e_re%f2( 2, i1-j+1, i2 ) )
            enddo
            b_re%f2( 3, i1, i2 ) = b_re%f2( 3, i1, i2 ) - dtdz * tmp_b &
                                + dtdr * ( e_re%f2( 1, i1, i2+1 ) - e_re%f2( 1, i1, i2 ))

            if (mode == 0) cycle

            b_re%f2( 1, i1, i2 ) = b_re%f2( 1, i1, i2 ) + dt_r * mode * e_im%f2(2, i1, i2)
            b_re%f2( 2, i1, i2 ) = b_re%f2( 2, i1, i2 ) - dt_r_2 * mode * e_im%f2(1, i1, i2)

            ! Now the Imaginary part
            !B1
            b_im%f2( 1, i1, i2 ) = b_im%f2(1, i1, i2) &
                        - tmp*( rp * e_im%f2( 3, i1, i2+1 ) - rm * e_im%f2( 3, i1, i2 )) &
                        - dt_r * mode * e_re%f2(2, i1, i2)

            !B2
            tmp_b = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_im%f2(3, i1+j, i2) - e_im%f2(3, i1-j+1, i2) )
            enddo
            b_im%f2( 2, i1, i2 ) = b_im%f2(2, i1, i2) + dtdz * tmp_b &
                                    + dt_r_2 * mode * e_re%f2(1, i1, i2)

            !B3
            tmp_b = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_im%f2(2, i1+j, i2) - e_im%f2(2, i1-j+1, i2) )
            enddo
            b_im%f2( 3, i1, i2 ) = b_im%f2( 3, i1, i2 ) - dtdz * tmp_b + &
                                   dtdr * ( e_im%f2( 1, i1, i2+1 ) - e_im%f2( 1, i1, i2 ))

          enddo ! i1

        enddo ! i2
        !$omp end parallel do

      endif

    endif ! end of doing tapered region

    ! do boundary
    ! note, the guard cells that reside below the axis do not hold true field values,
    ! but hold the field values appropriate for particle interpolation.
    ! i.e., we must have the following:
    ! B_z(-dr) = (-1)^m B_z(dr)
    ! B_r(-dr/2) = -(-1)^m B_r(dr/2)
    ! B_phi(-dr) = -(-1)^m B_phi(dr)
    ! Note, the B values below the axis aren't actually used in the field solver, so our
    ! choice here (unlike for the E values) is unimportant.
    ! if quartic interpolation is ever implemented, also handle cell -1 in r
    if ( gshift_i2 < 0 ) then

      if (mode == 0) then
        !$omp parallel do
        do i1 = i1_lb_ax, i1_ub_ax
          ! B1 can be nonzero on the axis
          b_re%f2( 1, i1, 1 ) = b_re%f2( 1, i1, 1 ) - 4 * dtdr * e_re%f2( 3, i1, 2 )
          b_re%f2( 1, i1, 0 ) = b_re%f2( 1, i1, 2 )

          ! B2 must be zero on the axis
          b_re%f2( 2, i1, 1 ) = - b_re%f2( 2, i1, 2 )
          b_re%f2( 2, i1, 0 ) = - b_re%f2( 2, i1, 3 )

          ! B3 must be zero on the axis
          b_re%f2( 3, i1, 1 ) = 0.0_p_k_fld
          b_re%f2( 3, i1, 0 ) = - b_re%f2( 3, i1, 2 )
        enddo
        !$omp end parallel do

      elseif (mode == 1) then
        ! First handle the regular regions (everywhere but tapered boundary)
        !$omp parallel do private(tmp_b,j)
        do i1 = i1_lb, i1_ub
          ! B1 must be zero on the axis
          b_re%f2( 1, i1, 1 ) = 0.0_p_k_fld
          b_re%f2( 1, i1, 0 ) = - b_re%f2( 1, i1, 2 )
          b_im%f2( 1, i1, 1 ) = 0.0_p_k_fld
          b_im%f2( 1, i1, 0 ) = - b_im%f2( 1, i1, 2 )

          ! B2 can be nonzero on the axis
          b_re%f2( 2, i1, 1 ) = b_re%f2( 2, i1, 2 )
          b_re%f2( 2, i1, 0 ) = b_re%f2( 2, i1, 3 )
          b_im%f2( 2, i1, 1 ) = b_im%f2( 2, i1, 2 )
          b_im%f2( 2, i1, 0 ) = b_im%f2( 2, i1, 3 )

          ! B3 can be nonzero on the axis
          ! Old method
          ! Instead of dtdr*(e_re%f2( 1, i1, 2 ) - e_re%f2( 1, i1, 1 )) or
          ! dtdr*(e_im%f2( 1, i1, 2 ) - e_im%f2( 1, i1, 1 ))
          ! we now use dtdr*2*e_re%f2( 1, i1, 2 ). They are equivalent, but it is
          ! more clear this way because e_re%f2( 1, i1, 1 ) and e_im%f2( 1, i1, 1 ) have
          ! a (-1)^m term multiplying them that can be confusing.
          ! tmp_b = 0.0_p_k_fld
          ! do j = 1, this%n_coef
          ! tmp_b = tmp_b + this%coef_e(j) * ( e_re%f2(2, i1+j, 1) - e_re%f2(2, i1-j+1, 1) )
          ! enddo
          ! b_re%f2( 3, i1, 1 ) = b_re%f2( 3, i1, 1 ) - dtdz * tmp_b + &
          ! dtdr * 2 * e_re%f2( 1, i1, 2 )
          ! tmp_b = 0.0_p_k_fld
          ! do j = 1, this%n_coef
          ! tmp_b = tmp_b + this%coef_e(j) * ( e_im%f2(2, i1+j, 1) - e_im%f2(2, i1-j+1, 1) )
          ! enddo
          ! b_im%f2( 3, i1, 1 ) = b_im%f2( 3, i1, 1 ) - dtdz * tmp_b + &
          ! dtdr * 2 * e_im%f2( 1, i1, 2 )

          ! New method
          ! Rather than solving for b_phi(r=0) in the normal way as listed in the paper,
          ! we instead solve ensuring that db_r/dr(r=0) = 0 and that the relation
          ! b_r = i b_phi is satisfied. See documentation for Smilei about this. For our
          ! grid this translates to b_phi_re = b_r_im and b_phi_im = -b_r_re at r = 0.
          b_re%f2( 3, i1, 1 ) = n_e * b_im%f2( 2, i1, 2 ) - o_e * b_im%f2( 2, i1, 3 )
          b_re%f2( 3, i1, 0 ) = b_re%f2( 3, i1, 2 )
          b_im%f2( 3, i1, 1 ) =-n_e * b_re%f2( 2, i1, 2 ) + o_e * b_re%f2( 2, i1, 3 )
          b_im%f2( 3, i1, 0 ) = b_im%f2( 3, i1, 2 )
        enddo
        !$omp end parallel do

        ! tapered region
        if ( this%taper_bnd ) then

          ! lower boundary
          if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
               bnd_type(p_lower, 1) /= p_bc_other_node ) then

            ord_hf = 0
            do i1 = 0, i1_lb-1
              ord_hf = ord_hf + 1

              ! B1 must be zero on the axis
              b_re%f2( 1, i1, 1 ) = 0.0_p_k_fld
              b_re%f2( 1, i1, 0 ) = - b_re%f2( 1, i1, 2 )
              b_im%f2( 1, i1, 1 ) = 0.0_p_k_fld
              b_im%f2( 1, i1, 0 ) = - b_im%f2( 1, i1, 2 )

              ! B2 can be nonzero on the axis
              b_re%f2( 2, i1, 1 ) = b_re%f2( 2, i1, 2 )
              b_re%f2( 2, i1, 0 ) = b_re%f2( 2, i1, 3 )
              b_im%f2( 2, i1, 1 ) = b_im%f2( 2, i1, 2 )
              b_im%f2( 2, i1, 0 ) = b_im%f2( 2, i1, 3 )

              ! B3 can be nonzero on the axis
              ! Old method
              ! tmp_b = 0.0_p_k_fld
              ! do j = 1, ord_hf
              ! tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_re%f2(2, i1+j, 1) - e_re%f2(2, i1-j+1, 1) )
              ! enddo
              ! b_re%f2( 3, i1, 1 ) = b_re%f2( 3, i1, 1 ) - dtdz * tmp_b + &
              ! dtdr * 2 * e_re%f2( 1, i1, 2 )
              ! tmp_b = 0.0_p_k_fld
              ! do j = 1, ord_hf
              ! tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_im%f2(2, i1+j, 1) - e_im%f2(2, i1-j+1, 1) )
              ! enddo
              ! b_im%f2( 3, i1, 1 ) = b_im%f2( 3, i1, 1 ) - dtdz * tmp_b + &
              ! dtdr * 2 * e_im%f2( 1, i1, 2 )

              ! New method
              b_re%f2( 3, i1, 1 ) = n_e * b_im%f2( 2, i1, 2 ) - o_e * b_im%f2( 2, i1, 3 )
              b_re%f2( 3, i1, 0 ) = b_re%f2( 3, i1, 2 )
              b_im%f2( 3, i1, 1 ) =-n_e * b_re%f2( 2, i1, 2 ) + o_e * b_re%f2( 2, i1, 3 )
              b_im%f2( 3, i1, 0 ) = b_im%f2( 3, i1, 2 )

            enddo

          endif

          ! upper boundary
          if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
               bnd_type(p_upper, 1) /= p_bc_other_node ) then

            ord_hf = offset + 1
            do i1 = i1_ub+1, b_re%nx_(1)+1
              ord_hf = ord_hf - 1

              ! B1 must be zero on the axis
              b_re%f2( 1, i1, 1 ) = 0.0_p_k_fld
              b_re%f2( 1, i1, 0 ) = - b_re%f2( 1, i1, 2 )
              b_im%f2( 1, i1, 1 ) = 0.0_p_k_fld
              b_im%f2( 1, i1, 0 ) = - b_im%f2( 1, i1, 2 )

              ! B2 can be nonzero on the axis
              b_re%f2( 2, i1, 1 ) = b_re%f2( 2, i1, 2 )
              b_re%f2( 2, i1, 0 ) = b_re%f2( 2, i1, 3 )
              b_im%f2( 2, i1, 1 ) = b_im%f2( 2, i1, 2 )
              b_im%f2( 2, i1, 0 ) = b_im%f2( 2, i1, 3 )

              ! B3 can be nonzero on the axis
              ! Old method
              ! tmp_b = 0.0_p_k_fld
              ! do j = 1, ord_hf
              ! tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_re%f2(2, i1+j, 1) - e_re%f2(2, i1-j+1, 1) )
              ! enddo
              ! b_re%f2( 3, i1, 1 ) = b_re%f2( 3, i1, 1 ) - dtdz * tmp_b + &
              ! dtdr * 2 * e_re%f2( 1, i1, 2 )
              ! tmp_b = 0.0_p_k_fld
              ! do j = 1, ord_hf
              ! tmp_b = tmp_b + this%taper_coef_e(j,ord_hf) * ( e_im%f2(2, i1+j, 1) - e_im%f2(2, i1-j+1, 1) )
              ! enddo
              ! b_im%f2( 3, i1, 1 ) = b_im%f2( 3, i1, 1 ) - dtdz * tmp_b + &
              ! dtdr * 2 * e_im%f2( 1, i1, 2 )

              ! New method
              b_re%f2( 3, i1, 1 ) = n_e * b_im%f2( 2, i1, 2 ) - o_e * b_im%f2( 2, i1, 3 )
              b_re%f2( 3, i1, 0 ) = b_re%f2( 3, i1, 2 )
              b_im%f2( 3, i1, 1 ) =-n_e * b_re%f2( 2, i1, 2 ) + o_e * b_re%f2( 2, i1, 3 )
              b_im%f2( 3, i1, 0 ) = b_im%f2( 3, i1, 2 )

            enddo
          endif

        endif

      else ! mode >= 2

        if ( mod(mode,2) == 0 ) then
          phase_factor = 1.0_p_k_fld
        else
          phase_factor = -1.0_p_k_fld
        endif

        ! We need to set the on-axis field value and the first (mode - 2) derivatives
        ! equal to 0 for Bz. For mode m, this can be done with the following:
        !
        ! Bz(i * dr) = i^m / m^m * Bz(m * dr) for 1 <= i < m - 1

        mm = mode - 1

        !$omp parallel do
        do i1 = i1_lb_ax, i1_ub_ax
          do i2 = 1, mode - 2
            ! Zero out first few derivatives for B1
            b_re%f2( 1, i1, i2+1 ) = real(i2**mm, p_k_fld) / real(mm**mm, p_k_fld) * b_re%f2( 1, i1, mm+1 )
            b_im%f2( 1, i1, i2+1 ) = real(i2**mm, p_k_fld) / real(mm**mm, p_k_fld) * b_im%f2( 1, i1, mm+1 )
          enddo

          ! B1 must be zero on the axis
          b_re%f2( 1, i1, 1 ) = 0.0_p_k_fld
          b_re%f2( 1, i1, 0 ) = phase_factor * b_re%f2( 1, i1, 2 )
          b_im%f2( 1, i1, 1 ) = 0.0_p_k_fld
          b_im%f2( 1, i1, 0 ) = phase_factor * b_im%f2( 1, i1, 2 )

          ! B2 must be zero on the axis
          b_re%f2( 2, i1, 1 ) = -phase_factor * b_re%f2( 2, i1, 2 )
          b_re%f2( 2, i1, 0 ) = -phase_factor * b_re%f2( 2, i1, 3 )
          b_im%f2( 2, i1, 1 ) = -phase_factor * b_im%f2( 2, i1, 2 )
          b_im%f2( 2, i1, 0 ) = -phase_factor * b_im%f2( 2, i1, 3 )

          ! B3 must be zero on the axis
          b_re%f2( 3, i1, 1 ) = 0.0_p_k_fld
          b_re%f2( 3, i1, 0 ) = -phase_factor * b_re%f2( 3, i1, 2 )
          b_im%f2( 3, i1, 1 ) = 0.0_p_k_fld
          b_im%f2( 3, i1, 0 ) = -phase_factor * b_im%f2( 3, i1, 2 )
        enddo
        !$omp end parallel do

      endif

    endif ! ( gshift_i2 < 0 )

  enddo ! mode

end subroutine dbdt_2d_cyl_modes_fei
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! E-field advance
!-----------------------------------------------------------------------------------------
subroutine dedt_2d_cyl_modes_fei( this, e_cyl_m, b_cyl_m, jay_cyl_m, dt, bnd_type )

  implicit none

  class ( t_emf_solver_cyl_modes_fei ), intent(in) :: this
  type( t_cyl_modes ), intent(inout) :: e_cyl_m
  type( t_cyl_modes ), intent(in) :: b_cyl_m, jay_cyl_m
  real(p_double), intent(in) :: dt
  integer, intent(in), dimension(2,p_max_dim) :: bnd_type

  real(p_k_fld) :: dtdz, dtdr, dtif, dt_r, dt_r_2
  real(p_k_fld) :: tmp, rcp, rcm, tmp_e, f_t, o_t, phase_factor
  integer :: i1, i2, i2_0, gshift_i2, mode, mm
  integer :: j, i1_lb, i1_ub, offset, ord_hf, i1_lb_ax, i1_ub_ax

  type( t_vdf ), pointer :: b_re, b_im, e_re, e_im, jay_re, jay_im

  integer :: n_modes

  n_modes = ubound(e_cyl_m%pf_re,1)

  dtdz = real(dt/e_cyl_m%pf_re(0)%dx_(1), p_k_fld)
  dtdr = real(dt/e_cyl_m%pf_re(0)%dx_(2), p_k_fld)
  dtif = real( dt, p_k_fld )
  f_t = real(4.0_p_k_fld/3.0_p_k_fld, p_k_fld)
  o_t = real(1.0_p_k_fld/3.0_p_k_fld, p_k_fld)

  i1_lb = 2 - this%n_coef
  i1_ub = e_cyl_m%pf_re(0)%nx_(1) + this%n_coef

  i1_lb_ax = i1_lb
  i1_ub_ax = i1_ub

  if ( this%taper_bnd ) then

    offset = this%taper_order / 2

    if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
         bnd_type(p_lower, 1) /= p_bc_other_node ) then
      i1_lb = offset + 1
    endif

    if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
        bnd_type(p_upper, 1) /= p_bc_other_node ) then
      i1_ub = e_cyl_m%pf_re(0)%nx_(1) - offset + 1
    endif

  endif

  gshift_i2 = this%gir_pos - 2

  if ( gshift_i2 < 0 ) then
    ! This node contains the cylindrical axis so start the solver in cell 2
    i2_0 = 2
  else
    i2_0 = 1
  endif


  do mode = 0, n_modes
    b_re => b_cyl_m%pf_re(mode)
    e_re => e_cyl_m%pf_re(mode)
    jay_re => jay_cyl_m%pf_re(mode)
    if (mode > 0) then
      b_im => b_cyl_m%pf_im(mode)
      e_im => e_cyl_m%pf_im(mode)
      jay_im => jay_cyl_m%pf_im(mode)
    endif

    !$omp parallel do private(tmp,dt_r,dt_r_2,rcp,rcm,i1,tmp_e,j)
    do i2 = i2_0, e_re%nx_(2)+2

      tmp = dtdr / (( i2 + gshift_i2 ) - 0.5) ! shifted because of staggering
      dt_r = tmp
      dt_r_2 = dtdr / (( i2 + gshift_i2 ))
      ! dt_r = dtif / (( i2 + gshift_i2 ) - 0.5)
      ! dt_r_2 = dtif / (( i2 + gshift_i2 )) ! at i2 = 2 this is 1/1
      rcp = ( i2 + gshift_i2 ) ! this is 1 at the axis - huh? if i2 = 1 this is 0
      rcm = ( i2 + gshift_i2 - 1 ) ! this must become zero on axis. - if i2 = 1 this is -1

      ! do i1 = 1, e_re%nx_(1)+2
      do i1 = i1_lb, i1_ub
        ! Real Part first

        ! E1
        e_re%f2(1, i1, i2) = e_re%f2(1, i1, i2) &
                        - dtif * jay_re%f2(1, i1, i2) &
                        + tmp * ( rcp * b_re%f2(3, i1, i2) - rcm * b_re%f2(3, i1, i2-1) )
        ! E2
        tmp_e = 0.0_p_k_fld
        do j = 1, this%n_coef
          tmp_e = tmp_e + this%coef_b(j) * ( b_re%f2(3, i1+j-1, i2) - b_re%f2(3, i1-j, i2) )
        enddo
        e_re%f2(2, i1, i2) = e_re%f2(2, i1, i2) &
                            - dtif * jay_re%f2(2, i1, i2) &
                            - dtdz * tmp_e
        ! E3
        tmp_e = 0.0_p_k_fld
        do j = 1, this%n_coef
          tmp_e = tmp_e + this%coef_b(j) * ( b_re%f2(2, i1+j-1, i2) - b_re%f2(2, i1-j, i2) )
        enddo
        e_re%f2(3, i1, i2) = e_re%f2(3, i1, i2) &
                              - dtif * jay_re%f2(3, i1, i2) + dtdz * tmp_e &
                              - dtdr * ( b_re%f2(1, i1, i2) - b_re%f2(1, i1, i2-1))

        if (mode == 0) cycle

        e_re%f2(1, i1, i2) = e_re%f2(1, i1, i2) - dt_r * mode * b_im%f2(2, i1, i2)
        e_re%f2(2, i1, i2) = e_re%f2(2, i1, i2) + dt_r_2 * mode * b_im%f2(1, i1, i2)

        ! Imaginary part
        ! E1
        e_im%f2(1, i1, i2) = e_im%f2(1, i1, i2) &
                      - dtif * jay_im%f2(1, i1, i2) &
                      + tmp * ( rcp * b_im%f2(3, i1, i2) - rcm * b_im%f2(3, i1, i2-1) ) &
                      + dt_r * mode * b_re%f2(2, i1, i2)

        ! E2
        tmp_e = 0.0_p_k_fld
        do j = 1, this%n_coef
          tmp_e = tmp_e + this%coef_b(j) * ( b_im%f2(3, i1+j-1, i2) - b_im%f2(3, i1-j, i2) )
        enddo
        e_im%f2(2, i1, i2) = e_im%f2(2, i1, i2) &
                              - dtif * jay_im%f2(2, i1, i2) - dtdz * tmp_e &
                              - dt_r_2 * mode * b_re%f2(1, i1, i2)

        ! E3
        tmp_e = 0.0_p_k_fld
        do j = 1, this%n_coef
          tmp_e = tmp_e + this%coef_b(j) * ( b_im%f2(2, i1+j-1, i2) - b_im%f2(2, i1-j, i2) )
        enddo
        e_im%f2(3, i1, i2) = e_im%f2(3, i1, i2) &
                              - dtif * jay_im%f2(3, i1, i2) + dtdz * tmp_e &
                              - dtdr * ( b_im%f2(1, i1, i2) - b_im%f2(1, i1, i2-1))

      enddo ! i2
    enddo ! i1
    !$omp end parallel do

    ! do the tapered region
    if ( this%taper_bnd ) then

      ! advance lower boundary
      if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
           bnd_type(p_lower, 1) /= p_bc_other_node ) then

        !$omp parallel do private(tmp,dt_r,dt_r_2,rcp,rcm,ord_hf,tmp_e,j)
        do i2 = i2_0, e_re%nx_(2)+2

          tmp = dtdr / (( i2 + gshift_i2 ) - 0.5) ! shifted because of staggering
          dt_r = tmp
          dt_r_2 = dtdr / (( i2 + gshift_i2 ))
          ! dt_r = dtif / (( i2 + gshift_i2 ) - 0.5)
          ! dt_r_2 = dtif / (( i2 + gshift_i2 )) ! at i2 = 2 this is 1/1
          rcp = ( i2 + gshift_i2 ) ! this is 1 at the axis - huh? if i2 = 1 this is 0
          rcm = ( i2 + gshift_i2 - 1 ) ! this must become zero on axis. - if i2 = 1 this is -1

          ord_hf = 0
          do i1 = 1, i1_lb-1
            ! Real Part first
            ord_hf = ord_hf + 1

            ! E1
            e_re%f2(1, i1, i2) = e_re%f2(1, i1, i2) &
                        - dtif * jay_re%f2(1, i1, i2) &
                        + tmp * ( rcp * b_re%f2(3, i1, i2) - rcm * b_re%f2(3, i1, i2-1) )
            ! E2
            tmp_e = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_re%f2(3, i1+j-1, i2) - b_re%f2(3, i1-j, i2) )
            enddo
            e_re%f2(2, i1, i2) = e_re%f2(2, i1, i2) &
                                - dtif * jay_re%f2(2, i1, i2) &
                                - dtdz * tmp_e
            ! E3
            tmp_e = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_re%f2(2, i1+j-1, i2) - b_re%f2(2, i1-j, i2) )
            enddo
            e_re%f2(3, i1, i2) = e_re%f2(3, i1, i2) &
                                - dtif * jay_re%f2(3, i1, i2) + dtdz * tmp_e &
                                - dtdr * ( b_re%f2(1, i1, i2) - b_re%f2(1, i1, i2-1))

            if (mode == 0) cycle

            e_re%f2(1, i1, i2) = e_re%f2(1, i1, i2) - dt_r * mode * b_im%f2(2, i1, i2)
            e_re%f2(2, i1, i2) = e_re%f2(2, i1, i2) + dt_r_2 * mode * b_im%f2(1, i1, i2)

            ! Imaginary part
            ! E1
            e_im%f2(1, i1, i2) = e_im%f2(1, i1, i2) &
                      - dtif * jay_im%f2(1, i1, i2) &
                      + tmp * ( rcp * b_im%f2(3, i1, i2) - rcm * b_im%f2(3, i1, i2-1) ) &
                      + dt_r * mode * b_re%f2(2, i1, i2)

            ! E2
            tmp_e = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_im%f2(3, i1+j-1, i2) - b_im%f2(3, i1-j, i2) )
            enddo
            e_im%f2(2, i1, i2) = e_im%f2(2, i1, i2) &
                                - dtif * jay_im%f2(2, i1, i2) - dtdz * tmp_e &
                                - dt_r_2 * mode * b_re%f2(1, i1, i2)

            ! E3
            tmp_e = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_im%f2(2, i1+j-1, i2) - b_im%f2(2, i1-j, i2) )
            enddo
            e_im%f2(3, i1, i2) = e_im%f2(3, i1, i2) &
                                - dtif * jay_im%f2(3, i1, i2) + dtdz * tmp_e &
                                - dtdr * ( b_im%f2(1, i1, i2) - b_im%f2(1, i1, i2-1))

          enddo ! i2
        enddo ! i1
        !$omp end parallel do

      endif

      ! advance upper boundary
      if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
           bnd_type(p_upper, 1) /= p_bc_other_node ) then

        !$omp parallel do private(tmp,dt_r,dt_r_2,rcp,rcm,ord_hf,tmp_e,j)
        do i2 = i2_0, e_re%nx_(2)+2

          tmp = dtdr / (( i2 + gshift_i2 ) - 0.5) ! shifted because of staggering
          dt_r = tmp
          dt_r_2 = dtdr / (( i2 + gshift_i2 ))
          ! dt_r = dtif / (( i2 + gshift_i2 ) - 0.5)
          ! dt_r_2 = dtif / (( i2 + gshift_i2 )) ! at i2 = 2 this is 1/1
          rcp = ( i2 + gshift_i2 ) ! this is 1 at the axis - huh? if i2 = 1 this is 0
          rcm = ( i2 + gshift_i2 - 1 ) ! this must become zero on axis. - if i2 = 1 this is -1

          ord_hf = offset + 1
          do i1 = i1_ub+1, e_re%nx_(1)+1
            ! Real Part first
            ord_hf = ord_hf - 1

            ! E1
            e_re%f2(1, i1, i2) = e_re%f2(1, i1, i2) &
                        - dtif * jay_re%f2(1, i1, i2) &
                        + tmp * ( rcp * b_re%f2(3, i1, i2) - rcm * b_re%f2(3, i1, i2-1) )
            ! E2
            tmp_e = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_re%f2(3, i1+j-1, i2) - b_re%f2(3, i1-j, i2) )
            enddo
            e_re%f2(2, i1, i2) = e_re%f2(2, i1, i2) &
                                - dtif * jay_re%f2(2, i1, i2) &
                                - dtdz * tmp_e
            ! E3
            tmp_e = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_re%f2(2, i1+j-1, i2) - b_re%f2(2, i1-j, i2) )
            enddo
            e_re%f2(3, i1, i2) = e_re%f2(3, i1, i2) &
                                - dtif * jay_re%f2(3, i1, i2) + dtdz * tmp_e &
                                - dtdr * ( b_re%f2(1, i1, i2) - b_re%f2(1, i1, i2-1))

            if (mode == 0) cycle

            e_re%f2(1, i1, i2) = e_re%f2(1, i1, i2) - dt_r * mode * b_im%f2(2, i1, i2)
            e_re%f2(2, i1, i2) = e_re%f2(2, i1, i2) + dt_r_2 * mode * b_im%f2(1, i1, i2)

            ! Imaginary part
            ! E1
            e_im%f2(1, i1, i2) = e_im%f2(1, i1, i2) &
                      - dtif * jay_im%f2(1, i1, i2) &
                      + tmp * ( rcp * b_im%f2(3, i1, i2) - rcm * b_im%f2(3, i1, i2-1) ) &
                      + dt_r * mode * b_re%f2(2, i1, i2)

            ! E2
            tmp_e = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_im%f2(3, i1+j-1, i2) - b_im%f2(3, i1-j, i2) )
            enddo
            e_im%f2(2, i1, i2) = e_im%f2(2, i1, i2) &
                                  - dtif * jay_im%f2(2, i1, i2) - dtdz * tmp_e &
                                  - dt_r_2 * mode * b_re%f2(1, i1, i2)

            ! E3
            tmp_e = 0.0_p_k_fld
            do j = 1, ord_hf
              tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_im%f2(2, i1+j-1, i2) - b_im%f2(2, i1-j, i2) )
            enddo
            e_im%f2(3, i1, i2) = e_im%f2(3, i1, i2) &
                                - dtif * jay_im%f2(3, i1, i2) + dtdz * tmp_e &
                                - dtdr * ( b_im%f2(1, i1, i2) - b_im%f2(1, i1, i2-1))

          enddo ! i2
        enddo ! i1
        !$omp end parallel do

      endif

    endif ! end of tapered region

    ! do boundary
    ! note, the guard cells that reside below the axis do not hold true field values,
    ! but hold the field values appropriate for particle interpolation.
    ! i.e., we must have the following:
    ! E_z(-dr/2) = (-1)^m E_z(dr/2)
    ! E_r(-dr) = -(-1)^m E_r(dr)
    ! E_phi(-dr/2) = -(-1)^m E_phi(dr/2)
    ! if quartic interpolation is ever implemented, also handle cell -1 in r
    if ( gshift_i2 < 0 ) then

      if (mode == 0) then
        !$omp parallel do
        do i1 = i1_lb_ax, i1_ub_ax
          ! E1 can be nonzero on the axis
          e_re%f2( 1, i1, 1 ) = e_re%f2( 1, i1, 2 )
          e_re%f2( 1, i1, 0 ) = e_re%f2( 1, i1, 3 )

          ! E2 must be zero on the axis
          e_re%f2( 2, i1, 1 ) = 0.0_p_k_fld
          e_re%f2( 2, i1, 0 ) = - e_re%f2( 2, i1, 2 )

          ! E3 must be zero on the axis
          e_re%f2( 3, i1, 1 ) = - e_re%f2( 3, i1, 2 )
          e_re%f2( 3, i1, 0 ) = - e_re%f2( 3, i1, 3 )
        enddo
        !$omp end parallel do

      elseif (mode == 1) then
        ! First handle the regular regions (everywhere but tapered boundary)
        !$omp parallel do private(tmp_e,j)
        do i1 = i1_lb, i1_ub
          ! E1 must be zero on the axis
          e_re%f2( 1, i1, 1 ) = - e_re%f2( 1, i1, 2 )
          e_re%f2( 1, i1, 0 ) = - e_re%f2( 1, i1, 3 )
          e_im%f2( 1, i1, 1 ) = - e_im%f2( 1, i1, 2 )
          e_im%f2( 1, i1, 0 ) = - e_im%f2( 1, i1, 3 )

          ! E2 can be nonzero on the axis
          ! Rather than solving for e_r(r=0) in the normal way as listed in the paper, we
          ! instead ensure that de_r/dr(r=0) = 0. Setting the derivative to zero reduces
          ! spurious field growth in e_r(r=0) due to the jay term.
          ! tmp_e = 0.0_p_k_fld
          ! do j = 1, this%n_coef
          ! tmp_e = tmp_e + this%coef_b(j) * ( b_re%f2(3, i1+j-1, 1) - b_re%f2(3, i1-j, 1) )
          ! enddo
          ! e_re%f2( 2, i1, 1 ) = e_re%f2( 2, i1, 1 ) &
          ! - dtif * jay_re%f2(2,i1,1) - dtdz * tmp_e &
          ! + dtdr * b_im%f2(1, i1, 2)
          e_re%f2( 2, i1, 1 ) = f_t * e_re%f2( 2, i1, 2 ) - o_t * e_re%f2( 2, i1, 3 )

          e_re%f2( 2, i1, 0 ) = e_re%f2( 2, i1, 2 )

          ! tmp_e = 0.0_p_k_fld
          ! do j = 1, this%n_coef
          ! tmp_e = tmp_e + this%coef_b(j) * ( b_im%f2(3, i1+j-1, 1) - b_im%f2(3, i1-j, 1) )
          ! enddo
          ! e_im%f2( 2, i1, 1 ) = e_im%f2( 2, i1, 1 ) &
          ! - dtif * jay_im%f2(2,i1,1) - dtdz * tmp_e &
          ! - dtdr * b_re%f2(1, i1, 2)
          e_im%f2( 2, i1, 1 ) = f_t * e_im%f2( 2, i1, 2 ) - o_t * e_im%f2( 2, i1, 3 )

          e_im%f2( 2, i1, 0 ) = e_im%f2( 2, i1, 2 )

          ! E3 can be nonzero on the axis
          ! Rather than setting the guard cells equal to the simulation space values, we
          ! instead ensure that e_r = i e_phi at r = 0. For our grid this translates to
          ! e_r_re = -e_phi_im and e_r_im = e_phi_re at r = 0.
          ! e_re%f2( 3, i1, 1 ) = e_re%f2( 3, i1, 2 )
          ! e_re%f2( 3, i1, 0 ) = e_re%f2( 3, i1, 3 )
          ! e_im%f2( 3, i1, 1 ) = e_im%f2( 3, i1, 2 )
          ! e_im%f2( 3, i1, 0 ) = e_im%f2( 3, i1, 3 )
          e_re%f2( 3, i1, 1 ) = 2.0_p_k_fld*e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 2 )
          e_re%f2( 3, i1, 0 ) = 2.0_p_k_fld*e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 3 )
          e_im%f2( 3, i1, 1 ) =-2.0_p_k_fld*e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 2 )
          e_im%f2( 3, i1, 0 ) =-2.0_p_k_fld*e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 3 )
        enddo
        !$omp end parallel do

        ! tapered region
        if ( this%taper_bnd ) then

          ! lower boundary
          if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
               bnd_type(p_lower, 1) /= p_bc_other_node ) then

            ord_hf = 0
            do i1 = 1, i1_lb-1
              ord_hf = ord_hf + 1

              ! E1 must be zero on the axis
              e_re%f2( 1, i1, 1 ) = - e_re%f2( 1, i1, 2 )
              e_re%f2( 1, i1, 0 ) = - e_re%f2( 1, i1, 3 )
              e_im%f2( 1, i1, 1 ) = - e_im%f2( 1, i1, 2 )
              e_im%f2( 1, i1, 0 ) = - e_im%f2( 1, i1, 3 )

              ! E2 can be nonzero on the axis
              ! tmp_e = 0.0_p_k_fld
              ! do j = 1, ord_hf
              ! tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_re%f2(3, i1+j-1, 1) - b_re%f2(3, i1-j, 1) )
              ! enddo
              ! e_re%f2( 2, i1, 1 ) = e_re%f2( 2, i1, 1 ) &
              ! - dtif * jay_re%f2(2,i1,1) - dtdz * tmp_e &
              ! + dtdr * b_im%f2(1, i1, 2)
              e_re%f2( 2, i1, 1 ) = f_t * e_re%f2( 2, i1, 2 ) - o_t * e_re%f2( 2, i1, 3 )

              e_re%f2( 2, i1, 0 ) = e_re%f2( 2, i1, 2 )

              ! tmp_e = 0.0_p_k_fld
              ! do j = 1, ord_hf
              ! tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_im%f2(3, i1+j-1, 1) - b_im%f2(3, i1-j, 1) )
              ! enddo
              ! e_im%f2( 2, i1, 1 ) = e_im%f2( 2, i1, 1 ) &
              ! - dtif * jay_im%f2(2,i1,1) - dtdz * tmp_e &
              ! - dtdr * b_re%f2(1, i1, 2)
              e_im%f2( 2, i1, 1 ) = f_t * e_im%f2( 2, i1, 2 ) - o_t * e_im%f2( 2, i1, 3 )

              e_im%f2( 2, i1, 0 ) = e_im%f2( 2, i1, 2 )

              ! E3 can be nonzero on the axis
              ! e_re%f2( 3, i1, 1 ) = e_re%f2( 3, i1, 2 )
              ! e_re%f2( 3, i1, 0 ) = e_re%f2( 3, i1, 3 )
              ! e_im%f2( 3, i1, 1 ) = e_im%f2( 3, i1, 2 )
              ! e_im%f2( 3, i1, 0 ) = e_im%f2( 3, i1, 3 )
              e_re%f2( 3, i1, 1 ) = 2.0_p_k_fld*e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 2 )
              e_re%f2( 3, i1, 0 ) = 2.0_p_k_fld*e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 3 )
              e_im%f2( 3, i1, 1 ) =-2.0_p_k_fld*e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 2 )
              e_im%f2( 3, i1, 0 ) =-2.0_p_k_fld*e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 3 )
            enddo

          endif

          ! upper boundary
          if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
               bnd_type(p_upper, 1) /= p_bc_other_node ) then

            ord_hf = offset + 1
            do i1 = i1_ub+1, e_re%nx_(1)+1
              ord_hf = ord_hf - 1

              ! E1 must be zero on the axis
              e_re%f2( 1, i1, 1 ) = - e_re%f2( 1, i1, 2 )
              e_re%f2( 1, i1, 0 ) = - e_re%f2( 1, i1, 3 )
              e_im%f2( 1, i1, 1 ) = - e_im%f2( 1, i1, 2 )
              e_im%f2( 1, i1, 0 ) = - e_im%f2( 1, i1, 3 )

              ! E2 can be nonzero on the axis
              ! tmp_e = 0.0_p_k_fld
              ! do j = 1, ord_hf
              ! tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_re%f2(3, i1+j-1, 1) - b_re%f2(3, i1-j, 1) )
              ! enddo
              ! e_re%f2( 2, i1, 1 ) = e_re%f2( 2, i1, 1 ) &
              ! - dtif * jay_re%f2(2,i1,1) - dtdz * tmp_e &
              ! + dtdr * b_im%f2(1, i1, 2)
              e_re%f2( 2, i1, 1 ) = f_t * e_re%f2( 2, i1, 2 ) - o_t * e_re%f2( 2, i1, 3 )

              e_re%f2( 2, i1, 0 ) = e_re%f2( 2, i1, 2 )

              ! tmp_e = 0.0_p_k_fld
              ! do j = 1, ord_hf
              ! tmp_e = tmp_e + this%taper_coef_b(j,ord_hf) * ( b_im%f2(3, i1+j-1, 1) - b_im%f2(3, i1-j, 1) )
              ! enddo
              ! e_im%f2( 2, i1, 1 ) = e_im%f2( 2, i1, 1 ) &
              ! - dtif * jay_im%f2(2,i1,1) - dtdz * tmp_e &
              ! - dtdr * b_re%f2(1, i1, 2)
              e_im%f2( 2, i1, 1 ) = f_t * e_im%f2( 2, i1, 2 ) - o_t * e_im%f2( 2, i1, 3 )

              e_im%f2( 2, i1, 0 ) = e_im%f2( 2, i1, 2 )

              ! E3 can be nonzero on the axis
              ! e_re%f2( 3, i1, 1 ) = e_re%f2( 3, i1, 2 )
              ! e_re%f2( 3, i1, 0 ) = e_re%f2( 3, i1, 3 )
              ! e_im%f2( 3, i1, 1 ) = e_im%f2( 3, i1, 2 )
              ! e_im%f2( 3, i1, 0 ) = e_im%f2( 3, i1, 3 )
              e_re%f2( 3, i1, 1 ) = 2.0_p_k_fld*e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 2 )
              e_re%f2( 3, i1, 0 ) = 2.0_p_k_fld*e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 3 )
              e_im%f2( 3, i1, 1 ) =-2.0_p_k_fld*e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 2 )
              e_im%f2( 3, i1, 0 ) =-2.0_p_k_fld*e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 3 )
            enddo

          endif
        endif

      else ! mode >= 2

        if ( mod(mode,2) == 0 ) then
          phase_factor = 1.0_p_k_fld
        else
          phase_factor = -1.0_p_k_fld
        endif

        ! We need to set the on-axis field value and the first (mode - 2) derivatives
        ! equal to 0 for Ez. For mode m, this can be done with the following:
        !
        ! Ez((2*i-1) / 2 * dr) = (2*i - 1)^m / (2*m-1)^m * Ephi((2*m-1) / 2 * dr)
        ! for 1 <= i < m - 1

        mm = mode - 1

        !$omp parallel do
        do i1 = i1_lb_ax, i1_ub_ax
          do i2 = 1, mode - 2
            ! Zero out first few derivatives for E1
            e_re%f2( 1, i1, i2+1 ) = real((2*i2-1)**mm, p_k_fld) / real((2*mm-1)**mm, p_k_fld) * e_re%f2( 1, i1, mm+1 )
            e_im%f2( 1, i1, i2+1 ) = real((2*i2-1)**mm, p_k_fld) / real((2*mm-1)**mm, p_k_fld) * e_im%f2( 1, i1, mm+1 )
          enddo

          ! E1 must be zero on the axis
          e_re%f2( 1, i1, 1 ) = phase_factor * e_re%f2( 1, i1, 2 )
          e_re%f2( 1, i1, 0 ) = phase_factor * e_re%f2( 1, i1, 3 )
          e_im%f2( 1, i1, 1 ) = phase_factor * e_im%f2( 1, i1, 2 )
          e_im%f2( 1, i1, 0 ) = phase_factor * e_im%f2( 1, i1, 3 )

          ! E2 must be zero on the axis
          e_re%f2( 2, i1, 1 ) = 0.0_p_k_fld
          e_re%f2( 2, i1, 0 ) = -phase_factor * e_re%f2( 2, i1, 2 )
          e_im%f2( 2, i1, 1 ) = 0.0_p_k_fld
          e_im%f2( 2, i1, 0 ) = -phase_factor * e_im%f2( 2, i1, 2 )

          ! E3 must be zero on the axis
          e_re%f2( 3, i1, 1 ) = -phase_factor * e_re%f2( 3, i1, 2 )
          e_re%f2( 3, i1, 0 ) = -phase_factor * e_re%f2( 3, i1, 3 )
          e_im%f2( 3, i1, 1 ) = -phase_factor * e_im%f2( 3, i1, 2 )
          e_im%f2( 3, i1, 0 ) = -phase_factor * e_im%f2( 3, i1, 3 )
        enddo
        !$omp end parallel do

      endif

    endif ! ( gshift_i2 < 0 )

  enddo ! mode

end subroutine dedt_2d_cyl_modes_fei
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Test algorithm stability for resolution and time-step
!-------------------------------------------------------------------------------
subroutine test_stability_cyl_modes_fei( this, dt, dx )

    implicit none

    class( t_emf_solver_cyl_modes_fei ), intent(in) :: this
    real(p_double), intent(in) :: dt
    real(p_double), dimension(:), intent(in) :: dx

    real(p_double) :: cour_new, cour_old, tmp, rdtheta
    real(p_k_fld), dimension(p_max_n_coef) :: coef_e, coef_b
    procedure(kernel), pointer :: kernel_ptr
    integer :: i, j, cnt
    integer, parameter :: cnt_max = 100

    if (mpi_node()==0) print *,'Calculating Courant stability condition...'

    cnt = 0
    cour_new = 0.0_p_double
    cour_old = 0.0_p_double
    tmp = 0.0_p_double

    coef_e(1:this%n_coef) = this%coef_e(1:this%n_coef)
    coef_b(1:this%n_coef) = this%coef_b(1:this%n_coef)

    ! Add another effective "cell" of size dr * pi / n_cyl_modes in theta
    ! The extra term in the sum is then 2 / (dr * dtheta)
    if ( this%n_cyl_modes > 0 ) then
        rdtheta = 1.0_p_double / (dx(p_x_dim)*pi/(2.0_p_double*this%n_cyl_modes))**2
    else
        rdtheta = 1.0_p_double / (dx(p_x_dim)*pi)**2
    endif

    select case ( this%type )

    case ( p_emf_fei_std, p_emf_fei_bump, p_emf_fei_coef )

        do j = 1, this%n_coef
            cour_new = cour_new + abs( coef_e(j) )
            tmp = tmp + abs( coef_b(j) )
        enddo

        cour_new = cour_new * tmp / (dx(1))**2

        do i = 2, p_x_dim
            cour_new = cour_new + 1.0_p_double / (dx(i))**2
        enddo
        cour_new = cour_new + rdtheta

        cour_new = sqrt( 1.0_p_double / cour_new )

    case ( p_emf_fei_xu )

        ! set the CFL of Yee as the initial guess
        do i = 1, p_x_dim
            cour_old = cour_old + 1.0_p_double / (dx(i))**2
        enddo
        cour_old = cour_old + rdtheta
        cour_old = sqrt( 1.0_p_double / cour_old )

        ! iterate to get CFL
        kernel_ptr => kernel_xu
        cnt = 0
        do

            ! calculate the solver stencil coefficients
            call gen_solver_coef( kernel_ptr, (/ cour_old/dx(1) /), this%solver_ord, &
                                    this%n_coef, coef_e, this%weight_w, this%weight_n )
            coef_b(1:this%n_coef) = coef_e(1:this%n_coef)

            ! calculate the updated CFL
            tmp = 0.0_p_double
            cour_new = 0.0_p_double
            do j = 1, this%n_coef
                cour_new = cour_new + abs( coef_e(j) )
                tmp = tmp + abs( coef_b(j) )
            enddo

            cour_new = cour_new * tmp / (dx(1))**2

            do i = 2, p_x_dim
                cour_new = cour_new + 1.0_p_double / (dx(i))**2
            enddo
            cour_new = cour_new + rdtheta

            cour_new = sqrt( 1.0_p_double / cour_new )
            cnt = cnt + 1

            if ( abs( cour_new - cour_old ) / cour_old < 1.0d-6 ) then
                exit
            elseif ( cnt > cnt_max ) then
                if (mpi_node()==0) print *,'(*error*) Calculation of (Fei) EMF solver stability condition may not converge, aborting'
                call abort_program(p_err_invalid)
            else
                cour_old = cour_new
                cour_new = 0.0_p_double
            endif

        enddo

    case ( p_emf_fei_dual )

        ! set the CFL of Yee as the initial guess
        do i = 1, p_x_dim
            cour_old = cour_old + 1.0_p_double / (dx(i))**2
        enddo
        cour_old = cour_old + rdtheta
        cour_old = sqrt( 1.0_p_double / cour_old )

        ! iterate to get CFL
        cnt = 0
        do

            ! calculate the solver stencil coefficients
            kernel_ptr => kernel_dual_e
            call gen_solver_coef( kernel_ptr, (/ cour_old/dx(1) /), this%solver_ord, &
                                    this%n_coef, coef_e, this%weight_w, this%weight_n )
            kernel_ptr => kernel_dual_b
            call gen_solver_coef( kernel_ptr, (/ cour_old/dx(1) /), this%solver_ord, &
                                    this%n_coef, coef_b, this%weight_w, this%weight_n )

            ! calculate the updated CFL
            tmp = 0.0_p_double
            cour_new = 0.0_p_double
            do j = 1, this%n_coef
                cour_new = cour_new + abs( coef_e(j) )
                tmp = tmp + abs( coef_b(j) )
            enddo

            cour_new = cour_new * tmp / (dx(1))**2

            do i = 2, p_x_dim
                cour_new = cour_new + 1.0_p_double / (dx(i))**2
            enddo
            cour_new = cour_new + rdtheta

            cour_new = sqrt( 1.0_p_double / cour_new )
            cnt = cnt + 1

            if ( abs( cour_new - cour_old ) / cour_old < 1.0d-6 ) then
                exit
            elseif ( cnt > cnt_max ) then
                if (mpi_node()==0) print *,'(*error*) Calculation of (Fei) EMF solver stability condition may not converge, aborting'
                call abort_program(p_err_invalid)
            else
                cour_old = cour_new
                cour_new = 0.0_p_double
            endif

        enddo

    end select

    if (dt > cour_new) then
        if ( mpi_node() == 0 ) then
            print *, '(*error*) Customized (Fei) EMF solver stability condition violated, aborting'
            print *, 'Please reduce the time step and re-try.'
            print *, 'Note that CFL will change as the time step changes.'
            print *, 'dx   = ', dx(1:p_x_dim)
            print *, 'dt   = ', dt, ' dt(CFL) = ', cour_new
        endif
        call abort_program(p_err_invalid)
    endif

end subroutine test_stability_cyl_modes_fei

!-------------------------------------------------------------------------------
! Cleanup solver buffers if necessary
!-------------------------------------------------------------------------------
subroutine cleanup_cyl_modes_fei( this )

    implicit none

    class( t_emf_solver_cyl_modes_fei ), intent(inout) :: this

    call this % t_emf_solver_fei % cleanup()

    if ( associated(this%f_buffer_re) ) call freemem(this%f_buffer_re,"cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03",1301)
    if ( associated(this%f_buffer_im) ) call freemem(this%f_buffer_im,"cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03",1302)

end subroutine cleanup_cyl_modes_fei

!-----------------------------------------------------------------------------------------
! Filter / correct current
!-----------------------------------------------------------------------------------------
subroutine filter_correct_current_cyl_modes_2d( this, jay_cyl_m )
!-----------------------------------------------------------------------------------------

  use, intrinsic :: iso_c_binding

  implicit none

  class( t_emf_solver_cyl_modes_fei ), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m

  real(p_k_fld) :: filter_limit, filter_width
  real(p_k_fld) :: corr1, corr2
  real(p_k_fld) :: dm1, k1, k1_half
  integer :: gc_l1, gc_u1, gc_l2, gc_u2, n1, n2, m1, m2, dc, fft_m1
  integer, dimension(2,2) :: gc_num
  integer, dimension(3) :: bshape
  integer :: i1, i2, j, mode
  type( t_vdf ), pointer :: jay_re, jay_im

  filter_width = this%filter_width
  filter_limit = this%filter_limit
  gc_num = jay_cyl_m%pf_re(0)%gc_num()
  gc_l1 = gc_num(p_lower,1)
  gc_u1 = gc_num(p_upper,1)
  gc_l2 = gc_num(p_lower,2)
  gc_u2 = gc_num(p_upper,2)
  n1 = jay_cyl_m%pf_re(0)%nx_(1)
  n2 = jay_cyl_m%pf_re(0)%nx_(2)
  m1 = n1 + gc_l1 + gc_u1
  m2 = n2 + gc_l2 + gc_u2
  dm1 = 2.0_p_k_fld * pi / real(m1, p_k_fld)
  dc = this%n_damp_cell


  fft_m1 = 2 * ( m1/2 + 1 )

  do mode = 0, ubound(jay_cyl_m%pf_re, 1)

    jay_re => jay_cyl_m%pf_re(mode)
    if ( .not. associated( this%f_buffer_re ) ) then
      !$omp critical
      call alloc(this%f_buffer_re, (/ 3, fft_m1, m2 /),"cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03",1350)
      call this % fft % init( m1 )
      !$omp end critical
    else
      ! check if buffer should be resized (e.g., from dlb)
      bshape = shape(this%f_buffer_re)
      if ( bshape(2)/=fft_m1 .or. bshape(3)/=m2 ) then
        !$omp critical
        call freemem(this%f_buffer_re,"cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03",1358)
        call alloc(this%f_buffer_re, (/ 3, fft_m1, m2 /),"cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03",1359)
        ! re-init fft if necessary
        if ( bshape(2)/=fft_m1 ) then
          call this % fft % cleanup()
          call this % fft % init( m1 )
        endif
        !$omp end critical
      endif
    endif
    if ( mode > 0 ) then
      jay_im => jay_cyl_m%pf_im(mode)
      if ( .not. associated( this%f_buffer_im ) ) then
        !$omp critical
        call alloc(this%f_buffer_im, (/ 3, fft_m1, m2 /),"cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03",1372)
        !$omp end critical
      else
        ! check if buffer should be resized (e.g., from dlb)
        bshape = shape(this%f_buffer_im)
        if ( bshape(2)/=fft_m1 .or. bshape(3)/=m2 ) then
          !$omp critical
          call freemem(this%f_buffer_im,"cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03",1379)
          call alloc(this%f_buffer_im, (/ 3, fft_m1, m2 /),"cyl_modes/emf_solver/os-emf-solver-cyl-modes-fei.f03",1380)
          !$omp end critical
        endif
      endif
    endif

    this%f_buffer_re = 0.0_p_k_fld
    if ( mode > 0 ) this%f_buffer_im = 0.0_p_k_fld

    ! damp the current
    do i2 = 1, m2
      do i1 = 1, m1

        if ( i1 <= dc ) then
          corr1 = cos( 0.5_p_k_fld * pi * (i1-dc-1) / dc )**2
        elseif ( i1 > m1-dc ) then
          corr1 = cos( 0.5_p_k_fld * pi * (i1-m1+dc) / dc )**2
        else
          corr1 = 1.0_p_k_fld
        endif

        this%f_buffer_re(:,i1,i2) = jay_re%f2( :, i1-gc_l1, i2-gc_l2 ) * corr1

        if ( mode == 0 ) cycle

        this%f_buffer_im(:,i1,i2) = jay_im%f2( :, i1-gc_l1, i2-gc_l2 ) * corr1

      enddo
    enddo

    call this % fft % fft( this%f_buffer_re, m1, m2 )
    if ( mode > 0 ) call this % fft % fft( this%f_buffer_im, m1, m2 )

    ! filter and corret current
    do i2 = 1, m2
      do i1 = 2, fft_m1/2
        k1 = dm1 * real(i1-1, p_k_fld)
        k1_half = k1 * 0.5_p_k_fld

        ! current filter factor
        corr1 = 1.0_p_k_fld
        if ( ( k1 > filter_limit*pi ) .and. &
              ( this%filter_current .eqv. .true. ) ) then

          if ( k1 < ( filter_limit+filter_width) * pi ) then
            corr1 = cos( 0.5_p_k_fld * pi * (k1 - filter_limit*pi) / (filter_width*pi) )**2
          else
            corr1 = 0.0_p_k_fld
          endif

        endif

        ! current corretion factor
        if ( this%correct_current .eqv. .true. ) then

          corr2 = 0.0_p_k_fld
          do j = 1, this%n_coef
            corr2 = corr2 + this%coef_b(j) * sin( real(2*j-1, p_k_fld) * k1_half )
          enddo
          corr2 = sin(k1_half) / corr2

        else

          corr2 = 1.0_p_k_fld

        endif

        this%f_buffer_re(1,2*i1-1,i2) = this%f_buffer_re(1,2*i1-1,i2) * corr1 * corr2
        this%f_buffer_re(2,2*i1-1,i2) = this%f_buffer_re(2,2*i1-1,i2) * corr1 * corr2
        this%f_buffer_re(3,2*i1-1,i2) = this%f_buffer_re(3,2*i1-1,i2) * corr1
        this%f_buffer_re(1,2*i1 ,i2) = this%f_buffer_re(1,2*i1 ,i2) * corr1
        this%f_buffer_re(2,2*i1 ,i2) = this%f_buffer_re(2,2*i1 ,i2) * corr1
        this%f_buffer_re(3,2*i1 ,i2) = this%f_buffer_re(3,2*i1 ,i2) * corr1

        if ( mode == 0 ) cycle

        this%f_buffer_im(1,2*i1-1,i2) = this%f_buffer_im(1,2*i1-1,i2) * corr1 * corr2
        this%f_buffer_im(2,2*i1-1,i2) = this%f_buffer_im(2,2*i1-1,i2) * corr1 * corr2
        this%f_buffer_im(3,2*i1-1,i2) = this%f_buffer_im(3,2*i1-1,i2) * corr1
        this%f_buffer_im(1,2*i1 ,i2) = this%f_buffer_im(1,2*i1 ,i2) * corr1
        this%f_buffer_im(2,2*i1 ,i2) = this%f_buffer_im(2,2*i1 ,i2) * corr1
        this%f_buffer_im(3,2*i1 ,i2) = this%f_buffer_im(3,2*i1 ,i2) * corr1

      enddo
    enddo


    call this % fft % ifft( this%f_buffer_re, m1, m2 )
    if ( mode > 0 ) call this % fft % ifft( this%f_buffer_im, m1, m2 )

    ! direct copy
    do i2 = 1, m2
      do i1 = 1,m1
        jay_re%f2(:, i1-gc_l1, i2-gc_l2) = this%f_buffer_re(:,i1,i2)
        if ( mode == 0 ) cycle
        jay_im%f2(:, i1-gc_l1, i2-gc_l2) = this%f_buffer_im(:,i1,i2)
      enddo
    enddo

  enddo

end subroutine filter_correct_current_cyl_modes_2d
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Initializes the solver
!-------------------------------------------------------------------------------
subroutine init_cyl_modes_fei( this, gix_pos, coords )

    implicit none

    class( t_emf_solver_cyl_modes_fei ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords

    call this % t_emf_solver_fei % init( gix_pos, coords )

    ! Only handle cylindrical (r-z) coordinates solvers
    this % dedt_cyl_modes => dedt_2d_cyl_modes_fei
    this % dbdt_cyl_modes => dbdt_2d_cyl_modes_fei
    this % filter_corr_current_cyl_modes => filter_correct_current_cyl_modes_2d

end subroutine init_cyl_modes_fei

subroutine advance_cyl_modes_fei( this, e_cyl_m, b_cyl_m, jay_cyl_m, dt, bnd_con )

    implicit none

    class( t_emf_solver_cyl_modes_fei ), intent(inout) :: this
    type( t_cyl_modes ), intent( inout ) :: e_cyl_m, b_cyl_m, jay_cyl_m
    class( t_emf_bound_cyl_modes ), intent(inout) :: bnd_con
    real(p_double), intent(in) :: dt

    real(p_double) :: dt_b, dt_e

    ! Filter current if necessary
    if ( this%filter_current .or. this%correct_current ) then
        call this % filter_corr_current_cyl_modes( jay_cyl_m )
    endif

    ! Advance fields
    dt_b = dt / 2.0_p_double
    dt_e = dt






    ! Advance B half time step
    call this % dbdt_cyl_modes( b_cyl_m, e_cyl_m, dt_b, bnd_con%type )
    call bnd_con%update_boundary_b_cm( e_cyl_m, b_cyl_m, step = 1 )






    ! Advance E one full time step
    call this % dedt_cyl_modes( e_cyl_m, b_cyl_m, jay_cyl_m, dt_e, bnd_con%type )
    call bnd_con%update_boundary_e_cm( e_cyl_m, b_cyl_m )






    ! Advance B another half time step
    call this % dbdt_cyl_modes( b_cyl_m, e_cyl_m, dt_b, bnd_con%type )
    call bnd_con%update_boundary_b_cm( e_cyl_m, b_cyl_m, step = 2 )






end subroutine advance_cyl_modes_fei

end module m_emf_solver_cyl_modes_fei
