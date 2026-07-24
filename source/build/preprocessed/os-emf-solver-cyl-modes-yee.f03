# 1 "cyl_modes/emf_solver/os-emf-solver-cyl-modes-yee.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/emf_solver/os-emf-solver-cyl-modes-yee.f03" 2
!-----------------------------------------------------------------------------------------
! Yee solver
! K. YEE, NUMERICAL SOLUTION OF INITIAL BOUNDARY VALUE PROBLEMS INVOLVING MAXWELLS
! EQUATIONS IN ISOTROPIC MEDIA, IEEE Transactions on Antenna Propagation, vol. 14,
! no. 3, pp. 302-307, 1966.
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
# 9 "cyl_modes/emf_solver/os-emf-solver-cyl-modes-yee.f03" 2
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
# 10 "cyl_modes/emf_solver/os-emf-solver-cyl-modes-yee.f03" 2

module m_emf_solver_cyl_modes_yee

use m_parameters
use m_emf_solver_yee
use m_vdf_define, only : t_vdf
use m_cyl_modes, only : t_cyl_modes, check_nan
use m_emf_bound_cyl_modes, only : t_emf_bound_cyl_modes
use m_math, only : pi

private

type, extends( t_emf_solver_yee ), public :: t_emf_solver_cyl_modes_yee

    integer :: n_cyl_modes = -1
    procedure(dedt_cyl_modes_yee), pointer :: dedt_cyl_modes => null()
    procedure(dbdt_cyl_modes_yee), pointer :: dbdt_cyl_modes => null()

contains

    procedure :: init => init_cyl_modes_yee
    procedure :: advance_cyl_modes => advance_cyl_modes_yee
    procedure :: test_stability => test_stability_cyl_modes_yee

end type t_emf_solver_cyl_modes_yee

interface
    subroutine dedt_cyl_modes_yee( this, e_cyl_m, b_cyl_m, jay_cyl_m, dt )
        import t_emf_solver_cyl_modes_yee, t_cyl_modes, p_double
        class( t_emf_solver_cyl_modes_yee ), intent(in) :: this
        type( t_cyl_modes ), intent(inout) :: e_cyl_m
        type( t_cyl_modes ), intent(in) :: b_cyl_m, jay_cyl_m
        real(p_double), intent(in) :: dt
    end subroutine
end interface

interface
    subroutine dbdt_cyl_modes_yee( this, b_cyl_m, e_cyl_m, dt )
        import t_emf_solver_cyl_modes_yee, t_cyl_modes, p_double
        class( t_emf_solver_cyl_modes_yee ), intent(in) :: this
        type( t_cyl_modes ), intent(inout) :: b_cyl_m
        type( t_cyl_modes ), intent(in) :: e_cyl_m
        real(p_double), intent(in) :: dt
    end subroutine
end interface


contains

!-----------------------------------------------------------------------------------------
! B-field advance
!-----------------------------------------------------------------------------------------
subroutine dbdt_2d_cyl_modes_yee( this, b_cyl_m, e_cyl_m, dt )

  implicit none

  class ( t_emf_solver_cyl_modes_yee ), intent(in) :: this
  type( t_cyl_modes ), intent(inout) :: b_cyl_m
  type( t_cyl_modes ), intent(in) :: e_cyl_m
  real(p_double), intent(in) :: dt

  real(p_k_fld) :: dtdz, dtdr, dtif, dt_r, dt_r_2
  real(p_k_fld) :: tmp, rm, rp, n_e, o_e, phase_factor
  integer :: i1, i2, i2_0, gshift_i2, mode, mm
  type( t_vdf ), pointer :: b_re, b_im, e_re, e_im
  integer :: n_modes

  n_modes = ubound(e_cyl_m%pf_re,1)

  dtdz = real(dt/b_cyl_m%pf_re(0)%dx_(1), p_k_fld)
  dtdr = real(dt/b_cyl_m%pf_re(0)%dx_(2), p_k_fld)
  dtif = real( dt, p_k_fld )
  n_e = real(9.0_p_k_fld/8.0_p_k_fld, p_k_fld)
  o_e = real(1.0_p_k_fld/8.0_p_k_fld, p_k_fld)

  gshift_i2 = this % gir_pos - 2

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

    !$omp parallel do private(rp,rm,tmp,dt_r,dt_r_2,i1)
    do i2 = i2_0, b_re%nx_(2)+1
      rp = i2 + gshift_i2 + 0.5 ! position of the upper edge of the cell (normalized to dr) ! at the axis this is 1.5
      rm = i2 + gshift_i2 - 0.5 ! position of the lower edge of the cell (normalized to dr) ! at the axis this is 0.5
      tmp = dtdr/(i2 + gshift_i2) ! (dt/dr) / position of the middle of the cell (normalized to dr)
      dt_r = tmp
      dt_r_2 = dtdr/(i2 + gshift_i2 - 0.5)
      !dt_r = dtif/(i2 + gshift_i2) ! (dt) / position of the middle of the cell (normalized to dr) ! at the axis this is 1
      !dt_r_2 = dtif/(i2 + gshift_i2 - 0.5) ! (dt) / position of the middle of the cell (normalized to dr)

      do i1 = 0, b_re%nx_(1)+1
        ! First the Real part
        !B1
        b_re%f2( 1, i1, i2 ) = b_re%f2(1, i1, i2) &
                               - tmp * ( rp * e_re%f2( 3, i1, i2+1 ) - rm * e_re%f2( 3, i1, i2 ))
        !B2
        b_re%f2( 2, i1, i2 ) = b_re%f2(2, i1, i2) &
                               + dtdz * ( e_re%f2( 3, i1+1, i2 ) - e_re%f2( 3, i1, i2))
        !B3
        b_re%f2( 3, i1, i2 ) = b_re%f2( 3, i1, i2 ) &
                        - dtdz * ( e_re%f2( 2, i1+1, i2 ) - e_re%f2( 2, i1, i2 )) &
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
        b_im%f2( 2, i1, i2 ) = b_im%f2(2, i1, i2) &
                               + dtdz * ( e_im%f2(3, i1+1, i2) - e_im%f2(3, i1, i2)) &
                               + dt_r_2 * mode * e_re%f2(1, i1, i2)

        !B3
        b_im%f2( 3, i1, i2 ) = b_im%f2( 3, i1, i2 ) &
                               - dtdz * ( e_im%f2( 2, i1+1, i2 ) - e_im%f2( 2, i1, i2 )) + &
                         dtdr * ( e_im%f2( 1, i1, i2+1 ) - e_im%f2( 1, i1, i2 ))

      enddo ! i1

    enddo ! i2
    !$omp end parallel do

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
        do i1 = 0, b_re%nx_(1)+1
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
        !$omp parallel do
        do i1 = 0, b_re%nx_(1)+1
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
          ! b_re%f2( 3, i1, 1 ) = b_re%f2( 3, i1, 1 ) - &
          ! dtdz * ( e_re%f2( 2, i1+1, 1 ) - e_re%f2( 2, i1, 1 )) + &
          ! dtdr * 2 * e_re%f2( 1, i1, 2 )
          ! b_im%f2( 3, i1, 1 ) = b_im%f2( 3, i1, 1 ) - &
          ! dtdz * ( e_im%f2( 2, i1+1, 1 ) - e_im%f2( 2, i1, 1 )) + &
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
        do i1 = 0, b_re%nx_(1)+1
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

end subroutine dbdt_2d_cyl_modes_yee
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! E-field advance
!-----------------------------------------------------------------------------------------
subroutine dedt_2d_cyl_modes_yee( this, e_cyl_m, b_cyl_m, jay_cyl_m, dt )

  implicit none

  class ( t_emf_solver_cyl_modes_yee ), intent(in) :: this
  type( t_cyl_modes ), intent(inout) :: e_cyl_m
  type( t_cyl_modes ), intent(in) :: b_cyl_m, jay_cyl_m
  real(p_double), intent(in) :: dt

  real(p_k_fld) :: dtdz, dtdr, dtif, dt_r, dt_r_2
  real(p_k_fld) :: tmp, rcp, rcm, f_t, o_t, phase_factor
  integer :: i1, i2, i2_0, gshift_i2, mode, mm

  type( t_vdf ), pointer :: b_re, b_im, e_re, e_im, jay_re, jay_im

  integer :: n_modes

  n_modes = ubound(e_cyl_m%pf_re,1)

  dtdz = real(dt/e_cyl_m%pf_re(0)%dx_(1), p_k_fld)
  dtdr = real(dt/e_cyl_m%pf_re(0)%dx_(2), p_k_fld)
  dtif = real( dt, p_k_fld )
  f_t = real(4.0_p_k_fld/3.0_p_k_fld, p_k_fld)
  o_t = real(1.0_p_k_fld/3.0_p_k_fld, p_k_fld)

  gshift_i2 = this % gir_pos - 2

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

    !$omp parallel do private(tmp,dt_r,dt_r_2,rcp,rcm,i1)
    do i2 = i2_0, e_re%nx_(2)+2

      tmp = dtdr / (( i2 + gshift_i2 ) - 0.5) ! shifted because of staggering
      dt_r = tmp
      dt_r_2 = dtdr / (( i2 + gshift_i2 ))
      ! dt_r = dtif / (( i2 + gshift_i2 ) - 0.5)
      ! dt_r_2 = dtif / (( i2 + gshift_i2 )) ! at i2 = 2 this is 1/1
      rcp = ( i2 + gshift_i2 ) ! this is 1 at the axis - huh? if i2 = 1 this is 0
      rcm = ( i2 + gshift_i2 - 1 ) ! this must become zero on axis. - if i2 = 1 this is -1

      do i1 = 1, e_re%nx_(1)+1
        ! Real Part first

        ! E1
        e_re%f2(1, i1, i2) = e_re%f2(1, i1, i2) &
                          - dtif * jay_re%f2(1, i1, i2) &
                          + tmp * ( rcp * b_re%f2(3, i1, i2) - rcm * b_re%f2(3, i1, i2-1) )
        ! E2
        e_re%f2(2, i1, i2) = e_re%f2(2, i1, i2) &
                          - dtif * jay_re%f2(2, i1, i2) &
                          - dtdz * ( b_re%f2(3, i1, i2) - b_re%f2(3, i1-1, i2) )
        ! E3
        e_re%f2(3, i1, i2) = e_re%f2(3, i1, i2) &
                          - dtif * jay_re%f2(3, i1, i2) &
                          + dtdz * ( b_re%f2(2, i1, i2) - b_re%f2(2, i1-1, i2)) &
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
        e_im%f2(2, i1, i2) = e_im%f2(2, i1, i2) &
                          - dtif * jay_im%f2(2, i1, i2) &
                          - dtdz * ( b_im%f2(3, i1, i2) - b_im%f2(3, i1-1, i2) ) &
                          - dt_r_2 * mode * b_re%f2(1, i1, i2)

        ! E3
        e_im%f2(3, i1, i2) = e_im%f2(3, i1, i2) &
                        - dtif * jay_im%f2(3, i1, i2) &
                        + dtdz * ( b_im%f2(2, i1, i2) - b_im%f2(2, i1-1, i2)) &
                        - dtdr * ( b_im%f2(1, i1, i2) - b_im%f2(1, i1, i2-1))

      enddo ! i2
    enddo ! i1
    !$omp end parallel do

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
        do i1 = 1, e_re%nx_(1)+1
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
        !$omp parallel do
        do i1 = 1, e_re%nx_(1)+1
          ! E1 must be zero on the axis
          e_re%f2( 1, i1, 1 ) = - e_re%f2( 1, i1, 2 )
          e_re%f2( 1, i1, 0 ) = - e_re%f2( 1, i1, 3 )
          e_im%f2( 1, i1, 1 ) = - e_im%f2( 1, i1, 2 )
          e_im%f2( 1, i1, 0 ) = - e_im%f2( 1, i1, 3 )

          ! E2 can be nonzero on the axis
          ! Rather than solving for e_r(r=0) in the normal way as listed in the paper, we
          ! instead ensure that de_r/dr(r=0) = 0. Setting the derivative to zero reduces
          ! spurious field growth in e_r(r=0) due to the jay term.
          ! e_re%f2( 2, i1, 1 ) = e_re%f2( 2, i1, 1 ) &
          ! - dtif * jay_re%f2(2,i1,1) &
          ! - dtdz * ( b_re%f2(3, i1, 1) - b_re%f2(3, i1-1, 1) ) &
          ! + dtdr * b_im%f2(1, i1, 2)
          e_re%f2( 2, i1, 1 ) = f_t * e_re%f2( 2, i1, 2 ) - o_t * e_re%f2( 2, i1, 3 )
          e_re%f2( 2, i1, 0 ) = e_re%f2( 2, i1, 2 )
          ! e_im%f2( 2, i1, 1 ) = e_im%f2( 2, i1, 1 ) &
          ! - dtif * jay_im%f2(2,i1,1) &
          ! - dtdz * ( b_im%f2(3, i1, 1) - b_im%f2(3, i1-1, 1) ) &
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
          e_re%f2( 3, i1, 1 ) = 2.0_p_k_fld * e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 2 )
          e_re%f2( 3, i1, 0 ) = 2.0_p_k_fld * e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 3 )
          e_im%f2( 3, i1, 1 ) =-2.0_p_k_fld * e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 2 )
          e_im%f2( 3, i1, 0 ) =-2.0_p_k_fld * e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 3 )
        enddo
        !$omp end parallel do

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
        do i1 = 1, e_re%nx_(1)+1
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

end subroutine dedt_2d_cyl_modes_yee
!-----------------------------------------------------------------------------------------

subroutine test_stability_cyl_modes_yee( this, dt, dx )

    implicit none

    class( t_emf_solver_cyl_modes_yee ), intent(in) :: this
    real(p_double), intent(in) :: dt
    real(p_double), dimension(:), intent(in) :: dx

    real(p_double) :: cour
    integer :: i

    cour = 0.0
    do i = 1, p_x_dim
        cour = cour + 1.0/(dx(i))**2
    enddo
    ! Add another effective "cell" of size dr * pi / n_cyl_modes
    ! The extra term in the sum is then 2 / (dr * dtheta)
    if ( this%n_cyl_modes > 0 ) then
        cour = cour + 1.0_p_double / (dx(p_x_dim)*pi/(2.0_p_double*this%n_cyl_modes))**2
    else
        cour = cour + 1.0_p_double / (dx(p_x_dim)*pi)**2
    endif
    cour = sqrt(1.0/cour)

    if (dt > cour) then
        if ( mpi_node() == 0 ) then
           print *, '(*error*) Yee EMF solver stability condition violated, aborting'
           print *, 'dx   = ', dx(1:p_x_dim)
           print *, 'dt   = ', dt, ' max(dt) = ', cour
        endif
        call abort_program(p_err_invalid)
    endif

end subroutine test_stability_cyl_modes_yee

subroutine init_cyl_modes_yee( this, gix_pos, coords )

    implicit none

    class( t_emf_solver_cyl_modes_yee ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords

    ! Only handle cylindrical (r-z) coordinates solvers
    this % dedt_cyl_modes => dedt_2d_cyl_modes_yee
    this % dbdt_cyl_modes => dbdt_2d_cyl_modes_yee
    ! Radial position of lower cell
    this % gir_pos = gix_pos(2)

end subroutine init_cyl_modes_yee

subroutine advance_cyl_modes_yee( this, e_cyl_m, b_cyl_m, jay_cyl_m, dt, bnd_con )

    implicit none

    class( t_emf_solver_cyl_modes_yee ), intent(inout) :: this
    type( t_cyl_modes ), intent( inout ) :: e_cyl_m, b_cyl_m, jay_cyl_m
    class( t_emf_bound_cyl_modes ), intent(inout) :: bnd_con
    real(p_double), intent(in) :: dt

    real(p_double) :: dt_b, dt_e

    dt_b = dt / 2.0_p_double
    dt_e = dt






    ! Advance B half time step
    call this % dbdt_cyl_modes( b_cyl_m, e_cyl_m, dt_b )
    call bnd_con%update_boundary_b_cm( e_cyl_m, b_cyl_m, step = 1 )






    ! Advance E one full time step
    call this % dedt_cyl_modes( e_cyl_m, b_cyl_m, jay_cyl_m, dt_e )
    call bnd_con%update_boundary_e_cm( e_cyl_m, b_cyl_m )






    ! Advance B another half time step
    call this % dbdt_cyl_modes( b_cyl_m, e_cyl_m, dt_b )
    call bnd_con%update_boundary_b_cm( e_cyl_m, b_cyl_m, step = 2 )






end subroutine advance_cyl_modes_yee

end module m_emf_solver_cyl_modes_yee
