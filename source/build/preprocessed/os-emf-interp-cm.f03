# 1 "cyl_modes/os-emf-interp-cm.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-emf-interp-cm.f03" 2
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
# 2 "cyl_modes/os-emf-interp-cm.f03" 2
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
# 3 "cyl_modes/os-emf-interp-cm.f03" 2

module m_emf_interpolate_cm

use m_emf_define, only: t_emf
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_vdf_define, only: t_vdf

use m_parameters

implicit none

private


interface get_emf_cyl_modes
  module procedure get_emf_cyl_modes
end interface

public :: get_emf_cyl_modes


contains


!--------------------------------------------------------------------------------------------------
subroutine get_emf_cyl_modes( this, mode, bp_re, ep_re, bp_im, ep_im, &
                              ix, x, np, interpolation )
!--------------------------------------------------------------------------------------------------
! Interpolate fields at particle positions for cell based positions. The fields pointed to by
! this%e_part and this%b_part already include smoothed and/or external fields
!--------------------------------------------------------------------------------------------------

  implicit none

  class( t_emf_cyl_modes ), intent(in), target :: this
  integer, intent(in) :: mode
  real(p_k_part), dimension(:,:), intent( out ) :: bp_re, ep_re
  real(p_k_part), dimension(:,:), intent( out ) :: bp_im, ep_im
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np
  integer, intent(in) :: interpolation

  type( t_vdf ), pointer :: e_re, e_im
  type( t_vdf ), pointer :: b_re, b_im

  e_re => this%e_part_cyl_m%pf_re( mode )
  e_im => this%e_part_cyl_m%pf_im( mode )

  b_re => this%b_part_cyl_m%pf_re( mode )
  b_im => this%b_part_cyl_m%pf_im( mode )

  select case (interpolation)
  case(p_linear)
    call get_emf_cyl_m_s1( bp_re, ep_re, bp_im, ep_im, &
                           b_re, e_re, b_im, e_im, ix, x, np )

  case(p_quadratic)
    call get_emf_cyl_m_s2( bp_re, ep_re, bp_im, ep_im, &
                           b_re, e_re, b_im, e_im, ix, x, np )

  case(p_cubic)
    call get_emf_cyl_m_s3( bp_re, ep_re, bp_im, ep_im, &
                           b_re, e_re, b_im, e_im, ix, x, np )

  case(p_quartic)
    call get_emf_cyl_m_s4( bp_re, ep_re, bp_im, ep_im, &
                           b_re, e_re, b_im, e_im, ix, x, np )

  case default
    write(err_buf__,*) 'Not implemented yet';call err__("cyl_modes/os-emf-interp-cm.f03",73)
    call abort_program( p_err_notimplemented )
  end select


end subroutine get_emf_cyl_modes
!---------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine get_emf_cyl_m_s1( bp_re, ep_re, bp_im, ep_im, &
                           b_re, e_re, b_im, e_im, ix, x, np )

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp_re, ep_re
  real(p_k_part), dimension(:,:), intent( out ) :: bp_im, ep_im

  type( t_vdf ), intent(in) :: b_re, e_re
  type( t_vdf ), intent(in) :: b_im, e_im

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  integer :: i, j, ih, jh, l
  real(p_k_fld) :: dx1, dx2, dx1h, dx2h
  real(p_k_fld), dimension(0:1) :: w1, w1h, w2, w2h

  do l = 1, np

    i = ix(1,l)
    j = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    ih = i - signbit(dx1)
    jh = j - signbit(dx2)
    dx1h = dx1 - 0.5_p_k_fld + (i-ih)
    dx2h = dx2 - 0.5_p_k_fld + (j-jh)

    ! get spline weitghts for x and y
    w1(0) = 0.5 - dx1
    w1(1) = 0.5 + dx1

    w1h(0) = 0.5 - dx1h
    w1h(1) = 0.5 + dx1h

    w2(0) = 0.5 - dx2
    w2(1) = 0.5 + dx2

    w2h(0) = 0.5 - dx2h
    w2h(1) = 0.5 + dx2h

    ! Interpolate Fields
    ep_re(1,l) = ( e_re%f2(1,ih ,j ) * w1h(0) + &
                    e_re%f2(1,ih+1,j ) * w1h(1) ) * w2(0) + &
                  ( e_re%f2(1,ih ,j+1) * w1h(0) + &
                    e_re%f2(1,ih+1,j+1) * w1h(1) ) * w2(1)

    ep_re(2,l) = ( e_re%f2(2,i ,jh ) * w1(0) + &
                    e_re%f2(2,i+1,jh ) * w1(1) ) * w2h(0) + &
                  ( e_re%f2(2,i ,jh+1) * w1(0) + &
                    e_re%f2(2,i+1,jh+1) * w1(1) ) * w2h(1)

    ep_re(3,l) = ( e_re%f2(3,i ,j ) * w1(0) + &
                    e_re%f2(3,i+1,j ) * w1(1) ) * w2(0) + &
                  ( e_re%f2(3,i ,j+1) * w1(0) + &
                    e_re%f2(3,i+1,j+1) * w1(1) ) * w2(1)


    bp_re(1,l) = ( b_re%f2(1,i ,jh ) * w1(0) + &
                    b_re%f2(1,i+1,jh ) * w1(1) ) * w2h(0) + &
                  ( b_re%f2(1,i ,jh+1) * w1(0) + &
                    b_re%f2(1,i+1,jh+1) * w1(1) ) * w2h(1)

    bp_re(2,l) = ( b_re%f2(2,ih ,j ) * w1h(0) + &
                    b_re%f2(2,ih+1,j ) * w1h(1) ) * w2(0) + &
                  ( b_re%f2(2,ih ,j+1) * w1h(0) + &
                    b_re%f2(2,ih+1,j+1) * w1h(1) ) * w2(1)

    bp_re(3,l) = ( b_re%f2(3,ih ,jh ) * w1h(0) + &
                    b_re%f2(3,ih+1,jh ) * w1h(1) ) * w2h(0) + &
                  ( b_re%f2(3,ih ,jh+1) * w1h(0) + &
                    b_re%f2(3,ih+1,jh+1) * w1h(1) ) * w2h(1)



    ep_im(1,l) = ( e_im%f2(1,ih ,j ) * w1h(0) + &
                    e_im%f2(1,ih+1,j ) * w1h(1) ) * w2(0) + &
                  ( e_im%f2(1,ih ,j+1) * w1h(0) + &
                    e_im%f2(1,ih+1,j+1) * w1h(1) ) * w2(1)

    ep_im(2,l) = ( e_im%f2(2,i ,jh ) * w1(0) + &
                    e_im%f2(2,i+1,jh ) * w1(1) ) * w2h(0) + &
                  ( e_im%f2(2,i ,jh+1) * w1(0) + &
                    e_im%f2(2,i+1,jh+1) * w1(1) ) * w2h(1)

    ep_im(3,l) = ( e_im%f2(3,i ,j ) * w1(0) + &
                    e_im%f2(3,i+1,j ) * w1(1) ) * w2(0) + &
                  ( e_im%f2(3,i ,j+1) * w1(0) + &
                    e_im%f2(3,i+1,j+1) * w1(1) ) * w2(1)


    bp_im(1,l) = ( b_im%f2(1,i ,jh ) * w1(0) + &
                    b_im%f2(1,i+1,jh ) * w1(1) ) * w2h(0) + &
                  ( b_im%f2(1,i ,jh+1) * w1(0) + &
                    b_im%f2(1,i+1,jh+1) * w1(1) ) * w2h(1)

    bp_im(2,l) = ( b_im%f2(2,ih ,j ) * w1h(0) + &
                    b_im%f2(2,ih+1,j ) * w1h(1) ) * w2(0) + &
                  ( b_im%f2(2,ih ,j+1) * w1h(0) + &
                    b_im%f2(2,ih+1,j+1) * w1h(1) ) * w2(1)

    bp_im(3,l) = ( b_im%f2(3,ih ,jh ) * w1h(0) + &
                    b_im%f2(3,ih+1,jh ) * w1h(1) ) * w2h(0) + &
                  ( b_im%f2(3,ih ,jh+1) * w1h(0) + &
                    b_im%f2(3,ih+1,jh+1) * w1h(1) ) * w2h(1)


  enddo


end subroutine get_emf_cyl_m_s1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine get_emf_cyl_m_s2( bp_re, ep_re, bp_im, ep_im, &
                           b_re, e_re, b_im, e_im, ix, x, np )

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp_re, ep_re
  real(p_k_part), dimension(:,:), intent( out ) :: bp_im, ep_im

  type( t_vdf ), intent(in) :: b_re, e_re
  type( t_vdf ), intent(in) :: b_im, e_im

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np


  real(p_k_fld) :: dx1, dx2, dx1_h, dx2_h
  real(p_k_fld), dimension(3) :: w1, w1_h, w2, w2_h

  integer :: i, j, ih, jh, l

  do l = 1, np

    ! particle cell is the nearest grid point
    i = ix(1,l)
    dx1 = real( x(1,l), p_k_fld )

    ! Note that these weights are only valid for -0.5 <= dx1 <= 0.5
    w1(1) = 0.5_p_k_fld*(0.5_p_k_fld - dx1)**2
    w1(2) = 0.75_p_k_fld - dx1**2
    w1(3) = 0.5_p_k_fld*(0.5_p_k_fld + dx1)**2

    ! safe
    ih = i - signbit(dx1)

    dx1_h = dx1 - 0.5_p_k_fld + (i-ih)

    w1_h(1) = 0.5_p_k_fld*(0.5_p_k_fld - dx1_h)**2
    w1_h(2) = 0.75_p_k_fld - dx1_h**2
    w1_h(3) = 0.5_p_k_fld*(0.5_p_k_fld + dx1_h)**2

    ! get interpolation indexes and weights for y
    j = ix(2,l)
    dx2 = real( x(2,l), p_k_fld )

    w2(1) = 0.5_p_k_fld*(0.5_p_k_fld - dx2)**2
    w2(2) = 0.75_p_k_fld - dx2**2
    w2(3) = 0.5_p_k_fld*(0.5_p_k_fld + dx2)**2

    jh = j - signbit(dx2)

    dx2_h = dx2 - 0.5_p_k_fld + (j-jh)

    w2_h(1) = 0.5_p_k_fld*(0.5_p_k_fld - dx2_h)**2
    w2_h(2) = 0.75_p_k_fld - dx2_h**2
    w2_h(3) = 0.5_p_k_fld*(0.5_p_k_fld + dx2_h)**2

    ! interpolate fields
    ep_re(1,l) = ( e_re%f2(1,ih-1,j-1) * w1_h(1) + &
                    e_re%f2(1,ih ,j-1) * w1_h(2) + &
                    e_re%f2(1,ih+1,j-1) * w1_h(3) ) * w2(1) + &
                  ( e_re%f2(1,ih-1,j ) * w1_h(1) + &
                    e_re%f2(1,ih ,j ) * w1_h(2) + &
                    e_re%f2(1,ih+1,j ) * w1_h(3) ) * w2(2) + &
                  ( e_re%f2(1,ih-1,j+1) * w1_h(1) + &
                    e_re%f2(1,ih ,j+1) * w1_h(2) + &
                    e_re%f2(1,ih+1,j+1) * w1_h(3) ) * w2(3)

    ep_re(2,l) = ( e_re%f2(2,i-1,jh-1) * w1(1) + &
                    e_re%f2(2,i ,jh-1) * w1(2) + &
                    e_re%f2(2,i+1,jh-1) * w1(3) ) * w2_h(1) + &
                  ( e_re%f2(2,i-1,jh ) * w1(1) + &
                    e_re%f2(2,i ,jh ) * w1(2) + &
                    e_re%f2(2,i+1,jh ) * w1(3) ) * w2_h(2) + &
                  ( e_re%f2(2,i-1,jh+1) * w1(1) + &
                    e_re%f2(2,i ,jh+1) * w1(2) + &
                    e_re%f2(2,i+1,jh+1) * w1(3) ) * w2_h(3)

    ep_re(3,l) = ( e_re%f2(3,i-1,j-1) * w1(1) + &
                    e_re%f2(3,i ,j-1) * w1(2) + &
                    e_re%f2(3,i+1,j-1) * w1(3) ) * w2(1) + &
                  ( e_re%f2(3,i-1,j ) * w1(1) + &
                    e_re%f2(3,i ,j ) * w1(2) + &
                    e_re%f2(3,i+1,j ) * w1(3) ) * w2(2) + &
                  ( e_re%f2(3,i-1,j+1) * w1(1) + &
                    e_re%f2(3,i ,j+1) * w1(2) + &
                    e_re%f2(3,i+1,j+1) * w1(3) ) * w2(3)

    bp_re(1,l) = ( b_re%f2(1,i-1,jh-1) * w1(1) + &
                    b_re%f2(1,i ,jh-1) * w1(2) + &
                    b_re%f2(1,i+1,jh-1) * w1(3) ) * w2_h(1) + &
                  ( b_re%f2(1,i-1,jh ) * w1(1) + &
                    b_re%f2(1,i ,jh ) * w1(2) + &
                    b_re%f2(1,i+1,jh ) * w1(3) ) * w2_h(2) + &
                  ( b_re%f2(1,i-1,jh+1) * w1(1) + &
                    b_re%f2(1,i ,jh+1) * w1(2) + &
                    b_re%f2(1,i+1,jh+1) * w1(3) ) * w2_h(3)

    bp_re(2,l) = ( b_re%f2(2,ih-1,j-1) * w1_h(1) + &
                    b_re%f2(2,ih ,j-1) * w1_h(2) + &
                    b_re%f2(2,ih+1,j-1) * w1_h(3) ) * w2(1) + &
                  ( b_re%f2(2,ih-1,j ) * w1_h(1) + &
                    b_re%f2(2,ih ,j ) * w1_h(2) + &
                    b_re%f2(2,ih+1,j ) * w1_h(3) ) * w2(2) + &
                  ( b_re%f2(2,ih-1,j+1) * w1_h(1) + &
                    b_re%f2(2,ih ,j+1) * w1_h(2) + &
                    b_re%f2(2,ih+1,j+1) * w1_h(3) ) * w2(3)

    bp_re(3,l) = ( b_re%f2(3,ih-1,jh-1) * w1_h(1) + &
                    b_re%f2(3,ih ,jh-1) * w1_h(2) + &
                    b_re%f2(3,ih+1,jh-1) * w1_h(3) ) * w2_h(1) + &
                  ( b_re%f2(3,ih-1,jh ) * w1_h(1) + &
                    b_re%f2(3,ih ,jh ) * w1_h(2) + &
                    b_re%f2(3,ih+1,jh ) * w1_h(3) ) * w2_h(2) + &
                  ( b_re%f2(3,ih-1,jh+1) * w1_h(1) + &
                    b_re%f2(3,ih ,jh+1) * w1_h(2) + &
                    b_re%f2(3,ih+1,jh+1) * w1_h(3) ) * w2_h(3)

    ep_im(1,l) = ( e_im%f2(1,ih-1,j-1) * w1_h(1) + &
                    e_im%f2(1,ih ,j-1) * w1_h(2) + &
                    e_im%f2(1,ih+1,j-1) * w1_h(3) ) * w2(1) + &
                  ( e_im%f2(1,ih-1,j ) * w1_h(1) + &
                    e_im%f2(1,ih ,j ) * w1_h(2) + &
                    e_im%f2(1,ih+1,j ) * w1_h(3) ) * w2(2) + &
                  ( e_im%f2(1,ih-1,j+1) * w1_h(1) + &
                    e_im%f2(1,ih ,j+1) * w1_h(2) + &
                    e_im%f2(1,ih+1,j+1) * w1_h(3) ) * w2(3)

    ep_im(2,l) = ( e_im%f2(2,i-1,jh-1) * w1(1) + &
                    e_im%f2(2,i ,jh-1) * w1(2) + &
                    e_im%f2(2,i+1,jh-1) * w1(3) ) * w2_h(1) + &
                  ( e_im%f2(2,i-1,jh ) * w1(1) + &
                    e_im%f2(2,i ,jh ) * w1(2) + &
                    e_im%f2(2,i+1,jh ) * w1(3) ) * w2_h(2) + &
                  ( e_im%f2(2,i-1,jh+1) * w1(1) + &
                    e_im%f2(2,i ,jh+1) * w1(2) + &
                    e_im%f2(2,i+1,jh+1) * w1(3) ) * w2_h(3)

    ep_im(3,l) = ( e_im%f2(3,i-1,j-1) * w1(1) + &
                    e_im%f2(3,i ,j-1) * w1(2) + &
                    e_im%f2(3,i+1,j-1) * w1(3) ) * w2(1) + &
                  ( e_im%f2(3,i-1,j ) * w1(1) + &
                    e_im%f2(3,i ,j ) * w1(2) + &
                    e_im%f2(3,i+1,j ) * w1(3) ) * w2(2) + &
                  ( e_im%f2(3,i-1,j+1) * w1(1) + &
                    e_im%f2(3,i ,j+1) * w1(2) + &
                    e_im%f2(3,i+1,j+1) * w1(3) ) * w2(3)

    bp_im(1,l) = ( b_im%f2(1,i-1,jh-1) * w1(1) + &
                    b_im%f2(1,i ,jh-1) * w1(2) + &
                    b_im%f2(1,i+1,jh-1) * w1(3) ) * w2_h(1) + &
                  ( b_im%f2(1,i-1,jh ) * w1(1) + &
                    b_im%f2(1,i ,jh ) * w1(2) + &
                    b_im%f2(1,i+1,jh ) * w1(3) ) * w2_h(2) + &
                  ( b_im%f2(1,i-1,jh+1) * w1(1) + &
                    b_im%f2(1,i ,jh+1) * w1(2) + &
                    b_im%f2(1,i+1,jh+1) * w1(3) ) * w2_h(3)

    bp_im(2,l) = ( b_im%f2(2,ih-1,j-1) * w1_h(1) + &
                    b_im%f2(2,ih ,j-1) * w1_h(2) + &
                    b_im%f2(2,ih+1,j-1) * w1_h(3) ) * w2(1) + &
                  ( b_im%f2(2,ih-1,j ) * w1_h(1) + &
                    b_im%f2(2,ih ,j ) * w1_h(2) + &
                    b_im%f2(2,ih+1,j ) * w1_h(3) ) * w2(2) + &
                  ( b_im%f2(2,ih-1,j+1) * w1_h(1) + &
                    b_im%f2(2,ih ,j+1) * w1_h(2) + &
                    b_im%f2(2,ih+1,j+1) * w1_h(3) ) * w2(3)

    bp_im(3,l) = ( b_im%f2(3,ih-1,jh-1) * w1_h(1) + &
                    b_im%f2(3,ih ,jh-1) * w1_h(2) + &
                    b_im%f2(3,ih+1,jh-1) * w1_h(3) ) * w2_h(1) + &
                  ( b_im%f2(3,ih-1,jh ) * w1_h(1) + &
                    b_im%f2(3,ih ,jh ) * w1_h(2) + &
                    b_im%f2(3,ih+1,jh ) * w1_h(3) ) * w2_h(2) + &
                  ( b_im%f2(3,ih-1,jh+1) * w1_h(1) + &
                    b_im%f2(3,ih ,jh+1) * w1_h(2) + &
                    b_im%f2(3,ih+1,jh+1) * w1_h(3) ) * w2_h(3)

  enddo

end subroutine get_emf_cyl_m_s2
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine get_emf_cyl_m_s3( bp_re, ep_re, bp_im, ep_im, &
                           b_re, e_re, b_im, e_im, ix, x, np )

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp_re, ep_re
  real(p_k_part), dimension(:,:), intent( out ) :: bp_im, ep_im

  type( t_vdf ), intent(in) :: b_re, e_re
  type( t_vdf ), intent(in) :: b_im, e_im

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np


  integer :: i, j, ih, jh, l
  real(p_k_fld) :: dx1, dx2, dx1h, dx2h
  real(p_k_fld), dimension(-1:2) :: w1, w1h, w2, w2h

  do l = 1, np

    ! order 3 field interpolation
    ! generated automatically by z-2.1

    i = ix(1,l)
    j = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    ih = i - signbit(dx1)
    jh = j - signbit(dx2)
    dx1h = dx1 - 0.5_p_k_fld + (i-ih)
    dx2h = dx2 - 0.5_p_k_fld + (j-jh)

    ! get spline weitghts for x and y
    w1(-1) = -(-0.5 + dx1)**3/6.
    w1(0) = (4 - 6*(0.5 + dx1)**2 + 3*(0.5 + dx1)**3)/6.
    w1(1) = (23 + 30*dx1 - 12*dx1**2 - 24*dx1**3)/48.
    w1(2) = (0.5 + dx1)**3/6.

    w1h(-1) = -(-0.5 + dx1h)**3/6.
    w1h(0) = (4 - 6*(0.5 + dx1h)**2 + 3*(0.5 + dx1h)**3)/6.
    w1h(1) = (23 + 30*dx1h - 12*dx1h**2 - 24*dx1h**3)/48.
    w1h(2) = (0.5 + dx1h)**3/6.

    w2(-1) = -(-0.5 + dx2)**3/6.
    w2(0) = (4 - 6*(0.5 + dx2)**2 + 3*(0.5 + dx2)**3)/6.
    w2(1) = (23 + 30*dx2 - 12*dx2**2 - 24*dx2**3)/48.
    w2(2) = (0.5 + dx2)**3/6.

    w2h(-1) = -(-0.5 + dx2h)**3/6.
    w2h(0) = (4 - 6*(0.5 + dx2h)**2 + 3*(0.5 + dx2h)**3)/6.
    w2h(1) = (23 + 30*dx2h - 12*dx2h**2 - 24*dx2h**3)/48.
    w2h(2) = (0.5 + dx2h)**3/6.

    ! Interpolate Fields
    ep_re(1,l) = ( e_re%f2(1,ih-1,j-1) * w1h(-1) + &
                    e_re%f2(1,ih ,j-1) * w1h(0) + &
                    e_re%f2(1,ih+1,j-1) * w1h(1) + &
                    e_re%f2(1,ih+2,j-1) * w1h(2) ) * w2(-1) + &
                  ( e_re%f2(1,ih-1,j ) * w1h(-1) + &
                    e_re%f2(1,ih ,j ) * w1h(0) + &
                    e_re%f2(1,ih+1,j ) * w1h(1) + &
                    e_re%f2(1,ih+2,j ) * w1h(2) ) * w2(0) + &
                  ( e_re%f2(1,ih-1,j+1) * w1h(-1) + &
                    e_re%f2(1,ih ,j+1) * w1h(0) + &
                    e_re%f2(1,ih+1,j+1) * w1h(1) + &
                    e_re%f2(1,ih+2,j+1) * w1h(2) ) * w2(1) + &
                  ( e_re%f2(1,ih-1,j+2) * w1h(-1) + &
                    e_re%f2(1,ih ,j+2) * w1h(0) + &
                    e_re%f2(1,ih+1,j+2) * w1h(1) + &
                    e_re%f2(1,ih+2,j+2) * w1h(2) ) * w2(2)

    ep_re(2,l) = ( e_re%f2(2,i-1,jh-1) * w1(-1) + &
                    e_re%f2(2,i ,jh-1) * w1(0) + &
                    e_re%f2(2,i+1,jh-1) * w1(1) + &
                    e_re%f2(2,i+2,jh-1) * w1(2) ) * w2h(-1) + &
                  ( e_re%f2(2,i-1,jh ) * w1(-1) + &
                    e_re%f2(2,i ,jh ) * w1(0) + &
                    e_re%f2(2,i+1,jh ) * w1(1) + &
                    e_re%f2(2,i+2,jh ) * w1(2) ) * w2h(0) + &
                  ( e_re%f2(2,i-1,jh+1) * w1(-1) + &
                    e_re%f2(2,i ,jh+1) * w1(0) + &
                    e_re%f2(2,i+1,jh+1) * w1(1) + &
                    e_re%f2(2,i+2,jh+1) * w1(2) ) * w2h(1) + &
                  ( e_re%f2(2,i-1,jh+2) * w1(-1) + &
                    e_re%f2(2,i ,jh+2) * w1(0) + &
                    e_re%f2(2,i+1,jh+2) * w1(1) + &
                    e_re%f2(2,i+2,jh+2) * w1(2) ) * w2h(2)

    ep_re(3,l) = ( e_re%f2(3,i-1,j-1) * w1(-1) + &
                    e_re%f2(3,i ,j-1) * w1(0) + &
                    e_re%f2(3,i+1,j-1) * w1(1) + &
                    e_re%f2(3,i+2,j-1) * w1(2) ) * w2(-1) + &
                  ( e_re%f2(3,i-1,j ) * w1(-1) + &
                    e_re%f2(3,i ,j ) * w1(0) + &
                    e_re%f2(3,i+1,j ) * w1(1) + &
                    e_re%f2(3,i+2,j ) * w1(2) ) * w2(0) + &
                  ( e_re%f2(3,i-1,j+1) * w1(-1) + &
                    e_re%f2(3,i ,j+1) * w1(0) + &
                    e_re%f2(3,i+1,j+1) * w1(1) + &
                    e_re%f2(3,i+2,j+1) * w1(2) ) * w2(1) + &
                  ( e_re%f2(3,i-1,j+2) * w1(-1) + &
                    e_re%f2(3,i ,j+2) * w1(0) + &
                    e_re%f2(3,i+1,j+2) * w1(1) + &
                    e_re%f2(3,i+2,j+2) * w1(2) ) * w2(2)

    bp_re(1,l) = ( b_re%f2(1,i-1,jh-1) * w1(-1) + &
                    b_re%f2(1,i ,jh-1) * w1(0) + &
                    b_re%f2(1,i+1,jh-1) * w1(1) + &
                    b_re%f2(1,i+2,jh-1) * w1(2) ) * w2h(-1) + &
                  ( b_re%f2(1,i-1,jh ) * w1(-1) + &
                    b_re%f2(1,i ,jh ) * w1(0) + &
                    b_re%f2(1,i+1,jh ) * w1(1) + &
                    b_re%f2(1,i+2,jh ) * w1(2) ) * w2h(0) + &
                  ( b_re%f2(1,i-1,jh+1) * w1(-1) + &
                    b_re%f2(1,i ,jh+1) * w1(0) + &
                    b_re%f2(1,i+1,jh+1) * w1(1) + &
                    b_re%f2(1,i+2,jh+1) * w1(2) ) * w2h(1) + &
                  ( b_re%f2(1,i-1,jh+2) * w1(-1) + &
                    b_re%f2(1,i ,jh+2) * w1(0) + &
                    b_re%f2(1,i+1,jh+2) * w1(1) + &
                    b_re%f2(1,i+2,jh+2) * w1(2) ) * w2h(2)

    bp_re(2,l) = ( b_re%f2(2,ih-1,j-1) * w1h(-1) + &
                    b_re%f2(2,ih ,j-1) * w1h(0) + &
                    b_re%f2(2,ih+1,j-1) * w1h(1) + &
                    b_re%f2(2,ih+2,j-1) * w1h(2) ) * w2(-1) + &
                  ( b_re%f2(2,ih-1,j ) * w1h(-1) + &
                    b_re%f2(2,ih ,j ) * w1h(0) + &
                    b_re%f2(2,ih+1,j ) * w1h(1) + &
                    b_re%f2(2,ih+2,j ) * w1h(2) ) * w2(0) + &
                  ( b_re%f2(2,ih-1,j+1) * w1h(-1) + &
                    b_re%f2(2,ih ,j+1) * w1h(0) + &
                    b_re%f2(2,ih+1,j+1) * w1h(1) + &
                    b_re%f2(2,ih+2,j+1) * w1h(2) ) * w2(1) + &
                  ( b_re%f2(2,ih-1,j+2) * w1h(-1) + &
                    b_re%f2(2,ih ,j+2) * w1h(0) + &
                    b_re%f2(2,ih+1,j+2) * w1h(1) + &
                    b_re%f2(2,ih+2,j+2) * w1h(2) ) * w2(2)

    bp_re(3,l) = ( b_re%f2(3,ih-1,jh-1) * w1h(-1) + &
                    b_re%f2(3,ih ,jh-1) * w1h(0) + &
                    b_re%f2(3,ih+1,jh-1) * w1h(1) + &
                    b_re%f2(3,ih+2,jh-1) * w1h(2) ) * w2h(-1) + &
                  ( b_re%f2(3,ih-1,jh ) * w1h(-1) + &
                    b_re%f2(3,ih ,jh ) * w1h(0) + &
                    b_re%f2(3,ih+1,jh ) * w1h(1) + &
                    b_re%f2(3,ih+2,jh ) * w1h(2) ) * w2h(0) + &
                  ( b_re%f2(3,ih-1,jh+1) * w1h(-1) + &
                    b_re%f2(3,ih ,jh+1) * w1h(0) + &
                    b_re%f2(3,ih+1,jh+1) * w1h(1) + &
                    b_re%f2(3,ih+2,jh+1) * w1h(2) ) * w2h(1) + &
                  ( b_re%f2(3,ih-1,jh+2) * w1h(-1) + &
                    b_re%f2(3,ih ,jh+2) * w1h(0) + &
                    b_re%f2(3,ih+1,jh+2) * w1h(1) + &
                    b_re%f2(3,ih+2,jh+2) * w1h(2) ) * w2h(2)

    ep_im(1,l) = ( e_im%f2(1,ih-1,j-1) * w1h(-1) + &
                    e_im%f2(1,ih ,j-1) * w1h(0) + &
                    e_im%f2(1,ih+1,j-1) * w1h(1) + &
                    e_im%f2(1,ih+2,j-1) * w1h(2) ) * w2(-1) + &
                  ( e_im%f2(1,ih-1,j ) * w1h(-1) + &
                    e_im%f2(1,ih ,j ) * w1h(0) + &
                    e_im%f2(1,ih+1,j ) * w1h(1) + &
                    e_im%f2(1,ih+2,j ) * w1h(2) ) * w2(0) + &
                  ( e_im%f2(1,ih-1,j+1) * w1h(-1) + &
                    e_im%f2(1,ih ,j+1) * w1h(0) + &
                    e_im%f2(1,ih+1,j+1) * w1h(1) + &
                    e_im%f2(1,ih+2,j+1) * w1h(2) ) * w2(1) + &
                  ( e_im%f2(1,ih-1,j+2) * w1h(-1) + &
                    e_im%f2(1,ih ,j+2) * w1h(0) + &
                    e_im%f2(1,ih+1,j+2) * w1h(1) + &
                    e_im%f2(1,ih+2,j+2) * w1h(2) ) * w2(2)

    ep_im(2,l) = ( e_im%f2(2,i-1,jh-1) * w1(-1) + &
                    e_im%f2(2,i ,jh-1) * w1(0) + &
                    e_im%f2(2,i+1,jh-1) * w1(1) + &
                    e_im%f2(2,i+2,jh-1) * w1(2) ) * w2h(-1) + &
                  ( e_im%f2(2,i-1,jh ) * w1(-1) + &
                    e_im%f2(2,i ,jh ) * w1(0) + &
                    e_im%f2(2,i+1,jh ) * w1(1) + &
                    e_im%f2(2,i+2,jh ) * w1(2) ) * w2h(0) + &
                  ( e_im%f2(2,i-1,jh+1) * w1(-1) + &
                    e_im%f2(2,i ,jh+1) * w1(0) + &
                    e_im%f2(2,i+1,jh+1) * w1(1) + &
                    e_im%f2(2,i+2,jh+1) * w1(2) ) * w2h(1) + &
                  ( e_im%f2(2,i-1,jh+2) * w1(-1) + &
                    e_im%f2(2,i ,jh+2) * w1(0) + &
                    e_im%f2(2,i+1,jh+2) * w1(1) + &
                    e_im%f2(2,i+2,jh+2) * w1(2) ) * w2h(2)

    ep_im(3,l) = ( e_im%f2(3,i-1,j-1) * w1(-1) + &
                    e_im%f2(3,i ,j-1) * w1(0) + &
                    e_im%f2(3,i+1,j-1) * w1(1) + &
                    e_im%f2(3,i+2,j-1) * w1(2) ) * w2(-1) + &
                  ( e_im%f2(3,i-1,j ) * w1(-1) + &
                    e_im%f2(3,i ,j ) * w1(0) + &
                    e_im%f2(3,i+1,j ) * w1(1) + &
                    e_im%f2(3,i+2,j ) * w1(2) ) * w2(0) + &
                  ( e_im%f2(3,i-1,j+1) * w1(-1) + &
                    e_im%f2(3,i ,j+1) * w1(0) + &
                    e_im%f2(3,i+1,j+1) * w1(1) + &
                    e_im%f2(3,i+2,j+1) * w1(2) ) * w2(1) + &
                  ( e_im%f2(3,i-1,j+2) * w1(-1) + &
                    e_im%f2(3,i ,j+2) * w1(0) + &
                    e_im%f2(3,i+1,j+2) * w1(1) + &
                    e_im%f2(3,i+2,j+2) * w1(2) ) * w2(2)

    bp_im(1,l) = ( b_im%f2(1,i-1,jh-1) * w1(-1) + &
                    b_im%f2(1,i ,jh-1) * w1(0) + &
                    b_im%f2(1,i+1,jh-1) * w1(1) + &
                    b_im%f2(1,i+2,jh-1) * w1(2) ) * w2h(-1) + &
                  ( b_im%f2(1,i-1,jh ) * w1(-1) + &
                    b_im%f2(1,i ,jh ) * w1(0) + &
                    b_im%f2(1,i+1,jh ) * w1(1) + &
                    b_im%f2(1,i+2,jh ) * w1(2) ) * w2h(0) + &
                  ( b_im%f2(1,i-1,jh+1) * w1(-1) + &
                    b_im%f2(1,i ,jh+1) * w1(0) + &
                    b_im%f2(1,i+1,jh+1) * w1(1) + &
                    b_im%f2(1,i+2,jh+1) * w1(2) ) * w2h(1) + &
                  ( b_im%f2(1,i-1,jh+2) * w1(-1) + &
                    b_im%f2(1,i ,jh+2) * w1(0) + &
                    b_im%f2(1,i+1,jh+2) * w1(1) + &
                    b_im%f2(1,i+2,jh+2) * w1(2) ) * w2h(2)

    bp_im(2,l) = ( b_im%f2(2,ih-1,j-1) * w1h(-1) + &
                    b_im%f2(2,ih ,j-1) * w1h(0) + &
                    b_im%f2(2,ih+1,j-1) * w1h(1) + &
                    b_im%f2(2,ih+2,j-1) * w1h(2) ) * w2(-1) + &
                  ( b_im%f2(2,ih-1,j ) * w1h(-1) + &
                    b_im%f2(2,ih ,j ) * w1h(0) + &
                    b_im%f2(2,ih+1,j ) * w1h(1) + &
                    b_im%f2(2,ih+2,j ) * w1h(2) ) * w2(0) + &
                  ( b_im%f2(2,ih-1,j+1) * w1h(-1) + &
                    b_im%f2(2,ih ,j+1) * w1h(0) + &
                    b_im%f2(2,ih+1,j+1) * w1h(1) + &
                    b_im%f2(2,ih+2,j+1) * w1h(2) ) * w2(1) + &
                  ( b_im%f2(2,ih-1,j+2) * w1h(-1) + &
                    b_im%f2(2,ih ,j+2) * w1h(0) + &
                    b_im%f2(2,ih+1,j+2) * w1h(1) + &
                    b_im%f2(2,ih+2,j+2) * w1h(2) ) * w2(2)

    bp_im(3,l) = ( b_im%f2(3,ih-1,jh-1) * w1h(-1) + &
                    b_im%f2(3,ih ,jh-1) * w1h(0) + &
                    b_im%f2(3,ih+1,jh-1) * w1h(1) + &
                    b_im%f2(3,ih+2,jh-1) * w1h(2) ) * w2h(-1) + &
                  ( b_im%f2(3,ih-1,jh ) * w1h(-1) + &
                    b_im%f2(3,ih ,jh ) * w1h(0) + &
                    b_im%f2(3,ih+1,jh ) * w1h(1) + &
                    b_im%f2(3,ih+2,jh ) * w1h(2) ) * w2h(0) + &
                  ( b_im%f2(3,ih-1,jh+1) * w1h(-1) + &
                    b_im%f2(3,ih ,jh+1) * w1h(0) + &
                    b_im%f2(3,ih+1,jh+1) * w1h(1) + &
                    b_im%f2(3,ih+2,jh+1) * w1h(2) ) * w2h(1) + &
                  ( b_im%f2(3,ih-1,jh+2) * w1h(-1) + &
                    b_im%f2(3,ih ,jh+2) * w1h(0) + &
                    b_im%f2(3,ih+1,jh+2) * w1h(1) + &
                    b_im%f2(3,ih+2,jh+2) * w1h(2) ) * w2h(2)

    ! end of automatic code
  enddo


end subroutine get_emf_cyl_m_s3
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine get_emf_cyl_m_s4( bp_re, ep_re, bp_im, ep_im, &
                           b_re, e_re, b_im, e_im, ix, x, np )

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp_re, ep_re
  real(p_k_part), dimension(:,:), intent( out ) :: bp_im, ep_im

  type( t_vdf ), intent(in) :: b_re, e_re
  type( t_vdf ), intent(in) :: b_im, e_im

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np


  integer :: i, j, ih, jh, l
  real(p_k_fld) :: dx1, dx2, dx1h, dx2h
  real(p_k_fld), dimension(-2:2) :: w1, w1h, w2, w2h

  do l = 1, np

    ! order 4 field interpolation
    ! generated automatically by z-2.0

    i = ix(1,l)
    j = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    ih = i - signbit(dx1)
    jh = j - signbit(dx2)
    dx1h = dx1 - 0.5_p_k_fld + (i-ih)
    dx2h = dx2 - 0.5_p_k_fld + (j-jh)

    ! get spline weitghts for x and y
    w1(-2) = (1 - 2*dx1)**4/384.
    w1(-1) = (19 - 44*dx1 + 24*dx1**2 + 16*dx1**3 - 16*dx1**4)/96.
    w1(0) = 0.5989583333333334 - (5*dx1**2)/8. + dx1**4/4.
    w1(1) = (19 + 44*dx1 + 24*dx1**2 - 16*dx1**3 - 16*dx1**4)/96.
    w1(2) = (1 + 2*dx1)**4/384.

    w1h(-2) = (1 - 2*dx1h)**4/384.
    w1h(-1) = (19 - 44*dx1h + 24*dx1h**2 + 16*dx1h**3 - 16*dx1h**4)/96.
    w1h(0) = 0.5989583333333334 - (5*dx1h**2)/8. + dx1h**4/4.
    w1h(1) = (19 + 44*dx1h + 24*dx1h**2 - 16*dx1h**3 - 16*dx1h**4)/96.
    w1h(2) = (1 + 2*dx1h)**4/384.

    w2(-2) = (1 - 2*dx2)**4/384.
    w2(-1) = (19 - 44*dx2 + 24*dx2**2 + 16*dx2**3 - 16*dx2**4)/96.
    w2(0) = 0.5989583333333334 - (5*dx2**2)/8. + dx2**4/4.
    w2(1) = (19 + 44*dx2 + 24*dx2**2 - 16*dx2**3 - 16*dx2**4)/96.
    w2(2) = (1 + 2*dx2)**4/384.

    w2h(-2) = (1 - 2*dx2h)**4/384.
    w2h(-1) = (19 - 44*dx2h + 24*dx2h**2 + 16*dx2h**3 - 16*dx2h**4)/96.
    w2h(0) = 0.5989583333333334 - (5*dx2h**2)/8. + dx2h**4/4.
    w2h(1) = (19 + 44*dx2h + 24*dx2h**2 - 16*dx2h**3 - 16*dx2h**4)/96.
    w2h(2) = (1 + 2*dx2h)**4/384.

    ! Interpolate Fields
    ep_re(1,l) = ( e_re%f2(1,ih-2,j-2) * w1h(-2) + &
                    e_re%f2(1,ih-1,j-2) * w1h(-1) + &
                    e_re%f2(1,ih ,j-2) * w1h(0) + &
                    e_re%f2(1,ih+1,j-2) * w1h(1) + &
                    e_re%f2(1,ih+2,j-2) * w1h(2) ) * w2(-2) + &
                  ( e_re%f2(1,ih-2,j-1) * w1h(-2) + &
                    e_re%f2(1,ih-1,j-1) * w1h(-1) + &
                    e_re%f2(1,ih ,j-1) * w1h(0) + &
                    e_re%f2(1,ih+1,j-1) * w1h(1) + &
                    e_re%f2(1,ih+2,j-1) * w1h(2) ) * w2(-1) + &
                  ( e_re%f2(1,ih-2,j ) * w1h(-2) + &
                    e_re%f2(1,ih-1,j ) * w1h(-1) + &
                    e_re%f2(1,ih ,j ) * w1h(0) + &
                    e_re%f2(1,ih+1,j ) * w1h(1) + &
                    e_re%f2(1,ih+2,j ) * w1h(2) ) * w2(0) + &
                  ( e_re%f2(1,ih-2,j+1) * w1h(-2) + &
                    e_re%f2(1,ih-1,j+1) * w1h(-1) + &
                    e_re%f2(1,ih ,j+1) * w1h(0) + &
                    e_re%f2(1,ih+1,j+1) * w1h(1) + &
                    e_re%f2(1,ih+2,j+1) * w1h(2) ) * w2(1) + &
                  ( e_re%f2(1,ih-2,j+2) * w1h(-2) + &
                    e_re%f2(1,ih-1,j+2) * w1h(-1) + &
                    e_re%f2(1,ih ,j+2) * w1h(0) + &
                    e_re%f2(1,ih+1,j+2) * w1h(1) + &
                    e_re%f2(1,ih+2,j+2) * w1h(2) ) * w2(2)

    ep_re(2,l) = ( e_re%f2(2,i-2,jh-2) * w1(-2) + &
                    e_re%f2(2,i-1,jh-2) * w1(-1) + &
                    e_re%f2(2,i ,jh-2) * w1(0) + &
                    e_re%f2(2,i+1,jh-2) * w1(1) + &
                    e_re%f2(2,i+2,jh-2) * w1(2) ) * w2h(-2) + &
                  ( e_re%f2(2,i-2,jh-1) * w1(-2) + &
                    e_re%f2(2,i-1,jh-1) * w1(-1) + &
                    e_re%f2(2,i ,jh-1) * w1(0) + &
                    e_re%f2(2,i+1,jh-1) * w1(1) + &
                    e_re%f2(2,i+2,jh-1) * w1(2) ) * w2h(-1) + &
                  ( e_re%f2(2,i-2,jh ) * w1(-2) + &
                    e_re%f2(2,i-1,jh ) * w1(-1) + &
                    e_re%f2(2,i ,jh ) * w1(0) + &
                    e_re%f2(2,i+1,jh ) * w1(1) + &
                    e_re%f2(2,i+2,jh ) * w1(2) ) * w2h(0) + &
                  ( e_re%f2(2,i-2,jh+1) * w1(-2) + &
                    e_re%f2(2,i-1,jh+1) * w1(-1) + &
                    e_re%f2(2,i ,jh+1) * w1(0) + &
                    e_re%f2(2,i+1,jh+1) * w1(1) + &
                    e_re%f2(2,i+2,jh+1) * w1(2) ) * w2h(1) + &
                  ( e_re%f2(2,i-2,jh+2) * w1(-2) + &
                    e_re%f2(2,i-1,jh+2) * w1(-1) + &
                    e_re%f2(2,i ,jh+2) * w1(0) + &
                    e_re%f2(2,i+1,jh+2) * w1(1) + &
                    e_re%f2(2,i+2,jh+2) * w1(2) ) * w2h(2)

    ep_re(3,l) = ( e_re%f2(3,i-2,j-2) * w1(-2) + &
                    e_re%f2(3,i-1,j-2) * w1(-1) + &
                    e_re%f2(3,i ,j-2) * w1(0) + &
                    e_re%f2(3,i+1,j-2) * w1(1) + &
                    e_re%f2(3,i+2,j-2) * w1(2) ) * w2(-2) + &
                  ( e_re%f2(3,i-2,j-1) * w1(-2) + &
                    e_re%f2(3,i-1,j-1) * w1(-1) + &
                    e_re%f2(3,i ,j-1) * w1(0) + &
                    e_re%f2(3,i+1,j-1) * w1(1) + &
                    e_re%f2(3,i+2,j-1) * w1(2) ) * w2(-1) + &
                  ( e_re%f2(3,i-2,j ) * w1(-2) + &
                    e_re%f2(3,i-1,j ) * w1(-1) + &
                    e_re%f2(3,i ,j ) * w1(0) + &
                    e_re%f2(3,i+1,j ) * w1(1) + &
                    e_re%f2(3,i+2,j ) * w1(2) ) * w2(0) + &
                  ( e_re%f2(3,i-2,j+1) * w1(-2) + &
                    e_re%f2(3,i-1,j+1) * w1(-1) + &
                    e_re%f2(3,i ,j+1) * w1(0) + &
                    e_re%f2(3,i+1,j+1) * w1(1) + &
                    e_re%f2(3,i+2,j+1) * w1(2) ) * w2(1) + &
                  ( e_re%f2(3,i-2,j+2) * w1(-2) + &
                    e_re%f2(3,i-1,j+2) * w1(-1) + &
                    e_re%f2(3,i ,j+2) * w1(0) + &
                    e_re%f2(3,i+1,j+2) * w1(1) + &
                    e_re%f2(3,i+2,j+2) * w1(2) ) * w2(2)

    bp_re(1,l) = ( b_re%f2(1,i-2,jh-2) * w1(-2) + &
                    b_re%f2(1,i-1,jh-2) * w1(-1) + &
                    b_re%f2(1,i ,jh-2) * w1(0) + &
                    b_re%f2(1,i+1,jh-2) * w1(1) + &
                    b_re%f2(1,i+2,jh-2) * w1(2) ) * w2h(-2) + &
                  ( b_re%f2(1,i-2,jh-1) * w1(-2) + &
                    b_re%f2(1,i-1,jh-1) * w1(-1) + &
                    b_re%f2(1,i ,jh-1) * w1(0) + &
                    b_re%f2(1,i+1,jh-1) * w1(1) + &
                    b_re%f2(1,i+2,jh-1) * w1(2) ) * w2h(-1) + &
                  ( b_re%f2(1,i-2,jh ) * w1(-2) + &
                    b_re%f2(1,i-1,jh ) * w1(-1) + &
                    b_re%f2(1,i ,jh ) * w1(0) + &
                    b_re%f2(1,i+1,jh ) * w1(1) + &
                    b_re%f2(1,i+2,jh ) * w1(2) ) * w2h(0) + &
                  ( b_re%f2(1,i-2,jh+1) * w1(-2) + &
                    b_re%f2(1,i-1,jh+1) * w1(-1) + &
                    b_re%f2(1,i ,jh+1) * w1(0) + &
                    b_re%f2(1,i+1,jh+1) * w1(1) + &
                    b_re%f2(1,i+2,jh+1) * w1(2) ) * w2h(1) + &
                  ( b_re%f2(1,i-2,jh+2) * w1(-2) + &
                    b_re%f2(1,i-1,jh+2) * w1(-1) + &
                    b_re%f2(1,i ,jh+2) * w1(0) + &
                    b_re%f2(1,i+1,jh+2) * w1(1) + &
                    b_re%f2(1,i+2,jh+2) * w1(2) ) * w2h(2)

    bp_re(2,l) = ( b_re%f2(2,ih-2,j-2) * w1h(-2) + &
                    b_re%f2(2,ih-1,j-2) * w1h(-1) + &
                    b_re%f2(2,ih ,j-2) * w1h(0) + &
                    b_re%f2(2,ih+1,j-2) * w1h(1) + &
                    b_re%f2(2,ih+2,j-2) * w1h(2) ) * w2(-2) + &
                  ( b_re%f2(2,ih-2,j-1) * w1h(-2) + &
                    b_re%f2(2,ih-1,j-1) * w1h(-1) + &
                    b_re%f2(2,ih ,j-1) * w1h(0) + &
                    b_re%f2(2,ih+1,j-1) * w1h(1) + &
                    b_re%f2(2,ih+2,j-1) * w1h(2) ) * w2(-1) + &
                  ( b_re%f2(2,ih-2,j ) * w1h(-2) + &
                    b_re%f2(2,ih-1,j ) * w1h(-1) + &
                    b_re%f2(2,ih ,j ) * w1h(0) + &
                    b_re%f2(2,ih+1,j ) * w1h(1) + &
                    b_re%f2(2,ih+2,j ) * w1h(2) ) * w2(0) + &
                  ( b_re%f2(2,ih-2,j+1) * w1h(-2) + &
                    b_re%f2(2,ih-1,j+1) * w1h(-1) + &
                    b_re%f2(2,ih ,j+1) * w1h(0) + &
                    b_re%f2(2,ih+1,j+1) * w1h(1) + &
                    b_re%f2(2,ih+2,j+1) * w1h(2) ) * w2(1) + &
                  ( b_re%f2(2,ih-2,j+2) * w1h(-2) + &
                    b_re%f2(2,ih-1,j+2) * w1h(-1) + &
                    b_re%f2(2,ih ,j+2) * w1h(0) + &
                    b_re%f2(2,ih+1,j+2) * w1h(1) + &
                    b_re%f2(2,ih+2,j+2) * w1h(2) ) * w2(2)

    bp_re(3,l) = ( b_re%f2(3,ih-2,jh-2) * w1h(-2) + &
                    b_re%f2(3,ih-1,jh-2) * w1h(-1) + &
                    b_re%f2(3,ih ,jh-2) * w1h(0) + &
                    b_re%f2(3,ih+1,jh-2) * w1h(1) + &
                    b_re%f2(3,ih+2,jh-2) * w1h(2) ) * w2h(-2) + &
                  ( b_re%f2(3,ih-2,jh-1) * w1h(-2) + &
                    b_re%f2(3,ih-1,jh-1) * w1h(-1) + &
                    b_re%f2(3,ih ,jh-1) * w1h(0) + &
                    b_re%f2(3,ih+1,jh-1) * w1h(1) + &
                    b_re%f2(3,ih+2,jh-1) * w1h(2) ) * w2h(-1) + &
                  ( b_re%f2(3,ih-2,jh ) * w1h(-2) + &
                    b_re%f2(3,ih-1,jh ) * w1h(-1) + &
                    b_re%f2(3,ih ,jh ) * w1h(0) + &
                    b_re%f2(3,ih+1,jh ) * w1h(1) + &
                    b_re%f2(3,ih+2,jh ) * w1h(2) ) * w2h(0) + &
                  ( b_re%f2(3,ih-2,jh+1) * w1h(-2) + &
                    b_re%f2(3,ih-1,jh+1) * w1h(-1) + &
                    b_re%f2(3,ih ,jh+1) * w1h(0) + &
                    b_re%f2(3,ih+1,jh+1) * w1h(1) + &
                    b_re%f2(3,ih+2,jh+1) * w1h(2) ) * w2h(1) + &
                  ( b_re%f2(3,ih-2,jh+2) * w1h(-2) + &
                    b_re%f2(3,ih-1,jh+2) * w1h(-1) + &
                    b_re%f2(3,ih ,jh+2) * w1h(0) + &
                    b_re%f2(3,ih+1,jh+2) * w1h(1) + &
                    b_re%f2(3,ih+2,jh+2) * w1h(2) ) * w2h(2)

    ep_im(1,l) = ( e_im%f2(1,ih-2,j-2) * w1h(-2) + &
                    e_im%f2(1,ih-1,j-2) * w1h(-1) + &
                    e_im%f2(1,ih ,j-2) * w1h(0) + &
                    e_im%f2(1,ih+1,j-2) * w1h(1) + &
                    e_im%f2(1,ih+2,j-2) * w1h(2) ) * w2(-2) + &
                  ( e_im%f2(1,ih-2,j-1) * w1h(-2) + &
                    e_im%f2(1,ih-1,j-1) * w1h(-1) + &
                    e_im%f2(1,ih ,j-1) * w1h(0) + &
                    e_im%f2(1,ih+1,j-1) * w1h(1) + &
                    e_im%f2(1,ih+2,j-1) * w1h(2) ) * w2(-1) + &
                  ( e_im%f2(1,ih-2,j ) * w1h(-2) + &
                    e_im%f2(1,ih-1,j ) * w1h(-1) + &
                    e_im%f2(1,ih ,j ) * w1h(0) + &
                    e_im%f2(1,ih+1,j ) * w1h(1) + &
                    e_im%f2(1,ih+2,j ) * w1h(2) ) * w2(0) + &
                  ( e_im%f2(1,ih-2,j+1) * w1h(-2) + &
                    e_im%f2(1,ih-1,j+1) * w1h(-1) + &
                    e_im%f2(1,ih ,j+1) * w1h(0) + &
                    e_im%f2(1,ih+1,j+1) * w1h(1) + &
                    e_im%f2(1,ih+2,j+1) * w1h(2) ) * w2(1) + &
                  ( e_im%f2(1,ih-2,j+2) * w1h(-2) + &
                    e_im%f2(1,ih-1,j+2) * w1h(-1) + &
                    e_im%f2(1,ih ,j+2) * w1h(0) + &
                    e_im%f2(1,ih+1,j+2) * w1h(1) + &
                    e_im%f2(1,ih+2,j+2) * w1h(2) ) * w2(2)

    ep_im(2,l) = ( e_im%f2(2,i-2,jh-2) * w1(-2) + &
                    e_im%f2(2,i-1,jh-2) * w1(-1) + &
                    e_im%f2(2,i ,jh-2) * w1(0) + &
                    e_im%f2(2,i+1,jh-2) * w1(1) + &
                    e_im%f2(2,i+2,jh-2) * w1(2) ) * w2h(-2) + &
                  ( e_im%f2(2,i-2,jh-1) * w1(-2) + &
                    e_im%f2(2,i-1,jh-1) * w1(-1) + &
                    e_im%f2(2,i ,jh-1) * w1(0) + &
                    e_im%f2(2,i+1,jh-1) * w1(1) + &
                    e_im%f2(2,i+2,jh-1) * w1(2) ) * w2h(-1) + &
                  ( e_im%f2(2,i-2,jh ) * w1(-2) + &
                    e_im%f2(2,i-1,jh ) * w1(-1) + &
                    e_im%f2(2,i ,jh ) * w1(0) + &
                    e_im%f2(2,i+1,jh ) * w1(1) + &
                    e_im%f2(2,i+2,jh ) * w1(2) ) * w2h(0) + &
                  ( e_im%f2(2,i-2,jh+1) * w1(-2) + &
                    e_im%f2(2,i-1,jh+1) * w1(-1) + &
                    e_im%f2(2,i ,jh+1) * w1(0) + &
                    e_im%f2(2,i+1,jh+1) * w1(1) + &
                    e_im%f2(2,i+2,jh+1) * w1(2) ) * w2h(1) + &
                  ( e_im%f2(2,i-2,jh+2) * w1(-2) + &
                    e_im%f2(2,i-1,jh+2) * w1(-1) + &
                    e_im%f2(2,i ,jh+2) * w1(0) + &
                    e_im%f2(2,i+1,jh+2) * w1(1) + &
                    e_im%f2(2,i+2,jh+2) * w1(2) ) * w2h(2)

    ep_im(3,l) = ( e_im%f2(3,i-2,j-2) * w1(-2) + &
                    e_im%f2(3,i-1,j-2) * w1(-1) + &
                    e_im%f2(3,i ,j-2) * w1(0) + &
                    e_im%f2(3,i+1,j-2) * w1(1) + &
                    e_im%f2(3,i+2,j-2) * w1(2) ) * w2(-2) + &
                  ( e_im%f2(3,i-2,j-1) * w1(-2) + &
                    e_im%f2(3,i-1,j-1) * w1(-1) + &
                    e_im%f2(3,i ,j-1) * w1(0) + &
                    e_im%f2(3,i+1,j-1) * w1(1) + &
                    e_im%f2(3,i+2,j-1) * w1(2) ) * w2(-1) + &
                  ( e_im%f2(3,i-2,j ) * w1(-2) + &
                    e_im%f2(3,i-1,j ) * w1(-1) + &
                    e_im%f2(3,i ,j ) * w1(0) + &
                    e_im%f2(3,i+1,j ) * w1(1) + &
                    e_im%f2(3,i+2,j ) * w1(2) ) * w2(0) + &
                  ( e_im%f2(3,i-2,j+1) * w1(-2) + &
                    e_im%f2(3,i-1,j+1) * w1(-1) + &
                    e_im%f2(3,i ,j+1) * w1(0) + &
                    e_im%f2(3,i+1,j+1) * w1(1) + &
                    e_im%f2(3,i+2,j+1) * w1(2) ) * w2(1) + &
                  ( e_im%f2(3,i-2,j+2) * w1(-2) + &
                    e_im%f2(3,i-1,j+2) * w1(-1) + &
                    e_im%f2(3,i ,j+2) * w1(0) + &
                    e_im%f2(3,i+1,j+2) * w1(1) + &
                    e_im%f2(3,i+2,j+2) * w1(2) ) * w2(2)

    bp_im(1,l) = ( b_im%f2(1,i-2,jh-2) * w1(-2) + &
                    b_im%f2(1,i-1,jh-2) * w1(-1) + &
                    b_im%f2(1,i ,jh-2) * w1(0) + &
                    b_im%f2(1,i+1,jh-2) * w1(1) + &
                    b_im%f2(1,i+2,jh-2) * w1(2) ) * w2h(-2) + &
                  ( b_im%f2(1,i-2,jh-1) * w1(-2) + &
                    b_im%f2(1,i-1,jh-1) * w1(-1) + &
                    b_im%f2(1,i ,jh-1) * w1(0) + &
                    b_im%f2(1,i+1,jh-1) * w1(1) + &
                    b_im%f2(1,i+2,jh-1) * w1(2) ) * w2h(-1) + &
                  ( b_im%f2(1,i-2,jh ) * w1(-2) + &
                    b_im%f2(1,i-1,jh ) * w1(-1) + &
                    b_im%f2(1,i ,jh ) * w1(0) + &
                    b_im%f2(1,i+1,jh ) * w1(1) + &
                    b_im%f2(1,i+2,jh ) * w1(2) ) * w2h(0) + &
                  ( b_im%f2(1,i-2,jh+1) * w1(-2) + &
                    b_im%f2(1,i-1,jh+1) * w1(-1) + &
                    b_im%f2(1,i ,jh+1) * w1(0) + &
                    b_im%f2(1,i+1,jh+1) * w1(1) + &
                    b_im%f2(1,i+2,jh+1) * w1(2) ) * w2h(1) + &
                  ( b_im%f2(1,i-2,jh+2) * w1(-2) + &
                    b_im%f2(1,i-1,jh+2) * w1(-1) + &
                    b_im%f2(1,i ,jh+2) * w1(0) + &
                    b_im%f2(1,i+1,jh+2) * w1(1) + &
                    b_im%f2(1,i+2,jh+2) * w1(2) ) * w2h(2)

    bp_im(2,l) = ( b_im%f2(2,ih-2,j-2) * w1h(-2) + &
                    b_im%f2(2,ih-1,j-2) * w1h(-1) + &
                    b_im%f2(2,ih ,j-2) * w1h(0) + &
                    b_im%f2(2,ih+1,j-2) * w1h(1) + &
                    b_im%f2(2,ih+2,j-2) * w1h(2) ) * w2(-2) + &
                  ( b_im%f2(2,ih-2,j-1) * w1h(-2) + &
                    b_im%f2(2,ih-1,j-1) * w1h(-1) + &
                    b_im%f2(2,ih ,j-1) * w1h(0) + &
                    b_im%f2(2,ih+1,j-1) * w1h(1) + &
                    b_im%f2(2,ih+2,j-1) * w1h(2) ) * w2(-1) + &
                  ( b_im%f2(2,ih-2,j ) * w1h(-2) + &
                    b_im%f2(2,ih-1,j ) * w1h(-1) + &
                    b_im%f2(2,ih ,j ) * w1h(0) + &
                    b_im%f2(2,ih+1,j ) * w1h(1) + &
                    b_im%f2(2,ih+2,j ) * w1h(2) ) * w2(0) + &
                  ( b_im%f2(2,ih-2,j+1) * w1h(-2) + &
                    b_im%f2(2,ih-1,j+1) * w1h(-1) + &
                    b_im%f2(2,ih ,j+1) * w1h(0) + &
                    b_im%f2(2,ih+1,j+1) * w1h(1) + &
                    b_im%f2(2,ih+2,j+1) * w1h(2) ) * w2(1) + &
                  ( b_im%f2(2,ih-2,j+2) * w1h(-2) + &
                    b_im%f2(2,ih-1,j+2) * w1h(-1) + &
                    b_im%f2(2,ih ,j+2) * w1h(0) + &
                    b_im%f2(2,ih+1,j+2) * w1h(1) + &
                    b_im%f2(2,ih+2,j+2) * w1h(2) ) * w2(2)

    bp_im(3,l) = ( b_im%f2(3,ih-2,jh-2) * w1h(-2) + &
                    b_im%f2(3,ih-1,jh-2) * w1h(-1) + &
                    b_im%f2(3,ih ,jh-2) * w1h(0) + &
                    b_im%f2(3,ih+1,jh-2) * w1h(1) + &
                    b_im%f2(3,ih+2,jh-2) * w1h(2) ) * w2h(-2) + &
                  ( b_im%f2(3,ih-2,jh-1) * w1h(-2) + &
                    b_im%f2(3,ih-1,jh-1) * w1h(-1) + &
                    b_im%f2(3,ih ,jh-1) * w1h(0) + &
                    b_im%f2(3,ih+1,jh-1) * w1h(1) + &
                    b_im%f2(3,ih+2,jh-1) * w1h(2) ) * w2h(-1) + &
                  ( b_im%f2(3,ih-2,jh ) * w1h(-2) + &
                    b_im%f2(3,ih-1,jh ) * w1h(-1) + &
                    b_im%f2(3,ih ,jh ) * w1h(0) + &
                    b_im%f2(3,ih+1,jh ) * w1h(1) + &
                    b_im%f2(3,ih+2,jh ) * w1h(2) ) * w2h(0) + &
                  ( b_im%f2(3,ih-2,jh+1) * w1h(-2) + &
                    b_im%f2(3,ih-1,jh+1) * w1h(-1) + &
                    b_im%f2(3,ih ,jh+1) * w1h(0) + &
                    b_im%f2(3,ih+1,jh+1) * w1h(1) + &
                    b_im%f2(3,ih+2,jh+1) * w1h(2) ) * w2h(1) + &
                  ( b_im%f2(3,ih-2,jh+2) * w1h(-2) + &
                    b_im%f2(3,ih-1,jh+2) * w1h(-1) + &
                    b_im%f2(3,ih ,jh+2) * w1h(0) + &
                    b_im%f2(3,ih+1,jh+2) * w1h(1) + &
                    b_im%f2(3,ih+2,jh+2) * w1h(2) ) * w2h(2)
  enddo


end subroutine get_emf_cyl_m_s4
!-------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function signbit(x)
!---------------------------------------------------------------------------------------------------
  implicit none
  real(p_k_fld), intent(in) :: x

  integer :: signbit

  if ( x < 0 ) then
    signbit = 1
  else
    signbit = 0
  endif

end function signbit
!---------------------------------------------------------------------------------------------------

end module
