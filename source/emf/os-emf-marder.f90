! ----------------------------------------------------------------------------------------
! This module implements a Marder-Langdon correction for charge conservation following
! the improved method devised by Langdon in [1], as described in section 3, eq. (8).
!
! It is meant to be used when field solvers other than the Yee solver are being used,
! since we have exact charge conservation in that case.
!
! [1] LANGDON, A. (1992). On Enforcing Gauss Law in Electromagnetic Particle-in-Cell Codes.
!      Computer Physics Communications, 70(3), 447450. doi:10.1016/0010-4655(92)90105-8
! ----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_marder

use m_emf_define
use m_vdf_define
use m_parameters
use m_node_conf

use m_system

implicit none

private

interface update_charge_cons
  module procedure update_charge_cons
end interface

interface setup_marder_langdon
  module procedure setup_marder_langdon
end interface

interface marder_langdon
  module procedure marder_langdon
end interface

public :: update_charge_cons, setup_marder_langdon, marder_langdon

contains

! ----------------------------------------------------------------------------------------
!

! ----------------------------------------------------------------------------------------
subroutine setup_marder_langdon( emf, dx, dt )

  use m_vdf

  implicit none

  class( t_emf ), intent(inout) :: emf
  real( p_double ), dimension(:), intent(in) :: dx
  real( p_double ), intent(in) :: dt

  real( p_double ) :: dmax

  if ( emf%marder_d > 0 ) then

    ! This is for 2D
    dmax = (dx(1)**2 * dx(2)**2 / (dx(1)**2 + dx(2)**2)) / 2 / dt

    SCR_ROOT( 'dmax = ', dmax )

    call emf % f % new( emf%e, f_dim = 1 )

  endif

end subroutine setup_marder_langdon

! ----------------------------------------------------------------------------------------
! Update the charge conservation grid f = div.E - \rho
! ----------------------------------------------------------------------------------------
subroutine update_charge_cons( emf, rho, send_msg, recv_msg, no_co )

  use m_vdf
  use m_vdf_comm

  implicit none

  class( t_emf ), intent(inout) :: emf
  type( t_vdf ),    intent(in) :: rho
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  class( t_node_conf ), optional, intent(in) :: no_co

  ! Allocate charge conservation grid if required
  if ( emf%f%x_dim_ < 1 ) call emf % f % new( emf%e, f_dim = 1 )

  if ( emf%coordinates ==  p_cylindrical_b .or. p_x_dim /= 2) then
     ERROR('Not implemented yet')
     call abort_program( p_err_notimplemented )
  endif

  select case ( emf%solver )

    case ( p_emf_stencil )
      call charge_cons_2d_stencil( emf%e, rho, emf%f, &
                                   emf%stencil_k1, emf%stencil_k2 )
    case ( p_emf_4order )
      call charge_cons_2d_4order( emf%e, rho, emf%f )

    case ( p_emf_yee )
      call charge_cons_2d_yee( emf%e, rho, emf%f )

    case default
       ERROR('Not implemented yet')
       call abort_program( p_err_notimplemented )

  end select

  ! Update guard cell values if requested
  if ( present(no_co) ) then
    call update_boundary( emf%f, p_vdf_replace, no_co, send_msg, recv_msg )
  endif

end subroutine update_charge_cons
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! Charge conservation using Yee type derivatives
! ----------------------------------------------------------------------------------------
subroutine charge_cons_2d_yee( e, rho, f )

  implicit none

  type( t_vdf ), intent(in) :: e, rho
  type( t_vdf ), intent(inout) :: f

  integer :: i, i1, i2
  integer, dimension(2) :: lb, ub
  real( p_k_fld ), dimension(2) :: rdx

  do i = 1, 2
    lb(i) = lbound( e%f2, i+1 ) + 1
    ub(i) = ubound( e%f2, i+1 )

    rdx(i) = 1.0d0 / e%dx_(i)
  enddo

  ! Calculate f using the appropriate spatial derivatives
  do i2 = lb(2), ub(2)
    do i1 = lb(1), ub(1)

      f%f2(1, i1, i2 ) = rdx(1) * ( e%f2(1, i1,i2 ) - e%f2(1, i1-1,i2 ) ) + &
                         rdx(2) * ( e%f2(2, i1,i2 ) - e%f2(2, i1,i2-1 ) ) - &
                         rho%f2( 1, i1, i2 )

    enddo
  enddo

end subroutine charge_cons_2d_yee
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! Charge conservation using 4th order accurate derivatives
! ----------------------------------------------------------------------------------------
subroutine charge_cons_2d_4order( e, rho, f )

  implicit none

  type( t_vdf ), intent(in) :: e, rho
  type( t_vdf ), intent(inout) :: f

  integer :: i, i1, i2
  integer, dimension(2) :: lb, ub
  real( p_k_fld ), dimension(2) :: rdx

  do i = 1, 2
    lb(i) = lbound( e%f2, i+1 ) + 2
    ub(i) = ubound( e%f2, i+1 ) - 1

    rdx(i) = 1.0d0 / ( 24 * e%dx_(i) )
  enddo

  ! Calculate f using the appropriate spatial derivatives
  do i2 = lb(2), ub(2)
    do i1 = lb(1), ub(1)

      f%f2(1, i1, i2 ) = rdx(1) * ( e%f2(1, i1-2,i2 ) - 27*e%f2(1, i1-1,i2 ) + 27*e%f2(1, i1,i2 ) - e%f2(1, i1+1,i2 ) ) + &
                         rdx(2) * ( e%f2(2, i1,i2-2 ) - 27*e%f2(2, i1,i2-1 ) + 27*e%f2(2, i1,i2 ) - e%f2(2, i1, i2+1) ) - &
                         rho%f2( 1, i1, i2 )

    enddo
  enddo

end subroutine charge_cons_2d_4order
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! Charge conservation using stencil derivatives
! ----------------------------------------------------------------------------------------
subroutine charge_cons_2d_stencil( e, rho, f, k1, k2 )

  implicit none

  type( t_vdf ), intent(in) :: e, rho
  type( t_vdf ), intent(inout) :: f

  real(p_double),  intent(in) :: k1, k2

  integer :: i, i1, i2
  integer, dimension(2) :: lb, ub
  real( p_k_fld ), dimension(2) :: rdx
  real( p_k_fld ) :: a0, a1, a2

  ! Spatial difference coefficients
  a0 = 1 - k1 - k2
  a1 = k1 / 3
  a2 = k2 / 6

  do i = 1, 2
    lb(i) = lbound( e%f2, i+1 ) + 2
    ub(i) = ubound( e%f2, i+1 ) - 1

    rdx(i) = 1.0d0 / e%dx_(i)
  enddo

  ! Calculate f using the appropriate spatial derivatives
  do i2 = lb(2), ub(2)
    do i1 = lb(1), ub(1)

      ! Yee calculation
      ! f(i1,i2) = ( rdx(1) * ( e%f2(1, i1,i2 ) - e%f2(1, i1-1,i2 ) ) + &
      !              rdx(2) * ( e%f2(2, i1,i2 ) - e%f2(1, i1,i2-1 ) ) ) - &

      !            rho%f2( 1, i1, i2 ) )

      ! Stencil calculation
      f%f2(1,i1,i2) = rdx(1) * ( a0 * (   e%f2(1, i1,  i2   ) - e%f2(1, i1-1,i2   ) ) + &
                                 a1 * (   e%f2(1, i1+1,i2   ) - e%f2(1, i1-2,i2   ) ) + &
                                 a2 * ( ( e%f2(1, i1+1,i2+1 ) - e%f2(1, i1-2,i2+1 ) ) + &
                                        ( e%f2(1, i1+1,i2-1 ) - e%f2(1, i1-2,i2-1 ) ) ) ) + &

                   rdx(2) * ( a0 * (   e%f2(2, i1,  i2   ) - e%f2(2, i1,  i2-1 ) ) + &
                              a1 * (   e%f2(2, i1,  i2+1 ) - e%f2(2, i1,  i2-2 ) ) + &
                              a2 * ( ( e%f2(2, i1+1,i2+1 ) - e%f2(2, i1+1,i2-2 ) ) + &
                                     ( e%f2(2, i1-1,i2+1 ) - e%f2(2, i1-1,i2-2 ) ) ) ) - &

                   rho%f2(1,i1,i2)

    enddo
  enddo

end subroutine charge_cons_2d_stencil
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
!
! ----------------------------------------------------------------------------------------
subroutine ml_pass_2d_yee( e, f, ddt )

  implicit none

  type( t_vdf ), intent(inout) :: e
  type( t_vdf ), intent(in) :: f
  real( p_double ), intent(in) :: ddt ! d * dt

  real( p_k_fld ), dimension(2) :: rc
  integer :: i1, i2

  rc(1) = ddt / e%dx_(1)
  rc(2) = ddt / e%dx_(2)

  ! Add the correction to the electric field
  do i2 = lbound(e%f2,3), ubound(e%f2,3) - 1
    do i1 = lbound(e%f2,2), ubound(e%f2,2) - 1

      ! Yee calculation
      e%f2( 1, i1, i2 ) =  e%f2( 1, i1, i2 ) + rc(1) * ( f%f2(1,i1+1, i2) - f%f2(1,i1, i2) )

      e%f2( 2, i1, i2 ) =  e%f2( 2, i1, i2 ) + rc(2) * ( f%f2(1,i1, i2+1) - f%f2(1,i1, i2) )

    enddo
  enddo

end subroutine ml_pass_2d_yee
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
!
! ----------------------------------------------------------------------------------------
subroutine ml_pass_2d_4order( e, f, ddt )

  implicit none

  type( t_vdf ), intent(inout) :: e
  type( t_vdf ), intent(in) :: f
  real( p_double ), intent(in) :: ddt ! d * dt

  real( p_k_fld ), dimension(2) :: rc
  integer :: i1, i2

  rc(1) = ddt /( 24 * e%dx_(1) )
  rc(2) = ddt /( 24 * e%dx_(2) )

  ! Add the correction to the electric field
  do i2 = lbound(e%f2,3), ubound(e%f2,3) - 1
    do i1 = lbound(e%f2,2), ubound(e%f2,2) - 1

      ! 4th order
      e%f2( 1, i1, i2 ) =  e%f2( 1, i1, i2 ) + rc(1) * ( f%f2(1,i1-1, i2) - 27 * f%f2(1,i1, i2) + 27 * f%f2(1,i1+1,i2) - f%f2(1,i1+2,i2) )

      e%f2( 2, i1, i2 ) =  e%f2( 2, i1, i2 ) + rc(2) * ( f%f2(1,i1, i2-1) - 27 * f%f2(1,i1, i2) + 27 * f%f2(1,i1,i2+1) - f%f2(1,i1,i2+2) )

    enddo
  enddo

end subroutine ml_pass_2d_4order
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! 2D correction using 'stencil' spatial derivatives. See Greenwood, A., et. al. (2004),
! Journal of Computational Physics, 201(2), 665684. doi:10.1016/j.jcp.2004.06.021 for
! details.
!
! Special cases:
!   k1 = 0, k2 = 0      : Yee scheme
!   k1 = -1/8, k2 = 0   : 4th order spatially accurate scheme
! ----------------------------------------------------------------------------------------
subroutine ml_pass_2d_stencil( e, f, ddt, k1, k2 )

  implicit none

  type( t_vdf ), intent(inout) :: e
  type( t_vdf ), intent(in) :: f
  real( p_double ), intent(in) :: ddt ! d * dt
  real( p_double ), intent(in) :: k1, k2

  real( p_k_fld ), dimension(2) :: rc
  real( p_k_fld ) :: a0, a1, a2
  integer :: i1, i2

  ! Spatial difference coefficients
  a0 = 1 - k1 - k2
  a1 = k1 / 3
  a2 = k2 / 6

  rc(1) = ddt / e%dx_(1)
  rc(2) = ddt / e%dx_(2)

   ! Add the correction to the electric field
   do i2 = lbound(e%f2,3)+1, ubound(e%f2,3) - 2
     do i1 = lbound(e%f2,2)+1, ubound(e%f2,2) - 2

        ! Stencil calculation
        e%f2( 1, i1, i2 ) =  e%f2( 1, i1, i2 ) + &
                             rc(1) * ( a0 * (   f%f2(1,i1+1, i2  ) - f%f2(1,i1,   i2  ) ) + &
                                       a1 * (   f%f2(1,i1+2, i2  ) - f%f2(1,i1-1, i2  ) ) + &
                                       a2 * ( ( f%f2(1,i1+2, i2+1) - f%f2(1,i1-1, i2+1) ) + &
                                              ( f%f2(1,i1+2, i2-1) - f%f2(1,i1-1, i2-1) ) ) )

        e%f2( 2, i1, i2 ) =  e%f2( 2, i1, i2 ) + &
                             rc(2) * ( a0 * (   f%f2(1,i1  , i2+1) - f%f2(1,i1,   i2  ) ) + &
                                       a1 * (   f%f2(1,i1  , i2+2) - f%f2(1,i1,   i2-1) ) + &
                                       a2 * ( ( f%f2(1,i1+1, i2+2) - f%f2(1,i1+1, i2-1) ) + &
                                              ( f%f2(1,i1-1, i2+2) - f%f2(1,i1-1, i2-1) ) ) )

     enddo
   enddo

end subroutine ml_pass_2d_stencil
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
!
! ----------------------------------------------------------------------------------------
subroutine ml_pass( emf, ddt, send_msg, recv_msg, no_co )

  use m_vdf_comm

  implicit none

  class( t_emf ), intent(inout) :: emf
  real( p_double ), intent(in) :: ddt
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  class( t_node_conf ), optional, intent(in) :: no_co

  if ( emf%coordinates ==  p_cylindrical_b .or. p_x_dim /= 2) then
     ERROR('Not implemented yet')
     call abort_program( p_err_notimplemented )
  endif

  select case ( emf%solver )

    case ( p_emf_stencil )
      call ml_pass_2d_stencil( emf%e, emf%f, ddt, &
                               emf%stencil_k1, emf%stencil_k2 )

    case ( p_emf_yee )
      call ml_pass_2d_yee( emf%e, emf%f, ddt )

    case ( p_emf_4order )
      call ml_pass_2d_4order( emf%e, emf%f, ddt )

    case default
       ERROR('Not implemented yet')
       call abort_program( p_err_notimplemented )

  end select

  ! Update guard cell values if requested
  if ( present(no_co) ) then
    call update_boundary( emf%e, p_vdf_replace, no_co, send_msg, recv_msg )
  endif

end subroutine ml_pass

! ----------------------------------------------------------------------------------------
!
! ----------------------------------------------------------------------------------------
subroutine marder_langdon( emf, rho, dt, no_co, send_msg, recv_msg )

  use m_vdf_math
  use m_vdf_comm, only : t_vdf_msg

  implicit none

  class( t_emf ), intent(inout) :: emf
  type( t_vdf ),    intent(in) :: rho
  real(p_double),   intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  real( p_double ) :: ddt
  real(p_k_fld), dimension(1) :: max_Fa, max_Fb
  integer :: i

  ddt = emf%marder_d * dt

  do i = 1, emf%marder_n
     ! Update the charge conservation error value
     call update_charge_cons( emf, rho, send_msg, recv_msg, no_co )

     ! (* debug *)

     max_Fa = maxval_abs( emf%f, no_co )
     SCR_ROOT( i, ' - max(|f|) = ',max_Fa )

     ! Do a Marder-Langdon pass
     call ml_pass( emf, ddt, send_msg, recv_msg, no_co )
  enddo

  ! (* debug *)

  call update_charge_cons( emf, rho, send_msg, recv_msg, no_co )
  max_Fb = maxval_abs( emf%f, no_co )
  SCR_ROOT( ' Final - max(|f|) = ',max_Fb )

end subroutine marder_langdon
! ----------------------------------------------------------------------------------------

end module m_emf_marder
