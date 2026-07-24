# 1 "gr/os-species-push-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-species-push-gr.f03" 2
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
# 2 "gr/os-species-push-gr.f03" 2
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
# 3 "gr/os-species-push-gr.f03" 2

module m_species_push_gr

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
# 7 "gr/os-species-push-gr.f03" 2

use m_system
use m_parameters

use m_emf_define, only : t_emf
use m_vdf_define, only : t_vdf
use m_time_step, only : t_time_step
use m_species_define_gr, only : t_species_gr
use m_math, only : pi, pi_2

implicit none

private

interface advance_deposit_gr
    module procedure advance_deposit_gr
end interface advance_deposit_gr

interface dudt_boris_gr
  module procedure dudt_boris_gr
end interface

interface get_kvec_lvec_rk_gr
  module procedure get_kvec_lvec_rk_gr
end interface

public :: advance_deposit_gr
public :: dudt_boris_gr
public :: get_kvec_lvec_rk_gr

contains

!-------------------------------------------------------------------------------
! Push particles and deposit electric current
!-------------------------------------------------------------------------------
subroutine advance_deposit_gr( this, emf, current, t, tstep, tid, &
  n_threads )

  use m_time_step, only : dt
  use m_geometry_gr

  implicit none

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: current

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep

  integer, intent(in) :: tid ! local thread id
  integer, intent(in) :: n_threads ! total number of threads

  ! local variables

  real(p_double) :: dtcycle, energy
  integer :: chunk, i0, i1

  ! executable statements

  ! if before push start time return silently
  if ( t < this%push_start_time ) return

  ! initialize time centered energy diagnostic
  energy = 0.0_p_double

  dtcycle = real( dt(tstep), p_k_part )

  ! range of particles for each thread
  chunk = ( this%num_par + n_threads - 1 ) / n_threads
  i0 = tid * chunk + 1
  i1 = min( (tid+1) * chunk, this%num_par )

  ! Push particles. Boundary crossings will be checked at update_boundary
  select case ( p_x_dim )
  case (2)
    select case ( this%geometry%metric )
    case ( p_geometry_minkowski )
      select case (this%push_type)
      case(p_std)
        call advance_deposit_2d_mink_gr( this, emf, current, energy, dtcycle, i0, i1, t )
      case default
        call advance_deposit_2d_rk_gr( this, emf, current, energy, dtcycle, i0, i1 )
      end select

    case ( p_geometry_schwarzschild )
      call advance_deposit_2d_rk_gr( this, emf, current, energy, dtcycle, i0, i1 )

    case ( p_geometry_kerr_slow )
      call advance_deposit_2d_rk_gr( this, emf, current, energy, dtcycle, i0, i1 )

    case ( p_geometry_kerr )
      call advance_deposit_2d_rk_gr( this, emf, current, energy, dtcycle, i0, i1 )

    end select
  case default
    write(err_buf__,*) 'Not implemented for x_dim = ',p_x_dim;call err__("gr/os-species-push-gr.f03",103)
    call abort_program(p_err_invalid)
  end select

  this%energy(tid+1) = energy

end subroutine advance_deposit_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine advance_deposit_2d_mink_gr( this, emf, jay, energy, gdt, i0, i1, time )
!-----------------------------------------------------------------------------

  use m_species_push, only : ntrim
  use m_geometry_gr
  use m_species_current_gr, only : dep_current_2d_srect_mink_gr,dep_current_2d_sint_mink_gr
  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: i0, i1
  real(p_double), intent(in) :: time

  integer :: i, pp, np, ptrcur
  real(p_k_part), dimension(rank+3 ,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi
  real(p_k_part), dimension(p_cache_size) :: rgamma

  class( t_geometry_gr ), pointer :: g
  real(p_k_part) :: r, t, dt

  g => this%geometry

  ! update momenta
  if ( .not. this%free_stream ) call this % dudt( emf, gdt, i0, i1, energy, time )

  dt = real(g%dt, p_k_part)

  ! advance position of particles i0 to i1 in chunks of p_cache_size
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
        np = i1 - ptrcur + 1
    else
        np = p_cache_size
    endif

    call this % emit( emf, gdt, ptrcur, np )

    ! this type of loop is actually faster than
    ! using a forall construct
    pp = ptrcur

    do i=1,np

      rgamma(i) = 1.0_p_k_part / &
        sqrt( ( ( 1.0_p_k_part + this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

      pp = pp + 1
    end do

    pp = ptrcur
    do i=1,np
      xbuf(3,i) = this%x(3,pp) + this%p(1,pp) * rgamma(i) * gdt
      xbuf(4,i) = this%x(4,pp) + this%p(2,pp) * rgamma(i) * gdt
      xbuf(5,i) = this%x(5,pp) + this%p(3,pp) * rgamma(i) * gdt

      r = sqrt( xbuf(3,i)**2 + xbuf(4,i)**2 +xbuf(5,i)**2 )
      t = acos( xbuf(5,i) / r )

      dxi(1,i) = ntrim( (r - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
      xbuf(1,i) = (r - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

      dxi(2,i) = ntrim( (t - g%t%f1(2, this%ix(2,pp) )) / dt )
      xbuf(2,i) = (t - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

      pp = pp + 1
    end do

    ! Deposit current
    select case ( g%spline_rule )
        case ( p_spline_rectangular )
          call dep_current_2d_srect_mink_gr( this, jay, dxi, xbuf, &
                              this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                              this%q(ptrcur:), rgamma, &
                              this%p(:,ptrcur:), &
                              np, gdt )
        case default !p_spline_integral
          call dep_current_2d_sint_mink_gr( this, jay, dxi, xbuf, &
                              this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                              this%q(ptrcur:), rgamma, &
                              this%p(:,ptrcur:), &
                              np, gdt )
    end select

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np

      this%x(3,pp) = xbuf(3,i)
      this%x(4,pp) = xbuf(4,i)
      this%x(5,pp) = xbuf(5,i)

      this%x(1,pp) = xbuf(1,i)
      this%x(2,pp) = xbuf(2,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

      pp = pp + 1
    end do

  enddo

end subroutine advance_deposit_2d_mink_gr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_boris_gr( this, emf, dt, i0, i1, energy, time )

  use m_species_define, only : t_species
  use m_emf_interpolate, only : get_emf
  use m_geometry_gr, only : t_geometry_gr

  implicit none

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  real(p_double), intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: time

  real(p_k_part) :: tem
  integer :: i, ptrcur, np, pp
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp

  real(p_k_part) :: gamma

  real(p_k_part), dimension(p_cache_size) :: gam_tem, otsq
  real(p_k_part) :: u2
  real(p_double) :: loc_ene

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bptemp, eptemp
  real(p_double) :: r, rc, st, ct, sp, cp

  class( t_geometry_gr ), pointer :: geometry

  select type( this )
    class is ( t_species_gr )
      geometry => this%geometry
    class default
      write(err_buf__,*) 'dudt_boris_gr must be called with t_species_gr objects';call err__("gr/os-species-push-gr.f03",260)
      call abort_program( p_err_invalid )
  end select

  ! executable statements
  ! print *, '(*warn*) Using dudt_simd'

  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * (dt / this%rqm), p_k_part )

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! modify bp & ep to include timestep and charge-to-mass ratio
    ! and perform half the electric field acceleration.
    ! Result is stored in UTEMP.

    call get_emf( emf, bptemp, eptemp, this%ix(:,ptrcur:), this%x(:,ptrcur:), &
      np, this%interpolation )

    ! transform E and B to cartesian coordinates
    pp = ptrcur
    do i=1, np

        r = sqrt( this%x(3,pp)**2 + this%x(4,pp)**2 + this%x(5,pp)**2 )
        rc = sqrt( this%x(3,pp)**2 + this%x(4,pp)**2 )

        st = rc / r
        ct = this%x(5,pp) / r

        sp = this%x(4,pp) / rc
        cp = this%x(3,pp) / rc

        ep(1,i) = st * cp * eptemp(1,i) + ct * cp * eptemp(2,i) - sp * eptemp(3,i)
        ep(2,i) = st * sp * eptemp(1,i) + ct * sp * eptemp(2,i) + cp * eptemp(3,i)
        ep(3,i) = ct * eptemp(1,i) - st * eptemp(2,i)

        bp(1,i) = st * cp * bptemp(1,i) + ct * cp * bptemp(2,i) - sp * bptemp(3,i)
        bp(2,i) = st * sp * bptemp(1,i) + ct * sp * bptemp(2,i) + cp * bptemp(3,i)
        bp(3,i) = ct * bptemp(1,i) - st * bptemp(2,i)

        pp = pp + 1
    end do

    do i=1, np
      ep(1,i) = ep(1,i) * tem
      ep(2,i) = ep(2,i) * tem
      ep(3,i) = ep(3,i) * tem
    end do

    loc_ene = 0

    pp = ptrcur
    do i=1,np
      utemp(1,i) = this%p(1,pp) + ep(1,i)
      utemp(2,i) = this%p(2,pp) + ep(2,i)
      utemp(3,i) = this%p(3,pp) + ep(3,i)

      ! Get time centered gamma
      u2 = (utemp(1,i)**2 + utemp(2,i)**2) + utemp(3,i)**2

      gamma = sqrt( u2 + 1 )

      ! accumulate time centered energy
      ! this is done in double precicion always
      loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1.0_p_double)

      gam_tem(i)= tem / gamma

      pp = pp + 1
    enddo

    ! accumulate global energy
    energy = energy + loc_ene

    do i=1,np
      bp(1,i) = bp(1,i)*gam_tem(i)
      bp(2,i) = bp(2,i)*gam_tem(i)
      bp(3,i) = bp(3,i)*gam_tem(i)
    end do

    pp = ptrcur
    do i=1,np
      this%p(1,pp) = utemp(1,i) + utemp(2,i) * bp(3,i)
      this%p(2,pp) = utemp(2,i) + utemp(3,i) * bp(1,i)
      this%p(3,pp) = utemp(3,i) + utemp(1,i) * bp(2,i)
      pp = pp + 1
    end do

    pp = ptrcur
    do i=1,np
      this%p(1,pp) = this%p(1,pp) - utemp(3,i) * bp(2,i)
      this%p(2,pp) = this%p(2,pp) - utemp(1,i) * bp(3,i)
      this%p(3,pp) = this%p(3,pp) - utemp(2,i) * bp(1,i)
      pp = pp + 1
    end do

    do i=1,np
      otsq(i) = 2.0_p_k_part / ( ((1.0_p_k_part + bp(1,i)**2) + bp(2,i)**2) + bp(3,i)**2)
    end do

    do i=1,np
      bp(1,i) = bp(1,i) * otsq(i)
      bp(2,i) = bp(2,i) * otsq(i)
      bp(3,i) = bp(3,i) * otsq(i)
    end do

    pp = ptrcur
    do i=1,np
      utemp(1,i) = utemp(1,i) + this%p(2,pp) * bp(3,i)
      utemp(2,i) = utemp(2,i) + this%p(3,pp) * bp(1,i)
      utemp(3,i) = utemp(3,i) + this%p(1,pp) * bp(2,i)
      pp = pp + 1
    end do

    pp = ptrcur
    do i=1,np
      utemp(1,i) = utemp(1,i) - this%p(3,pp) * bp(2,i)
      utemp(2,i) = utemp(2,i) - this%p(1,pp) * bp(3,i)
      utemp(3,i) = utemp(3,i) - this%p(2,pp) * bp(1,i)
      pp = pp + 1
    end do

    ! Perform second half of electric field acceleration.
    pp = ptrcur
    do i=1,np
      this%p(1,pp) = utemp(1,i) + ep(1,i)
      this%p(2,pp) = utemp(2,i) + ep(2,i)
      this%p(3,pp) = utemp(3,i) + ep(3,i)
      pp = pp + 1
    end do

  enddo

end subroutine dudt_boris_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine advance_deposit_2d_rk_gr( this, emf, jay, energy, gdt, i0, i1 )
!-----------------------------------------------------------------------------

  use m_species_push, only : ntrim
  use m_geometry_gr
  use m_species_current_gr, only : dep_current_2d_srect_curved_gr,dep_current_2d_sint_curved_gr
  use m_species_define_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 4

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: i0, i1

  integer :: i, pp, np, ptrcur

  real(p_k_part), dimension(rank + 1 , p_cache_size) :: xbuf, pinit, xbufinit
  real(p_k_part), dimension(rank , p_cache_size) :: xbufg,xbufginit
  integer, dimension(rank, p_cache_size) :: dxi, dxig, ixinit
  real(p_k_part), dimension(p_cache_size) :: rgamma, p2init, p2new
  real(p_k_part) :: dt

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  !-------------------------------------------
  ! t = t0 inputs: (xp0,wp0)
  !-------------------------------------------
  ! xp0 is the particle's local position (t=t0)
  ! wp0 is the particle's local momenta  (t=t0)

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! loop to get particle's position
    ! xbuf(:,i) = (r,theta,phi) at t=t0
    pp = ptrcur
    do i=1,np
      xbufginit(1,i) = this%x ( 1, pp )
      xbufginit(2,i) = this%x ( 2, pp )

      ixinit(1,i) = this%ix( 1, pp )
      ixinit(2,i) = this%ix( 2, pp )

      pinit(1,i) = this%p(1,pp)
      pinit(2,i) = this%p(2,pp)
      pinit(3,i) = this%p(3,pp)

      p2init(i) = pinit(1,i)**2 + pinit(2,i)**2 + pinit(3,i)**2

      pp = pp + 1
    end do

    ! selection of spherical or cartesian pusher
    if (g%south_axis .and. g%north_axis) then
      select case (this%push_type)
        case(p_rk4)
          select case (this%geometry%classical_rr)
            case(.false.)
              call dudt_dxdt_rk4_cart_gr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
            case(.true.)
              call dudt_dxdt_rk4_cart_gr_rr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
          end select
        case(p_rk6)
          select case (this%geometry%classical_rr)
            case(.false.)
              call dudt_dxdt_rk6_cart_gr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
            case(.true.)
              call dudt_dxdt_rk6_cart_gr_rr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
          end select
      end select
    else if (g%south_axis) then
      select case (this%push_type)
        case(p_rk4)
          select case (this%geometry%classical_rr)
            case(.false.)
              print*,"should not be here!!!!!!------------------------"
              !call dudt_dxdt_rk4_cartsouth_gr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
            case(.true.)
              print*,"should not be here!!!!!!------------------------"
              !call dudt_dxdt_rk4_cartsouth_gr_rr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
          end select
        case(p_rk6)
          select case (this%geometry%classical_rr)
            case(.false.)
              call dudt_dxdt_rk6_cartsouth_gr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
            case(.true.)
              print*,"should not be here!!!!!!------------------------"
              !call dudt_dxdt_rk6_cartsouth_gr_rr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
          end select
      end select
    else if (g%north_axis) then
      select case (this%push_type)
        case(p_rk4)
          select case (this%geometry%classical_rr)
            case(.false.)
              print*,"should not be here!!!!!!------------------------"
              !call dudt_dxdt_rk4_cartnorth_gr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
            case(.true.)
              print*,"should not be here!!!!!!------------------------"
              !call dudt_dxdt_rk4_cartnorth_gr_rr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
          end select
        case(p_rk6)
          select case (this%geometry%classical_rr)
            case(.false.)
              call dudt_dxdt_rk6_cartnorth_gr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
            case(.true.)
              print*,"should not be here!!!!!!------------------------"
              !call dudt_dxdt_rk6_cartnorth_gr_rr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
          end select
      end select
    else
      select case (this%push_type)
        case(p_rk4)
          select case (this%geometry%classical_rr)
            case(.false.)
              call dudt_dxdt_rk4_gr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
            case(.true.)
              call dudt_dxdt_rk4_gr_rr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
          end select
        case(p_rk6)
          select case (this%geometry%classical_rr)
            case(.false.)
              call dudt_dxdt_rk6_gr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
            case(.true.)
              call dudt_dxdt_rk6_gr_rr(this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit)
          end select
      end select
    endif


    ! local gamma factor (required for current deposition scheme at t = n+1/2)
    pp = ptrcur
    do i=1,np
      p2new(i) = this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2
      rgamma(i) = 1.0_p_k_part / sqrt( 1.0_p_k_part + 0.5_p_k_part * (p2init(i) + p2new(i)) )

      pp = pp + 1
    end do

    ! loop to get particle's grid quantities at t= t0 + dt
    pp = ptrcur
    do i=1,np
      dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
      xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

      dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
      xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

      ! Global x variation
      dxig(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, ixinit(1,i) )) / g%dr%f1(1, ixinit(1,i) ) )
      dxig(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, ixinit(2,i) )) / dt )

      pp = pp + 1
    end do

    ! loop to get initial ix
    pp = ptrcur
    do i=1,np
      this%x(1,pp) = xbufginit(1,i)
      this%x(2,pp) = xbufginit(2,i)

      this%ix(1,pp) = ixinit(1,i)
      this%ix(2,pp) = ixinit(2,i)

      pp = pp + 1
    end do

    call this % emit( emf, gdt, ptrcur, np )

    ! Deposit current
    select case ( g%spline_rule )
        case ( p_spline_rectangular )
          call dep_current_2d_srect_curved_gr( this, jay, dxig, xbufg, xbuf, &
                              ixinit, xbufginit, xbufinit,&
                              this%q(ptrcur:), rgamma, &
                              this%p(:,ptrcur:), pinit, &
                              np, gdt )
        case default !p_spline_integral
          call dep_current_2d_sint_curved_gr( this, jay, dxig, xbufg, xbuf, &
                              ixinit, xbufginit, xbufinit,&
                              this%q(ptrcur:), rgamma, &
                              this%p(:,ptrcur:), pinit, &
                              np, gdt )
    end select

    ! gets grid quantities at t= t0 + dt
    pp = ptrcur
    do i=1,np
      this%x(1,pp) = xbufg(1,i)
      this%x(2,pp) = xbufg(2,i)
      this%x(3,pp) = xbuf(3,i)
      this%ix(1,pp) = this%ix(1,pp) + dxig(1,i)
      this%ix(2,pp) = this%ix(2,pp) + dxig(2,i)

      pp = pp + 1
    end do


  enddo


end subroutine advance_deposit_2d_rk_gr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk4_gr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit )
!-----------------------------------------------------------------------------

  use m_species_push, only : ntrim
  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 4

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbufinit

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt2, gdt6

  real(p_double), dimension(rank, p_cache_size) :: xbufg
  integer, dimension(rank, p_cache_size) :: dxi
  real(p_double), dimension(rank + 1, p_cache_size) :: ptempinit
  real(p_double), dimension(order, p_p_dim, p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt2 = real( 0.5_p_double * gdt, p_double )
  gdt6 = real( gdt / 6._p_double, p_double )

  !-------------------------------------------
  ! t = t0 inputs: (xp0,wp0)
  !-------------------------------------------
  ! xp0 is the particle's local position (t=t0)
  ! wp0 is the particle's local momenta  (t=t0)

  ! loop to get particle's position
  ! xbuf(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbuf(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbuf(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbuf(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! save initial position vector xp0 and wp0
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)

    pp = pp + 1
  end do

  ! first step
  step = 1

  ! gets the values of kvec and lvec for step = 1
  ! inputs: (xp0,wp0)
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt2
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt2
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt2

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt2
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt2
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt2

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt/2
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! second step
  step = 2

  ! gets the values of kvec and lvec for step = 2
  ! inputs: (xp1,wp1)
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt2
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt2
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt2

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt2
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt2
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt2

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt/2
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! third step
  step = 3

  ! gets the values of kvec and lvec for step = 3
  ! inputs: (xp2,wp2)
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! fourth step
  step = 4

  ! gets the values of kvec and lvec for step = 4
  ! inputs: (xp3,wp3)
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt
  !-------------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i)+2._p_double*kvec(2,1,i)+&
                                 2._p_double*kvec(3,1,i)+kvec(4,1,i)) * gdt6
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i)+2._p_double*kvec(2,2,i)+&
                                 2._p_double*kvec(3,2,i)+kvec(4,2,i)) * gdt6
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i)+2._p_double*kvec(2,3,i)+&
                                 2._p_double*kvec(3,3,i)+kvec(4,3,i)) * gdt6

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i)+2._p_double*lvec(2,1,i)+&
                                     2._p_double*lvec(3,1,i)+lvec(4,1,i)) * gdt6
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i)+2._p_double*lvec(2,2,i)+&
                                     2._p_double*lvec(3,2,i)+lvec(4,2,i)) * gdt6
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i)+2._p_double*lvec(2,3,i)+&
                                     2._p_double*lvec(3,3,i)+lvec(4,3,i)) * gdt6

    pp = pp + 1
  end do


end subroutine dudt_dxdt_rk4_gr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk4_gr_rr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit )
!-----------------------------------------------------------------------------

  use m_species_push, only : ntrim
  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 4

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbufinit

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt2, gdt6

  real(p_double), dimension(rank , p_cache_size) :: xbufg
  integer, dimension(rank , p_cache_size) :: dxi
  real(p_double), dimension(rank + 1 , p_cache_size) :: ptempinit
  real(p_double), dimension(order ,p_p_dim,p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt2 = real( 0.5_p_double * gdt, p_double )
  gdt6 = real( gdt / 6._p_double, p_double )

  !-------------------------------------------
  ! t = t0 inputs: (xp0,wp0)
  !-------------------------------------------
  ! xp0 is the particle's local position (t=t0)
  ! wp0 is the particle's local momenta  (t=t0)

  ! loop to get particle's position
  ! xbuf(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbuf(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbuf(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbuf(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! save initial position vector xp0 and wp0
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)

    pp = pp + 1
  end do

  ! first step
  step = 1

  ! gets the values of kvec and lvec for step = 1
  ! inputs: (xp0,wp0)
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt2
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt2
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt2

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt2
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt2
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt2

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt/2
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! second step
  step = 2

  ! gets the values of kvec and lvec for step = 2
  ! inputs: (xp1,wp1)
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt2
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt2
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt2

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt2
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt2
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt2

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt/2
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! third step
  step = 3

  ! gets the values of kvec and lvec for step = 3
  ! inputs: (xp2,wp2)
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! fourth step
  step = 4

  ! gets the values of kvec and lvec for step = 4
  ! inputs: (xp3,wp3)
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt
  !-------------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i)+2._p_double*kvec(2,1,i)+&
                                 2._p_double*kvec(3,1,i)+kvec(4,1,i)) * gdt6
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i)+2._p_double*kvec(2,2,i)+&
                                 2._p_double*kvec(3,2,i)+kvec(4,2,i)) * gdt6
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i)+2._p_double*kvec(2,3,i)+&
                                 2._p_double*kvec(3,3,i)+kvec(4,3,i)) * gdt6

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i)+2._p_double*lvec(2,1,i)+&
                                     2._p_double*lvec(3,1,i)+lvec(4,1,i)) * gdt6
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i)+2._p_double*lvec(2,2,i)+&
                                     2._p_double*lvec(3,2,i)+lvec(4,2,i)) * gdt6
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i)+2._p_double*lvec(2,3,i)+&
                                     2._p_double*lvec(3,3,i)+lvec(4,3,i)) * gdt6

    pp = pp + 1
  end do


end subroutine dudt_dxdt_rk4_gr_rr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk4_cart_gr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufsinit )
!-----------------------------------------------------------------------------

  use m_species_push, only : ntrim
  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 4

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf,xbufsinit
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt2, gdt6
  real(p_double) :: r, t, p, xpos, ypos, zpos
  real(p_double) :: st, ct, sp, cp
  real(p_double) :: pr, pt, pphi, px, py, pz, alpha

  real(p_double), dimension(rank , p_cache_size) :: xbufg
  integer, dimension(rank , p_cache_size) :: dxi
  real(p_k_part), dimension(rank + 1 , p_cache_size) :: xbufinit, ptempinit, xbufs
  real(p_double), dimension(order ,p_p_dim,p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt2 = real( 0.5_p_double * gdt, p_double )
  gdt6 = real( gdt / 6._p_double, p_double )

  !-------------------------------------------
  ! t = t0 inputs: (xp0,wp0)
  !-------------------------------------------
  ! xp0 is the particle's local position (t=t0)
  ! wp0 is the particle's local momenta  (t=t0)

  ! loop to get particle's position
  ! xbuf(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbufs(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbufs(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbufs(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! From this point on xbuf and this%p are cartesian
  pp = ptrcur
  do i=1,np
    r = xbufs(1, i)
    t = xbufs(2, i)
    p = xbufs(3, i)

    xbufsinit(1, i) = r
    xbufsinit(2, i) = t
    xbufsinit(3, i) = p

    st = sin(t)
    ct = cos(t)
    sp = sin(p)
    cp = cos(p)

    xbuf(1, i) = r * st * cp
    xbuf(2, i) = r * st * sp
    xbuf(3, i) = r * ct

    pr = this%p(1,pp)
    pt = this%p(2,pp)
    pphi = this%p(3,pp)

    alpha = sqrt(1._p_double - g%rs/r)

    this%p(1,pp) = st * cp * pr * alpha + ct * cp * pt - sp * pphi
    this%p(2,pp) = st * sp * pr * alpha + ct * sp * pt + cp * pphi
    this%p(3,pp) = ct * pr * alpha - st * pt

    pp = pp + 1
  end do

  ! save initial position vector xp0 and wp0
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)

    pp = pp + 1
  end do

  ! first step
  step = 1

  ! gets the values of kvec and lvec for step = 1
  ! inputs: (xp0,wp0)
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt2
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt2
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt2

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt2
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt2
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt2

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt/2
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! second step
  step = 2

  ! gets the values of kvec and lvec for step = 2
  ! inputs: (xp1,wp1)
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt2
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt2
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt2

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt2
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt2
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt2

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt/2
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! third step
  step = 3

  ! gets the values of kvec and lvec for step = 3
  ! inputs: (xp2,wp2)
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! fourth step
  step = 4

  ! gets the values of kvec and lvec for step = 4
  ! inputs: (xp3,wp3)
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt
  !-------------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i)+2._p_double*kvec(2,1,i)+&
                                 2._p_double*kvec(3,1,i)+kvec(4,1,i)) * gdt6
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i)+2._p_double*kvec(2,2,i)+&
                                 2._p_double*kvec(3,2,i)+kvec(4,2,i)) * gdt6
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i)+2._p_double*kvec(2,3,i)+&
                                 2._p_double*kvec(3,3,i)+kvec(4,3,i)) * gdt6

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i)+2._p_double*lvec(2,1,i)+&
                                     2._p_double*lvec(3,1,i)+lvec(4,1,i)) * gdt6
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i)+2._p_double*lvec(2,2,i)+&
                                     2._p_double*lvec(3,2,i)+lvec(4,2,i)) * gdt6
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i)+2._p_double*lvec(2,3,i)+&
                                     2._p_double*lvec(3,3,i)+lvec(4,3,i)) * gdt6

    pp = pp + 1
  end do

  ! convert xbuf and this%p back to spherical
  pp = ptrcur
  do i=1,np
    xpos = xbuf(1,i)
    ypos = xbuf(2,i)
    zpos = xbuf(3,i)

    xbuf(1,i) = sqrt(xpos**2 + ypos**2 + zpos**2)
    xbuf(2,i) = acos( zpos / xbuf(1,i) )
    xbuf(3,i) = atan2(ypos,xpos)

    alpha = sqrt(1._p_double-g%rs/xbuf(1,i))

    ct = cos(xbuf(2,i))
    st = sin(xbuf(2,i))

    cp = cos(xbuf(3,i))
    sp = sin(xbuf(3,i))

    px = this%p(1,pp)
    py = this%p(2,pp)
    pz = this%p(3,pp)

    this%p(1,pp) = (st * cp * px + st * sp * py + ct * pz) / alpha
    this%p(2,pp) = ct * cp * px + ct * sp * py - st * pz
    this%p(3,pp) = - sp * px + cp* py

    pp = pp + 1
  end do



end subroutine dudt_dxdt_rk4_cart_gr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk4_cart_gr_rr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufsinit )
!-----------------------------------------------------------------------------

  use m_species_push, only : ntrim
  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 4

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf,xbufsinit
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt2, gdt6
  real(p_double) :: r, t, p, xpos, ypos, zpos
  real(p_double) :: st, ct, sp, cp
  real(p_double) :: pr, pt, pphi, px, py, pz, alpha

  real(p_double), dimension(rank , p_cache_size) :: xbufg
  integer, dimension(rank , p_cache_size) :: dxi
  real(p_k_part), dimension(rank + 1 , p_cache_size) :: xbufinit, ptempinit, xbufs
  real(p_double), dimension(order ,p_p_dim,p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt2 = real( 0.5_p_double * gdt, p_double )
  gdt6 = real( gdt / 6._p_double, p_double )

  !-------------------------------------------
  ! t = t0 inputs: (xp0,wp0)
  !-------------------------------------------
  ! xp0 is the particle's local position (t=t0)
  ! wp0 is the particle's local momenta  (t=t0)

  ! loop to get particle's position
  ! xbuf(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbufs(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbufs(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbufs(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! From this point on xbuf and this%p are cartesian
  pp = ptrcur
  do i=1,np
    r = xbufs(1, i)
    t = xbufs(2, i)
    p = xbufs(3, i)

    xbufsinit(1, i) = r
    xbufsinit(2, i) = t
    xbufsinit(3, i) = p

    st = sin(t)
    ct = cos(t)
    sp = sin(p)
    cp = cos(p)

    xbuf(1, i) = r * st * cp
    xbuf(2, i) = r * st * sp
    xbuf(3, i) = r * ct

    pr = this%p(1,pp)
    pt = this%p(2,pp)
    pphi = this%p(3,pp)

    alpha = sqrt(1._p_double - g%rs/r)

    this%p(1,pp) = st * cp * pr * alpha + ct * cp * pt - sp * pphi
    this%p(2,pp) = st * sp * pr * alpha + ct * sp * pt + cp * pphi
    this%p(3,pp) = ct * pr * alpha - st * pt

    pp = pp + 1
  end do

  ! save initial position vector xp0 and wp0
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)

    pp = pp + 1
  end do

  ! first step
  step = 1

  ! gets the values of kvec and lvec for step = 1
  ! inputs: (xp0,wp0)
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt2
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt2
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt2

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt2
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt2
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt2

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt/2
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! second step
  step = 2

  ! gets the values of kvec and lvec for step = 2
  ! inputs: (xp1,wp1)
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt2
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt2
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt2

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt2
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt2
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt2

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt/2
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! third step
  step = 3

  ! gets the values of kvec and lvec for step = 3
  ! inputs: (xp2,wp2)
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt
  !-------------------------------------------
  ! evaluates the next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt

    pp = pp + 1
  end do

  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets grid quantities at t= t0 + dt
  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  ! fourth step
  step = 4

  ! gets the values of kvec and lvec for step = 4
  ! inputs: (xp3,wp3)
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------------
  ! t = t0 + dt
  !-------------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i)+2._p_double*kvec(2,1,i)+&
                                 2._p_double*kvec(3,1,i)+kvec(4,1,i)) * gdt6
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i)+2._p_double*kvec(2,2,i)+&
                                 2._p_double*kvec(3,2,i)+kvec(4,2,i)) * gdt6
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i)+2._p_double*kvec(2,3,i)+&
                                 2._p_double*kvec(3,3,i)+kvec(4,3,i)) * gdt6

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i)+2._p_double*lvec(2,1,i)+&
                                     2._p_double*lvec(3,1,i)+lvec(4,1,i)) * gdt6
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i)+2._p_double*lvec(2,2,i)+&
                                     2._p_double*lvec(3,2,i)+lvec(4,2,i)) * gdt6
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i)+2._p_double*lvec(2,3,i)+&
                                     2._p_double*lvec(3,3,i)+lvec(4,3,i)) * gdt6

    pp = pp + 1
  end do

  ! convert xbuf and this%p back to spherical
  pp = ptrcur
  do i=1,np
    xpos = xbuf(1,i)
    ypos = xbuf(2,i)
    zpos = xbuf(3,i)

    xbuf(1,i) = sqrt(xpos**2 + ypos**2 + zpos**2)
    xbuf(2,i) = acos( zpos / xbuf(1,i) )
    xbuf(3,i) = atan2(ypos,xpos)

    alpha = sqrt(1._p_double-g%rs/xbuf(1,i))

    ct = cos(xbuf(2,i))
    st = sin(xbuf(2,i))

    cp = cos(xbuf(3,i))
    sp = sin(xbuf(3,i))

    px = this%p(1,pp)
    py = this%p(2,pp)
    pz = this%p(3,pp)

    this%p(1,pp) = (st * cp * px + st * sp * py + ct * pz) / alpha
    this%p(2,pp) = ct * cp * px + ct * sp * py - st * pz
    this%p(3,pp) = - sp * px + cp* py

    pp = pp + 1
  end do



end subroutine dudt_dxdt_rk4_cart_gr_rr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk6_gr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit )
!-----------------------------------------------------------------------------

  use m_species_push, only : ntrim
  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 6 + 1 !6th order requires 7 saved quantities

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbufinit

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt3, gdt12, gdt16, gdt8, gdt44, gdt120

  real(p_double), dimension(rank , p_cache_size) :: xbufg
  integer, dimension(rank , p_cache_size) :: dxi
  real(p_double), dimension(rank + 1 , p_cache_size) :: ptempinit
  real(p_double), dimension(order ,p_p_dim,p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt3 = real( gdt / 3._p_double , p_k_part )
  gdt12 = real( gdt / 12._p_double , p_k_part )
  gdt16 = real( gdt / 16._p_double , p_k_part )
  gdt8 = real( gdt / 8._p_double , p_k_part )
  gdt44 = real( gdt / 44._p_double , p_k_part )
  gdt120 = real( gdt / 120._p_double, p_k_part )

  ! xp0 is the particle's local position
  ! wp0 is the particle's local momenta

  ! loop to get particle's position
  ! xbuf(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbuf(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbuf(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbuf(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! save initial position vector xp0 and momenta
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)
    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 (xp0,wp0)
  !-------------------------------------
  step = 1
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt3

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 2
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + 2._p_double * kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + 2._p_double * kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + 2._p_double * kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + 2._p_double * lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + 2._p_double * lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + 2._p_double * lvec(step,3,i) * gdt3

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + 2dt/3
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 3
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i) + 4._p_double * kvec(2,1,i) - kvec(step,1,i)) * gdt12
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i) + 4._p_double * kvec(2,2,i) - kvec(step,2,i)) * gdt12
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i) + 4._p_double * kvec(2,3,i) - kvec(step,3,i)) * gdt12

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i) + 4._p_double * lvec(2,1,i) - lvec(step,1,i)) * gdt12
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i) + 4._p_double * lvec(2,2,i) - lvec(step,2,i)) * gdt12
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i) + 4._p_double * lvec(2,3,i) - lvec(step,3,i)) * gdt12

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 4
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp4,wp4)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (-kvec(1,1,i) + 18._p_double * kvec(2,1,i) -&
         3._p_double * kvec(3,1,i) - 6._p_double * kvec(step,1,i)) * gdt16
    xbuf(2,i) = xbufinit(2,i) + (-kvec(1,2,i) + 18._p_double * kvec(2,2,i) -&
         3._p_double * kvec(3,2,i) - 6._p_double * kvec(step,2,i)) * gdt16
    xbuf(3,i) = xbufinit(3,i) + (-kvec(1,3,i) + 18._p_double * kvec(2,3,i) -&
         3._p_double * kvec(3,3,i) - 6._p_double * kvec(step,3,i)) * gdt16

    this%p(1,pp) = ptempinit(1,i) + (-lvec(1,1,i) + 18._p_double * lvec(2,1,i) -&
         3._p_double * lvec(3,1,i) - 6._p_double * lvec(step,1,i)) * gdt16
    this%p(2,pp) = ptempinit(2,i) + (-lvec(1,2,i) + 18._p_double * lvec(2,2,i) -&
         3._p_double * lvec(3,2,i) - 6._p_double * lvec(step,2,i)) * gdt16
    this%p(3,pp) = ptempinit(3,i) + (-lvec(1,3,i) + 18._p_double * lvec(2,3,i) -&
         3._p_double * lvec(3,3,i) - 6._p_double * lvec(step,3,i)) * gdt16

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 5
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp5,wp5)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(2,1,i) - 3._p_double * kvec(3,1,i) -&
         6._p_double * kvec(4,1,i) + 4._p_double * kvec(step,1,i)) * gdt8
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(2,2,i) - 3._p_double * kvec(3,2,i) -&
         6._p_double * kvec(4,2,i) + 4._p_double * kvec(step,2,i)) * gdt8
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(2,3,i) - 3._p_double * kvec(3,3,i) -&
         6._p_double * kvec(4,3,i) + 4._p_double * kvec(step,3,i)) * gdt8

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(2,1,i) - 3._p_double * lvec(3,1,i) -&
         6._p_double * lvec(4,1,i) + 4._p_double * lvec(step,1,i)) * gdt8
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(2,2,i) - 3._p_double * lvec(3,2,i) -&
         6._p_double * lvec(4,2,i) + 4._p_double * lvec(step,2,i)) * gdt8
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(2,3,i) - 3._p_double * lvec(3,3,i) -&
         6._p_double * lvec(4,3,i) + 4._p_double * lvec(step,3,i)) * gdt8

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 6
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp6,wp6)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(1,1,i) - 36._p_double * kvec(2,1,i) +&
         63._p_double * kvec(3,1,i) + 72._p_double * kvec(4,1,i) - 64._p_double * kvec(step,1,i)) * gdt44
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(1,2,i) - 36._p_double * kvec(2,2,i) +&
         63._p_double * kvec(3,2,i) + 72._p_double * kvec(4,2,i) - 64._p_double * kvec(step,2,i)) * gdt44
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(1,3,i) - 36._p_double * kvec(2,3,i) +&
         63._p_double * kvec(3,3,i) + 72._p_double * kvec(4,3,i) - 64._p_double * kvec(step,3,i)) * gdt44

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(1,1,i) - 36._p_double * lvec(2,1,i) +&
         63._p_double * lvec(3,1,i) + 72._p_double * lvec(4,1,i) - 64._p_double * lvec(step,1,i)) * gdt44
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(1,2,i) - 36._p_double * lvec(2,2,i) +&
         63._p_double * lvec(3,2,i) + 72._p_double * lvec(4,2,i) - 64._p_double * lvec(step,2,i)) * gdt44
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(1,3,i) - 36._p_double * lvec(2,3,i) +&
         63._p_double * lvec(3,3,i) + 72._p_double * lvec(4,3,i) - 64._p_double * lvec(step,3,i)) * gdt44

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 7
  call get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------
  ! t = t0 + dt (xp,wp)
  !-------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (11._p_double * kvec(1,1,i) + 81._p_double * (kvec(3,1,i) + kvec(4,1,i)) -&
          32._p_double * (kvec(5,1,i) + kvec(6,1,i)) + 11._p_double * kvec(step,1,i)) * gdt120
    xbuf(2,i) = xbufinit(2,i) + (11._p_double * kvec(1,2,i) + 81._p_double * (kvec(3,2,i) + kvec(4,2,i)) -&
          32._p_double * (kvec(5,2,i) + kvec(6,2,i)) + 11._p_double * kvec(step,2,i)) * gdt120
    xbuf(3,i) = xbufinit(3,i) + (11._p_double * kvec(1,3,i) + 81._p_double * (kvec(3,3,i) + kvec(4,3,i)) -&
          32._p_double * (kvec(5,3,i) + kvec(6,3,i)) + 11._p_double * kvec(step,3,i)) * gdt120

    this%p(1,pp) = ptempinit(1,i) + (11._p_double * lvec(1,1,i) + 81._p_double * (lvec(3,1,i) + lvec(4,1,i)) -&
          32._p_double * (lvec(5,1,i) + lvec(6,1,i)) + 11._p_double * lvec(step,1,i)) * gdt120
    this%p(2,pp) = ptempinit(2,i) + (11._p_double * lvec(1,2,i) + 81._p_double * (lvec(3,2,i) + lvec(4,2,i)) -&
          32._p_double * (lvec(5,2,i) + lvec(6,2,i)) + 11._p_double * lvec(step,2,i)) * gdt120
    this%p(3,pp) = ptempinit(3,i) + (11._p_double * lvec(1,3,i) + 81._p_double * (lvec(3,3,i) + lvec(4,3,i)) -&
          32._p_double * (lvec(5,3,i) + lvec(6,3,i)) + 11._p_double * lvec(step,3,i)) * gdt120

    pp = pp + 1
  end do

end subroutine dudt_dxdt_rk6_gr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk6_gr_rr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufinit )
!-----------------------------------------------------------------------------

  use m_species_push, only : ntrim
  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 6 + 1 !6th order requires 7 saved quantities

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbufinit

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt3, gdt12, gdt16, gdt8, gdt44, gdt120

  real(p_double), dimension(rank , p_cache_size) :: xbufg
  integer, dimension(rank , p_cache_size) :: dxi
  real(p_double), dimension(rank + 1 , p_cache_size) :: ptempinit
  real(p_double), dimension(order ,p_p_dim,p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt3 = real( gdt / 3._p_double , p_k_part )
  gdt12 = real( gdt / 12._p_double , p_k_part )
  gdt16 = real( gdt / 16._p_double , p_k_part )
  gdt8 = real( gdt / 8._p_double , p_k_part )
  gdt44 = real( gdt / 44._p_double , p_k_part )
  gdt120 = real( gdt / 120._p_double, p_k_part )

  ! xp0 is the particle's local position
  ! wp0 is the particle's local momenta

  ! loop to get particle's position
  ! xbuf(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbuf(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbuf(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbuf(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! save initial position vector xp0 and momenta
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)
    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 (xp0,wp0)
  !-------------------------------------
  step = 1
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt3

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 2
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + 2._p_double * kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + 2._p_double * kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + 2._p_double * kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + 2._p_double * lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + 2._p_double * lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + 2._p_double * lvec(step,3,i) * gdt3

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + 2dt/3
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 3
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i) + 4._p_double * kvec(2,1,i) - kvec(step,1,i)) * gdt12
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i) + 4._p_double * kvec(2,2,i) - kvec(step,2,i)) * gdt12
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i) + 4._p_double * kvec(2,3,i) - kvec(step,3,i)) * gdt12

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i) + 4._p_double * lvec(2,1,i) - lvec(step,1,i)) * gdt12
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i) + 4._p_double * lvec(2,2,i) - lvec(step,2,i)) * gdt12
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i) + 4._p_double * lvec(2,3,i) - lvec(step,3,i)) * gdt12

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 4
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp4,wp4)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (-kvec(1,1,i) + 18._p_double * kvec(2,1,i) -&
         3._p_double * kvec(3,1,i) - 6._p_double * kvec(step,1,i)) * gdt16
    xbuf(2,i) = xbufinit(2,i) + (-kvec(1,2,i) + 18._p_double * kvec(2,2,i) -&
         3._p_double * kvec(3,2,i) - 6._p_double * kvec(step,2,i)) * gdt16
    xbuf(3,i) = xbufinit(3,i) + (-kvec(1,3,i) + 18._p_double * kvec(2,3,i) -&
         3._p_double * kvec(3,3,i) - 6._p_double * kvec(step,3,i)) * gdt16

    this%p(1,pp) = ptempinit(1,i) + (-lvec(1,1,i) + 18._p_double * lvec(2,1,i) -&
         3._p_double * lvec(3,1,i) - 6._p_double * lvec(step,1,i)) * gdt16
    this%p(2,pp) = ptempinit(2,i) + (-lvec(1,2,i) + 18._p_double * lvec(2,2,i) -&
         3._p_double * lvec(3,2,i) - 6._p_double * lvec(step,2,i)) * gdt16
    this%p(3,pp) = ptempinit(3,i) + (-lvec(1,3,i) + 18._p_double * lvec(2,3,i) -&
         3._p_double * lvec(3,3,i) - 6._p_double * lvec(step,3,i)) * gdt16

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 5
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp5,wp5)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(2,1,i) - 3._p_double * kvec(3,1,i) -&
         6._p_double * kvec(4,1,i) + 4._p_double * kvec(step,1,i)) * gdt8
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(2,2,i) - 3._p_double * kvec(3,2,i) -&
         6._p_double * kvec(4,2,i) + 4._p_double * kvec(step,2,i)) * gdt8
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(2,3,i) - 3._p_double * kvec(3,3,i) -&
         6._p_double * kvec(4,3,i) + 4._p_double * kvec(step,3,i)) * gdt8

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(2,1,i) - 3._p_double * lvec(3,1,i) -&
         6._p_double * lvec(4,1,i) + 4._p_double * lvec(step,1,i)) * gdt8
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(2,2,i) - 3._p_double * lvec(3,2,i) -&
         6._p_double * lvec(4,2,i) + 4._p_double * lvec(step,2,i)) * gdt8
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(2,3,i) - 3._p_double * lvec(3,3,i) -&
         6._p_double * lvec(4,3,i) + 4._p_double * lvec(step,3,i)) * gdt8

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 6
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp6,wp6)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(1,1,i) - 36._p_double * kvec(2,1,i) +&
         63._p_double * kvec(3,1,i) + 72._p_double * kvec(4,1,i) - 64._p_double * kvec(step,1,i)) * gdt44
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(1,2,i) - 36._p_double * kvec(2,2,i) +&
         63._p_double * kvec(3,2,i) + 72._p_double * kvec(4,2,i) - 64._p_double * kvec(step,2,i)) * gdt44
    xbuf(3,i) = xbufinit(3,i) + (9 * kvec(1,3,i) - 36 * kvec(2,3,i) +&
         63._p_double * kvec(3,3,i) + 72._p_double * kvec(4,3,i) - 64._p_double * kvec(step,3,i)) * gdt44

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(1,1,i) - 36._p_double * lvec(2,1,i) +&
         63._p_double * lvec(3,1,i) + 72._p_double * lvec(4,1,i) - 64._p_double * lvec(step,1,i)) * gdt44
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(1,2,i) - 36._p_double * lvec(2,2,i) +&
         63._p_double * lvec(3,2,i) + 72._p_double * lvec(4,2,i) - 64._p_double * lvec(step,2,i)) * gdt44
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(1,3,i) - 36._p_double * lvec(2,3,i) +&
         63._p_double * lvec(3,3,i) + 72._p_double * lvec(4,3,i) - 64._p_double * lvec(step,3,i)) * gdt44

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    dxi(1,i) = ntrim( (xbuf(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbuf(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbuf(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbuf(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 7
  call get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------
  ! t = t0 + dt (xp,wp)
  !-------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (11._p_double * kvec(1,1,i) + 81._p_double * (kvec(3,1,i) + kvec(4,1,i)) -&
          32._p_double * (kvec(5,1,i) + kvec(6,1,i)) + 11._p_double * kvec(step,1,i)) * gdt120
    xbuf(2,i) = xbufinit(2,i) + (11._p_double * kvec(1,2,i) + 81._p_double * (kvec(3,2,i) + kvec(4,2,i)) -&
          32._p_double * (kvec(5,2,i) + kvec(6,2,i)) + 11._p_double * kvec(step,2,i)) * gdt120
    xbuf(3,i) = xbufinit(3,i) + (11._p_double * kvec(1,3,i) + 81._p_double * (kvec(3,3,i) + kvec(4,3,i)) -&
          32._p_double * (kvec(5,3,i) + kvec(6,3,i)) + 11._p_double * kvec(step,3,i)) * gdt120

    this%p(1,pp) = ptempinit(1,i) + (11._p_double * lvec(1,1,i) + 81._p_double * (lvec(3,1,i) + lvec(4,1,i)) -&
          32._p_double * (lvec(5,1,i) + lvec(6,1,i)) + 11._p_double * lvec(step,1,i)) * gdt120
    this%p(2,pp) = ptempinit(2,i) + (11._p_double * lvec(1,2,i) + 81._p_double * (lvec(3,2,i) + lvec(4,2,i)) -&
          32._p_double * (lvec(5,2,i) + lvec(6,2,i)) + 11._p_double * lvec(step,2,i)) * gdt120
    this%p(3,pp) = ptempinit(3,i) + (11._p_double * lvec(1,3,i) + 81._p_double * (lvec(3,3,i) + lvec(4,3,i)) -&
          32._p_double * (lvec(5,3,i) + lvec(6,3,i)) + 11._p_double * lvec(step,3,i)) * gdt120

    pp = pp + 1
  end do

end subroutine dudt_dxdt_rk6_gr_rr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk6_cart_gr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufsinit )
!-----------------------------------------------------------------------------

  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 6 + 1 !6th order requires 7 saved quantities

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf, xbufsinit
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt3, gdt12, gdt16, gdt8, gdt44, gdt120
  real(p_double) :: r, t, p, xpos, ypos, zpos
  real(p_double) :: st, ct, sp, cp
  real(p_double) :: pr, pt, pphi, px, py, pz, alpha

  real(p_double), dimension(rank , p_cache_size) :: xbufg
  real(p_k_part), dimension(rank + 1 , p_cache_size) :: xbufinit, ptempinit, xbufs
  real(p_double), dimension(order ,p_p_dim,p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt
  real(p_double) :: rho, auxr, auxrho

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt3 = real( gdt / 3._p_double , p_k_part )
  gdt12 = real( gdt / 12._p_double , p_k_part )
  gdt16 = real( gdt / 16._p_double , p_k_part )
  gdt8 = real( gdt / 8._p_double , p_k_part )
  gdt44 = real( gdt / 44._p_double , p_k_part )
  gdt120 = real( gdt / 120._p_double, p_k_part )

  ! xp0 is the particle's local position
  ! wp0 is the particle's local momenta

  ! loop to get particle's position
  ! xbufs(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbufs(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbufs(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbufs(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! From this point on xbuf and this%p are cartesian
  pp = ptrcur
  do i=1,np
    r = xbufs(1, i)
    t = xbufs(2, i)
    p = xbufs(3, i)

    xbufsinit(1, i) = r
    xbufsinit(2, i) = t
    xbufsinit(3, i) = p

    st = sin(t)
    ct = cos(t)
    sp = sin(p)
    cp = cos(p)

    xbuf(1, i) = r * st * cp
    xbuf(2, i) = r * st * sp
    xbuf(3, i) = r * ct

    pr = this%p(1,pp)
    pt = this%p(2,pp)
    pphi = this%p(3,pp)

    alpha = sqrt(1._p_double - g%rs/r)

    this%p(1,pp) = st * cp * pr * alpha + ct * cp * pt - sp * pphi
    this%p(2,pp) = st * sp * pr * alpha + ct * sp * pt + cp * pphi
    this%p(3,pp) = ct * pr * alpha - st * pt

    pp = pp + 1
  end do

  ! save initial position vector xp0 and momenta
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 (xp0,wp0)
  !-------------------------------------
  step = 1
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt3

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 2
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + 2._p_double * kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + 2._p_double * kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + 2._p_double * kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + 2._p_double * lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + 2._p_double * lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + 2._p_double * lvec(step,3,i) * gdt3

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + 2dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 3
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i) + 4._p_double * kvec(2,1,i) - kvec(step,1,i)) * gdt12
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i) + 4._p_double * kvec(2,2,i) - kvec(step,2,i)) * gdt12
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i) + 4._p_double * kvec(2,3,i) - kvec(step,3,i)) * gdt12

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i) + 4._p_double * lvec(2,1,i) - lvec(step,1,i)) * gdt12
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i) + 4._p_double * lvec(2,2,i) - lvec(step,2,i)) * gdt12
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i) + 4._p_double * lvec(2,3,i) - lvec(step,3,i)) * gdt12

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 4
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp4,wp4)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (-kvec(1,1,i) + 18._p_double * kvec(2,1,i) -&
         3._p_double * kvec(3,1,i) - 6._p_double * kvec(step,1,i)) * gdt16
    xbuf(2,i) = xbufinit(2,i) + (-kvec(1,2,i) + 18._p_double * kvec(2,2,i) -&
         3._p_double * kvec(3,2,i) - 6._p_double * kvec(step,2,i)) * gdt16
    xbuf(3,i) = xbufinit(3,i) + (-kvec(1,3,i) + 18._p_double * kvec(2,3,i) -&
         3._p_double * kvec(3,3,i) - 6._p_double * kvec(step,3,i)) * gdt16

    this%p(1,pp) = ptempinit(1,i) + (-lvec(1,1,i) + 18._p_double * lvec(2,1,i) -&
         3._p_double * lvec(3,1,i) - 6._p_double * lvec(step,1,i)) * gdt16
    this%p(2,pp) = ptempinit(2,i) + (-lvec(1,2,i) + 18._p_double * lvec(2,2,i) -&
         3._p_double * lvec(3,2,i) - 6._p_double * lvec(step,2,i)) * gdt16
    this%p(3,pp) = ptempinit(3,i) + (-lvec(1,3,i) + 18._p_double * lvec(2,3,i) -&
         3._p_double * lvec(3,3,i) - 6._p_double * lvec(step,3,i)) * gdt16

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 5
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp5,wp5)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(2,1,i) - 3._p_double * kvec(3,1,i) -&
         6._p_double * kvec(4,1,i) + 4._p_double * kvec(step,1,i)) * gdt8
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(2,2,i) - 3._p_double * kvec(3,2,i) -&
         6._p_double * kvec(4,2,i) + 4._p_double * kvec(step,2,i)) * gdt8
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(2,3,i) - 3._p_double * kvec(3,3,i) -&
         6._p_double * kvec(4,3,i) + 4._p_double * kvec(step,3,i)) * gdt8

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(2,1,i) - 3._p_double * lvec(3,1,i) -&
         6._p_double * lvec(4,1,i) + 4._p_double * lvec(step,1,i)) * gdt8
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(2,2,i) - 3._p_double * lvec(3,2,i) -&
         6._p_double * lvec(4,2,i) + 4._p_double * lvec(step,2,i)) * gdt8
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(2,3,i) - 3._p_double * lvec(3,3,i) -&
         6._p_double * lvec(4,3,i) + 4._p_double * lvec(step,3,i)) * gdt8

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 6
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp6,wp6)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(1,1,i) - 36._p_double * kvec(2,1,i) +&
         63._p_double * kvec(3,1,i) + 72._p_double * kvec(4,1,i) - 64._p_double * kvec(step,1,i)) * gdt44
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(1,2,i) - 36._p_double * kvec(2,2,i) +&
         63._p_double * kvec(3,2,i) + 72._p_double * kvec(4,2,i) - 64._p_double * kvec(step,2,i)) * gdt44
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(1,3,i) - 36._p_double * kvec(2,3,i) +&
         63._p_double * kvec(3,3,i) + 72._p_double * kvec(4,3,i) - 64._p_double * kvec(step,3,i)) * gdt44

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(1,1,i) - 36._p_double * lvec(2,1,i) +&
         63._p_double * lvec(3,1,i) + 72._p_double * lvec(4,1,i) - 64._p_double * lvec(step,1,i)) * gdt44
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(1,2,i) - 36._p_double * lvec(2,2,i) +&
         63._p_double * lvec(3,2,i) + 72._p_double * lvec(4,2,i) - 64._p_double * lvec(step,2,i)) * gdt44
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(1,3,i) - 36._p_double * lvec(2,3,i) +&
         63._p_double * lvec(3,3,i) + 72._p_double * lvec(4,3,i) - 64._p_double * lvec(step,3,i)) * gdt44

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 7
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------
  ! t = t0 + dt (xp,wp)
  !-------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (11._p_double * kvec(1,1,i) + 81._p_double * (kvec(3,1,i) + kvec(4,1,i)) -&
          32._p_double * (kvec(5,1,i) + kvec(6,1,i)) + 11._p_double * kvec(step,1,i)) * gdt120
    xbuf(2,i) = xbufinit(2,i) + (11._p_double * kvec(1,2,i) + 81._p_double * (kvec(3,2,i) + kvec(4,2,i)) -&
          32._p_double * (kvec(5,2,i) + kvec(6,2,i)) + 11._p_double * kvec(step,2,i)) * gdt120
    xbuf(3,i) = xbufinit(3,i) + (11._p_double * kvec(1,3,i) + 81._p_double * (kvec(3,3,i) + kvec(4,3,i)) -&
          32._p_double * (kvec(5,3,i) + kvec(6,3,i)) + 11._p_double * kvec(step,3,i)) * gdt120

    this%p(1,pp) = ptempinit(1,i) + (11._p_double * lvec(1,1,i) + 81._p_double * (lvec(3,1,i) + lvec(4,1,i)) -&
          32._p_double * (lvec(5,1,i) + lvec(6,1,i)) + 11._p_double * lvec(step,1,i)) * gdt120
    this%p(2,pp) = ptempinit(2,i) + (11._p_double * lvec(1,2,i) + 81._p_double * (lvec(3,2,i) + lvec(4,2,i)) -&
          32._p_double * (lvec(5,2,i) + lvec(6,2,i)) + 11._p_double * lvec(step,2,i)) * gdt120
    this%p(3,pp) = ptempinit(3,i) + (11._p_double * lvec(1,3,i) + 81._p_double * (lvec(3,3,i) + lvec(4,3,i)) -&
          32._p_double * (lvec(5,3,i) + lvec(6,3,i)) + 11._p_double * lvec(step,3,i)) * gdt120

    pp = pp + 1
  end do

  ! convert xbuf and this%p back to spherical
  pp = ptrcur
  do i=1,np
    xpos = xbuf(1,i)
    ypos = xbuf(2,i)
    zpos = xbuf(3,i)
    rho = sqrt(xpos * xpos + ypos * ypos)

    xbuf(1,i) = sqrt(xpos**2 + ypos**2 + zpos**2)

    ! evaluate sin and cos of theta
    auxr = 1.0_p_double / xbuf(1,i)
    st = rho * auxr
    ct = zpos * auxr

    ! evaluate sin and cos of phi
    auxrho = 1.0_p_double / rho
    sp = ypos * auxrho
    cp = xpos * auxrho

    xbuf(2,i) = acos(zpos * auxr)
    xbuf(3,i) = atan2(ypos,xpos)

    alpha = sqrt(1._p_double - g%rs * auxr)

    px = this%p(1,pp)
    py = this%p(2,pp)
    pz = this%p(3,pp)

    this%p(1,pp) = (st * cp * px + st * sp * py + ct * pz) / alpha
    this%p(2,pp) = ct * cp * px + ct * sp * py - st * pz
    this%p(3,pp) = - sp * px + cp* py

    pp = pp + 1
  end do

end subroutine dudt_dxdt_rk6_cart_gr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk6_cartnorth_gr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufsinit )
!-----------------------------------------------------------------------------

  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 6 + 1 !6th order requires 7 saved quantities

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf, xbufsinit
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt3, gdt12, gdt16, gdt8, gdt44, gdt120
  real(p_double) :: r, t, p, xpos, ypos, zpos
  real(p_double) :: st, ct, sp, cp
  real(p_double) :: pr, pt, pphi, px, py, pz, alpha

  real(p_double), dimension(rank , p_cache_size) :: xbufg
  real(p_k_part), dimension(rank + 1 , p_cache_size) :: xbufinit, ptempinit, xbufs
  real(p_double), dimension(order ,p_p_dim,p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt
  real(p_double) :: rho, auxr, auxrho
  real(p_double) :: xc, xc12, xc32, xc52, xc72, xc92
  real(p_double) :: t2, t3, t4, t5, t6, t7, t8, t9, t10
  real(p_double) :: p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15
  real(p_double) :: x, xa, xa2, xa3, xa5, xa7, xa9, xa11, xa13, xa15
  real(p_double) :: xac, xac2, xac3, xac5, xac6, xac7, xac9
  real(p_double) :: xac10, xac11, xac13, xac14, xac15
  real(p_double) :: pi2, pi4, atanxa, signx, signypos, signdif

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt3 = real( gdt / 3._p_double , p_k_part )
  gdt12 = real( gdt / 12._p_double , p_k_part )
  gdt16 = real( gdt / 16._p_double , p_k_part )
  gdt8 = real( gdt / 8._p_double , p_k_part )
  gdt44 = real( gdt / 44._p_double , p_k_part )
  gdt120 = real( gdt / 120._p_double, p_k_part )

  ! used to evaluate the polynomial approx of trig functions
  pi2 = pi_2
  pi4 = 0.500000000000000_p_double * pi2

  ! xp0 is the particle's local position
  ! wp0 is the particle's local momenta

  ! loop to get particle's position
  ! xbufs(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbufs(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbufs(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbufs(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! From this point on xbuf and this%p are cartesian
  pp = ptrcur
  do i=1,np
    r = xbufs(1, i)
    t = xbufs(2, i)
    p = xbufs(3, i)

    xbufsinit(1, i) = r
    xbufsinit(2, i) = t
    xbufsinit(3, i) = p

    ! taylor expansion of sin and cos around theta ~ 0
    t2 = t * t
    t3 = t2 * t
    t4 = t3 * t
    t5 = t4 * t
    t6 = t5 * t
    t7 = t6 * t
    t8 = t7 * t
    t9 = t8 * t
    t10 = t9 * t

    st = t - 0.166666666666667_p_double * t3 + 0.008333333333333_p_double * t5 -&
             0.000198412698413_p_double * t7 + 0.000002755731922_p_double * t9

    ct = 1.0_p_double - 0.500000000000000_p_double * t2 +&
         0.041666666666667_p_double * t4 - 0.001388888888889_p_double * t6 +&
         0.000024801587302_p_double * t8 - 0.000000275573192_p_double * t10

    ! taylor expansion of sin and cos around phi ~ 0
    p2 = p * p
    p3 = p2 * p
    p4 = p3 * p
    p5 = p4 * p
    p6 = p5 * p
    p7 = p6 * p
    p8 = p7 * p
    p9 = p8 * p
    p10 = p9 * p
    p11 = p10 * p
    p12 = p11 * p
    p13 = p12 * p
    p14 = p13 * p
    p15 = p14 * p

    sp = p - 0.166666666666667_p_double * p3 +&
             0.008333333333333_p_double * p5 - 0.000198412698413_p_double * p7 +&
             0.000002755731922_p_double * p9 - 0.000000025052108_p_double * p11 +&
             0.000000000160590_p_double * p13 - 0.000000000000765_p_double * p15

    cp = 1.0_p_double - 0.500000000000000_p_double * p2 +&
                        0.041666666666667_p_double * p4 - 0.001388888888889_p_double * p6 +&
                        0.000024801587302_p_double * p8 - 0.000000275573192_p_double * p10 +&
                        0.000000002087676_p_double * p12 - 0.000000000011471_p_double * p14

    xbuf(1, i) = r * st * cp
    xbuf(2, i) = r * st * sp
    xbuf(3, i) = r * ct

    pr = this%p(1,pp)
    pt = this%p(2,pp)
    pphi = this%p(3,pp)

    alpha = sqrt(1._p_double - g%rs/r)

    this%p(1,pp) = st * cp * pr * alpha + ct * cp * pt - sp * pphi
    this%p(2,pp) = st * sp * pr * alpha + ct * sp * pt + cp * pphi
    this%p(3,pp) = ct * pr * alpha - st * pt

    pp = pp + 1
  end do

  ! save initial position vector xp0 and momenta
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 (xp0,wp0)
  !-------------------------------------
  step = 1
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt3

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartnorth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 2
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + 2._p_double * kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + 2._p_double * kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + 2._p_double * kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + 2._p_double * lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + 2._p_double * lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + 2._p_double * lvec(step,3,i) * gdt3

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + 2dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartnorth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 3
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i) + 4._p_double * kvec(2,1,i) - kvec(step,1,i)) * gdt12
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i) + 4._p_double * kvec(2,2,i) - kvec(step,2,i)) * gdt12
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i) + 4._p_double * kvec(2,3,i) - kvec(step,3,i)) * gdt12

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i) + 4._p_double * lvec(2,1,i) - lvec(step,1,i)) * gdt12
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i) + 4._p_double * lvec(2,2,i) - lvec(step,2,i)) * gdt12
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i) + 4._p_double * lvec(2,3,i) - lvec(step,3,i)) * gdt12

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartnorth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 4
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp4,wp4)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (-kvec(1,1,i) + 18._p_double * kvec(2,1,i) -&
         3._p_double * kvec(3,1,i) - 6._p_double * kvec(step,1,i)) * gdt16
    xbuf(2,i) = xbufinit(2,i) + (-kvec(1,2,i) + 18._p_double * kvec(2,2,i) -&
         3._p_double * kvec(3,2,i) - 6._p_double * kvec(step,2,i)) * gdt16
    xbuf(3,i) = xbufinit(3,i) + (-kvec(1,3,i) + 18._p_double * kvec(2,3,i) -&
         3._p_double * kvec(3,3,i) - 6._p_double * kvec(step,3,i)) * gdt16

    this%p(1,pp) = ptempinit(1,i) + (-lvec(1,1,i) + 18._p_double * lvec(2,1,i) -&
         3._p_double * lvec(3,1,i) - 6._p_double * lvec(step,1,i)) * gdt16
    this%p(2,pp) = ptempinit(2,i) + (-lvec(1,2,i) + 18._p_double * lvec(2,2,i) -&
         3._p_double * lvec(3,2,i) - 6._p_double * lvec(step,2,i)) * gdt16
    this%p(3,pp) = ptempinit(3,i) + (-lvec(1,3,i) + 18._p_double * lvec(2,3,i) -&
         3._p_double * lvec(3,3,i) - 6._p_double * lvec(step,3,i)) * gdt16

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartnorth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 5
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp5,wp5)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(2,1,i) - 3._p_double * kvec(3,1,i) -&
         6._p_double * kvec(4,1,i) + 4._p_double * kvec(step,1,i)) * gdt8
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(2,2,i) - 3._p_double * kvec(3,2,i) -&
         6._p_double * kvec(4,2,i) + 4._p_double * kvec(step,2,i)) * gdt8
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(2,3,i) - 3._p_double * kvec(3,3,i) -&
         6._p_double * kvec(4,3,i) + 4._p_double * kvec(step,3,i)) * gdt8

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(2,1,i) - 3._p_double * lvec(3,1,i) -&
         6._p_double * lvec(4,1,i) + 4._p_double * lvec(step,1,i)) * gdt8
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(2,2,i) - 3._p_double * lvec(3,2,i) -&
         6._p_double * lvec(4,2,i) + 4._p_double * lvec(step,2,i)) * gdt8
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(2,3,i) - 3._p_double * lvec(3,3,i) -&
         6._p_double * lvec(4,3,i) + 4._p_double * lvec(step,3,i)) * gdt8

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartnorth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 6
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp6,wp6)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(1,1,i) - 36._p_double * kvec(2,1,i) +&
         63._p_double * kvec(3,1,i) + 72._p_double * kvec(4,1,i) - 64._p_double * kvec(step,1,i)) * gdt44
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(1,2,i) - 36._p_double * kvec(2,2,i) +&
         63._p_double * kvec(3,2,i) + 72._p_double * kvec(4,2,i) - 64._p_double * kvec(step,2,i)) * gdt44
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(1,3,i) - 36._p_double * kvec(2,3,i) +&
         63._p_double * kvec(3,3,i) + 72._p_double * kvec(4,3,i) - 64._p_double * kvec(step,3,i)) * gdt44

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(1,1,i) - 36._p_double * lvec(2,1,i) +&
         63._p_double * lvec(3,1,i) + 72._p_double * lvec(4,1,i) - 64._p_double * lvec(step,1,i)) * gdt44
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(1,2,i) - 36._p_double * lvec(2,2,i) +&
         63._p_double * lvec(3,2,i) + 72._p_double * lvec(4,2,i) - 64._p_double * lvec(step,2,i)) * gdt44
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(1,3,i) - 36._p_double * lvec(2,3,i) +&
         63._p_double * lvec(3,3,i) + 72._p_double * lvec(4,3,i) - 64._p_double * lvec(step,3,i)) * gdt44

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartnorth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 7
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------
  ! t = t0 + dt (xp,wp)
  !-------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (11._p_double * kvec(1,1,i) + 81._p_double * (kvec(3,1,i) + kvec(4,1,i)) -&
          32._p_double * (kvec(5,1,i) + kvec(6,1,i)) + 11._p_double * kvec(step,1,i)) * gdt120
    xbuf(2,i) = xbufinit(2,i) + (11._p_double * kvec(1,2,i) + 81._p_double * (kvec(3,2,i) + kvec(4,2,i)) -&
          32._p_double * (kvec(5,2,i) + kvec(6,2,i)) + 11._p_double * kvec(step,2,i)) * gdt120
    xbuf(3,i) = xbufinit(3,i) + (11._p_double * kvec(1,3,i) + 81._p_double * (kvec(3,3,i) + kvec(4,3,i)) -&
          32._p_double * (kvec(5,3,i) + kvec(6,3,i)) + 11._p_double * kvec(step,3,i)) * gdt120

    this%p(1,pp) = ptempinit(1,i) + (11._p_double * lvec(1,1,i) + 81._p_double * (lvec(3,1,i) + lvec(4,1,i)) -&
          32._p_double * (lvec(5,1,i) + lvec(6,1,i)) + 11._p_double * lvec(step,1,i)) * gdt120
    this%p(2,pp) = ptempinit(2,i) + (11._p_double * lvec(1,2,i) + 81._p_double * (lvec(3,2,i) + lvec(4,2,i)) -&
          32._p_double * (lvec(5,2,i) + lvec(6,2,i)) + 11._p_double * lvec(step,2,i)) * gdt120
    this%p(3,pp) = ptempinit(3,i) + (11._p_double * lvec(1,3,i) + 81._p_double * (lvec(3,3,i) + lvec(4,3,i)) -&
          32._p_double * (lvec(5,3,i) + lvec(6,3,i)) + 11._p_double * lvec(step,3,i)) * gdt120

    pp = pp + 1
  end do

  ! convert xbuf and this%p back to spherical
  pp = ptrcur
  do i=1,np
    xpos = xbuf(1,i)
    ypos = xbuf(2,i)
    zpos = xbuf(3,i)
    rho = sqrt(xpos * xpos + ypos * ypos)

    xbuf(1,i) = sqrt(xpos**2 + ypos**2 + zpos**2)

    ! evaluate sin and cos of theta
    auxr = 1.0_p_double / xbuf(1,i)
    st = rho * auxr
    ct = zpos * auxr

    ! evaluate sin and cos of phi
    auxrho = 1.0_p_double / rho
    sp = ypos * auxrho
    cp = xpos * auxrho

    ! evaluate the arccos of x around x~1
    xc = 1.0_p_double - zpos * auxr
    xc12 = sqrt(xc)
    xc32 = xc12 * xc
    xc52 = xc32 * xc
    xc72 = xc52 * xc
    xc92 = xc72 * xc

    ! theta = acos(z/r) , z>0
    xbuf(2,i) = 1.414213562373095_p_double * xc12 + 0.117851130197758_p_double * xc32 +&
                0.026516504294496_p_double * xc52 + 0.007891816754314_p_double * xc72 +&
                0.002685409867787_p_double * xc92

    ! evaluate atan2 of x
    ! argument that enters atan(x)
    x = ypos / xpos
    xa = abs( x )

    ! loop that defines the three expansion points (0,1,+inf)
    if (xa .le. 0.418597947318945_p_double) then
      !here goes the expansion around x ~ 0
      xa2 = xa * xa
      xa3 = xa2 * xa
      xa5 = xa3 * xa2
      xa7 = xa5 * xa2
      xa9 = xa7 * xa2
      xa11 = xa9 * xa2
      xa13 = xa11 * xa2
      xa15 = xa13 * xa2

      ! evaluate atan(abs(x)) around x ~ 0
      atanxa = xa - 0.333333333333333_p_double * xa3 +&
                    0.200000000000000_p_double * xa5 - 0.142857142857143_p_double * xa7 +&
                    0.111111111111111_p_double * xa9 - 0.090909090909091_p_double * xa11 +&
                    0.076923076923077_p_double * xa13 - 0.066666666666667_p_double * xa15

    else if (xa .le. 1.817313073051244_p_double) then
      !here goes the expansion around x ~ 1
      xac = xa - 1.0_p_double
      xac2 = xac * xac
      xac3 = xac2 * xac
      xac5 = xac3 * xac2
      xac7 = xac5 * xac2
      xac9 = xac7 * xac2
      xac11 = xac9 * xac2
      xac13 = xac11 * xac2
      xac15 = xac13 * xac2

      xac6 = xac5 * xac
      xac10 = xac9 * xac
      xac14 = xac13 * xac

      atanxa = pi4 + 0.500000000000000_p_double * xac - 0.250000000000000_p_double * xac2 +&
                     0.083333333333333_p_double * xac3 - 0.025000000000000_p_double * xac5 +&
                     0.020833333333333_p_double * xac6 - 0.008928571428571_p_double * xac7 +&
                     0.003472222222222_p_double * xac9 - 0.003125000000000_p_double * xac10 +&
                     0.001420454545455_p_double * xac11 - 0.000600961538462_p_double * xac13 +&
                     0.000558035714286_p_double * xac14 - 0.000260416666667_p_double * xac15

    else
      !here goes the expansion around x ~ +inf
      xac = 1.0_p_double / xa
      xac2 = xac * xac
      xac3 = xac2 * xac
      xac5 = xac3 * xac2
      xac7 = xac5 * xac2
      xac9 = xac7 * xac2
      xac11 = xac9 * xac2
      xac13 = xac11 * xac2
      xac15 = xac13 * xac2

      atanxa = pi2 - xac + 0.333333333333333_p_double * xac3 -&
                     0.200000000000000_p_double * xac5 + 0.142857142857143_p_double * xac7 -&
                     0.111111111111111_p_double * xac9 + 0.090909090909091_p_double * xac11 -&
                     0.076923076923077_p_double * xac13 + 0.066666666666667_p_double * xac15

    endif

    ! evaluate the signs of variables
    signx = sign(1.0_p_double, x)
    signypos = sign(1.0_p_double, ypos)
    signdif = signypos - signx

    ! evaluate atan2(ypos/xpos)
    xbuf(3,i) = signx * atanxa + pi2 * signdif

    alpha = sqrt(1._p_double - g%rs * auxr)

    px = this%p(1,pp)
    py = this%p(2,pp)
    pz = this%p(3,pp)

    this%p(1,pp) = (st * cp * px + st * sp * py + ct * pz) / alpha
    this%p(2,pp) = ct * cp * px + ct * sp * py - st * pz
    this%p(3,pp) = - sp * px + cp* py

    pp = pp + 1
  end do

end subroutine dudt_dxdt_rk6_cartnorth_gr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk6_cartsouth_gr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufsinit )
!-----------------------------------------------------------------------------

  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 6 + 1 !6th order requires 7 saved quantities

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf, xbufsinit
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt3, gdt12, gdt16, gdt8, gdt44, gdt120
  real(p_double) :: r, t, p, xpos, ypos, zpos
  real(p_double) :: st, ct, sp, cp
  real(p_double) :: pr, pt, pphi, px, py, pz, alpha

  real(p_double), dimension(rank , p_cache_size) :: xbufg
  real(p_k_part), dimension(rank + 1 , p_cache_size) :: xbufinit, ptempinit, xbufs
  real(p_double), dimension(order ,p_p_dim,p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt
  real(p_double) :: rho, auxr, auxrho
  real(p_double) :: shift, xc, xc12, xc32, xc52, xc72, xc92
  real(p_double) :: tc, tc2, tc3, tc4, tc5, tc6, tc7, tc8, tc9, tc10
  real(p_double) :: p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15
  real(p_double) :: x, xa, xa2, xa3, xa5, xa7, xa9, xa11, xa13, xa15
  real(p_double) :: xac, xac2, xac3, xac5, xac6, xac7, xac9
  real(p_double) :: xac10, xac11, xac13, xac14, xac15
  real(p_double) :: pi2, pi4, atanxa, signx, signypos, signdif

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt3 = real( gdt / 3._p_double , p_k_part )
  gdt12 = real( gdt / 12._p_double , p_k_part )
  gdt16 = real( gdt / 16._p_double , p_k_part )
  gdt8 = real( gdt / 8._p_double , p_k_part )
  gdt44 = real( gdt / 44._p_double , p_k_part )
  gdt120 = real( gdt / 120._p_double, p_k_part )

  ! used to evaluate the polynomial approx of trig functions
  shift = pi
  pi2 = pi_2
  pi4 = 0.500000000000000_p_double * pi2

  ! xp0 is the particle's local position
  ! wp0 is the particle's local momenta

  ! loop to get particle's position
  ! xbufs(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbufs(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbufs(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbufs(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! From this point on xbuf and this%p are cartesian
  pp = ptrcur
  do i=1,np
    r = xbufs(1, i)
    t = xbufs(2, i)
    p = xbufs(3, i)

    xbufsinit(1, i) = r
    xbufsinit(2, i) = t
    xbufsinit(3, i) = p

    ! taylor expansion of sin and cos around theta ~ pi
    tc = t - shift

    tc2 = tc * tc
    tc3 = tc2 * tc
    tc4 = tc3 * tc
    tc5 = tc4 * tc
    tc6 = tc5 * tc
    tc7 = tc6 * tc
    tc8 = tc7 * tc
    tc9 = tc8 * tc
    tc10 = tc9 * tc

    st = - tc + 0.166666666666667_p_double * tc3 - 0.008333333333333_p_double * tc5 +&
                0.000198412698413_p_double * tc7 - 0.000002755731922_p_double * tc9

    ct = -1.0_p_double + 0.500000000000000_p_double * tc2 -&
         0.041666666666667_p_double * tc4 + 0.001388888888889_p_double * tc6 -&
         0.000024801587302_p_double * tc8 + 0.000000275573192_p_double * tc10

    ! taylor expansion of sin and cos around phi ~ 0
    p2 = p * p
    p3 = p2 * p
    p4 = p3 * p
    p5 = p4 * p
    p6 = p5 * p
    p7 = p6 * p
    p8 = p7 * p
    p9 = p8 * p
    p10 = p9 * p
    p11 = p10 * p
    p12 = p11 * p
    p13 = p12 * p
    p14 = p13 * p
    p15 = p14 * p

    sp = p - 0.166666666666667_p_double * p3 +&
             0.008333333333333_p_double * p5 - 0.000198412698413_p_double * p7 +&
             0.000002755731922_p_double * p9 - 0.000000025052108_p_double * p11 +&
             0.000000000160590_p_double * p13 - 0.000000000000765_p_double * p15

    cp = 1.0_p_double - 0.500000000000000_p_double * p2 +&
                        0.041666666666667_p_double * p4 - 0.001388888888889_p_double * p6 +&
                        0.000024801587302_p_double * p8 - 0.000000275573192_p_double * p10 +&
                        0.000000002087676_p_double * p12 - 0.000000000011471_p_double * p14

    xbuf(1, i) = r * st * cp
    xbuf(2, i) = r * st * sp
    xbuf(3, i) = r * ct

    pr = this%p(1,pp)
    pt = this%p(2,pp)
    pphi = this%p(3,pp)

    alpha = sqrt(1._p_double - g%rs/r)

    this%p(1,pp) = st * cp * pr * alpha + ct * cp * pt - sp * pphi
    this%p(2,pp) = st * sp * pr * alpha + ct * sp * pt + cp * pphi
    this%p(3,pp) = ct * pr * alpha - st * pt

    pp = pp + 1
  end do

  ! save initial position vector xp0 and momenta
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 (xp0,wp0)
  !-------------------------------------
  step = 1
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt3

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartsouth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 2
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + 2._p_double * kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + 2._p_double * kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + 2._p_double * kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + 2._p_double * lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + 2._p_double * lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + 2._p_double * lvec(step,3,i) * gdt3

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + 2dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartsouth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 3
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i) + 4._p_double * kvec(2,1,i) - kvec(step,1,i)) * gdt12
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i) + 4._p_double * kvec(2,2,i) - kvec(step,2,i)) * gdt12
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i) + 4._p_double * kvec(2,3,i) - kvec(step,3,i)) * gdt12

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i) + 4._p_double * lvec(2,1,i) - lvec(step,1,i)) * gdt12
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i) + 4._p_double * lvec(2,2,i) - lvec(step,2,i)) * gdt12
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i) + 4._p_double * lvec(2,3,i) - lvec(step,3,i)) * gdt12

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartsouth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 4
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp4,wp4)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (-kvec(1,1,i) + 18._p_double * kvec(2,1,i) -&
         3._p_double * kvec(3,1,i) - 6._p_double * kvec(step,1,i)) * gdt16
    xbuf(2,i) = xbufinit(2,i) + (-kvec(1,2,i) + 18._p_double * kvec(2,2,i) -&
         3._p_double * kvec(3,2,i) - 6._p_double * kvec(step,2,i)) * gdt16
    xbuf(3,i) = xbufinit(3,i) + (-kvec(1,3,i) + 18._p_double * kvec(2,3,i) -&
         3._p_double * kvec(3,3,i) - 6._p_double * kvec(step,3,i)) * gdt16

    this%p(1,pp) = ptempinit(1,i) + (-lvec(1,1,i) + 18._p_double * lvec(2,1,i) -&
         3._p_double * lvec(3,1,i) - 6._p_double * lvec(step,1,i)) * gdt16
    this%p(2,pp) = ptempinit(2,i) + (-lvec(1,2,i) + 18._p_double * lvec(2,2,i) -&
         3._p_double * lvec(3,2,i) - 6._p_double * lvec(step,2,i)) * gdt16
    this%p(3,pp) = ptempinit(3,i) + (-lvec(1,3,i) + 18._p_double * lvec(2,3,i) -&
         3._p_double * lvec(3,3,i) - 6._p_double * lvec(step,3,i)) * gdt16

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartsouth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 5
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp5,wp5)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(2,1,i) - 3._p_double * kvec(3,1,i) -&
         6._p_double * kvec(4,1,i) + 4._p_double * kvec(step,1,i)) * gdt8
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(2,2,i) - 3._p_double * kvec(3,2,i) -&
         6._p_double * kvec(4,2,i) + 4._p_double * kvec(step,2,i)) * gdt8
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(2,3,i) - 3._p_double * kvec(3,3,i) -&
         6._p_double * kvec(4,3,i) + 4._p_double * kvec(step,3,i)) * gdt8

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(2,1,i) - 3._p_double * lvec(3,1,i) -&
         6._p_double * lvec(4,1,i) + 4._p_double * lvec(step,1,i)) * gdt8
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(2,2,i) - 3._p_double * lvec(3,2,i) -&
         6._p_double * lvec(4,2,i) + 4._p_double * lvec(step,2,i)) * gdt8
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(2,3,i) - 3._p_double * lvec(3,3,i) -&
         6._p_double * lvec(4,3,i) + 4._p_double * lvec(step,3,i)) * gdt8

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartsouth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 6
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp6,wp6)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(1,1,i) - 36._p_double * kvec(2,1,i) +&
         63._p_double * kvec(3,1,i) + 72._p_double * kvec(4,1,i) - 64._p_double * kvec(step,1,i)) * gdt44
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(1,2,i) - 36._p_double * kvec(2,2,i) +&
         63._p_double * kvec(3,2,i) + 72._p_double * kvec(4,2,i) - 64._p_double * kvec(step,2,i)) * gdt44
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(1,3,i) - 36._p_double * kvec(2,3,i) +&
         63._p_double * kvec(3,3,i) + 72._p_double * kvec(4,3,i) - 64._p_double * kvec(step,3,i)) * gdt44

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(1,1,i) - 36._p_double * lvec(2,1,i) +&
         63._p_double * lvec(3,1,i) + 72._p_double * lvec(4,1,i) - 64._p_double * lvec(step,1,i)) * gdt44
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(1,2,i) - 36._p_double * lvec(2,2,i) +&
         63._p_double * lvec(3,2,i) + 72._p_double * lvec(4,2,i) - 64._p_double * lvec(step,2,i)) * gdt44
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(1,3,i) - 36._p_double * lvec(2,3,i) +&
         63._p_double * lvec(3,3,i) + 72._p_double * lvec(4,3,i) - 64._p_double * lvec(step,3,i)) * gdt44

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_cartsouth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 7
  call get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------
  ! t = t0 + dt (xp,wp)
  !-------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (11._p_double * kvec(1,1,i) + 81._p_double * (kvec(3,1,i) + kvec(4,1,i)) -&
          32._p_double * (kvec(5,1,i) + kvec(6,1,i)) + 11._p_double * kvec(step,1,i)) * gdt120
    xbuf(2,i) = xbufinit(2,i) + (11._p_double * kvec(1,2,i) + 81._p_double * (kvec(3,2,i) + kvec(4,2,i)) -&
          32._p_double * (kvec(5,2,i) + kvec(6,2,i)) + 11._p_double * kvec(step,2,i)) * gdt120
    xbuf(3,i) = xbufinit(3,i) + (11._p_double * kvec(1,3,i) + 81._p_double * (kvec(3,3,i) + kvec(4,3,i)) -&
          32._p_double * (kvec(5,3,i) + kvec(6,3,i)) + 11._p_double * kvec(step,3,i)) * gdt120

    this%p(1,pp) = ptempinit(1,i) + (11._p_double * lvec(1,1,i) + 81._p_double * (lvec(3,1,i) + lvec(4,1,i)) -&
          32._p_double * (lvec(5,1,i) + lvec(6,1,i)) + 11._p_double * lvec(step,1,i)) * gdt120
    this%p(2,pp) = ptempinit(2,i) + (11._p_double * lvec(1,2,i) + 81._p_double * (lvec(3,2,i) + lvec(4,2,i)) -&
          32._p_double * (lvec(5,2,i) + lvec(6,2,i)) + 11._p_double * lvec(step,2,i)) * gdt120
    this%p(3,pp) = ptempinit(3,i) + (11._p_double * lvec(1,3,i) + 81._p_double * (lvec(3,3,i) + lvec(4,3,i)) -&
          32._p_double * (lvec(5,3,i) + lvec(6,3,i)) + 11._p_double * lvec(step,3,i)) * gdt120

    pp = pp + 1
  end do

  ! convert xbuf and this%p back to spherical
  pp = ptrcur
  do i=1,np
    xpos = xbuf(1,i)
    ypos = xbuf(2,i)
    zpos = xbuf(3,i)
    rho = sqrt(xpos * xpos + ypos * ypos)

    xbuf(1,i) = sqrt(xpos**2 + ypos**2 + zpos**2)

    ! evaluate sin and cos of theta
    auxr = 1.0_p_double / xbuf(1,i)
    st = rho * auxr
    ct = zpos * auxr

    ! evaluate sin and cos of phi
    auxrho = 1.0_p_double / rho
    sp = ypos * auxrho
    cp = xpos * auxrho

    xc = 1.0_p_double + zpos * auxr
    xc12 = - sqrt(xc)
    xc32 = xc12 * xc
    xc52 = xc32 * xc
    xc72 = xc52 * xc
    xc92 = xc72 * xc

    ! theta = acos(z/r) , z<0
    xbuf(2,i) = 1.414213562373095_p_double * xc12 + 0.117851130197758_p_double * xc32 +&
                0.026516504294496_p_double * xc52 + 0.007891816754314_p_double * xc72 +&
                0.002685409867787_p_double * xc92 + shift

    ! evaluate the atan2 of x
    ! argument that enters atan(x)
    x = ypos / xpos
    xa = abs( x )

    ! loop that selects the expansion point x~0,1,inf
    if (xa .le. 0.418597947318945_p_double) then
      !here goes the expansion around x ~ 0
      xa2 = xa * xa
      xa3 = xa2 * xa
      xa5 = xa3 * xa2
      xa7 = xa5 * xa2
      xa9 = xa7 * xa2
      xa11 = xa9 * xa2
      xa13 = xa11 * xa2
      xa15 = xa13 * xa2

      ! evaluate atan(abs(x)) around x ~ 0
      atanxa = xa - 0.333333333333333_p_double * xa3 +&
                    0.200000000000000_p_double * xa5 - 0.142857142857143_p_double * xa7 +&
                    0.111111111111111_p_double * xa9 - 0.090909090909091_p_double * xa11 +&
                    0.076923076923077_p_double * xa13 - 0.066666666666667_p_double * xa15

    else if (xa .le. 1.817313073051244_p_double) then
      !here goes the expansion around x ~ 1
      xac = xa - 1.0_p_double
      xac2 = xac * xac
      xac3 = xac2 * xac
      xac5 = xac3 * xac2
      xac7 = xac5 * xac2
      xac9 = xac7 * xac2
      xac11 = xac9 * xac2
      xac13 = xac11 * xac2
      xac15 = xac13 * xac2

      xac6 = xac5 * xac
      xac10 = xac9 * xac
      xac14 = xac13 * xac

      atanxa = pi4 + 0.500000000000000_p_double * xac - 0.250000000000000_p_double * xac2 +&
                     0.083333333333333_p_double * xac3 - 0.025000000000000_p_double * xac5 +&
                     0.020833333333333_p_double * xac6 - 0.008928571428571_p_double * xac7 +&
                     0.003472222222222_p_double * xac9 - 0.003125000000000_p_double * xac10 +&
                     0.001420454545455_p_double * xac11 - 0.000600961538462_p_double * xac13 +&
                     0.000558035714286_p_double * xac14 - 0.000260416666667_p_double * xac15

    else
      !here goes the expansion around x ~ +inf
      xac = 1.0_p_double / xa
      xac2 = xac * xac
      xac3 = xac2 * xac
      xac5 = xac3 * xac2
      xac7 = xac5 * xac2
      xac9 = xac7 * xac2
      xac11 = xac9 * xac2
      xac13 = xac11 * xac2
      xac15 = xac13 * xac2

      atanxa = pi2 - xac + 0.333333333333333_p_double * xac3 -&
                     0.200000000000000_p_double * xac5 + 0.142857142857143_p_double * xac7 -&
                     0.111111111111111_p_double * xac9 + 0.090909090909091_p_double * xac11 -&
                     0.076923076923077_p_double * xac13 + 0.066666666666667_p_double * xac15

    endif

    ! signs of variables to evaluate atan2
    signx = sign(1.0_p_double, x)
    signypos = sign(1.0_p_double, ypos)
    signdif = signypos - signx

    ! evaluate atan2(ypos/xpos)
    xbuf(3,i) = signx * atanxa + pi2 * signdif

    alpha = sqrt(1._p_double - g%rs * auxr)

    px = this%p(1,pp)
    py = this%p(2,pp)
    pz = this%p(3,pp)

    this%p(1,pp) = (st * cp * px + st * sp * py + ct * pz) / alpha
    this%p(2,pp) = ct * cp * px + ct * sp * py - st * pz
    this%p(3,pp) = - sp * px + cp* py

    pp = pp + 1
  end do

end subroutine dudt_dxdt_rk6_cartsouth_gr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dudt_dxdt_rk6_cart_gr_rr( this, emf, jay, energy, xbuf, gdt, ptrcur, np, xbufsinit )
!-----------------------------------------------------------------------------

  use m_species_push, only : ntrim
  use m_geometry_gr, only : t_geometry_gr

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: order = 6 + 1 !6th order requires 7 saved quantities

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(out) :: xbuf, xbufsinit
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: ptrcur, np

  integer :: i, pp
  integer :: ir, it, step

  real(p_double) :: gdt3, gdt12, gdt16, gdt8, gdt44, gdt120
  real(p_double) :: r, t, p, xpos, ypos, zpos
  real(p_double) :: st, ct, sp, cp
  real(p_double) :: pr, pt, pphi, px, py, pz, alpha

  real(p_double), dimension(rank , p_cache_size) :: xbufg
  integer, dimension(rank , p_cache_size) :: dxi
  real(p_k_part), dimension(rank + 1 , p_cache_size) :: xbufinit, ptempinit, xbufs
  real(p_double), dimension(order ,p_p_dim,p_cache_size) :: kvec, lvec
  real(p_k_part) :: dt

  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  dt = real(g%dt, p_k_part)

  ! required for fractional step
  gdt3 = real( gdt / 3._p_double , p_k_part )
  gdt12 = real( gdt / 12._p_double , p_k_part )
  gdt16 = real( gdt / 16._p_double , p_k_part )
  gdt8 = real( gdt / 8._p_double , p_k_part )
  gdt44 = real( gdt / 44._p_double , p_k_part )
  gdt120 = real( gdt / 120._p_double, p_k_part )

  ! xp0 is the particle's local position
  ! wp0 is the particle's local momenta

  ! loop to get particle's position
  ! xbufs(:,i) = (r,theta,phi) at t=t0
  pp = ptrcur
  do i=1,np
    xbufg(1,i) = this%x ( 1, pp )
    xbufg(2,i) = this%x ( 2, pp )

    ir = this%ix( 1, pp )
    xbufs(1, i) = g%r%f1(2, ir) + g%dr%f1(1, ir) * xbufg(1,i)

    it = this%ix( 2, pp )
    xbufs(2, i) = g%t%f1(2, it) + dt * xbufg(2,i)

    xbufs(3, i) = this%x ( 3, pp )

    pp = pp + 1
  end do

  ! From this point on xbuf and this%p are cartesian
  pp = ptrcur
  do i=1,np
    r = xbufs(1, i)
    t = xbufs(2, i)
    p = xbufs(3, i)

    xbufsinit(1, i) = r
    xbufsinit(2, i) = t
    xbufsinit(3, i) = p

    st = sin(t)
    ct = cos(t)
    sp = sin(p)
    cp = cos(p)

    xbuf(1, i) = r * st * cp
    xbuf(2, i) = r * st * sp
    xbuf(3, i) = r * ct

    pr = this%p(1,pp)
    pt = this%p(2,pp)
    pphi = this%p(3,pp)

    alpha = sqrt(1._p_double - g%rs/r)

    this%p(1,pp) = st * cp * pr * alpha + ct * cp * pt - sp * pphi
    this%p(2,pp) = st * sp * pr * alpha + ct * sp * pt + cp * pphi
    this%p(3,pp) = ct * pr * alpha - st * pt

    pp = pp + 1
  end do

  ! save initial position vector xp0 and momenta
  pp = ptrcur
  do i=1,np
    xbufinit(1, i) = xbuf(1, i)
    xbufinit(2, i) = xbuf(2, i)
    xbufinit(3, i) = xbuf(3, i)
    ptempinit(1,i) = this%p(1,pp)
    ptempinit(2,i) = this%p(2,pp)
    ptempinit(3,i) = this%p(3,pp)

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 (xp0,wp0)
  !-------------------------------------
  step = 1
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp1,wp1)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + lvec(step,3,i) * gdt3

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! get particle's grid quantities for emf interpolation
  call get_sph_grid_pos_gr( this, xbuf, xbufs, ptrcur, np, dt)

  step = 2
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp2,wp2)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + 2._p_double * kvec(step,1,i) * gdt3
    xbuf(2,i) = xbufinit(2,i) + 2._p_double * kvec(step,2,i) * gdt3
    xbuf(3,i) = xbufinit(3,i) + 2._p_double * kvec(step,3,i) * gdt3

    this%p(1,pp) = ptempinit(1,i) + 2._p_double * lvec(step,1,i) * gdt3
    this%p(2,pp) = ptempinit(2,i) + 2._p_double * lvec(step,2,i) * gdt3
    this%p(3,pp) = ptempinit(3,i) + 2._p_double * lvec(step,3,i) * gdt3

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + 2dt/3
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 3
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp3,wp3)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (kvec(1,1,i) + 4._p_double * kvec(2,1,i) - kvec(step,1,i)) * gdt12
    xbuf(2,i) = xbufinit(2,i) + (kvec(1,2,i) + 4._p_double * kvec(2,2,i) - kvec(step,2,i)) * gdt12
    xbuf(3,i) = xbufinit(3,i) + (kvec(1,3,i) + 4._p_double * kvec(2,3,i) - kvec(step,3,i)) * gdt12

    this%p(1,pp) = ptempinit(1,i) + (lvec(1,1,i) + 4._p_double * lvec(2,1,i) - lvec(step,1,i)) * gdt12
    this%p(2,pp) = ptempinit(2,i) + (lvec(1,2,i) + 4._p_double * lvec(2,2,i) - lvec(step,2,i)) * gdt12
    this%p(3,pp) = ptempinit(3,i) + (lvec(1,3,i) + 4._p_double * lvec(2,3,i) - lvec(step,3,i)) * gdt12

    pp = pp + 1
  end do
  !-------------------------------------
  ! t = t0 + dt/3
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 4
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp4,wp4)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (-kvec(1,1,i) + 18._p_double * kvec(2,1,i) -&
         3._p_double * kvec(3,1,i) - 6._p_double * kvec(step,1,i)) * gdt16
    xbuf(2,i) = xbufinit(2,i) + (-kvec(1,2,i) + 18._p_double * kvec(2,2,i) -&
         3._p_double * kvec(3,2,i) - 6._p_double * kvec(step,2,i)) * gdt16
    xbuf(3,i) = xbufinit(3,i) + (-kvec(1,3,i) + 18._p_double * kvec(2,3,i) -&
         3._p_double * kvec(3,3,i) - 6._p_double * kvec(step,3,i)) * gdt16

    this%p(1,pp) = ptempinit(1,i) + (-lvec(1,1,i) + 18._p_double * lvec(2,1,i) -&
         3._p_double * lvec(3,1,i) - 6._p_double * lvec(step,1,i)) * gdt16
    this%p(2,pp) = ptempinit(2,i) + (-lvec(1,2,i) + 18._p_double * lvec(2,2,i) -&
         3._p_double * lvec(3,2,i) - 6._p_double * lvec(step,2,i)) * gdt16
    this%p(3,pp) = ptempinit(3,i) + (-lvec(1,3,i) + 18._p_double * lvec(2,3,i) -&
         3._p_double * lvec(3,3,i) - 6._p_double * lvec(step,3,i)) * gdt16

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 5
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp5,wp5)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(2,1,i) - 3._p_double * kvec(3,1,i) -&
         6._p_double * kvec(4,1,i) + 4._p_double * kvec(step,1,i)) * gdt8
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(2,2,i) - 3._p_double * kvec(3,2,i) -&
         6._p_double * kvec(4,2,i) + 4._p_double * kvec(step,2,i)) * gdt8
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(2,3,i) - 3._p_double * kvec(3,3,i) -&
         6._p_double * kvec(4,3,i) + 4._p_double * kvec(step,3,i)) * gdt8

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(2,1,i) - 3._p_double * lvec(3,1,i) -&
         6._p_double * lvec(4,1,i) + 4._p_double * lvec(step,1,i)) * gdt8
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(2,2,i) - 3._p_double * lvec(3,2,i) -&
         6._p_double * lvec(4,2,i) + 4._p_double * lvec(step,2,i)) * gdt8
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(2,3,i) - 3._p_double * lvec(3,3,i) -&
         6._p_double * lvec(4,3,i) + 4._p_double * lvec(step,3,i)) * gdt8

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt/2
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 6
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  ! next step's position and momenta (xp6,wp6)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (9._p_double * kvec(1,1,i) - 36._p_double * kvec(2,1,i) +&
         63._p_double * kvec(3,1,i) + 72._p_double * kvec(4,1,i) - 64._p_double * kvec(step,1,i)) * gdt44
    xbuf(2,i) = xbufinit(2,i) + (9._p_double * kvec(1,2,i) - 36._p_double * kvec(2,2,i) +&
         63._p_double * kvec(3,2,i) + 72._p_double * kvec(4,2,i) - 64._p_double * kvec(step,2,i)) * gdt44
    xbuf(3,i) = xbufinit(3,i) + (9._p_double * kvec(1,3,i) - 36._p_double * kvec(2,3,i) +&
         63._p_double * kvec(3,3,i) + 72._p_double * kvec(4,3,i) - 64._p_double * kvec(step,3,i)) * gdt44

    this%p(1,pp) = ptempinit(1,i) + (9._p_double * lvec(1,1,i) - 36._p_double * lvec(2,1,i) +&
         63._p_double * lvec(3,1,i) + 72._p_double * lvec(4,1,i) - 64._p_double * lvec(step,1,i)) * gdt44
    this%p(2,pp) = ptempinit(2,i) + (9._p_double * lvec(1,2,i) - 36._p_double * lvec(2,2,i) +&
         63._p_double * lvec(3,2,i) + 72._p_double * lvec(4,2,i) - 64._p_double * lvec(step,2,i)) * gdt44
    this%p(3,pp) = ptempinit(3,i) + (9._p_double * lvec(1,3,i) - 36._p_double * lvec(2,3,i) +&
         63._p_double * lvec(3,3,i) + 72._p_double * lvec(4,3,i) - 64._p_double * lvec(step,3,i)) * gdt44

    pp = pp + 1
  end do

  !-------------------------------------
  ! t = t0 + dt
  !-------------------------------------
  ! loop to get particle's grid quantities for emf interpolation
  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  pp = ptrcur
  do i=1,np
    this%x(1,pp) = xbufg(1,i)
    this%x(2,pp) = xbufg(2,i)
    this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
    this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

    pp = pp + 1
  end do

  step = 7
  call get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  !-------------------------------------
  ! t = t0 + dt (xp,wp)
  !-------------------------------------
  ! final position and momenta (xp,wp)
  pp = ptrcur
  do i=1,np
    xbuf(1,i) = xbufinit(1,i) + (11._p_double * kvec(1,1,i) + 81._p_double * (kvec(3,1,i) + kvec(4,1,i)) -&
          32._p_double * (kvec(5,1,i) + kvec(6,1,i)) + 11._p_double * kvec(step,1,i)) * gdt120
    xbuf(2,i) = xbufinit(2,i) + (11._p_double * kvec(1,2,i) + 81._p_double * (kvec(3,2,i) + kvec(4,2,i)) -&
          32._p_double * (kvec(5,2,i) + kvec(6,2,i)) + 11._p_double * kvec(step,2,i)) * gdt120
    xbuf(3,i) = xbufinit(3,i) + (11._p_double * kvec(1,3,i) + 81._p_double * (kvec(3,3,i) + kvec(4,3,i)) -&
          32._p_double * (kvec(5,3,i) + kvec(6,3,i)) + 11._p_double * kvec(step,3,i)) * gdt120

    this%p(1,pp) = ptempinit(1,i) + (11._p_double * lvec(1,1,i) + 81._p_double * (lvec(3,1,i) + lvec(4,1,i)) -&
          32._p_double * (lvec(5,1,i) + lvec(6,1,i)) + 11._p_double * lvec(step,1,i)) * gdt120
    this%p(2,pp) = ptempinit(2,i) + (11._p_double * lvec(1,2,i) + 81._p_double * (lvec(3,2,i) + lvec(4,2,i)) -&
          32._p_double * (lvec(5,2,i) + lvec(6,2,i)) + 11._p_double * lvec(step,2,i)) * gdt120
    this%p(3,pp) = ptempinit(3,i) + (11._p_double * lvec(1,3,i) + 81._p_double * (lvec(3,3,i) + lvec(4,3,i)) -&
          32._p_double * (lvec(5,3,i) + lvec(6,3,i)) + 11._p_double * lvec(step,3,i)) * gdt120

    pp = pp + 1
  end do

  ! convert xbuf and this%p back to spherical
  pp = ptrcur
  do i=1,np
    xpos = xbuf(1,i)
    ypos = xbuf(2,i)
    zpos = xbuf(3,i)

    xbuf(1,i) = sqrt(xpos**2 + ypos**2 + zpos**2)
    xbuf(2,i) = acos( zpos / xbuf(1,i) )
    xbuf(3,i) = atan2(ypos,xpos)

    alpha = sqrt(1._p_double-g%rs/xbuf(1,i))

    ct = cos(xbuf(2,i))
    st = sin(xbuf(2,i))

    cp = cos(xbuf(3,i))
    sp = sin(xbuf(3,i))

    px = this%p(1,pp)
    py = this%p(2,pp)
    pz = this%p(3,pp)

    this%p(1,pp) = (st * cp * px + st * sp * py + ct * pz) / alpha
    this%p(2,pp) = ct * cp * px + ct * sp * py - st * pz
    this%p(3,pp) = - sp * px + cp* py

    pp = pp + 1
  end do

end subroutine dudt_dxdt_rk6_cart_gr_rr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine get_kvec_lvec_rk_gr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  use m_emf_interpolate, only : get_emf
  use m_geometry_gr

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(inout) :: xbuf
  integer, intent(in) :: ptrcur, np, step, order
  real(p_double), dimension(order ,p_p_dim,p_cache_size), intent(inout) :: kvec, lvec


  integer :: i, pp

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bptemp, eptemp
  real(p_double), dimension(p_cache_size) :: alphaptemp, gpart
  real(p_double), dimension(p_p_dim,p_cache_size) :: betaptemp, hvecptemp
  real(p_double), dimension(p_p_dim,p_cache_size) :: gravptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfsrptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfstptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfspptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: htensorptemp


  ! executable statements

  ! get local gr quantaties (analitical expressions)
  select case (this%geometry%metric)
  case (p_geometry_minkowski)
    call get_geometric_quant_sph(this, alphaptemp, betaptemp, hvecptemp, cfsrptemp, &
    cfstptemp, cfspptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_schwarzschild)
    call get_geometric_quant_sch(this, alphaptemp, betaptemp, hvecptemp, cfsrptemp, &
    cfstptemp, cfspptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_kerr_slow)
    call get_geometric_quant_kslow(this, alphaptemp, betaptemp, hvecptemp, cfsrptemp, &
    cfstptemp, cfspptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_kerr)
    !not implemented yet
  end select

  ! get em field in particle's local position (interpolation)
  call get_emf( emf, bptemp, eptemp, this%ix(:,ptrcur:), this%x(:,ptrcur:), &
    np, this%interpolation )

  ! local gamma factor
  pp = ptrcur
  do i=1,np
    gpart(i) = sqrt( ( ( 1.0_p_k_part + this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

    pp = pp + 1
  enddo

  ! K is the value of dx/dt (used to update positions)
  ! K(step) -> kvec(step,:,:)
  pp = ptrcur
  do i=1,np
    kvec(step,1,i) = (alphaptemp(i) * this%p(1,pp) / gpart(i) &
              - betaptemp(1, i)) / hvecptemp(1, i)
    kvec(step,2,i) = (alphaptemp(i) * this%p(2,pp) / gpart(i) &
              - betaptemp(2, i)) / hvecptemp(2, i)
    kvec(step,3,i) = (alphaptemp(i) * this%p(3,pp) / gpart(i) &
              - betaptemp(3, i)) / hvecptemp(3, i)

    pp = pp + 1
  end do

  ! L is the value of dw/dt (used to update momenta)
  ! L(step) -> lvec(step,:,:)
  pp = ptrcur
  do i=1,np
    lvec(step,1,i) = alphaptemp(i)*( (eptemp(1,i) + &
           (this%p(2,pp)*bptemp(3,i)-this%p(3,pp)*bptemp(2,i))/ gpart(i) )/this%rqm +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(1,i) +&
           htensorptemp(1,1,i) * this%p(1,pp) +&
           htensorptemp(1,2,i) * this%p(2,pp) +&
           htensorptemp(1,3,i) * this%p(3,pp) )) -&
           (cfsrptemp(1,2,i) * kvec(step,1,i)*this%p(2,pp) +&
            cfsrptemp(1,3,i) * kvec(step,1,i)*this%p(3,pp) +&
            cfsrptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfsrptemp(2,3,i) * kvec(step,2,i)*this%p(3,pp) +&
            cfsrptemp(3,2,i) * kvec(step,3,i)*this%p(2,pp) +&
            cfsrptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )

    lvec(step,2,i) = alphaptemp(i)*( (eptemp(2,i) + &
           (this%p(3,pp)*bptemp(1,i)-this%p(1,pp)*bptemp(3,i))/ gpart(i) )/this%rqm +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(2,i) +&
           htensorptemp(2,1,i) * this%p(1,pp) +&
           htensorptemp(2,2,i) * this%p(2,pp) +&
           htensorptemp(2,3,i) * this%p(3,pp) )) -&
           (cfstptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfstptemp(1,3,i) * kvec(step,1,i)*this%p(3,pp) +&
            cfstptemp(2,1,i) * kvec(step,2,i)*this%p(1,pp) +&
            cfstptemp(2,3,i) * kvec(step,2,i)*this%p(3,pp) +&
            cfstptemp(3,1,i) * kvec(step,3,i)*this%p(1,pp) +&
            cfstptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )

    lvec(step,3,i) = alphaptemp(i)*( (eptemp(3,i) + &
           (this%p(1,pp)*bptemp(2,i)-this%p(2,pp)*bptemp(1,i))/ gpart(i) )/this%rqm +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(3,i) +&
           htensorptemp(3,1,i) * this%p(1,pp) +&
           htensorptemp(3,2,i) * this%p(2,pp) +&
           htensorptemp(3,3,i) * this%p(3,pp) )) -&
           (cfspptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfspptemp(1,2,i) * kvec(step,1,i)*this%p(2,pp) +&
            cfspptemp(2,1,i) * kvec(step,2,i)*this%p(1,pp) +&
            cfspptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfspptemp(3,1,i) * kvec(step,3,i)*this%p(1,pp) +&
            cfspptemp(3,2,i) * kvec(step,3,i)*this%p(2,pp) )

    pp = pp + 1
  end do

end subroutine get_kvec_lvec_rk_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine get_kvec_lvec_rk_cart_gr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  use m_emf_interpolate, only : get_emf
  use m_geometry_gr

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(inout) :: xbuf
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(in) :: xbufs
  integer, intent(in) :: ptrcur, np, step, order
  real(p_double), dimension(order ,p_p_dim,p_cache_size), intent(inout) :: kvec, lvec


  integer :: i, pp
  real(p_double) :: r, st, ct, sp, cp, alpha
  real(p_double) :: rho, auxr, auxrho, rsloc

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bptemp, eptemp, ep, bp
  real(p_double), dimension(p_cache_size) :: alphaptemp, gpart, smdet
  real(p_double), dimension(p_p_dim,p_cache_size) :: betaptemp
  real(p_double), dimension(p_p_dim,p_cache_size) :: gravptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: smcov, smctr
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfsxptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfsyptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfszptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: htensorptemp

  ! executable statements

  ! get local gr quantaties (analitical expressions)
  select case (this%geometry%metric)
  case (p_geometry_minkowski)
    call get_geometric_quant_sph_cart(this, alphaptemp, betaptemp, smcov, smctr, smdet,&
    cfsxptemp, cfsyptemp, cfszptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_schwarzschild)
    call get_geometric_quant_sch_cart(this, alphaptemp, betaptemp, smcov, smctr, smdet,&
    cfsxptemp, cfsyptemp, cfszptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_kerr_slow)
    call get_geometric_quant_kslow_cart(this, alphaptemp, betaptemp, smcov, smctr, smdet,&
    cfsxptemp, cfsyptemp, cfszptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_kerr)
    !not implemented yet
  end select

  ! get em field in particle's local position (interpolation)
  call get_emf( emf, bptemp, eptemp, this%ix(:,ptrcur:), this%x(:,ptrcur:), &
    np, this%interpolation )

  rsloc = this%geometry%rs

  ! convert spherical fields to cartesian
  pp = ptrcur
  do i=1,np
    r = xbufs(1,i)
    rho = sqrt(xbuf(1,i) * xbuf(1,i) + xbuf(2,i) * xbuf(2,i))

    ! trigonometric functions of theta
    auxr = 1.0_p_double / r
    st = rho * auxr
    ct = xbuf(3,i) * auxr

    ! trigonometric functions of phi
    auxrho = 1.0_p_double / rho
    sp = xbuf(2,i) * auxrho
    cp = xbuf(1,i) * auxrho

    alpha = sqrt(1._p_double - rsloc * auxr)

    ep(1,i) = st * cp * eptemp(1,i) * alpha + ct * cp * eptemp(2,i) - sp * eptemp(3,i)
    ep(2,i) = st * sp * eptemp(1,i) * alpha + ct * sp * eptemp(2,i) + cp * eptemp(3,i)
    ep(3,i) = ct * eptemp(1,i) * alpha - st * eptemp(2,i)

    bp(1,i) = st * cp * bptemp(1,i) * alpha + ct * cp * bptemp(2,i) - sp * bptemp(3,i)
    bp(2,i) = st * sp * bptemp(1,i) * alpha + ct * sp * bptemp(2,i) + cp * bptemp(3,i)
    bp(3,i) = ct * bptemp(1,i) * alpha - st * bptemp(2,i)

    pp = pp + 1
  enddo

  ! local gamma factor
  pp = ptrcur
  do i=1,np
    gpart(i) = sqrt( ( ( 1.0_p_k_part + smcov(1,1,i) * this%p(1,pp)**2 ) +&
          smcov(2,2,i) * this%p(2,pp)**2 ) +&
          smcov(3,3,i) * this%p(3,pp)**2 +&
          2._p_double * smcov(1,2,i) * this%p(1,pp) * this%p(2,pp) +&
          2._p_double * smcov(1,3,i) * this%p(1,pp) * this%p(3,pp) +&
          2._p_double * smcov(2,3,i) * this%p(2,pp) * this%p(3,pp) )

    pp = pp + 1
  enddo

  ! K is the value of dx/dt (used to update positions)
  ! K(step) -> kvec(step,:,:)
  pp = ptrcur
  do i=1,np
    kvec(step,1,i) = alphaptemp(i) * this%p(1,pp) / gpart(i) &
              - betaptemp(1, i)
    kvec(step,2,i) = alphaptemp(i) * this%p(2,pp) / gpart(i) &
              - betaptemp(2, i)
    kvec(step,3,i) = alphaptemp(i) * this%p(3,pp) / gpart(i) &
              - betaptemp(3, i)

    pp = pp + 1
  end do

  ! L is the value of dw/dt (used to update momenta)
  ! L(step) -> lvec(step,:,:)
  pp = ptrcur
  do i=1,np
    lvec(step,1,i) = alphaptemp(i)*( (ep(1,i) + sqrt(smdet(i)) * &
           (smctr(1,1,i)*(this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)) +&
            smctr(1,2,i)*(this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)) +&
            smctr(1,3,i)*(this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)))/ gpart(i) )/this%rqm +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(1,i) +&
           this%p(1,pp) * (htensorptemp(1,1,i) * smctr(1,1,i) +&
                           htensorptemp(2,1,i) * smctr(1,2,i) +&
                           htensorptemp(3,1,i) * smctr(1,3,i)) +&
           this%p(2,pp) * (htensorptemp(1,2,i) * smctr(1,1,i) +&
                           htensorptemp(2,2,i) * smctr(1,2,i) +&
                           htensorptemp(3,2,i) * smctr(1,3,i)) +&
           this%p(3,pp) * (htensorptemp(1,3,i) * smctr(1,1,i) +&
                           htensorptemp(2,3,i) * smctr(1,2,i) +&
                           htensorptemp(3,3,i) * smctr(1,3,i)) )) -&
           (cfsxptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfsxptemp(1,2,i) * (kvec(step,1,i)*this%p(2,pp) + kvec(step,2,i)*this%p(1,pp)) +&
            cfsxptemp(1,3,i) * (kvec(step,1,i)*this%p(3,pp) + kvec(step,3,i)*this%p(1,pp)) +&
            cfsxptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfsxptemp(2,3,i) * (kvec(step,2,i)*this%p(3,pp) + kvec(step,3,i)*this%p(2,pp)) +&
            cfsxptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )

    lvec(step,2,i) = alphaptemp(i)*( (ep(2,i) + sqrt(smdet(i)) * &
           (smctr(2,1,i)*(this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)) +&
            smctr(2,2,i)*(this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)) +&
            smctr(2,3,i)*(this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)))/ gpart(i) )/this%rqm +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(2,i) +&
           this%p(1,pp) * (htensorptemp(1,1,i) * smctr(2,1,i) +&
                           htensorptemp(2,1,i) * smctr(2,2,i) +&
                           htensorptemp(3,1,i) * smctr(2,3,i)) +&
           this%p(2,pp) * (htensorptemp(1,2,i) * smctr(2,1,i) +&
                           htensorptemp(2,2,i) * smctr(2,2,i) +&
                           htensorptemp(3,2,i) * smctr(2,3,i)) +&
           this%p(3,pp) * (htensorptemp(1,3,i) * smctr(2,1,i) +&
                           htensorptemp(2,3,i) * smctr(2,2,i) +&
                           htensorptemp(3,3,i) * smctr(2,3,i)) )) -&
           (cfsyptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfsyptemp(1,2,i) * (kvec(step,1,i)*this%p(2,pp) + kvec(step,2,i)*this%p(1,pp)) +&
            cfsyptemp(1,3,i) * (kvec(step,1,i)*this%p(3,pp) + kvec(step,3,i)*this%p(1,pp)) +&
            cfsyptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfsyptemp(2,3,i) * (kvec(step,2,i)*this%p(3,pp) + kvec(step,3,i)*this%p(2,pp)) +&
            cfsyptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )

    lvec(step,3,i) = alphaptemp(i)*( (ep(3,i) + sqrt(smdet(i)) * &
           (smctr(3,1,i)*(this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)) +&
            smctr(3,2,i)*(this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)) +&
            smctr(3,3,i)*(this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)))/ gpart(i) )/this%rqm +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(3,i) +&
           this%p(1,pp) * (htensorptemp(1,1,i) * smctr(3,1,i) +&
                           htensorptemp(2,1,i) * smctr(3,2,i) +&
                           htensorptemp(3,1,i) * smctr(3,3,i)) +&
           this%p(2,pp) * (htensorptemp(1,2,i) * smctr(3,1,i) +&
                           htensorptemp(2,2,i) * smctr(3,2,i) +&
                           htensorptemp(3,2,i) * smctr(3,3,i)) +&
           this%p(3,pp) * (htensorptemp(1,3,i) * smctr(3,1,i) +&
                           htensorptemp(2,3,i) * smctr(3,2,i) +&
                           htensorptemp(3,3,i) * smctr(3,3,i)) )) -&
           (cfszptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfszptemp(1,2,i) * (kvec(step,1,i)*this%p(2,pp) + kvec(step,2,i)*this%p(1,pp)) +&
            cfszptemp(1,3,i) * (kvec(step,1,i)*this%p(3,pp) + kvec(step,3,i)*this%p(1,pp)) +&
            cfszptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfszptemp(2,3,i) * (kvec(step,2,i)*this%p(3,pp) + kvec(step,3,i)*this%p(2,pp)) +&
            cfszptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )

    pp = pp + 1
  end do

end subroutine get_kvec_lvec_rk_cart_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine get_kvec_lvec_rk_gr_rr( this, emf, xbuf, kvec, lvec, ptrcur, np, order, step )

  use m_emf_interpolate, only : get_emf
  use m_geometry_gr

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(inout) :: xbuf
  integer, intent(in) :: ptrcur, np, step, order
  real(p_double), dimension(order ,p_p_dim,p_cache_size), intent(inout) :: kvec, lvec


  integer :: i, pp

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bptemp, eptemp
  real(p_double), dimension(p_cache_size) :: alphaptemp, gpart
  real(p_double), dimension(p_cache_size) :: FLr,FLt,FLp,uE,auxRR
  real(p_double), dimension(p_cache_size) :: RRr,RRt,RRp
  real(p_double), dimension(p_p_dim,p_cache_size) :: betaptemp, hvecptemp
  real(p_double), dimension(p_p_dim,p_cache_size) :: gravptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfsrptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfstptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfspptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: htensorptemp


  ! executable statements

  ! get local gr quantaties (analitical expressions)
  select case (this%geometry%metric)
  case (p_geometry_minkowski)
    call get_geometric_quant_sph(this, alphaptemp, betaptemp, hvecptemp, cfsrptemp, &
    cfstptemp, cfspptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_schwarzschild)
    call get_geometric_quant_sch(this, alphaptemp, betaptemp, hvecptemp, cfsrptemp, &
    cfstptemp, cfspptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_kerr_slow)
    call get_geometric_quant_kslow(this, alphaptemp, betaptemp, hvecptemp, cfsrptemp, &
    cfstptemp, cfspptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_kerr)
    !not implemented yet
  end select

  ! get em field in particle's local position (interpolation)
  call get_emf( emf, bptemp, eptemp, this%ix(:,ptrcur:), this%x(:,ptrcur:), &
    np, this%interpolation )

  ! local gamma factor
  pp = ptrcur
  do i=1,np
    gpart(i) = sqrt( ( ( 1.0_p_k_part + this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

    pp = pp + 1
  enddo

  ! K is the value of dx/dt (used to update positions)
  ! K(step) -> kvec(step,:,:)
  pp = ptrcur
  do i=1,np
    kvec(step,1,i) = (alphaptemp(i) * this%p(1,pp) / gpart(i) &
              - betaptemp(1, i)) / hvecptemp(1, i)
    kvec(step,2,i) = (alphaptemp(i) * this%p(2,pp) / gpart(i) &
              - betaptemp(2, i)) / hvecptemp(2, i)
    kvec(step,3,i) = (alphaptemp(i) * this%p(3,pp) / gpart(i) &
              - betaptemp(3, i)) / hvecptemp(3, i)

    pp = pp + 1
  end do

  ! Lorentz Force contribution
  pp = ptrcur
  do i=1,np
    FLr(i) = (eptemp(1,i) + &
           (this%p(2,pp)*bptemp(3,i)-this%p(3,pp)*bptemp(2,i))/ gpart(i) )
    FLt(i) = (eptemp(2,i) + &
           (this%p(3,pp)*bptemp(1,i)-this%p(1,pp)*bptemp(3,i))/ gpart(i) )
    FLp(i) = (eptemp(3,i) + &
           (this%p(1,pp)*bptemp(2,i)-this%p(2,pp)*bptemp(1,i))/ gpart(i) )
    uE(i) = this%p(1,pp)*eptemp(1,i) +&
             this%p(2,pp)*eptemp(2,i) +&
             this%p(3,pp)*eptemp(3,i)
    auxRR(i) = uE(i)*uE(i)/(gpart(i)*gpart(i)) - &
               FLr(i)*FLr(i) - FLt(i)*FLt(i) - FLp(i)*FLp(i)
    pp = pp + 1
  end do

  ! Radiation Reaction contribution
  pp = ptrcur
  do i=1,np
    RRr(i) = eptemp(2,i)*bptemp(3,i) - eptemp(3,i)*bptemp(2,i) + &
             (bptemp(1,i)*(this%p(3,pp)*bptemp(3,i) + this%p(2,pp)*bptemp(2,i)) - &
             this%p(1,pp) * (bptemp(2,i)*bptemp(2,i) + bptemp(3,i)*bptemp(3,i)) + &
             eptemp(1,i)*uE(i) )/gpart(i) + &
             gpart(i)*this%p(1,pp)*auxRR(i)

    RRt(i) = eptemp(3,i)*bptemp(1,i) - eptemp(1,i)*bptemp(3,i) + &
             (bptemp(2,i)*(this%p(3,pp)*bptemp(3,i) + this%p(1,pp)*bptemp(1,i)) - &
             this%p(2,pp) * (bptemp(1,i)*bptemp(1,i) + bptemp(3,i)*bptemp(3,i)) + &
             eptemp(2,i)*uE(i) )/gpart(i) + &
             gpart(i)*this%p(2,pp)*auxRR(i)

    RRp(i) = eptemp(1,i)*bptemp(2,i) - eptemp(2,i)*bptemp(1,i) + &
             (bptemp(3,i)*(this%p(2,pp)*bptemp(2,i) + this%p(1,pp)*bptemp(1,i)) - &
             this%p(3,pp) * (bptemp(1,i)*bptemp(1,i) + bptemp(2,i)*bptemp(2,i)) + &
             eptemp(3,i)*uE(i) )/gpart(i) + &
             gpart(i)*this%p(3,pp)*auxRR(i)
    pp = pp + 1
  end do

  ! L is the value of dw/dt (used to update momenta)
  ! L(step) -> lvec(step,:,:)
  pp = ptrcur
  do i=1,np
    lvec(step,1,i) = alphaptemp(i)*((FLr(i)/this%rqm+RRr(i)*this%geometry%rr_norm) +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(1,i) +&
           htensorptemp(1,1,i) * this%p(1,pp) +&
           htensorptemp(1,2,i) * this%p(2,pp) +&
           htensorptemp(1,3,i) * this%p(3,pp) )) -&
           (cfsrptemp(1,2,i) * kvec(step,1,i)*this%p(2,pp) +&
            cfsrptemp(1,3,i) * kvec(step,1,i)*this%p(3,pp) +&
            cfsrptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfsrptemp(2,3,i) * kvec(step,2,i)*this%p(3,pp) +&
            cfsrptemp(3,2,i) * kvec(step,3,i)*this%p(2,pp) +&
            cfsrptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )

    lvec(step,2,i) = alphaptemp(i)*((FLt(i)/this%rqm+RRt(i)*this%geometry%rr_norm) +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(2,i) +&
           htensorptemp(2,1,i) * this%p(1,pp) +&
           htensorptemp(2,2,i) * this%p(2,pp) +&
           htensorptemp(2,3,i) * this%p(3,pp) )) -&
           (cfstptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfstptemp(1,3,i) * kvec(step,1,i)*this%p(3,pp) +&
            cfstptemp(2,1,i) * kvec(step,2,i)*this%p(1,pp) +&
            cfstptemp(2,3,i) * kvec(step,2,i)*this%p(3,pp) +&
            cfstptemp(3,1,i) * kvec(step,3,i)*this%p(1,pp) +&
            cfstptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )

    lvec(step,3,i) = alphaptemp(i)*((FLp(i)/this%rqm+RRp(i)*this%geometry%rr_norm) +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(3,i) +&
           htensorptemp(3,1,i) * this%p(1,pp) +&
           htensorptemp(3,2,i) * this%p(2,pp) +&
           htensorptemp(3,3,i) * this%p(3,pp) )) -&
           (cfspptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfspptemp(1,2,i) * kvec(step,1,i)*this%p(2,pp) +&
            cfspptemp(2,1,i) * kvec(step,2,i)*this%p(1,pp) +&
            cfspptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfspptemp(3,1,i) * kvec(step,3,i)*this%p(1,pp) +&
            cfspptemp(3,2,i) * kvec(step,3,i)*this%p(2,pp) )

    pp = pp + 1
  end do

end subroutine get_kvec_lvec_rk_gr_rr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine get_kvec_lvec_rk_cart_gr_rr( this, emf, xbuf, xbufs, kvec, lvec, ptrcur, np, order, step )

  use m_emf_interpolate, only : get_emf
  use m_geometry_gr

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(inout) :: xbuf
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(in) :: xbufs
  integer, intent(in) :: ptrcur, np, step, order
  real(p_double), dimension(order ,p_p_dim,p_cache_size), intent(inout) :: kvec, lvec


  integer :: i, pp
  real(p_double) :: r, st, ct, sp, cp, alpha

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bptemp, eptemp, ep, bp
  real(p_double), dimension(p_cache_size) :: alphaptemp, gpart, smdet
  real(p_double), dimension(p_p_dim,p_cache_size) :: betaptemp
  real(p_double), dimension(p_p_dim,p_cache_size) :: gravptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: smcov, smctr
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfsxptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfsyptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: cfszptemp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size) :: htensorptemp
  real(p_double), dimension(p_cache_size) :: K1,K2,K3,K4,K5,K6
  real(p_double), dimension(p_cache_size) :: FLx,FLy,FLz,uE
  real(p_double), dimension(p_cache_size) :: auxEB1,auxEB2,auxEB3,auxRR
  real(p_double), dimension(p_cache_size) :: RRx,RRy,RRz

  ! executable statements

  ! get local gr quantaties (analitical expressions)
  select case (this%geometry%metric)
  case (p_geometry_minkowski)
    call get_geometric_quant_sph_cart(this, alphaptemp, betaptemp, smcov, smctr, smdet,&
    cfsxptemp, cfsyptemp, cfszptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_schwarzschild)
    call get_geometric_quant_sch_cart(this, alphaptemp, betaptemp, smcov, smctr, smdet,&
    cfsxptemp, cfsyptemp, cfszptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_kerr_slow)
    call get_geometric_quant_kslow_cart(this, alphaptemp, betaptemp, smcov, smctr, smdet,&
    cfsxptemp, cfsyptemp, cfszptemp, gravptemp, htensorptemp, xbuf, np)!, eptemp, bptemp )
  case (p_geometry_kerr)
    !not implemented yet
  end select

  ! get em field in particle's local position (interpolation)
  call get_emf( emf, bptemp, eptemp, this%ix(:,ptrcur:), this%x(:,ptrcur:), &
    np, this%interpolation )

  ! convert spherical fields to cartesian
  pp = ptrcur
  do i=1,np
    r = xbufs(1,i)

    st = sin(xbufs(2,i))
    ct = cos(xbufs(2,i))

    sp = sin(xbufs(3,i))
    cp = cos(xbufs(3,i))

    alpha = sqrt(1._p_double - this%geometry%rs/r)

    ep(1,i) = st * cp * eptemp(1,i) * alpha + ct * cp * eptemp(2,i) - sp * eptemp(3,i)
    ep(2,i) = st * sp * eptemp(1,i) * alpha + ct * sp * eptemp(2,i) + cp * eptemp(3,i)
    ep(3,i) = ct * eptemp(1,i) * alpha - st * eptemp(2,i)

    bp(1,i) = st * cp * bptemp(1,i) * alpha + ct * cp * bptemp(2,i) - sp * bptemp(3,i)
    bp(2,i) = st * sp * bptemp(1,i) * alpha + ct * sp * bptemp(2,i) + cp * bptemp(3,i)
    bp(3,i) = ct * bptemp(1,i) * alpha - st * bptemp(2,i)

    pp = pp + 1
  enddo

  ! local gamma factor
  pp = ptrcur
  do i=1,np
    gpart(i) = sqrt( ( ( 1.0_p_k_part + smcov(1,1,i) * this%p(1,pp)**2 ) +&
          smcov(2,2,i) * this%p(2,pp)**2 ) +&
          smcov(3,3,i) * this%p(3,pp)**2 +&
          2._p_double * smcov(1,2,i) * this%p(1,pp) * this%p(2,pp) +&
          2._p_double * smcov(1,3,i) * this%p(1,pp) * this%p(3,pp) +&
          2._p_double * smcov(2,3,i) * this%p(2,pp) * this%p(3,pp) )

    pp = pp + 1
  enddo

  ! K is the value of dx/dt (used to update positions)
  ! K(step) -> kvec(step,:,:)
  pp = ptrcur
  do i=1,np
    kvec(step,1,i) = alphaptemp(i) * this%p(1,pp) / gpart(i) &
              - betaptemp(1, i)
    kvec(step,2,i) = alphaptemp(i) * this%p(2,pp) / gpart(i) &
              - betaptemp(2, i)
    kvec(step,3,i) = alphaptemp(i) * this%p(3,pp) / gpart(i) &
              - betaptemp(3, i)

    pp = pp + 1
  end do

  ! Lorentz Force contribution (w/o rqm)
  pp = ptrcur
  do i=1,np
    FLx(i) = (ep(1,i) + sqrt(smdet(i)) * &
           (smctr(1,1,i)*(this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)) +&
            smctr(1,2,i)*(this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)) +&
            smctr(1,3,i)*(this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)))/ gpart(i) )

    FLy(i) = (ep(2,i) + sqrt(smdet(i)) * &
           (smctr(2,1,i)*(this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)) +&
            smctr(2,2,i)*(this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)) +&
            smctr(2,3,i)*(this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)))/ gpart(i) )

    FLz(i) = (ep(3,i) + sqrt(smdet(i)) * &
           (smctr(3,1,i)*(this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)) +&
            smctr(3,2,i)*(this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)) +&
            smctr(3,3,i)*(this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)))/ gpart(i) )

    uE(i) = ep(1,i) * (smcov(1,1,i)*this%p(1,pp) +&
                       smcov(1,2,i)*this%p(2,pp) +&
                       smcov(1,3,i)*this%p(3,pp)) + &
            ep(2,i) * (smcov(2,1,i)*this%p(1,pp) +&
                       smcov(2,2,i)*this%p(2,pp) +&
                       smcov(2,3,i)*this%p(3,pp)) + &
            ep(3,i) * (smcov(3,1,i)*this%p(1,pp) +&
                       smcov(3,2,i)*this%p(2,pp) +&
                       smcov(3,3,i)*this%p(3,pp))
    pp = pp + 1
  end do

  ! evaluate auxiliary quantities
  pp = ptrcur
  do i=1,np
    K1(i) = smctr(1,1,i)*smctr(2,2,i) - smctr(1,2,i)*smctr(1,2,i)
    K2(i) = smctr(1,1,i)*smctr(2,3,i) - smctr(1,2,i)*smctr(1,3,i)
    K3(i) = smctr(1,1,i)*smctr(3,3,i) - smctr(1,3,i)*smctr(1,3,i)
    K4(i) = smctr(1,2,i)*smctr(3,3,i) - smctr(1,3,i)*smctr(2,3,i)
    K5(i) = smctr(1,2,i)*smctr(2,3,i) - smctr(1,3,i)*smctr(2,2,i)
    K6(i) = smctr(2,2,i)*smctr(3,3,i) - smctr(2,3,i)*smctr(2,3,i)

    auxEB1(i) = ep(2,i)*bp(3,i)-ep(3,i)*bp(2,i)
    auxEB2(i) = ep(3,i)*bp(1,i)-ep(1,i)*bp(3,i)
    auxEB3(i) = ep(1,i)*bp(2,i)-ep(2,i)*bp(1,i)

    auxRR(i) = uE(i)*uE(i)/(gpart(i)*gpart(i)) -&
                  smcov(1,1,i)*FLx(i)*FLx(i) - smcov(2,2,i)*FLy(i)*FLy(i) -&
                  smcov(3,3,i)*FLz(i)*FLz(i) - 2._p_double*smcov(1,2,i)*FLx(i)*FLy(i) -&
                  2._p_double*smcov(1,3,i)*FLx(i)*FLz(i) - 2._p_double*smcov(2,3,i)*FLy(i)*FLz(i)
    pp = pp + 1
  end do

  ! Radiation Reaction contribution
  pp = ptrcur
  do i=1,np
    RRx(i) = sqrt(smdet(i)) * &
              (smctr(1,1,i)*auxEB1(i) + &
               smctr(1,2,i)*auxEB2(i) + &
               smctr(1,3,i)*auxEB3(i) ) + &
             (smdet(i) * &
                (bp(1,i) * (bp(1,i)*(-this%p(2,pp)*K4(i)+this%p(3,pp)*K5(i)) +&
                            bp(2,i)*( this%p(1,pp)*K4(i)+this%p(2,pp)*K3(i)-&
                                                        this%p(3,pp)*K2(i)) +&
                            bp(3,i)*(-this%p(1,pp)*K5(i)-this%p(2,pp)*K2(i)+&
                                                        this%p(3,pp)*K1(i))) +&
                 bp(2,i) * (bp(2,i)*(-this%p(1,pp)*K3(i)) +&
                            bp(3,i)*(2*this%p(1,pp)*K2(i))) +&
                 bp(3,i) * (bp(3,i)*(-this%p(1,pp)*K1(i)))) +&
              ep(1,i)*uE(i) ) / gpart(i) +&
             gpart(i)*this%p(1,pp)*auxRR(i)


    RRy(i) = sqrt(smdet(i)) * &
              (smctr(2,1,i)*auxEB1(i) + &
               smctr(2,2,i)*auxEB2(i) + &
               smctr(2,3,i)*auxEB3(i) ) + &
              (smdet(i) * &
                (bp(1,i) * (bp(1,i)*(-this%p(2,pp)*K6(i)) +&
                            bp(2,i)*(this%p(1,pp)*K6(i)+this%p(2,pp)*K4(i)+&
                                                        this%p(3,pp)*K5(i)) +&
                            bp(3,i)*(-2*this%p(2,pp)*K5(i))) +&
                  bp(2,i) * (bp(2,i)*(-this%p(1,pp)*K4(i)-this%p(3,pp)*K2(i)) +&
                            bp(3,i)*(this%p(1,pp)*K5(i)+this%p(2,pp)*K2(i)+&
                                                        this%p(3,pp)*K1(i))) +&
                  bp(3,i) * (bp(3,i)*(-this%p(2,pp)*K1(i)))) +&
              ep(2,i)*uE(i) ) / gpart(i) +&
              gpart(i)*this%p(2,pp)*auxRR(i)

    RRz(i) = sqrt(smdet(i)) * &
              (smctr(3,1,i)*auxEB1(i) + &
               smctr(3,2,i)*auxEB2(i) + &
               smctr(3,3,i)*auxEB3(i) ) + &
              (smdet(i) * &
                (bp(1,i) * (bp(1,i)*(-this%p(3,pp)*K6(i)) +&
                            bp(2,i)*(2*this%p(3,pp)*K4(i)) +&
                            bp(3,i)*(this%p(1,pp)*K6(i)-this%p(2,pp)*K4(i)-&
                                                        this%p(3,pp)*K5(i))) +&
                  bp(2,i) * (bp(2,i)*(-this%p(3,pp)*K3(i)) +&
                            bp(3,i)*(-this%p(1,pp)*K4(i)+this%p(2,pp)*K3(i)+&
                                                        this%p(3,pp)*K2(i))) +&
                  bp(3,i) * (bp(3,i)*(this%p(1,pp)*K5(i)-this%p(2,pp)*K2(i)))) +&
              ep(3,i)*uE(i) ) / gpart(i) +&
              gpart(i)*this%p(3,pp)*auxRR(i)
    pp = pp + 1
  end do

  ! L is the value of dw/dt (used to update momenta)
  ! L(step) -> lvec(step,:,:)
  pp = ptrcur
  do i=1,np
    lvec(step,1,i) = alphaptemp(i)*( (FLx(i)/this%rqm + RRx(i) * this%geometry%rr_norm) +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(1,i) +&
           this%p(1,pp) * (htensorptemp(1,1,i) * smctr(1,1,i) +&
                           htensorptemp(2,1,i) * smctr(1,2,i) +&
                           htensorptemp(3,1,i) * smctr(1,3,i)) +&
           this%p(2,pp) * (htensorptemp(1,2,i) * smctr(1,1,i) +&
                           htensorptemp(2,2,i) * smctr(1,2,i) +&
                           htensorptemp(3,2,i) * smctr(1,3,i)) +&
           this%p(3,pp) * (htensorptemp(1,3,i) * smctr(1,1,i) +&
                           htensorptemp(2,3,i) * smctr(1,2,i) +&
                           htensorptemp(3,3,i) * smctr(1,3,i)) )) -&
           (cfsxptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfsxptemp(1,2,i) * (kvec(step,1,i)*this%p(2,pp) + kvec(step,2,i)*this%p(1,pp)) +&
            cfsxptemp(1,3,i) * (kvec(step,1,i)*this%p(3,pp) + kvec(step,3,i)*this%p(1,pp)) +&
            cfsxptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfsxptemp(2,3,i) * (kvec(step,2,i)*this%p(3,pp) + kvec(step,3,i)*this%p(2,pp)) +&
            cfsxptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )

    lvec(step,2,i) = alphaptemp(i)*( (FLy(i)/this%rqm + RRy(i) * this%geometry%rr_norm) +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(2,i) +&
           this%p(1,pp) * (htensorptemp(1,1,i) * smctr(2,1,i) +&
                           htensorptemp(2,1,i) * smctr(2,2,i) +&
                           htensorptemp(3,1,i) * smctr(2,3,i)) +&
           this%p(2,pp) * (htensorptemp(1,2,i) * smctr(2,1,i) +&
                           htensorptemp(2,2,i) * smctr(2,2,i) +&
                           htensorptemp(3,2,i) * smctr(2,3,i)) +&
           this%p(3,pp) * (htensorptemp(1,3,i) * smctr(2,1,i) +&
                           htensorptemp(2,3,i) * smctr(2,2,i) +&
                           htensorptemp(3,3,i) * smctr(2,3,i)) )) -&
           (cfsyptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfsyptemp(1,2,i) * (kvec(step,1,i)*this%p(2,pp) + kvec(step,2,i)*this%p(1,pp)) +&
            cfsyptemp(1,3,i) * (kvec(step,1,i)*this%p(3,pp) + kvec(step,3,i)*this%p(1,pp)) +&
            cfsyptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfsyptemp(2,3,i) * (kvec(step,2,i)*this%p(3,pp) + kvec(step,3,i)*this%p(2,pp)) +&
            cfsyptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )

    lvec(step,3,i) = alphaptemp(i)*( (FLz(i)/this%rqm + RRz(i) * this%geometry%rr_norm) +&
           this%geometry%grav_switch * ( gpart(i) * gravptemp(3,i) +&
           this%p(1,pp) * (htensorptemp(1,1,i) * smctr(3,1,i) +&
                           htensorptemp(2,1,i) * smctr(3,2,i) +&
                           htensorptemp(3,1,i) * smctr(3,3,i)) +&
           this%p(2,pp) * (htensorptemp(1,2,i) * smctr(3,1,i) +&
                           htensorptemp(2,2,i) * smctr(3,2,i) +&
                           htensorptemp(3,2,i) * smctr(3,3,i)) +&
           this%p(3,pp) * (htensorptemp(1,3,i) * smctr(3,1,i) +&
                           htensorptemp(2,3,i) * smctr(3,2,i) +&
                           htensorptemp(3,3,i) * smctr(3,3,i)) )) -&
           (cfszptemp(1,1,i) * kvec(step,1,i)*this%p(1,pp) +&
            cfszptemp(1,2,i) * (kvec(step,1,i)*this%p(2,pp) + kvec(step,2,i)*this%p(1,pp)) +&
            cfszptemp(1,3,i) * (kvec(step,1,i)*this%p(3,pp) + kvec(step,3,i)*this%p(1,pp)) +&
            cfszptemp(2,2,i) * kvec(step,2,i)*this%p(2,pp) +&
            cfszptemp(2,3,i) * (kvec(step,2,i)*this%p(3,pp) + kvec(step,3,i)*this%p(2,pp)) +&
            cfszptemp(3,3,i) * kvec(step,3,i)*this%p(3,pp) )
    pp = pp + 1
  end do

end subroutine get_kvec_lvec_rk_cart_gr_rr
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Gets position dependent quantities for GR particle push - Spherical
!-------------------------------------------------------------------------------
subroutine get_geometric_quant_sph(this, alpha, beta, hvec, cfsr, cfst, cfsp, grav, htensor, xbuf, np)!, ep, bp)

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(in) :: this
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(in) :: xbuf
  integer, intent(in) :: np
  real(p_double), dimension(p_cache_size), intent(inout) :: alpha
  real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: beta, hvec, grav
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsr
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfst
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: htensor
  !real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: bp, ep

  integer :: i
  real(p_double) :: shift, tc, tc2, tc3, tc4, tc5, tc6, tc7, tc8, tc9, tc10
  real(p_double) :: st, ct, auxr

  ! used to Taylor expand trig functions
  shift = pi_2

  do i=1, np
    ! temporary functions
    ! evaluate the sin and cos of theta around pi/2
    tc = xbuf(2,i) - shift
    tc2 = tc * tc
    tc3 = tc2 * tc
    tc4 = tc3 * tc
    tc5 = tc4 * tc
    tc6 = tc5 * tc
    tc7 = tc6 * tc
    tc8 = tc7 * tc
    tc9 = tc8 * tc
    tc10 = tc9 * tc

    st = 1.0_p_double - 0.500000000000000_p_double * tc2 +&
                        0.041666666666667_p_double * tc4 - 0.001388888888889_p_double * tc6 +&
                        0.000024801587302_p_double * tc8 - 0.000000275573192_p_double * tc10

    ct = - tc + 0.166666666666667_p_double * tc3 - 0.008333333333333_p_double * tc5 +&
                0.000198412698413_p_double * tc7 - 0.000002755731922_p_double * tc9

    ! auxiliar 1/r function
    auxr = 1._p_double / xbuf(1,i)

    ! lapse function
    alpha(i) = 1._p_double

    ! shift vector
    beta(1,i) = 0._p_double
    beta(2,i) = 0._p_double
    beta(3,i) = 0._p_double

    ! normalizing vector
    hvec(1,i) = 1._p_double
    hvec(2,i) = xbuf(1,i)
    hvec(3,i) = xbuf(1,i) * st !(singularity)

    ! radial component Christoffel symb (normalized)
    cfsr(1,1,i) = 0._p_double
    cfsr(1,2,i) = 0._p_double
    cfsr(1,3,i) = 0._p_double
    cfsr(2,1,i) = 0._p_double
    cfsr(2,2,i) = -1._p_double
    cfsr(2,3,i) = 0._p_double
    cfsr(3,1,i) = 0._p_double
    cfsr(3,2,i) = 0._p_double
    cfsr(3,3,i) = -st

    ! poloidal component Christoffel symb (normalized)
    cfst(1,1,i) = 0._p_double
    cfst(1,2,i) = 1._p_double * auxr
    cfst(1,3,i) = 0._p_double
    cfst(2,1,i) = 1._p_double
    cfst(2,2,i) = 0._p_double
    cfst(2,3,i) = 0._p_double
    cfst(3,1,i) = 0._p_double
    cfst(3,2,i) = 0._p_double
    cfst(3,3,i) = -ct

    ! azimuthal component Christoffel symb (normalized)
    cfsp(1,1,i) = 0._p_double
    cfsp(1,2,i) = 0._p_double
    cfsp(1,3,i) = 1._p_double * auxr
    cfsp(2,1,i) = 0._p_double
    cfsp(2,2,i) = 0._p_double
    cfsp(2,3,i) = ct / st !(singularity)
    cfsp(3,1,i) = st
    cfsp(3,2,i) = ct
    cfsp(3,3,i) = 0._p_double

    ! gravitational acceleration
    grav(1,i) = 0._p_double
    grav(2,i) = 0._p_double
    grav(3,i) = 0._p_double

    ! gravitomagnetic tensor
    htensor(1,1,i) = 0._p_double
    htensor(1,2,i) = 0._p_double
    htensor(1,3,i) = 0._p_double
    htensor(2,1,i) = 0._p_double
    htensor(2,2,i) = 0._p_double
    htensor(2,3,i) = 0._p_double
    htensor(3,1,i) = 0._p_double
    htensor(3,2,i) = 0._p_double
    htensor(3,3,i) = 0._p_double

    !ep(1,i) = 0._p_double
    !ep(2,i) = 0._p_double
    !ep(3,i) = 0._p_double
    !bp(1,i) = ct
    !bp(2,i) = -st
    !bp(3,i) = 0._p_double
  end do


end subroutine get_geometric_quant_sph
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Gets position dependent quantities for GR particle push - Schwarzschild
!-------------------------------------------------------------------------------
subroutine get_geometric_quant_sch(this, alpha, beta, hvec, cfsr, cfst, cfsp, grav, htensor, xbuf, np)!, ep, bp)

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(in) :: this
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(in) :: xbuf
  integer, intent(in) :: np
  real(p_double), dimension(p_cache_size), intent(inout) :: alpha
  real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: beta, hvec, grav
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsr
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfst
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: htensor
  !real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: bp, ep

  integer :: i
  real(p_double) :: shift, tc, tc2, tc3, tc4, tc5, tc6, tc7, tc8, tc9, tc10
  real(p_double) :: st, ct, auxr, alpha2, rsloc

  ! used to Taylor expand trig functions
  shift = pi_2

  ! avoid read memory in loop
  rsloc = this%geometry%rs

  do i=1, np
    ! temporary functions
    ! evaluate the sin and cos of theta around pi/2
    tc = xbuf(2,i) - shift
    tc2 = tc * tc
    tc3 = tc2 * tc
    tc4 = tc3 * tc
    tc5 = tc4 * tc
    tc6 = tc5 * tc
    tc7 = tc6 * tc
    tc8 = tc7 * tc
    tc9 = tc8 * tc
    tc10 = tc9 * tc

    st = 1.0_p_double - 0.500000000000000_p_double * tc2 +&
                        0.041666666666667_p_double * tc4 - 0.001388888888889_p_double * tc6 +&
                        0.000024801587302_p_double * tc8 - 0.000000275573192_p_double * tc10

    ct = - tc + 0.166666666666667_p_double * tc3 - 0.008333333333333_p_double * tc5 +&
                0.000198412698413_p_double * tc7 - 0.000002755731922_p_double * tc9

    ! auxiliar 1/r function
    auxr = 1._p_double / xbuf(1,i)

    alpha2 = 1._p_double - rsloc * auxr

    ! lapse function
    alpha(i) = sqrt(alpha2)

    ! shift vector
    beta(1,i) = 0._p_double
    beta(2,i) = 0._p_double
    beta(3,i) = 0._p_double

    ! normalizing vector
    hvec(1,i) = 1._p_double / alpha(i)
    hvec(2,i) = xbuf(1,i)
    hvec(3,i) = xbuf(1,i) * st !(singularity)

    ! radial component Christoffel symb (normalized)
    cfsr(1,1,i) = - rsloc / (2._p_double * xbuf(1,i)**2 * alpha2)
    cfsr(1,2,i) = 0._p_double
    cfsr(1,3,i) = 0._p_double
    cfsr(2,1,i) = 0._p_double
    cfsr(2,2,i) = -alpha(i)
    cfsr(2,3,i) = 0._p_double
    cfsr(3,1,i) = 0._p_double
    cfsr(3,2,i) = 0._p_double
    cfsr(3,3,i) = -alpha(i) * st

    ! poloidal component Christoffel symb (normalized)
    cfst(1,1,i) = 0._p_double
    cfst(1,2,i) = 1._p_double * auxr
    cfst(1,3,i) = 0._p_double
    cfst(2,1,i) = alpha(i)
    cfst(2,2,i) = 0._p_double
    cfst(2,3,i) = 0._p_double
    cfst(3,1,i) = 0._p_double
    cfst(3,2,i) = 0._p_double
    cfst(3,3,i) = -ct

    ! azimuthal component Christoffel symb (normalized)
    cfsp(1,1,i) = 0._p_double
    cfsp(1,2,i) = 0._p_double
    cfsp(1,3,i) = 1._p_double * auxr
    cfsp(2,1,i) = 0._p_double
    cfsp(2,2,i) = 0._p_double
    cfsp(2,3,i) = ct / st !(singularity)
    cfsp(3,1,i) = alpha(i) * st
    cfsp(3,2,i) = ct
    cfsp(3,3,i) = 0._p_double

    ! gravitational acceleration
    grav(1,i) = - rsloc / (2._p_double * xbuf(1,i)**2 * alpha(i))
    grav(2,i) = 0._p_double
    grav(3,i) = 0._p_double

    ! gravitomagnetic tensor
    htensor(1,1,i) = 0._p_double
    htensor(1,2,i) = 0._p_double
    htensor(1,3,i) = 0._p_double
    htensor(2,1,i) = 0._p_double
    htensor(2,2,i) = 0._p_double
    htensor(2,3,i) = 0._p_double
    htensor(3,1,i) = 0._p_double
    htensor(3,2,i) = 0._p_double
    htensor(3,3,i) = 0._p_double

    !ep(1,i) = 0._p_double
    !ep(2,i) = 0._p_double
    !ep(3,i) = 0._p_double
    !bp(1,i) = 10._p_double*ct
    !bp(2,i) = -10._p_double*st
    !bp(3,i) = 0._p_double
  end do


end subroutine get_geometric_quant_sch
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Gets position dependent quantities for GR particle push - Kerr-Slow
!-------------------------------------------------------------------------------
subroutine get_geometric_quant_kslow(this, alpha, beta, hvec, cfsr, cfst, cfsp, grav, htensor, xbuf, np)!, ep, bp)

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(in) :: this
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(in) :: xbuf
  integer, intent(in) :: np
  real(p_double), dimension(p_cache_size), intent(inout) :: alpha
  real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: beta, hvec, grav
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsr
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfst
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsp
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: htensor
  !real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: bp, ep

  integer :: i
  real(p_double) :: shift, tc, tc2, tc3, tc4, tc5, tc6, tc7, tc8, tc9, tc10
  real(p_double) :: st, ct, auxr, auxr2, auxr3, alpha2, rsloc

  ! used to Taylor expand trig functions
  shift = pi_2

  ! avoid read memory in loop
  rsloc = this%geometry%rs

  do i=1, np
    ! temporary functions
    ! evaluate the sin and cos of theta around pi/2
    tc = xbuf(2,i) - shift
    tc2 = tc * tc
    tc3 = tc2 * tc
    tc4 = tc3 * tc
    tc5 = tc4 * tc
    tc6 = tc5 * tc
    tc7 = tc6 * tc
    tc8 = tc7 * tc
    tc9 = tc8 * tc
    tc10 = tc9 * tc

    st = 1.0_p_double - 0.500000000000000_p_double * tc2 +&
                        0.041666666666667_p_double * tc4 - 0.001388888888889_p_double * tc6 +&
                        0.000024801587302_p_double * tc8 - 0.000000275573192_p_double * tc10

    ct = - tc + 0.166666666666667_p_double * tc3 - 0.008333333333333_p_double * tc5 +&
                0.000198412698413_p_double * tc7 - 0.000002755731922_p_double * tc9

    ! auxiliar 1/r function
    auxr = 1._p_double / xbuf(1,i)
    auxr2 = auxr * auxr
    auxr3 = auxr * auxr2

    alpha2 = 1._p_double - rsloc * auxr

    ! lapse function
    alpha(i) = sqrt(alpha2)

    ! shift vector
    beta(1,i) = 0._p_double
    beta(2,i) = 0._p_double
    beta(3,i) = - this%geometry%beta0 * st * auxr2

    ! normalizing vector
    hvec(1,i) = 1._p_double / alpha(i)
    hvec(2,i) = xbuf(1,i)
    hvec(3,i) = xbuf(1,i) * st !(singularity)

    ! radial component Christoffel symb (normalized)
    cfsr(1,1,i) = - rsloc * auxr2 / (2._p_double * alpha2)
    cfsr(1,2,i) = 0._p_double
    cfsr(1,3,i) = 0._p_double
    cfsr(2,1,i) = 0._p_double
    cfsr(2,2,i) = -alpha(i)
    cfsr(2,3,i) = 0._p_double
    cfsr(3,1,i) = 0._p_double
    cfsr(3,2,i) = 0._p_double
    cfsr(3,3,i) = -alpha(i) * st

    ! poloidal component Christoffel symb (normalized)
    cfst(1,1,i) = 0._p_double
    cfst(1,2,i) = 1._p_double * auxr
    cfst(1,3,i) = 0._p_double
    cfst(2,1,i) = alpha(i)
    cfst(2,2,i) = 0._p_double
    cfst(2,3,i) = 0._p_double
    cfst(3,1,i) = 0._p_double
    cfst(3,2,i) = 0._p_double
    cfst(3,3,i) = -ct

    ! azimuthal component Christoffel symb (normalized)
    cfsp(1,1,i) = 0._p_double
    cfsp(1,2,i) = 0._p_double
    cfsp(1,3,i) = 1._p_double * auxr
    cfsp(2,1,i) = 0._p_double
    cfsp(2,2,i) = 0._p_double
    cfsp(2,3,i) = ct / st !(singularity)
    cfsp(3,1,i) = alpha(i) * st
    cfsp(3,2,i) = ct
    cfsp(3,3,i) = 0._p_double

    ! gravitational acceleration
    grav(1,i) = - rsloc * auxr2 / (2._p_double * alpha(i))
    grav(2,i) = 0._p_double
    grav(3,i) = 0._p_double

    ! gravitomagnetic tensor
    htensor(1,1,i) = 0._p_double
    htensor(1,2,i) = 0._p_double
    htensor(1,3,i) = 2._p_double * this%geometry%beta0 * st * auxr3
    htensor(2,1,i) = 0._p_double
    htensor(2,2,i) = 0._p_double
    htensor(2,3,i) = - this%geometry%beta0 * ct * auxr3 / (alpha(i))
    htensor(3,1,i) = this%geometry%beta0 * st * auxr3
    htensor(3,2,i) = this%geometry%beta0 * ct * auxr3 / (alpha(i))
    htensor(3,3,i) = 0._p_double

    !ep(1,i) = 0._p_double
    !ep(2,i) = 0._p_double
    !ep(3,i) = 0._p_double
    !bp(1,i) = 10._p_double*ct
    !bp(2,i) = -10._p_double*st
    !bp(3,i) = 0._p_double
  end do


end subroutine get_geometric_quant_kslow
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Gets position dependent quantities for GR particle push - Spherical (cart)
!-------------------------------------------------------------------------------
subroutine get_geometric_quant_sph_cart(this, alpha, beta, smcov, smctr, smdet, cfsx, cfsy, cfsz, grav, htensor, xbuf, np)!, ep, bp)

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(in) :: this
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(in) :: xbuf
  integer, intent(in) :: np
  real(p_double), dimension(p_cache_size), intent(inout) :: alpha, smdet
  real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: beta, grav
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: smcov, smctr
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsx
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsy
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsz
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: htensor
  !real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: bp, ep

  ! lapse function
  alpha(:) = 1._p_double

  ! shift vector
  beta(1,:) = 0._p_double
  beta(2,:) = 0._p_double
  beta(3,:) = 0._p_double

  ! spatial metric tensor (covariant)
  smcov(1,1,:) = 1._p_double
  smcov(1,2,:) = 0._p_double
  smcov(1,3,:) = 0._p_double
  smcov(2,1,:) = 0._p_double
  smcov(2,2,:) = 1._p_double
  smcov(2,3,:) = 0._p_double
  smcov(3,1,:) = 0._p_double
  smcov(3,2,:) = 0._p_double
  smcov(3,3,:) = 1._p_double

  ! spatial metric determinant
  smdet(:) = 1._p_double

  ! spatial metric tensor (contravariant)
  smctr(1,1,:) = 1._p_double
  smctr(1,2,:) = 0._p_double
  smctr(1,3,:) = 0._p_double
  smctr(2,1,:) = 0._p_double
  smctr(2,2,:) = 1._p_double
  smctr(2,3,:) = 0._p_double
  smctr(3,1,:) = 0._p_double
  smctr(3,2,:) = 0._p_double
  smctr(3,3,:) = 1._p_double

  ! radial component Christoffel symb (normalized)
  cfsx(1,1,:) = 0._p_double
  cfsx(1,2,:) = 0._p_double
  cfsx(1,3,:) = 0._p_double
  cfsx(2,1,:) = 0._p_double
  cfsx(2,2,:) = 0._p_double
  cfsx(2,3,:) = 0._p_double
  cfsx(3,1,:) = 0._p_double
  cfsx(3,2,:) = 0._p_double
  cfsx(3,3,:) = 0._p_double

  ! poloidal component Christoffel symb (normalized)
  cfsy(1,1,:) = 0._p_double
  cfsy(1,2,:) = 0._p_double
  cfsy(1,3,:) = 0._p_double
  cfsy(2,1,:) = 0._p_double
  cfsy(2,2,:) = 0._p_double
  cfsy(2,3,:) = 0._p_double
  cfsy(3,1,:) = 0._p_double
  cfsy(3,2,:) = 0._p_double
  cfsy(3,3,:) = 0._p_double

  ! azimuthal component Christoffel symb (normalized)
  cfsz(1,1,:) = 0._p_double
  cfsz(1,2,:) = 0._p_double
  cfsz(1,3,:) = 0._p_double
  cfsz(2,1,:) = 0._p_double
  cfsz(2,2,:) = 0._p_double
  cfsz(2,3,:) = 0._p_double
  cfsz(3,1,:) = 0._p_double
  cfsz(3,2,:) = 0._p_double
  cfsz(3,3,:) = 0._p_double

  ! gravitational acceleration
  grav(1,:) = 0._p_double
  grav(2,:) = 0._p_double
  grav(3,:) = 0._p_double

  ! gravitomagnetic tensor
  htensor(1,1,:) = 0._p_double
  htensor(1,2,:) = 0._p_double
  htensor(1,3,:) = 0._p_double
  htensor(2,1,:) = 0._p_double
  htensor(2,2,:) = 0._p_double
  htensor(2,3,:) = 0._p_double
  htensor(3,1,:) = 0._p_double
  htensor(3,2,:) = 0._p_double
  htensor(3,3,:) = 0._p_double

  ! e and b fields analytical (physical - spherical)
  !ep(1,:) = 0._p_double
  !ep(2,:) = 0._p_double
  !ep(3,:) = 0._p_double
  !bp(1,:) = 0._p_double
  !bp(2,:) = 0._p_double
  !bp(3,:) = 1._p_double


end subroutine get_geometric_quant_sph_cart
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Gets position dependent quantities for GR particle push - Schwarzschild (cart)
!-------------------------------------------------------------------------------
subroutine get_geometric_quant_sch_cart(this, alpha, beta, smcov, smctr, smdet, cfsx, cfsy, cfsz, grav, htensor, xbuf, np)!, ep, bp)

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(in) :: this
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(in) :: xbuf
  integer, intent(in) :: np
  real(p_double), dimension(p_cache_size), intent(inout) :: alpha, smdet
  real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: beta, grav
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: smcov, smctr
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsx
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsy
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsz
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: htensor
  !real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: bp, ep

  integer :: i
  real(p_double) :: r, fac1, fac2, fac3, alpha2
  real(p_double), dimension(p_p_dim,p_p_dim) :: M, cfsb


  do i=1, np
    ! temporary functions
    r = sqrt(xbuf(1,i)**2 + xbuf(2,i)**2 + xbuf(3,i)**2)
    alpha2 = 1 - this%geometry%rs / r

    M(1,1) = - (xbuf(2,i)**2 + xbuf(3,i)**2)
    M(1,2) = xbuf(1,i) * xbuf(2,i)
    M(1,3) = xbuf(1,i) * xbuf(3,i)
    M(2,1) = xbuf(2,i) * xbuf(1,i)
    M(2,2) = - (xbuf(1,i)**2 + xbuf(3,i)**2)
    M(2,3) = xbuf(2,i) * xbuf(3,i)
    M(3,1) = xbuf(3,i) * xbuf(1,i)
    M(3,2) = xbuf(3,i) * xbuf(2,i)
    M(3,3) = - (xbuf(1,i)**2 + xbuf(2,i)**2)

    fac1 = - this%geometry%rs * (3._p_double*r - 2._p_double*this%geometry%rs) / (2._p_double* r**5 * (r-this%geometry%rs))

    cfsb(1,1) = fac1 * (M(1,1) + r**3 / (3._p_double*r - 2._p_double*this%geometry%rs))
    cfsb(1,2) = fac1 * M(1,2)
    cfsb(1,3) = fac1 * M(1,3)
    cfsb(2,1) = fac1 * M(2,1)
    cfsb(2,2) = fac1 * (M(2,2) + r**3 / (3._p_double*r - 2._p_double*this%geometry%rs))
    cfsb(2,3) = fac1 * M(2,3)
    cfsb(3,1) = fac1 * M(3,1)
    cfsb(3,2) = fac1 * M(3,2)
    cfsb(3,3) = fac1 * (M(3,3) + r**3 / (3._p_double*r - 2._p_double*this%geometry%rs))

    ! lapse function
    alpha(i) = sqrt(alpha2)

    ! shift vector
    beta(1,i) = 0._p_double
    beta(2,i) = 0._p_double
    beta(3,i) = 0._p_double

    ! spatial metric tensor (covariant)
    fac2 = this%geometry%rs / (r**2 * (r - this%geometry%rs))

    smcov(1,1,i) = fac2 * (M(1,1) + r**3 / this%geometry%rs)
    smcov(1,2,i) = fac2 * M(1,2)
    smcov(1,3,i) = fac2 * M(1,3)
    smcov(2,1,i) = fac2 * M(2,1)
    smcov(2,2,i) = fac2 * (M(2,2) + r**3 / this%geometry%rs)
    smcov(2,3,i) = fac2 * M(2,3)
    smcov(3,1,i) = fac2 * M(3,1)
    smcov(3,2,i) = fac2 * M(3,2)
    smcov(3,3,i) = fac2 * (M(3,3) + r**3 / this%geometry%rs)

    ! spatial metric determinant
    smdet(i) = 1._p_double / alpha2

    ! spatial metric tensor (contravariant)
    smctr(1,1,i) = - this%geometry%rs * M(1,1) / r**3 + alpha2
    smctr(1,2,i) = - this%geometry%rs * M(1,2) / r**3
    smctr(1,3,i) = - this%geometry%rs * M(1,3) / r**3
    smctr(2,1,i) = - this%geometry%rs * M(2,1) / r**3
    smctr(2,2,i) = - this%geometry%rs * M(2,2) / r**3 + alpha2
    smctr(2,3,i) = - this%geometry%rs * M(2,3) / r**3
    smctr(3,1,i) = - this%geometry%rs * M(3,1) / r**3
    smctr(3,2,i) = - this%geometry%rs * M(3,2) / r**3
    smctr(3,3,i) = - this%geometry%rs * M(3,3) / r**3 + alpha2

    ! radial component Christoffel symb (normalized)
    cfsx(1,1,i) = cfsb(1,1) * xbuf(1,i)
    cfsx(1,2,i) = cfsb(1,2) * xbuf(1,i)
    cfsx(1,3,i) = cfsb(1,3) * xbuf(1,i)
    cfsx(2,1,i) = cfsb(2,1) * xbuf(1,i)
    cfsx(2,2,i) = cfsb(2,2) * xbuf(1,i)
    cfsx(2,3,i) = cfsb(2,3) * xbuf(1,i)
    cfsx(3,1,i) = cfsb(3,1) * xbuf(1,i)
    cfsx(3,2,i) = cfsb(3,2) * xbuf(1,i)
    cfsx(3,3,i) = cfsb(3,3) * xbuf(1,i)

    ! poloidal component Christoffel symb (normalized)
    cfsy(1,1,i) = cfsb(1,1) * xbuf(2,i)
    cfsy(1,2,i) = cfsb(1,2) * xbuf(2,i)
    cfsy(1,3,i) = cfsb(1,3) * xbuf(2,i)
    cfsy(2,1,i) = cfsb(2,1) * xbuf(2,i)
    cfsy(2,2,i) = cfsb(2,2) * xbuf(2,i)
    cfsy(2,3,i) = cfsb(2,3) * xbuf(2,i)
    cfsy(3,1,i) = cfsb(3,1) * xbuf(2,i)
    cfsy(3,2,i) = cfsb(3,2) * xbuf(2,i)
    cfsy(3,3,i) = cfsb(3,3) * xbuf(2,i)

    ! azimuthal component Christoffel symb (normalized)
    cfsz(1,1,i) = cfsb(1,1) * xbuf(3,i)
    cfsz(1,2,i) = cfsb(1,2) * xbuf(3,i)
    cfsz(1,3,i) = cfsb(1,3) * xbuf(3,i)
    cfsz(2,1,i) = cfsb(2,1) * xbuf(3,i)
    cfsz(2,2,i) = cfsb(2,2) * xbuf(3,i)
    cfsz(2,3,i) = cfsb(2,3) * xbuf(3,i)
    cfsz(3,1,i) = cfsb(3,1) * xbuf(3,i)
    cfsz(3,2,i) = cfsb(3,2) * xbuf(3,i)
    cfsz(3,3,i) = cfsb(3,3) * xbuf(3,i)

    ! gravitational acceleration
    fac3 = - this%geometry%rs / (2._p_double * r**3)
    grav(1,i) = fac3 * xbuf(1,i)
    grav(2,i) = fac3 * xbuf(2,i)
    grav(3,i) = fac3 * xbuf(3,i)

    ! gravitomagnetic tensor
    htensor(1,1,i) = 0._p_double
    htensor(1,2,i) = 0._p_double
    htensor(1,3,i) = 0._p_double
    htensor(2,1,i) = 0._p_double
    htensor(2,2,i) = 0._p_double
    htensor(2,3,i) = 0._p_double
    htensor(3,1,i) = 0._p_double
    htensor(3,2,i) = 0._p_double
    htensor(3,3,i) = 0._p_double

    ! e and b fields analytical (physical - spherical)
    !ep(1,i) = 0._p_double
    !ep(2,i) = 0._p_double
    !ep(3,i) = 0._p_double
    !bp(1,i) = 10._p_double*xbuf(3,i)/r
    !bp(2,i) = -10._p_double*sqrt(xbuf(1,i)**2+xbuf(2,i)**2)/r
    !bp(3,i) = 0._p_double
  end do


end subroutine get_geometric_quant_sch_cart
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Gets position dependent quantities for GR particle push - Kerr-Slow (cart)
!-------------------------------------------------------------------------------
subroutine get_geometric_quant_kslow_cart(this, alpha, beta, smcov, smctr, smdet, cfsx, cfsy, cfsz, grav, htensor, xbuf, np)!, ep, bp)

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(in) :: this
  real(p_k_part), dimension(rank + 1 , p_cache_size), intent(in) :: xbuf
  integer, intent(in) :: np
  real(p_double), dimension(p_cache_size), intent(inout) :: alpha, smdet
  real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: beta, grav
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: smcov, smctr
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsx
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsy
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: cfsz
  real(p_double), dimension(p_p_dim,p_p_dim,p_cache_size), intent(inout) :: htensor
  !real(p_double), dimension(p_p_dim,p_cache_size), intent(inout) :: bp, ep

  integer :: i
  real(p_double) :: r, fac1, fac2, fac3, fac4, alpha2
  real(p_double), dimension(p_p_dim,p_p_dim) :: M, cfsb


  do i=1, np
    ! temporary functions
    r = sqrt(xbuf(1,i)**2 + xbuf(2,i)**2 + xbuf(3,i)**2)
    alpha2 = 1._p_double - this%geometry%rs / r

    M(1,1) = - (xbuf(2,i)**2 + xbuf(3,i)**2)
    M(1,2) = xbuf(1,i) * xbuf(2,i)
    M(1,3) = xbuf(1,i) * xbuf(3,i)
    M(2,1) = xbuf(2,i) * xbuf(1,i)
    M(2,2) = - (xbuf(1,i)**2 + xbuf(3,i)**2)
    M(2,3) = xbuf(2,i) * xbuf(3,i)
    M(3,1) = xbuf(3,i) * xbuf(1,i)
    M(3,2) = xbuf(3,i) * xbuf(2,i)
    M(3,3) = - (xbuf(1,i)**2 + xbuf(2,i)**2)

    fac1 = - this%geometry%rs * (3._p_double*r - 2._p_double*this%geometry%rs) / (2._p_double* r**5 * (r-this%geometry%rs))

    cfsb(1,1) = fac1 * (M(1,1) + r**3 / (3._p_double*r - 2._p_double*this%geometry%rs))
    cfsb(1,2) = fac1 * M(1,2)
    cfsb(1,3) = fac1 * M(1,3)
    cfsb(2,1) = fac1 * M(2,1)
    cfsb(2,2) = fac1 * (M(2,2) + r**3 / (3._p_double*r - 2._p_double*this%geometry%rs))
    cfsb(2,3) = fac1 * M(2,3)
    cfsb(3,1) = fac1 * M(3,1)
    cfsb(3,2) = fac1 * M(3,2)
    cfsb(3,3) = fac1 * (M(3,3) + r**3 / (3._p_double*r - 2._p_double*this%geometry%rs))

    ! lapse function
    alpha(i) = sqrt(alpha2)

    ! shift vector
    beta(1,i) = this%geometry%beta0 * xbuf(2,i) / r**3
    beta(2,i) = -this%geometry%beta0 * xbuf(1,i) / r**3
    beta(3,i) = 0._p_double

    ! spatial metric tensor (covariant)
    fac2 = this%geometry%rs / (r**2 * (r - this%geometry%rs))

    smcov(1,1,i) = fac2 * (M(1,1) + r**3 / this%geometry%rs)
    smcov(1,2,i) = fac2 * M(1,2)
    smcov(1,3,i) = fac2 * M(1,3)
    smcov(2,1,i) = fac2 * M(2,1)
    smcov(2,2,i) = fac2 * (M(2,2) + r**3 / this%geometry%rs)
    smcov(2,3,i) = fac2 * M(2,3)
    smcov(3,1,i) = fac2 * M(3,1)
    smcov(3,2,i) = fac2 * M(3,2)
    smcov(3,3,i) = fac2 * (M(3,3) + r**3 / this%geometry%rs)

    ! spatial metric determinant
    smdet(i) = 1._p_double / alpha2

    ! spatial metric tensor (contravariant)
    smctr(1,1,i) = - this%geometry%rs * M(1,1) / r**3 + alpha2
    smctr(1,2,i) = - this%geometry%rs * M(1,2) / r**3
    smctr(1,3,i) = - this%geometry%rs * M(1,3) / r**3
    smctr(2,1,i) = - this%geometry%rs * M(2,1) / r**3
    smctr(2,2,i) = - this%geometry%rs * M(2,2) / r**3 + alpha2
    smctr(2,3,i) = - this%geometry%rs * M(2,3) / r**3
    smctr(3,1,i) = - this%geometry%rs * M(3,1) / r**3
    smctr(3,2,i) = - this%geometry%rs * M(3,2) / r**3
    smctr(3,3,i) = - this%geometry%rs * M(3,3) / r**3 + alpha2

    ! radial component Christoffel symb (normalized)
    cfsx(1,1,i) = cfsb(1,1) * xbuf(1,i)
    cfsx(1,2,i) = cfsb(1,2) * xbuf(1,i)
    cfsx(1,3,i) = cfsb(1,3) * xbuf(1,i)
    cfsx(2,1,i) = cfsb(2,1) * xbuf(1,i)
    cfsx(2,2,i) = cfsb(2,2) * xbuf(1,i)
    cfsx(2,3,i) = cfsb(2,3) * xbuf(1,i)
    cfsx(3,1,i) = cfsb(3,1) * xbuf(1,i)
    cfsx(3,2,i) = cfsb(3,2) * xbuf(1,i)
    cfsx(3,3,i) = cfsb(3,3) * xbuf(1,i)

    ! poloidal component Christoffel symb (normalized)
    cfsy(1,1,i) = cfsb(1,1) * xbuf(2,i)
    cfsy(1,2,i) = cfsb(1,2) * xbuf(2,i)
    cfsy(1,3,i) = cfsb(1,3) * xbuf(2,i)
    cfsy(2,1,i) = cfsb(2,1) * xbuf(2,i)
    cfsy(2,2,i) = cfsb(2,2) * xbuf(2,i)
    cfsy(2,3,i) = cfsb(2,3) * xbuf(2,i)
    cfsy(3,1,i) = cfsb(3,1) * xbuf(2,i)
    cfsy(3,2,i) = cfsb(3,2) * xbuf(2,i)
    cfsy(3,3,i) = cfsb(3,3) * xbuf(2,i)

    ! azimuthal component Christoffel symb (normalized)
    cfsz(1,1,i) = cfsb(1,1) * xbuf(3,i)
    cfsz(1,2,i) = cfsb(1,2) * xbuf(3,i)
    cfsz(1,3,i) = cfsb(1,3) * xbuf(3,i)
    cfsz(2,1,i) = cfsb(2,1) * xbuf(3,i)
    cfsz(2,2,i) = cfsb(2,2) * xbuf(3,i)
    cfsz(2,3,i) = cfsb(2,3) * xbuf(3,i)
    cfsz(3,1,i) = cfsb(3,1) * xbuf(3,i)
    cfsz(3,2,i) = cfsb(3,2) * xbuf(3,i)
    cfsz(3,3,i) = cfsb(3,3) * xbuf(3,i)

    ! gravitational acceleration
    fac3 = - this%geometry%rs / (2._p_double * r**3)
    grav(1,i) = fac3 * xbuf(1,i)
    grav(2,i) = fac3 * xbuf(2,i)
    grav(3,i) = fac3 * xbuf(3,i)

    ! gravitomagnetic tensor
    fac4 = this%geometry%beta0 / (alpha(i) * r**5)
    htensor(1,1,i) = - 3._p_double * xbuf(1,i) * xbuf(2,i) * fac4
    htensor(1,2,i) = (3._p_double * xbuf(1,i)**2 - r**2) * fac4
    htensor(1,3,i) = 0._p_double
    htensor(2,1,i) = (r**2 - 3._p_double * xbuf(2,i)**2) * fac4
    htensor(2,2,i) = 3._p_double * xbuf(1,i) * xbuf(2,i) * fac4
    htensor(2,3,i) = 0._p_double
    htensor(3,1,i) = - 3._p_double * xbuf(2,i) * xbuf(3,i) * fac4
    htensor(3,2,i) = 3._p_double * xbuf(1,i) * xbuf(3,i) * fac4
    htensor(3,3,i) = 0._p_double

    ! e and b fields analytical (physical - spherical)
    !ep(1,i) = 0._p_double
    !ep(2,i) = 0._p_double
    !ep(3,i) = 0._p_double
    !bp(1,i) = 10._p_double*xbuf(3,i)/r
    !bp(2,i) = -10._p_double*sqrt(xbuf(1,i)**2+xbuf(2,i)**2)/r
    !bp(3,i) = 0._p_double
  end do


end subroutine get_geometric_quant_kslow_cart
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine get_sph_grid_pos_gr( this, xbuf, xbufs, ptrcur, np, dt)

  use m_species_push, only : ntrim
  use m_geometry_gr

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(inout) :: this
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(in) :: xbuf
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(inout) :: xbufs
  integer, intent(in) :: ptrcur, np
  real(p_k_part), intent(in) :: dt


  integer :: i, pp
  real(p_double), dimension(rank , p_cache_size) :: xbufg
  integer, dimension(rank , p_cache_size) :: dxi
  class( t_geometry_gr ), pointer :: g

  ! executable statements
  g => this%geometry

  ! loop to get particle's grid quantities for emf interpolation
  ! gets also the updated spherical positions - xbufs

  pp = ptrcur
  do i=1,np
    xbufs(1,i) = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    xbufs(2,i) = acos( xbuf(3,i) / xbufs(1,i) )
    xbufs(3,i) = atan2(xbuf(2,i),xbuf(1,i))

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets the updated spherical grid positions
    pp = ptrcur
    do i=1,np
      this%x(1,pp) = xbufg(1,i)
      this%x(2,pp) = xbufg(2,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

      pp = pp + 1
    end do

end subroutine get_sph_grid_pos_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine get_sph_grid_pos_cartnorth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  use m_species_push, only : ntrim
  use m_geometry_gr

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(inout) :: this
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(in) :: xbuf
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(inout) :: xbufs
  integer, intent(in) :: ptrcur, np
  real(p_k_part), intent(in) :: dt


  integer :: i, pp
  real(p_double), dimension(rank , p_cache_size) :: xbufg
  integer, dimension(rank , p_cache_size) :: dxi
  real(p_double) :: r, auxr, xc, xc12, xc32, xc52, xc72, xc92
  real(p_double) :: x, xa, xa2, xa3, xa5, xa7, xa9, xa11, xa13, xa15
  real(p_double) :: xac, xac2, xac3, xac5, xac6, xac7, xac9
  real(p_double) :: xac10, xac11, xac13, xac14, xac15
  real(p_double) :: pi2, pi4, atanxa, signx, signypos, signdif

  class( t_geometry_gr ), pointer :: g

  ! executable statements
  g => this%geometry

  ! loop to get particle's grid quantities for emf interpolation
  ! gets also the updated spherical positions - xbufs

  ! used to Taylor expand trig functions
  pi2 = pi_2
  pi4 = 0.500000000000000_p_double * pi2

  pp = ptrcur
  do i=1,np
    r = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    auxr = 1.0_p_double / r

    xc = 1.0_p_double - xbuf(3,i) * auxr
    xc12 = sqrt(xc)
    xc32 = xc12 * xc
    xc52 = xc32 * xc
    xc72 = xc52 * xc
    xc92 = xc72 * xc

    xbufs(1,i) = r

    ! evaluate acos(z/r)
    xbufs(2,i) = 1.414213562373095_p_double * xc12 + 0.117851130197758_p_double * xc32 +&
                 0.026516504294496_p_double * xc52 + 0.007891816754314_p_double * xc72 +&
                 0.002685409867787_p_double * xc92

    ! evaluate atan2 of x
    ! argument that enters atan(x)
    x = xbuf(2,i) / xbuf(1,i)
    xa = abs( x )

    ! loop to select the expansion point
    if (xa .le. 0.418597947318945_p_double) then
      !here goes the expansion around x ~ 0
      xa2 = xa * xa
      xa3 = xa2 * xa
      xa5 = xa3 * xa2
      xa7 = xa5 * xa2
      xa9 = xa7 * xa2
      xa11 = xa9 * xa2
      xa13 = xa11 * xa2
      xa15 = xa13 * xa2

      ! evaluate atan(abs(x)) around x ~ 0
      atanxa = xa - 0.333333333333333_p_double * xa3 +&
                    0.200000000000000_p_double * xa5 - 0.142857142857143_p_double * xa7 +&
                    0.111111111111111_p_double * xa9 - 0.090909090909091_p_double * xa11 +&
                    0.076923076923077_p_double * xa13 - 0.066666666666667_p_double * xa15

    else if (xa .le. 1.817313073051244_p_double) then
      !here goes the expansion around x ~ 1
      xac = xa - 1.0_p_double
      xac2 = xac * xac
      xac3 = xac2 * xac
      xac5 = xac3 * xac2
      xac7 = xac5 * xac2
      xac9 = xac7 * xac2
      xac11 = xac9 * xac2
      xac13 = xac11 * xac2
      xac15 = xac13 * xac2

      xac6 = xac5 * xac
      xac10 = xac9 * xac
      xac14 = xac13 * xac

      atanxa = pi4 + 0.500000000000000_p_double * xac - 0.250000000000000_p_double * xac2 +&
                     0.083333333333333_p_double * xac3 - 0.025000000000000_p_double * xac5 +&
                     0.020833333333333_p_double * xac6 - 0.008928571428571_p_double * xac7 +&
                     0.003472222222222_p_double * xac9 - 0.003125000000000_p_double * xac10 +&
                     0.001420454545455_p_double * xac11 - 0.000600961538462_p_double * xac13 +&
                     0.000558035714286_p_double * xac14 - 0.000260416666667_p_double * xac15

    else
      !here goes the expansion around x ~ +inf
      xac = 1.0_p_double / xa
      xac2 = xac * xac
      xac3 = xac2 * xac
      xac5 = xac3 * xac2
      xac7 = xac5 * xac2
      xac9 = xac7 * xac2
      xac11 = xac9 * xac2
      xac13 = xac11 * xac2
      xac15 = xac13 * xac2

      atanxa = pi2 - xac + 0.333333333333333_p_double * xac3 -&
                     0.200000000000000_p_double * xac5 + 0.142857142857143_p_double * xac7 -&
                     0.111111111111111_p_double * xac9 + 0.090909090909091_p_double * xac11 -&
                     0.076923076923077_p_double * xac13 + 0.066666666666667_p_double * xac15

    endif

    ! required signs to evaluate atan2(x)
    signx = sign(1.0_p_double, x)
    signypos = sign(1.0_p_k_part, xbuf(2,i))
    signdif = signypos - signx

    ! evaluate atan2(ypos/xpos)
    xbufs(3,i) = signx * atanxa + pi2 * signdif

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets the updated spherical grid positions
    pp = ptrcur
    do i=1,np
      this%x(1,pp) = xbufg(1,i)
      this%x(2,pp) = xbufg(2,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

      pp = pp + 1
    end do

end subroutine get_sph_grid_pos_cartnorth_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine get_sph_grid_pos_cartsouth_gr( this, xbuf, xbufs, ptrcur, np, dt)

  use m_species_push, only : ntrim
  use m_geometry_gr

  implicit none

  integer, parameter :: rank = 2

  class( t_species_gr ), intent(inout) :: this
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(in) :: xbuf
  real(p_k_part), dimension(rank +1 , p_cache_size), intent(inout) :: xbufs
  integer, intent(in) :: ptrcur, np
  real(p_k_part), intent(in) :: dt


  integer :: i, pp
  real(p_double), dimension(rank , p_cache_size) :: xbufg
  integer, dimension(rank , p_cache_size) :: dxi
  real(p_double) :: r, auxr, xc, xc12, xc32, xc52, xc72, xc92, shift
  real(p_double) :: x, xa, xa2, xa3, xa5, xa7, xa9, xa11, xa13, xa15
  real(p_double) :: xac, xac2, xac3, xac5, xac6, xac7, xac9
  real(p_double) :: xac10, xac11, xac13, xac14, xac15
  real(p_double) :: pi2, pi4, atanxa, signx, signypos, signdif

  class( t_geometry_gr ), pointer :: g

  ! executable statements
  g => this%geometry

  ! loop to get particle's grid quantities for emf interpolation
  ! gets also the updated spherical positions - xbufs

  ! used to Taylor expand trig functions
  shift = pi
  pi2 = pi_2
  pi4 = 0.500000000000000_p_double * pi2

  pp = ptrcur
  do i=1,np
    r = sqrt( xbuf(1,i)**2 + xbuf(2,i)**2 +xbuf(3,i)**2 )
    auxr = 1.0_p_double / r

    ! evaluate acos of (x = z/r)
    xc = 1.0_p_double + xbuf(3,i) * auxr
    xc12 = - sqrt(xc)
    xc32 = xc12 * xc
    xc52 = xc32 * xc
    xc72 = xc52 * xc
    xc92 = xc72 * xc

    xbufs(1,i) = r

    ! evaluate acos(z/r)
    xbufs(2,i) = 1.414213562373095_p_double * xc12 + 0.117851130197758_p_double * xc32 +&
                 0.026516504294496_p_double * xc52 + 0.007891816754314_p_double * xc72 +&
                 0.002685409867787_p_double * xc92 + shift

    ! evaluate atan2 of x
    ! argument that enters atan(x)
    x = xbuf(2,i) / xbuf(1,i)
    xa = abs( x )

    ! loop that selects the Taylor expansion points
    if (xa .le. 0.418597947318945_p_double) then
      !here goes the expansion around x ~ 0
      xa2 = xa * xa
      xa3 = xa2 * xa
      xa5 = xa3 * xa2
      xa7 = xa5 * xa2
      xa9 = xa7 * xa2
      xa11 = xa9 * xa2
      xa13 = xa11 * xa2
      xa15 = xa13 * xa2

      ! evaluate atan(abs(x)) around x ~ 0
      atanxa = xa - 0.333333333333333_p_double * xa3 +&
                    0.200000000000000_p_double * xa5 - 0.142857142857143_p_double * xa7 +&
                    0.111111111111111_p_double * xa9 - 0.090909090909091_p_double * xa11 +&
                    0.076923076923077_p_double * xa13 - 0.066666666666667_p_double * xa15

    else if (xa .le. 1.817313073051244_p_double) then
      !here goes the expansion around x ~ 1
      xac = xa - 1.0_p_double
      xac2 = xac * xac
      xac3 = xac2 * xac
      xac5 = xac3 * xac2
      xac7 = xac5 * xac2
      xac9 = xac7 * xac2
      xac11 = xac9 * xac2
      xac13 = xac11 * xac2
      xac15 = xac13 * xac2

      xac6 = xac5 * xac
      xac10 = xac9 * xac
      xac14 = xac13 * xac

      atanxa = pi4 + 0.500000000000000_p_double * xac - 0.250000000000000_p_double * xac2 +&
                     0.083333333333333_p_double * xac3 - 0.025000000000000_p_double * xac5 +&
                     0.020833333333333_p_double * xac6 - 0.008928571428571_p_double * xac7 +&
                     0.003472222222222_p_double * xac9 - 0.003125000000000_p_double * xac10 +&
                     0.001420454545455_p_double * xac11 - 0.000600961538462_p_double * xac13 +&
                     0.000558035714286_p_double * xac14 - 0.000260416666667_p_double * xac15

    else
      !here goes the expansion around x ~ +inf
      xac = 1.0_p_double / xa
      xac2 = xac * xac
      xac3 = xac2 * xac
      xac5 = xac3 * xac2
      xac7 = xac5 * xac2
      xac9 = xac7 * xac2
      xac11 = xac9 * xac2
      xac13 = xac11 * xac2
      xac15 = xac13 * xac2

      atanxa = pi2 - xac + 0.333333333333333_p_double * xac3 -&
                     0.200000000000000_p_double * xac5 + 0.142857142857143_p_double * xac7 -&
                     0.111111111111111_p_double * xac9 + 0.090909090909091_p_double * xac11 -&
                     0.076923076923077_p_double * xac13 + 0.066666666666667_p_double * xac15

    endif

    ! required signs to evaluate atan2
    signx = sign(1.0_p_double, x)
    signypos = sign(1.0_p_k_part, xbuf(2,i))
    signdif = signypos - signx

    ! evaluate atan2(ypos/xpos)
    xbufs(3,i) = signx * atanxa + pi2 * signdif

    dxi(1,i) = ntrim( (xbufs(1,i) - g%r%f1(2, this%ix(1,pp) )) / g%dr%f1(1, this%ix(1,pp) ) )
    xbufg(1,i) = (xbufs(1,i) - g%r%f1(2, this%ix(1,pp)+dxi(1,i) )) / g%dr%f1(1, this%ix(1,pp)+dxi(1,i) )

    dxi(2,i) = ntrim( (xbufs(2,i) - g%t%f1(2, this%ix(2,pp) )) / dt )
    xbufg(2,i) = (xbufs(2,i) - g%t%f1(2, this%ix(2,pp)+dxi(2,i) )) / dt

    pp = pp + 1
  end do

  ! gets the updated spherical grid positions
    pp = ptrcur
    do i=1,np
      this%x(1,pp) = xbufg(1,i)
      this%x(2,pp) = xbufg(2,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

      pp = pp + 1
    end do

end subroutine get_sph_grid_pos_cartsouth_gr
!-------------------------------------------------------------------------------

end module m_species_push_gr
