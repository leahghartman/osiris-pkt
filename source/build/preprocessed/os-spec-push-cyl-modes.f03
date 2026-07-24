# 1 "cyl_modes/os-spec-push-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-spec-push-cyl-modes.f03" 2
! m_species_push_cyl_modes module
!
! Handles the species advance/deposit for the quasi-3D algorithm

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
# 6 "cyl_modes/os-spec-push-cyl-modes.f03" 2
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
# 7 "cyl_modes/os-spec-push-cyl-modes.f03" 2

module m_species_push_cyl_modes

use m_parameters
use m_math
use m_species_cyl_modes_define, only: t_species_cyl_modes
use m_species_define, only: t_species, p_cell_low, p_cell_near, incr_linear, &
                                      incr_regressive, incr_unfreeze, p_ene_recalc
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_vdf_define, only: t_vdf
use m_time_step, only: t_time_step, dt, n
use m_species_current_cyl_modes,only: get_coef, getjr_cyl_m_s1, getjr_cyl_m_s2, getjr_cyl_m_s3
use m_cyl_modes, only: t_cyl_modes

implicit none

private

interface advance_deposit_cyl_modes
  module procedure advance_deposit_cyl_modes
end interface

interface accelerate_deposit_cyl_modes
  module procedure accelerate_deposit_cyl_modes
end interface

public :: advance_deposit_cyl_modes, accelerate_deposit_cyl_modes


contains
!-----------------------------------------------------------------------------------------
! Push particles and deposit electric current
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_cyl_modes( this, emf, jay_cyl_m, t, tstep, tid, n_threads )

  implicit none

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_emf_cyl_modes ), intent( in ) :: emf
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  integer, intent(in) :: tid ! local thread id
  integer, intent(in) :: n_threads ! total number of threads

  real(p_double) :: dtcycle, energy
  integer :: chunk, i0, i1

  ! call this%validate( "before advance deposit" )

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
  ! Only handle 2d with cylindrical coordinates
  call advance_deposit_2d_cyl_modes( this, emf, jay_cyl_m, energy, dtcycle, i0, i1, t )

  this%energy(tid+1) = energy

  ! call this%validate( "after advance deposit" )

end subroutine advance_deposit_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine advance_deposit_2d_cyl_modes( this, emf, jay_cyl_m, energy, gdt, i0, i1, t )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: p_cm_x = 3, p_cm_y = 4

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_emf_cyl_modes ), intent( in ) :: emf
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: i0, i1
  real(p_double), intent(in) :: t

  ! local variables
  real(p_k_part), dimension(p_x_dim) :: rdx
  integer :: i, pp, np, ptrcur

  real(p_k_part), dimension(rank+2,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi
  real(p_k_part), dimension(p_cache_size) :: rg, rgamma
  real(p_k_part), dimension(3,p_cache_size) :: p_cyl ! fourier decomp of q with theta

  integer :: gix2, shift_ix2, mode
  real(p_double) :: dr, rdr, r_shift
  real(p_double) :: r_new, r_half, x_new, y_new, x_half, y_half, sin_th, cos_th

  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )

  shift_ix2 = this%my_nx_p(p_lower, 2) - 2
  dr = this%dx(p_r_dim)
  rdr = 1.0_p_double/dr

  ! update momenta
  if ( (.not. this%free_stream) .or. this%push_type==p_gravity ) call this%dudt( emf, gdt, i0, i1, energy, t )

  ! Since r = 0 is at the center of cell ix2 = 1 the radial position of particles will be
  ! - Positions defined with regard to the center of the cell (odd interpolation)
  ! r = ( (gix2-1) + x2 ) * dr
  ! - Positions defined with regard to the corner of the cell (even interpolation)
  ! r = ( (gix2-1) + x2 - 0.5 ) * dr
  if ( this%pos_type == p_cell_near ) then
    r_shift = 0.5_p_double
  else
    r_shift = 0.0_p_double
  endif

  ! advance positions and deposit current
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! this type of loop is actually faster than
    ! using a forall construct
    pp = ptrcur

    do i=1,np
      rg(i) = 1.0_p_k_part / &
        sqrt( ( ( 1.0_p_k_part + this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

      rgamma(i) = gdt * rg(i)
      pp = pp + 1
    enddo

    ! advance particle position
    pp = ptrcur

    do i=1,np

      ! typical cartesian push
      xbuf(1,i) = this%x(1,pp) + ( this%p(1,pp) * rgamma(i) ) * rdx(1)

      ! Convert radial "cell" position to "box" position in double precision
      gix2 = this%ix(2,pp) + shift_ix2

      ! We no longer need r_old
      ! r_old = ( this%x(2,pp) + (gix2 - r_shift) ) * dr
      ! This could alternatively be calculated in the following manner:
      ! r_old = sqrt( x_old**2 + y_old**2 )
      ! gix2 = nint( r_old*rdx(2) + r_shift )

      ! Cartesian push in transverse plane
      ! Note: these coordinates are of "box" type
      x_new = this%x( p_cm_x, pp ) + this%p(2,pp) * rgamma(i)
      y_new = this%x( p_cm_y, pp ) + this%p(3,pp) * rgamma(i)

      ! Store cartesian coordinates (used for cyl_modes current deposition)
      xbuf( p_cm_x, i ) = x_new
      xbuf( p_cm_y, i ) = y_new
      r_new = sqrt( x_new**2 + y_new**2 )

      ! Convert new position to "cell" type position

      ! this is a protection against roundoff for cold plasmas
      if ( (this%x( p_cm_x, pp ) == x_new) .and. (this%x( p_cm_y, pp ) == y_new) ) then
        xbuf(2,i) = this%x(2,pp)
      else
        xbuf(2,i) = real( r_new * rdr - ( gix2 - r_shift ) , p_k_part )
      endif

      ! Get number of cells moved in each (z,r) direction
      dxi(1,i) = ntrim( xbuf(1,i) ) ! returns +1 , 0, or -1
      dxi(2,i) = ntrim( xbuf(2,i) )

      ! Half cartesian push in transverse plane (to be centered with momentum)
      x_half = this%x( p_cm_x, pp ) + 0.5_p_k_part * this%p(2,pp) * rgamma(i)
      y_half = this%x( p_cm_y, pp ) + 0.5_p_k_part * this%p(3,pp) * rgamma(i)
      r_half = sqrt( x_half**2 + y_half**2 )

      ! convert momentum into cylindrical coordinates before depositing current
      cos_th = x_half / r_half
      sin_th = y_half / r_half

      p_cyl(1,i) = this%p(1,pp)
      p_cyl(2,i) = this%p(2,pp)*cos_th + this%p(3,pp)*sin_th
      p_cyl(3,i) = -this%p(2,pp)*sin_th + this%p(3,pp)*cos_th

      pp = pp + 1
    enddo

    ! deposit current
    ! handle mode 0 first - for this part the original 2D deposit should work fine
    call this % dep_current_2d( jay_cyl_m%pf_re(0), dxi, xbuf(1:2,:), &
                                this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
                                this%q(ptrcur:), rg, p_cyl, &
                                np, gdt )

    ! handle high order modes - should use special cylindrical mode deposit to remain charge conserving
    do mode = 1, this%n_cyl_modes

      ! still need to port cylindrical mode current deposition functions
      select case (this%interpolation)
      case( p_linear )

        call getjr_cyl_m_s1( jay_cyl_m%pf_re(mode), jay_cyl_m%pf_im(mode), dxi, &
                              xbuf(:,:), this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                              this%q(ptrcur:), rg, p_cyl, &
                              np, gdt, shift_ix2, mode )

      case( p_quadratic )

        call getjr_cyl_m_s2( jay_cyl_m%pf_re(mode), jay_cyl_m%pf_im(mode), dxi, &
                              xbuf(:,:), this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                              this%q(ptrcur:), rg, p_cyl, &
                              np, gdt, shift_ix2, mode )

      case( p_cubic )

        call getjr_cyl_m_s3( jay_cyl_m%pf_re(mode), jay_cyl_m%pf_im(mode), dxi, &
                              xbuf(:,:), this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                              this%q(ptrcur:), rg, p_cyl, &
                              np, gdt, shift_ix2, mode )

      case default
        write(err_buf__,*) 'Not implemented yet';call err__("cyl_modes/os-spec-push-cyl-modes.f03",245)
        call abort_program( p_err_notimplemented )
      end select

    enddo ! mode

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np

      this%x(1,pp) = xbuf(1,i) - dxi(1,i)
      this%x(2,pp) = xbuf(2,i) - dxi(2,i)
      this%x(3,pp) = xbuf(3,i)
      this%x(4,pp) = xbuf(4,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

      pp = pp + 1
    enddo

  enddo ! ptrcur

end subroutine advance_deposit_2d_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit_cyl_modes( this, emf, jay_cyl_m, t, tstep, tid, n_threads )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_emf_cyl_modes ), intent(in) :: emf
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  integer, intent(in) :: tid ! local thread id
  integer, intent(in) :: n_threads ! total number of threads

  real(p_k_part) :: u_frac, q_frac
  integer :: chunk, i0, i1

  ! Time centered energy diagnostic is not yet available, force recalculation later
  this%energy = p_ene_recalc

  ! range of particles for each thread
  chunk = ( this%num_par + n_threads - 1 ) / n_threads
  i0 = tid * chunk + 1
  i1 = min( (tid+1) * chunk, this%num_par )

  ! ufl momentum acceleration in p1
  if ( this%udist%n_accelerate > 0 ) then
    ! calculate charge fraction range=(0., 1.]
    u_frac = 1.0
    select case ( this%udist%n_accelerate_type )
    case ( incr_linear )
      u_frac = real( grad_linear( n(tstep)+1, this%udist%n_accelerate ), p_k_part )

    case ( incr_regressive )
      u_frac = real( grad_regressive( n(tstep)+1, this%udist%n_accelerate ), p_k_part )

    case ( incr_unfreeze )
      ! No need to do anything, keep u_frac = 1.0

    case default
      write(err_buf__,*) 'Not implemented';call err__("cyl_modes/os-spec-push-cyl-modes.f03",312)
      call abort_program( p_err_notimplemented )
    end select

  else
    u_frac = 1.0
  endif

  ! increase charge over n steps
  if ( this%udist%n_q_incr > 0 ) then
    ! calculate charge fraction range=(0., 1.]
    select case ( this%udist%n_q_incr_type )
    case ( incr_linear )
      q_frac = real( grad_linear( n(tstep)+1, this%udist%n_q_incr ), p_k_part )

    case ( incr_regressive )
      q_frac = real( grad_regressive( n(tstep)+1, this%udist%n_q_incr ), p_k_part )

    case default
      write(err_buf__,*) 'Not implemented';call err__("cyl_modes/os-spec-push-cyl-modes.f03",331)
      call abort_program( p_err_notimplemented )
    end select
  else
    q_frac = 1.0
  endif

  if ( this%udist%n_accelerate > 0 .and. &
       this%udist%n_accelerate_type == incr_unfreeze ) then

    call unfreeze_deposit_2d_cyl_modes( this, emf, jay_cyl_m, dt(tstep), t, i0, i1, &
                                        u_frac, q_frac )

  else

    call accelerate_deposit_2d_cyl_modes( this, jay_cyl_m, dt(tstep), i0, i1, u_frac, &
                                          q_frac )

  endif




end subroutine accelerate_deposit_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit_2d_cyl_modes( this, jay_cyl_m, dt, i0, i1, u_frac, q_frac )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_species_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m
  real(p_double), intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac

  real(p_k_part) :: vdt_dx1, u_accel, dt_dx1
  integer :: i, mode
  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi
  real(p_k_part), dimension(p_cache_size) :: rg
  real(p_k_part), dimension(p_cache_size) :: q_reduced, q_cos, q_sin, coeff_re, coeff_im
  integer :: ptrcur, np, pp

  dt_dx1 = real( dt/this%dx(1), p_k_part )
  u_accel = u_frac * this%udist%ufl(1)
  vdt_dx1 = u_accel * dt_dx1 / sqrt( 1.0_p_k_part + u_accel**2 )

  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
    rg(i) = 0.0
    dxi(2,i) = 0
  enddo

  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    pp = ptrcur
    if(this%udist%use_particle_uacc) then
      do i=1,np
        u_accel = u_frac * this%p(1,pp)
        vdt_dx1 = u_accel * dt_dx1 / sqrt( 1.0_p_k_part + u_accel**2 )
        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      enddo
    else
      do i=1,np

        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      enddo
    endif

    q_reduced(1:np) = this%q(ptrcur:ptrcur+np-1) * q_frac

    ! No need to convert momentum to cylindrical coordinates since only movement in x1
    ! Use the original 2D deposit for the 0th mode
    call this % dep_current_2d( jay_cyl_m%pf_re(0), dxi, xbuf, &
                                this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
                                q_reduced, rg, &
                                this%p(:,ptrcur:), &
                                np, dt )

    ! In the accelerate_deposit stage particles are only accelerated in the x1 coordinate,
    ! which means there will not be any J_r or J_theta deposited. The original 2D
    ! deposition scheme is sufficient to deposit J_z for higher modes.
    do mode = 1, this%n_cyl_modes

      call get_coef( coeff_re, coeff_im, this%x(3:4,ptrcur:), np, mode )

      q_cos(1:np) = q_reduced(1:np) * coeff_re(1:np)
      q_sin(1:np) = q_reduced(1:np) * coeff_im(1:np)

      ! Deposit current
      call this % dep_current_2d( jay_cyl_m%pf_re(mode), dxi, xbuf, &
                                  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                                  q_cos, rg, this%p(:,ptrcur:), &
                                  np, dt )

      call this % dep_current_2d( jay_cyl_m%pf_im(mode), dxi, xbuf, &
                                  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                                  q_sin, rg, this%p(:,ptrcur:), &
                                  np, dt )

    enddo ! mode

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp) = xbuf(1,i) - dxi(1,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      pp = pp + 1
    enddo

  enddo

end subroutine accelerate_deposit_2d_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine unfreeze_deposit_2d_cyl_modes( this, emf, jay_cyl_m, dt, t, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_emf_cyl_modes ), intent(in) :: emf
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m
  real(p_double), intent(in) :: dt, t
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac

  real(p_k_part), dimension(rank+2,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi, ixbuf
  real(p_k_part), dimension(p_cache_size) :: rg
  real(p_k_part), dimension(p_cache_size) :: q_reduced, q_cos, q_sin, coeff_re, coeff_im
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: pbuf
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel, z_unfreeze
  integer :: i, ptrcur, np, pp, mode, n_reg_push
  real(p_double) :: energy

  dt_dx1 = real( dt/this%dx(1), p_k_part )
  u_accel = u_frac * this%udist%ufl(1)
  vdt_dx1 = u_accel * dt_dx1 / sqrt( 1.0_p_k_part + u_accel**2 )
  ! To the right (left) of this distance, perform the regular (accel) push
  z_unfreeze = real( this%udist%unfreeze_z0 + this%udist%unfreeze_vel * &
                     ( t - this%udist%tmin ), p_k_part )
  energy = 0.0_p_double ! Not used after calculation

  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
    rg(i) = 0.0
    dxi(2,i) = 0
  end do

  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! First get the global z position of the particles
    call this % get_position( 1, ptrcur, ptrcur + np - 1, xbuf(1,:) )

    ! Count the number of particles that need to be pushed normally
    n_reg_push = 0
    pp = ptrcur
    do i = 1, np

      if ( xbuf(1,i) >= z_unfreeze ) then
        n_reg_push = n_reg_push + 1
        dxi(1,i) = 0 ! Signal that this was a regular push
        q_reduced(i) = 0.0_p_k_part ! So no current is deposited in accel push
      else
        dxi(1,i) = 1 ! Signal that this is an accel push
        ixbuf(:,i) = this%ix(:,pp)
        xbuf(:,i) = this%x(:,pp)
        pbuf(:,i) = this%p(:,pp)
        q_reduced(i) = this%q(pp)
      endif

      pp = pp + 1
    enddo

    if ( n_reg_push > 0 ) then

      ! Set charge to zero for particles requiring accel push so no current is deposited
      if ( n_reg_push < np ) then
        pp = ptrcur
        do i = 1, np
          if ( dxi(1,i) == 1 ) then
            this%q(pp) = 0.0_p_k_part
          endif
          pp = pp + 1
        enddo
      endif

      ! Perform regular push
      call advance_deposit_2d_cyl_modes( this, emf, jay_cyl_m, energy, dt, ptrcur, &
                               ptrcur + np - 1, t )

      ! Reset requisite particle data for accel push
      if ( n_reg_push < np ) then
        pp = ptrcur
        do i = 1, np
          if ( dxi(1,i) == 1 ) then
            this%ix(:,pp) = ixbuf(:,i)
            this%x(:,pp) = xbuf(:,i)
            this%p(:,pp) = pbuf(:,i)
            this%q(pp) = q_reduced(i)
          endif
          pp = pp + 1
        enddo
      endif

    endif

    ! Perform accel push if necessary
    if ( n_reg_push < np ) then

      pp = ptrcur
      if(this%udist%use_particle_uacc) then
        do i=1,np

          if ( dxi(1,i) == 1 ) then
            ! Accel push
            u_accel = u_frac * this%p(1,pp)
            vdt_dx1 = u_accel * dt_dx1 / sqrt( 1.0_p_k_part + u_accel**2 )
            xbuf(1,i) = this%x(1,pp) + vdt_dx1
            xbuf(2,i) = this%x(2,pp)

            dxi(1,i) = ntrim( xbuf(1,i) )
          else
            ! Regular push (remember dxi(1,i) = 0)
            xbuf(1,i) = this%x(1,pp)
            xbuf(2,i) = this%x(2,pp)
          endif

          pp = pp + 1
        end do
      else
        do i=1,np

          if ( dxi(1,i) == 1 ) then

            ! Accel push
            xbuf(1,i) = this%x(1,pp) + vdt_dx1
            xbuf(2,i) = this%x(2,pp)

            dxi(1,i) = ntrim( xbuf(1,i) )

          else

            ! Regular push (remember dxi(1,i) = 0)
            xbuf(1,i) = this%x(1,pp)
            xbuf(2,i) = this%x(2,pp)

          endif

          pp = pp + 1
        end do
      endif

      ! Scale charge for accel push (particles from regular push should have no charge)
      q_reduced(1:np) = q_reduced(1:np) * q_frac

      ! No need to convert momentum to cylindrical coordinates since only movement in x1
      ! Use the original 2D deposit for the 0th mode
      call this % dep_current_2d( jay_cyl_m%pf_re(0), dxi, xbuf, &
                                  this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
                                  q_reduced, rg, &
                                  this%p(:,ptrcur:), &
                                  np, dt )

      ! In the accelerate_deposit stage particles are only accelerated in the x1 coordinate,
      ! which means there will not be any J_r or J_theta deposited. The original 2D
      ! deposition scheme is sufficient to deposit J_z for higher modes.
      do mode = 1, this%n_cyl_modes

        call get_coef( coeff_re, coeff_im, this%x(3:4,ptrcur:), np, mode )

        q_cos(1:np) = q_reduced(1:np) * coeff_re(1:np)
        q_sin(1:np) = q_reduced(1:np) * coeff_im(1:np)

        ! Deposit current
        call this % dep_current_2d( jay_cyl_m%pf_re(mode), dxi, xbuf, &
                                    this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                                    q_cos, rg, this%p(:,ptrcur:), &
                                    np, dt )

        call this % dep_current_2d( jay_cyl_m%pf_im(mode), dxi, xbuf, &
                                    this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                                    q_sin, rg, this%p(:,ptrcur:), &
                                    np, dt )

      enddo ! mode

      ! copy data from buffer to species data trimming positions
      pp = ptrcur
      do i = 1, np
        this%x(1,pp) = xbuf(1,i) - dxi(1,i)
        this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
        pp = pp + 1
      enddo

    endif

  end do

end subroutine unfreeze_deposit_2d_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function ntrim(x)
!-----------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
! This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
! routine that would break for x = nearest( 0.5, -1.0 )
!-----------------------------------------------------------------------------------------
  implicit none

  real(p_k_part), intent(in) :: x
  integer :: ntrim, a, b

  if ( x < -.5 ) then
    a = -1
  else
    a = 0
  endif

  if ( x >= .5 ) then
    b = +1
  else
    b = 0
  endif

  ntrim = a+b

end function ntrim
!-----------------------------------------------------------------------------------------

end module m_species_push_cyl_modes

!-----------------------------------------------------------------------------------------
! Particle pusher for cyl_modes simulations. Note that when using these simulations the
! push_type parameter is ignored
!-----------------------------------------------------------------------------------------
subroutine push_species_cyl_modes( this, emf, jay_cyl_m, t, tstep, tid, n_threads, options )

  use m_parameters
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_emf_cyl_modes, only : t_emf_cyl_modes
  use m_time_step, only : t_time_step, n
  use m_species_udist, only : use_accelerate, use_q_incr
  use m_species_push_cyl_modes, only : advance_deposit_cyl_modes, accelerate_deposit_cyl_modes
  use m_cyl_modes, only : t_cyl_modes

  implicit none

  class( t_species_cyl_modes ), intent(inout) :: this
  class( t_emf_cyl_modes ), intent( inout ) :: emf
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  type( t_options ), intent(in) :: options
  integer, intent(in) :: tid ! local thread id
  integer, intent(in) :: n_threads ! total number of threads

  integer :: push_type

  if ( t >= this%push_start_time ) then

    push_type = this%push_type
    if ( use_accelerate( this%udist, n(tstep) ) ) push_type = p_beam_accel
    if ( use_q_incr( this%udist, n(tstep) ) ) push_type = p_beam_accel

    select case ( push_type )

      case ( p_beam_accel )
        call accelerate_deposit_cyl_modes( this, emf, jay_cyl_m, t, tstep, tid, n_threads)

      case default
        call advance_deposit_cyl_modes( this, emf, jay_cyl_m, t, tstep, tid, n_threads )

    end select


  endif

end subroutine push_species_cyl_modes
!-----------------------------------------------------------------------------------------
