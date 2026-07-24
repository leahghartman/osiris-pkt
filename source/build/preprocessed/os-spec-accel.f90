# 1 "spec/os-spec-accel.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/os-spec-accel.f90" 2
!-----------------------------------------------------------------------------------------
! Particle accelerator (for beam initialization)
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
# 6 "spec/os-spec-accel.f90" 2
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
# 7 "spec/os-spec-accel.f90" 2

module m_species_accelerate

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
# 11 "spec/os-spec-accel.f90" 2

use m_parameters

use m_species_define
use m_species_current
use m_species_push, only : advance_deposit_1d, advance_deposit_2d, advance_deposit_3d, &
                           dudt_boris

use m_vdf_define

use m_emf_define
use m_emf
use m_emf_interpolate

use m_math

implicit none

private

interface accelerate_deposit
  module procedure accelerate_deposit
end interface

public :: accelerate_deposit

contains

!-----------------------------------------------------------------------------------------
! Accelerate particles using a fixed acceleration, advance positions only in the x1
! direction and deposit current
!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit( this, emf, jay, t, tstep, tid, n_threads )

  use m_time_step

  implicit none

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent(in) :: emf
  type( t_vdf ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep

  integer, intent(in) :: tid ! local thread id
  integer, intent(in) :: n_threads ! total number of threads

  ! local variables

  real(p_k_part) :: u_frac, q_frac
  integer :: chunk, i0, i1
  logical :: temp_dudt


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
      write(err_buf__,*) 'Not implemented';call err__("spec/os-spec-accel.f90",89)
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
      write(err_buf__,*) 'Not implemented';call err__("spec/os-spec-accel.f90",108)
      call abort_program( p_err_notimplemented )
    end select
  else
    q_frac = 1.0
  endif

  ! Push particles. Boundary crossings will be checked at update_boundary
  if ( this%udist%n_accelerate > 0 .and. &
       this%udist%n_accelerate_type == incr_unfreeze ) then

    ! If using simd pusher normally, use boris pusher for now
    if ( .not. associated(this%dudt) ) then
      temp_dudt = .true.
      this%dudt => dudt_boris
    else
      temp_dudt = .false.
    endif

    select case ( p_x_dim )
    case (1)
      call unfreeze_deposit_1d( this, emf, jay, dt(tstep), t, i0, i1, u_frac, q_frac )

    case (2)
      call unfreeze_deposit_2d( this, emf, jay, dt(tstep), t, i0, i1, u_frac, q_frac )

    case (3)
      call unfreeze_deposit_3d( this, emf, jay, dt(tstep), t, i0, i1, u_frac, q_frac )

    end select

    if ( temp_dudt ) then
      this%dudt => null()
    endif

  else

    select case ( p_x_dim )
    case (1)
      call accelerate_deposit_1d( this, jay, dt(tstep), i0, i1, u_frac, q_frac )

    case (2)
      call accelerate_deposit_2d( this, jay, dt(tstep), i0, i1, u_frac, q_frac )

    case (3)
      call accelerate_deposit_3d( this, jay, dt(tstep), i0, i1, u_frac, q_frac )

    end select

  endif


end subroutine accelerate_deposit
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit_1d( this, jay, dt, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  class( t_species ), intent(inout) :: this

  type( t_vdf ), intent(inout) :: jay
  real( p_double ), intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac
  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi
  real(p_k_part), dimension(p_cache_size) :: rg
  real(p_k_part), dimension(p_cache_size) :: q_reduced
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel
  integer :: i, ptrcur, np, pp

  dt_dx1 = real( dt/this%dx(1), p_k_part )
  u_accel = u_frac * this%udist%ufl(1)
  vdt_dx1 = u_accel * dt_dx1 / &
                    sqrt( 1.0_p_k_part + u_accel**2 )
  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
    rg(i) = 0.0
  end do

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

        vdt_dx1 = u_accel * dt_dx1 / &
                          sqrt( 1.0_p_k_part + u_accel**2 )

        xbuf(1,i) = this%x(1,pp) + vdt_dx1

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    else
      do i=1,np

        xbuf(1,i) = this%x(1,pp) + vdt_dx1

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    endif

    q_reduced(1:np) = this%q(ptrcur:ptrcur+np-1) * q_frac

    ! Deposit current
    call this % dep_current_1d( jay, dxi, xbuf, &
                this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                q_reduced, rg, &
                this%p(:,ptrcur:), &
                np, dt )

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp) = xbuf(1,i) - dxi(1,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      pp = pp + 1
    end do

  end do

end subroutine accelerate_deposit_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit_2d( this, jay, dt, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_species ), intent(inout) :: this

  type( t_vdf ), intent(inout) :: jay
  real( p_double ), intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi
  real(p_k_part), dimension(p_cache_size) :: rg
  real(p_k_part), dimension(p_cache_size) :: q_reduced
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel
  integer :: i, ptrcur, np, pp

  dt_dx1 = real( dt/this%dx(1), p_k_part )
  u_accel = u_frac * this%udist%ufl(1)
  vdt_dx1 = u_accel * dt_dx1 / &
                    sqrt( 1.0_p_k_part + u_accel**2 )

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

    pp = ptrcur
    if(this%udist%use_particle_uacc) then
      do i=1,np
        u_accel = u_frac * this%p(1,pp)
        vdt_dx1 = u_accel * dt_dx1 / &
                          sqrt( 1.0_p_k_part + u_accel**2 )
        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    else
      do i=1,np

        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    endif

    q_reduced(1:np) = this%q(ptrcur:ptrcur+np-1) * q_frac

    ! Deposit current
    call this % dep_current_2d( jay, dxi, xbuf, &
                this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                q_reduced, rg, &
                this%p(:,ptrcur:), &
                np, dt )

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp) = xbuf(1,i) - dxi(1,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      pp = pp + 1
    end do

  end do

end subroutine accelerate_deposit_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit_3d( this, jay, dt, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  class( t_species ), intent(inout) :: this

  type( t_vdf ), intent(inout) :: jay
  real( p_double ), intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi
  real(p_k_part), dimension(p_cache_size) :: q_reduced
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel
  integer :: i, ptrcur, np, pp


  dt_dx1 = real( dt/this%dx(1), p_k_part )
  u_accel = u_frac * this%udist%ufl(1)
  vdt_dx1 = u_accel * dt_dx1 / &
                    sqrt( 1.0_p_k_part + u_accel**2 )
  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
    dxi(2,i) = 0
    dxi(3,i) = 0
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

        vdt_dx1 = u_accel * dt_dx1 / &
                          sqrt( 1.0_p_k_part + u_accel**2 )

        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)
        xbuf(3,i) = this%x(3,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    else
      do i=1,np

        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)
        xbuf(3,i) = this%x(3,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    endif


    q_reduced(1:np) = this%q(ptrcur:ptrcur+np-1) * q_frac

    ! Deposit current
    call this % dep_current_3d( jay, dxi, xbuf, &
                this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                q_reduced, np, dt )

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp) = xbuf(1,i) - dxi(1,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      pp = pp + 1
    end do

  end do

end subroutine accelerate_deposit_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine unfreeze_deposit_1d( this, emf, jay, dt, t, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent(in) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(in) :: dt, t
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi, ixbuf
  real(p_k_part), dimension(p_cache_size) :: rg
  real(p_k_part), dimension(p_cache_size) :: q_reduced
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: pbuf
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel, z_unfreeze
  integer :: i, ptrcur, np, pp, n_reg_push
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
  end do

  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! Count the number of particles that need to be pushed normally
    n_reg_push = 0
    pp = ptrcur
    do i = 1, np

      if ( this%x(1,pp) >= z_unfreeze ) then
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
      call advance_deposit_1d( this, emf, jay, energy, dt, ptrcur, ptrcur + np - 1, t )

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

            dxi(1,i) = ntrim( xbuf(1,i) )
          else
            ! Regular push (remember dxi(1,i) = 0)
            xbuf(1,i) = this%x(1,pp)
          endif

          pp = pp + 1
        end do
      else
        do i=1,np

          if ( dxi(1,i) == 1 ) then

            ! Accel push
            xbuf(1,i) = this%x(1,pp) + vdt_dx1

            dxi(1,i) = ntrim( xbuf(1,i) )

          else

            ! Regular push (remember dxi(1,i) = 0)
            xbuf(1,i) = this%x(1,pp)

          endif

          pp = pp + 1
        end do
      endif

      ! Scale charge for accel push (particles from regular push should have no charge)
      q_reduced(1:np) = q_reduced(1:np) * q_frac

      ! Deposit current
      call this % dep_current_1d( jay, dxi, xbuf, &
                  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                  q_reduced, rg, &
                  this%p(:,ptrcur:), &
                  np, dt )

      ! copy data from buffer to species data trimming positions
      pp = ptrcur
      do i = 1, np
        this%x(1,pp) = xbuf(1,i) - dxi(1,i)
        this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
        pp = pp + 1
      end do

    endif

  end do

end subroutine unfreeze_deposit_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine unfreeze_deposit_2d( this, emf, jay, dt, t, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent(in) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(in) :: dt, t
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi, ixbuf
  real(p_k_part), dimension(p_cache_size) :: rg
  real(p_k_part), dimension(p_cache_size) :: q_reduced
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: pbuf
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel, z_unfreeze
  integer :: i, ptrcur, np, pp, n_reg_push
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

    ! Count the number of particles that need to be pushed normally
    n_reg_push = 0
    pp = ptrcur
    do i = 1, np

      if ( this%x(1,pp) >= z_unfreeze ) then
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
      call advance_deposit_2d( this, emf, jay, energy, dt, ptrcur, ptrcur + np - 1, t )

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

      ! Deposit current
      call this % dep_current_2d( jay, dxi, xbuf, &
                  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                  q_reduced, rg, &
                  this%p(:,ptrcur:), &
                  np, dt )

      ! copy data from buffer to species data trimming positions
      pp = ptrcur
      do i = 1, np
        this%x(1,pp) = xbuf(1,i) - dxi(1,i)
        this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
        pp = pp + 1
      end do

    endif

  end do

end subroutine unfreeze_deposit_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine unfreeze_deposit_3d( this, emf, jay, dt, t, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent(in) :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(in) :: dt, t
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi, ixbuf
  real(p_k_part), dimension(p_cache_size) :: q_reduced
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: pbuf
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel, z_unfreeze
  integer :: i, ptrcur, np, pp, n_reg_push
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
    dxi(2,i) = 0
    dxi(3,i) = 0
  end do

  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! Count the number of particles that need to be pushed normally
    n_reg_push = 0
    pp = ptrcur
    do i = 1, np

      if ( this%x(1,pp) >= z_unfreeze ) then
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
      call advance_deposit_3d( this, emf, jay, energy, dt, ptrcur, ptrcur + np - 1, t )

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
            xbuf(3,i) = this%x(3,pp)

            dxi(1,i) = ntrim( xbuf(1,i) )
          else
            ! Regular push (remember dxi(1,i) = 0)
            xbuf(1,i) = this%x(1,pp)
            xbuf(2,i) = this%x(2,pp)
            xbuf(3,i) = this%x(3,pp)
          endif

          pp = pp + 1
        end do
      else
        do i=1,np

          if ( dxi(1,i) == 1 ) then

            ! Accel push
            xbuf(1,i) = this%x(1,pp) + vdt_dx1
            xbuf(2,i) = this%x(2,pp)
            xbuf(3,i) = this%x(3,pp)

            dxi(1,i) = ntrim( xbuf(1,i) )

          else

            ! Regular push (remember dxi(1,i) = 0)
            xbuf(1,i) = this%x(1,pp)
            xbuf(2,i) = this%x(2,pp)
            xbuf(3,i) = this%x(3,pp)

          endif

          pp = pp + 1
        end do
      endif

      ! Scale charge for accel push (particles from regular push should have no charge)
      q_reduced(1:np) = q_reduced(1:np) * q_frac

      ! Deposit current
      call this % dep_current_3d( jay, dxi, xbuf, &
                  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                  q_reduced, np, dt )

      ! copy data from buffer to species data trimming positions
      pp = ptrcur
      do i = 1, np
        this%x(1,pp) = xbuf(1,i) - dxi(1,i)
        this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
        pp = pp + 1
      end do

    endif

  end do

end subroutine unfreeze_deposit_3d
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


end module m_species_accelerate
