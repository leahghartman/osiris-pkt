! m_species_shear module
!
! Handles particle advance / deposit for the Shear algorithm
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_shear

use m_emf_shear
use m_emf_define
use m_emf

use m_space
use m_grid_define
use m_node_conf
use m_restart
use m_time_step
use m_input_file

use m_species_define
use m_emf_define
use m_vdf_define

use m_system
use m_parameters

use m_species_current

implicit none

private

! parameter for timing event

type, extends( t_species ) :: t_species_shear

  ! specific shear parameter
  real( p_double )    :: alpha

contains

  procedure :: push => push_species_shear
  procedure :: list_algorithm => list_algorithm_spec_shear

end type t_species_shear


public :: t_species_shear

! ----------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_spec_shear( this )

  implicit none
  class( t_species_shear ), intent(in) :: this

  print *, ' '
  print *, trim(this%name),' :'
  print *, '- Shearing co-rotating frame pusher'

end subroutine list_algorithm_spec_shear
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Particle pusher for SHEAR simulations. Note that when using these simulations the
! push_type parameter is ignored
!-----------------------------------------------------------------------------------------
subroutine push_species_shear ( this, emf, current, t, tstep, tid, n_threads, options )

  use m_time_step
  use m_species_push
  use m_species_vpush
  use m_species_accelerate

  use m_logprof

  implicit none

  class( t_species_shear ), intent(inout) :: this

  ! since we are overloading the push method, we must have exactly the same interface i.e.
  ! we must declare t_emf here
  class( t_emf ), intent( inout )  ::  emf

  type( t_vdf ), intent(inout) :: current

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  type( t_options ), intent(in) :: options

  integer, intent(in) :: tid        ! local thread id
  integer, intent(in) :: n_threads  ! total number of threads

  select type( emf )

  class is ( t_emf_shear )

    ! advance_deposit_emf_shear only works with t_emf_shear objects, and emf was declared as
    ! t_emf we must use a select type construct to check that the object is indeed a
    ! t_emf_shear instance

    if ( t >= this%push_start_time ) then
        call advance_deposit_emf_shear( this, emf, current, t, tstep, tid, n_threads )
    endif

  class default

    ! This must never happen, push_species_shear must always be called with a
    ! t_emf_shear object
    call abort_program( p_err_invalid )

  end select
  ! TODO: REMOVE: this%n_current_spec = this%n_current_spec + 1
end subroutine push_species_shear
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Push particles and deposit electric current
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_emf_shear( this, emf_shear, jay, t, tstep, tid, n_threads )

  use m_time_step

  implicit none

  class( t_species_shear ), intent(inout) :: this
  class( t_emf_shear ), intent( inout )  ::  emf_shear
  type( t_vdf ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep

  integer, intent(in) :: tid        ! local thread id
  integer, intent(in) :: n_threads  ! total number of threads

  ! local variables

  real(p_double) :: dtcycle, energy
  integer :: chunk, i0, i1

  ! executable statements

  ! call validate( this, "before advance deposit" )

  ! if before push start time return silently
  if ( t < this%push_start_time ) return

  ! initialize time centered energy diagnostic
  energy = 0.0_p_double

  dtcycle = real( dt(tstep), p_k_part )

  ! range of particles for each thread
  chunk = ( this%num_par + n_threads - 1 ) / n_threads
  i0    = tid * chunk + 1
  i1    = min( (tid+1) * chunk, this%num_par )

  ! Push particles. Boundary crossings will be checked at update_boundary

  select case ( this%coordinates )

    case default
      select case ( p_x_dim )

      case( 1 )
        call advance_deposit_1d_shear( this, emf_shear, jay, energy, dtcycle, i0, i1 )
        
      case( 2 )
        call advance_deposit_2d_shear( this, emf_shear, jay, energy, dtcycle, i0, i1 )
        
      case( 3 )
        call advance_deposit_3d_shear( this, emf_shear, jay, energy, dtcycle, i0, i1 )

      case default
        ERROR('Not implemented for x_dim = ',p_x_dim)
        call abort_program(p_err_invalid)

      end select
      
  end select

end subroutine advance_deposit_emf_shear
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine advance_deposit_1d_shear( this, emf, jay, energy, gdt, i0, i1 )
!-----------------------------------------------------------------------------------------
!  In 1D the orientation of the shear is along x3 in order to have evolution of 
!  instabilities along the x1 direction
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1
  
  class( t_species_shear ), intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  type( t_vdf ),     intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_double),                     intent(in) :: gdt
  integer, intent(in) :: i0, i1

  ! local variables
  integer :: np, ptrcur, pp, i
  real(p_k_part), dimension(rank) :: rdx 

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg, rgamma
  real(p_k_part) :: dt
  
  real(p_k_part) :: e_coef, b_coef

  real(p_k_part) :: tem

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp, ptemp
  real(p_k_part), dimension(p_cache_size) :: gam_tem, otsq

  real( p_double ) :: loc_ene, gamma, u2

  ! executable statements
  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  dt = real( gdt, p_k_part )

  ! update pistons
#if 0
  havepiston = .false.
  do pst_id=1, this%num_pistons
    call update_piston( this%pistons(pst_id), t, dt, l_space )
    havepiston = havepiston .or. this%pistons(pst_id)%inbox
  enddo
#else
  if ( this%num_pistons > 0 ) then
    ERROR('Pistons are currently offline')
    call abort_program( p_err_invalid )
  endif
#endif

!--------------------------------advance momenta------------------------------------------

  ! if this is a free streaming species then return 
  if ( this%free_stream ) return
  
  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )
  e_coef = real(0.75_p_double * this%alpha * dt, p_k_part)
  b_coef = real(2.0_p_double * this%alpha * this%rqm, p_k_part)
  
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

    call this % get_emf( emf, bp, ep, np, ptrcur )
	  
	  do i=1, np  
		ep(1,i) = ep(1,i) * tem
		ep(2,i) = ep(2,i) * tem
		ep(3,i) = ep(3,i) * tem
	  end do

	  ! Perform first half of electric field acceleration.
	  ! and get time centered gamma
	  pp = ptrcur
	  
	  do i=1,np
	  	ptemp(2,i) = this%p(2,pp)
		utemp(1,i) = this%p(1,pp) + ep(1,i)
		utemp(2,i) = this%p(2,pp) + ep(2,i)
		utemp(3,i) = this%p(3,pp) + ep(3,i)

		! Get time centered gamma
        u2 = utemp(1,i)**2 + utemp(2,i)**2 + utemp(3,i)**2
        gamma = sqrt( u2 + 1 )
        
        ! accumulate time centered energy
        loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1)
  
        gam_tem(i)= tem / gamma

		pp = pp + 1
	  enddo

	  do i=1,np
		bp(1,i) = (bp(1,i)+b_coef*gamma)*gam_tem(i)
		bp(2,i) = bp(2,i)*gam_tem(i)
		bp(3,i) = bp(3,i)*gam_tem(i)
	  end do
	  
	  ! accumulate global energy
      energy = energy + loc_ene

      do i=1,np
        bp(1,i) = bp(1,i)*gam_tem(i)
        bp(2,i) = bp(2,i)*gam_tem(i)
        bp(3,i) = bp(3,i)*gam_tem(i)
      end do

	  ! Perform first half of the rotation and store in u.
	  pp = ptrcur
	  do i=1,np
		 this%p(1,pp)=utemp(1,i)+utemp(2,i)*bp(3,i)-utemp(3,i)*bp(2,i)
		 this%p(2,pp)=utemp(2,i)+utemp(3,i)*bp(1,i)-utemp(1,i)*bp(3,i)
		 this%p(3,pp)=utemp(3,i)+utemp(1,i)*bp(2,i)-utemp(2,i)*bp(1,i)
		 pp = pp + 1
	  end do
		
	  do i=1,np
		 otsq(i) = 2 / (1.0_p_k_part+bp(1,i)**2+bp(2,i)**2+bp(3,i)**2)
	  end do

	  do i=1,np
		 bp(1,i)= bp(1,i) * otsq(i)
		 bp(2,i)= bp(2,i) * otsq(i)
		 bp(3,i)= bp(3,i) * otsq(i)
	  end do

	  ! Perform second half of the rotation.
	  pp = ptrcur
	  do i=1,np
		 utemp(1,i) = utemp(1,i)+this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)
		 utemp(2,i) = utemp(2,i)+this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)
		 utemp(3,i) = utemp(3,i)+this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)
		 pp = pp + 1
	  end do

	  ! Perform second half of electric field acceleration.
	  pp = ptrcur
	  do i=1,np
		 this%p(1,pp) = utemp(1,i) + ep(1,i)
		 this%p(2,pp) = utemp(2,i) + ep(2,i)
		 this%p(3,pp) = utemp(3,i) + ep(3,i) + e_coef*(this%p(2,pp)+ptemp(2,i))
		 pp = pp + 1
	  end do
		  
  enddo
  
!--------------------------------advance positions----------------------------------------
  
  ! loop through all particles
  do ptrcur = i0, i1, p_cache_size
       
     ! check if last copy of table and set np
     if( ptrcur + p_cache_size > i1 ) then
         np = i1 - ptrcur + 1
     else
         np = p_cache_size
     endif

     pp = ptrcur
     do i=1,np
        rg(i) = 1.0_p_k_part / &
          sqrt( 1.0_p_k_part + &
          this%p(1,pp)**2 + &
          this%p(2,pp)**2 + &
          this%p(3,pp)**2 )
        rgamma(i) = dt * rg(i)
        pp = pp + 1
     end do
 
     pp = ptrcur
     do i=1,np
       xbuf (1,i) = this%x (1,pp) + this%p(1,pp) * rgamma(i) * rdx(1)
       dxi(1,i) = ntrim( xbuf(1,i) )

       pp = pp + 1
     end do

     ! Deposit current
     call this % dep_current_1d(  jay, dxi, xbuf, &
                          this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                          this%q(ptrcur:), rg,     &
                          this%p(:,ptrcur:),      &
                          np, gdt )

     ! copy data from buffer to species data
     pp = ptrcur
     do i = 1, np
       this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
       this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
       
       pp = pp + 1
     end do
     
  enddo

end subroutine advance_deposit_1d_shear
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine advance_deposit_2d_shear( this, emf, jay, energy, gdt, i0, i1 )
!-----------------------------------------------------------------------------------------
!  In 2D the orientation of the shear is along x3 in order to have evolution of 
!  the instability in the poloidal plane
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  ! dummy variables

  class( t_species_shear ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  type( t_vdf ),        intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_double),   intent(in) :: gdt
  integer, intent(in) :: i0, i1

  real(p_k_part), dimension(p_x_dim) :: rdx 
  integer :: i, pp, np, ptrcur

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg, rgamma
  real(p_k_part) :: dt
  
  real(p_k_part) :: tem

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp, ptemp
  real(p_k_part), dimension(p_cache_size) ::  gam_tem, otsq

  real( p_double ) :: loc_ene, gamma, u2
  
  real(p_k_part) :: e_coef, b_coef

      
  ! executable statements
  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
  dt = real( gdt, p_k_part )   
  
#if 0
  ! calculate piston position
  havepiston = .false.
  do pst_id=1, this%num_pistons
    call update_piston( this%pistons(pst_id), t, dt, l_space )
    havepiston = havepiston .or. this%pistons(pst_id)%inbox
  enddo
#else
  if ( this%num_pistons > 0 ) then
    ERROR('Pistons are currently offline')
    call abort_program( p_err_invalid )
  endif
#endif

!-------------------------------advance momenta-------------------------------------------
    ! executable statements

  ! if this is a free streaming species then return 
  if ( this%free_stream ) return
  
  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )
  e_coef = real(0.75_p_double * this%alpha * dt, p_k_part)
  b_coef = real(2.0_p_double * this%alpha * this%rqm, p_k_part)
  
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

     call this % get_emf( emf, bp, ep, np, ptrcur )

     do i=1, np
         ep(1,i) = ep(1,i) * tem
        ep(2,i) = ep(2,i) * tem
        ep(3,i) = ep(3,i) * tem
     end do

     ! Perform first half of electric field acceleration.
     ! and get time centered gamma

     loc_ene = 0

     pp = ptrcur
     do i=1,np
        ptemp(2,i) = this%p(2,pp)
        utemp(1,i) = this%p(1,pp) + ep(1,i)
        utemp(2,i) = this%p(2,pp) + ep(2,i)
        utemp(3,i) = this%p(3,pp) + ep(3,i)

        ! Get time centered gamma
        u2 = utemp(1,i)**2 + utemp(2,i)**2 + utemp(3,i)**2
        gamma = sqrt( u2 + 1 )

        ! accumulate time centered energy
        loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1)
  
        gam_tem(i)= tem / gamma

        pp = pp + 1
     enddo

     ! accumulate global energy
     energy = energy + loc_ene

     do i=1,np
       bp(1,i) = bp(1,i)*gam_tem(i)
       bp(2,i) = (bp(2,i)+b_coef*gamma)*gam_tem(i)
       bp(3,i) = bp(3,i)*gam_tem(i)
     end do

     ! Perform first half of the rotation and store in u.
     pp = ptrcur
     do i=1,np
        this%p(1,pp)=utemp(1,i)+utemp(2,i)*bp(3,i)-utemp(3,i)*bp(2,i)
        this%p(2,pp)=utemp(2,i)+utemp(3,i)*bp(1,i)-utemp(1,i)*bp(3,i)
        this%p(3,pp)=utemp(3,i)+utemp(1,i)*bp(2,i)-utemp(2,i)*bp(1,i)
        pp = pp + 1
     end do
  
     do i=1,np
        otsq(i) = 2. / (1.0_p_k_part+bp(1,i)**2+bp(2,i)**2+bp(3,i)**2)
     end do

     do i=1,np
        bp(1,i)= bp(1,i) * otsq(i)
        bp(2,i)= bp(2,i) * otsq(i)
        bp(3,i)= bp(3,i) * otsq(i)
     end do

     ! Perform second half of the rotation.
     pp = ptrcur
     do i=1,np
        utemp(1,i) = utemp(1,i)+this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)
        utemp(2,i) = utemp(2,i)+this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)
        utemp(3,i) = utemp(3,i)+this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)
        pp = pp + 1
     end do

     ! Perform second half of electric field acceleration.
     pp = ptrcur
     do i=1,np
        this%p(1,pp) = utemp(1,i) + ep(1,i)
        this%p(2,pp) = utemp(2,i) + ep(2,i)
        this%p(3,pp) = utemp(3,i) + ep(3,i)-e_coef*(this%p(1,pp)+ptemp(1,i))
        pp = pp + 1
     end do
               
  enddo
  
!--------advance position of particles i0 to i1 in chunks of p_cache_size-----------------
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
          sqrt( 1.0_p_k_part + &
            this%p(1,pp)**2 + &
            this%p(2,pp)**2 + &
            this%p(3,pp)**2 )                
        
        rgamma(i) = dt * rg(i)
        pp = pp + 1
     end do

     ! reference implementation
     pp = ptrcur
     do i=1,np
       xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rgamma(i) * rdx(1)
       xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rgamma(i) * rdx(2)

       dxi(1,i) = ntrim( xbuf(1,i) )
       dxi(2,i) = ntrim( xbuf(2,i) )
       
       pp = pp + 1
     end do
     
     ! Deposit current
     call this % dep_current_2d(  jay, dxi, xbuf, &
                         this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                         this%q(ptrcur:), rg,     &
                         this%p(:,ptrcur:),      &
                         np, gdt )

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
      this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
      
      pp = pp + 1
    end do

  enddo

end subroutine advance_deposit_2d_shear
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine advance_deposit_3d_shear( this, emf, jay, energy, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  ! dummy variables
  class( t_species_shear ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  type( t_vdf), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_double),                     intent(in) :: gdt
  integer, intent(in) :: i0, i1

  ! local variables
  real(p_k_part), dimension(p_max_dim) :: dt_dx 
  integer :: np, ptrcur, pp, i

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rgamma
  real(p_k_part) :: dt
  
  real(p_k_part) :: tem

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp, ptemp
  real(p_k_part), dimension(p_cache_size) ::  gam_tem, otsq

  real( p_double ) :: loc_ene, gamma, u2
  
  real(p_k_part) :: e_coef, b_coef

  ! executable statements


  dt_dx(1) = real( gdt/this%dx(1), p_k_part )
  dt_dx(2) = real( gdt/this%dx(2), p_k_part )
  dt_dx(3) = real( gdt/this%dx(3), p_k_part )
  dt = real( gdt, p_k_part ) 

#if 0
  havepiston = .false.
  do pst_id=1, this%num_pistons
    call update_piston( this%pistons(pst_id), t, dt, l_space )
    havepiston = havepiston .or. this%pistons(pst_id)%inbox
  enddo
#else
  if ( this%num_pistons > 0 ) then
    ERROR('Pistons are currently offline')
    call abort_program( p_err_invalid )
  endif
#endif

!---------------------------------advance momenta-----------------------------------------
  
! if this is a free streaming species then return 
  if ( this%free_stream ) return
  
  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )
  e_coef = real(3.0_p_double * this%alpha**2 * this%rq_real, p_k_part)
  b_coef = real(this%alpha * dt * this%rq_real, p_k_part)
  
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

     call this % get_emf( emf, bp, ep, np, ptrcur )
                  
     pp = ptrcur
     do i=1,np
        utemp(1,i) = this%p(1,pp)
        utemp(2,i) = this%p(2,pp)
        utemp(3,i) = this%p(3,pp)

        ! Get time centered gamma before rotation for E field correction
        u2 = utemp(1,i)**2 + utemp(2,i)**2 + utemp(3,i)**2
        gamma = sqrt( u2 + 1 )
     
        ! Variation of the E filed
        ep(1,i) = (ep(1,i) + e_coef*gamma*this%x(1,pp)) * tem
        ep(2,i) = ep(2,i) * tem
        ep(3,i) = ep(3,i) * tem
        
        pp = pp + 1
     end do

     ! Perform first half of electric field acceleration.
     ! and get time centered gamma

     loc_ene = 0

     pp = ptrcur
     do i=1,np
        ptemp(2,i) = this%p(2,pp)
        utemp(1,i) = this%p(1,pp) + ep(1,i)
        utemp(2,i) = this%p(2,pp) + ep(2,i)
        utemp(3,i) = this%p(3,pp) + ep(3,i)

        ! Get time centered gamma
        u2 = utemp(1,i)**2 + utemp(2,i)**2 + utemp(3,i)**2
        gamma = sqrt( u2 + 1 )

        ! accumulate time centered energy
        loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1)
  
        gam_tem= tem / gamma

        pp = pp + 1
     enddo

     ! accumulate global energy
     energy = energy + loc_ene

     do i=1,np
       bp(1,i) = bp(1,i)*gam_tem(i)
       bp(2,i) = bp(2,i)*gam_tem(i)
       bp(3,i) = (bp(3,i)+b_coef*gamma)*gam_tem(i)
     end do

     ! Perform first half of the rotation and store in u.
     pp = ptrcur
     do i=1,np
        this%p(1,pp)=utemp(1,i)+utemp(2,i)*bp(3,i)-utemp(3,i)*bp(2,i)
        this%p(2,pp)=utemp(2,i)+utemp(3,i)*bp(1,i)-utemp(1,i)*bp(3,i)
        this%p(3,pp)=utemp(3,i)+utemp(1,i)*bp(2,i)-utemp(2,i)*bp(1,i)
        pp = pp + 1
     end do
  
     do i=1,np
        otsq(i) = 2. / (1.0_p_k_part+bp(1,i)**2+bp(2,i)**2+bp(3,i)**2)
     end do

     do i=1,np
        bp(1,i)= bp(1,i) * otsq(i)
        bp(2,i)= bp(2,i) * otsq(i)
        bp(3,i)= bp(3,i) * otsq(i)
     end do

     ! Perform second half of the rotation.
     pp = ptrcur
     do i=1,np
        utemp(1,i) = utemp(1,i)+this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)
        utemp(2,i) = utemp(2,i)+this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)
        utemp(3,i) = utemp(3,i)+this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)
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
  
!--------advance position of particles i0 to i1 in chunks of p_cache_size-----------------

  ! loop through all the particles          
  do ptrcur = i0, i1, p_cache_size

     ! check if last copy of table and set np
     if( ptrcur + p_cache_size > i1 ) then
         np = i1 - ptrcur + 1
     else
         np = p_cache_size
     endif
     
     ! calculate dt/gamma
     pp = ptrcur
     do i=1, np
       rgamma(i) = 1.0_p_k_part / &
          sqrt( 1.0_p_k_part + &
          this%p(1,pp)**2 + &
          this%p(2,pp)**2 + &
          this%p(3,pp)**2 )
          pp = pp + 1
     end do

     ! advance particle position
     pp = ptrcur
     do i=1,np

       xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rgamma(i) * dt_dx(1)
       xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rgamma(i) * dt_dx(2)
       xbuf(3,i) = this%x(3,pp) + this%p(3,pp) * rgamma(i) * dt_dx(3)
       
       dxi(1,i) = ntrim( xbuf(1,i) )
       dxi(2,i) = ntrim( xbuf(2,i) )
       dxi(3,i) = ntrim( xbuf(3,i) )
                   
       pp = pp + 1

     end do

     ! Deposit current
     call this % dep_current_3d(  jay, dxi, xbuf, &
                       this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                       this%q(ptrcur:), np, gdt )

     pp = ptrcur
     do i = 1, np
       this%x(1,pp)  = xbuf(1,i)  - dxi(1,i) 
       this%x(2,pp)  = xbuf(2,i)  - dxi(2,i)
       this%x(3,pp)  = xbuf(3,i)  - dxi(3,i)
       this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
       this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
       this%ix(3,pp) = this%ix(3,pp) + dxi(3,i)
 
       pp = pp + 1
     end do

  enddo
  
end subroutine advance_deposit_3d_shear
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

end module m_species_shear
