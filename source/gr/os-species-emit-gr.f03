#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_emit_gr

#include "memory/memory.h"

use m_system
use m_parameters

use m_emf_define, only : t_emf
use m_vdf_define, only : t_vdf
use m_time_step, only : t_time_step
use m_species_define_gr, only : t_species_gr

implicit none

private

interface emit_qed_red
    module procedure emit_qed_red
end interface emit_qed_red

interface emit_qed_red_gr
    module procedure emit_qed_red_gr
end interface emit_qed_red_gr

public :: emit_qed_red, emit_qed_red_gr

contains

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine emit_qed_red( this, dt, ptrcur, np )

  implicit none

  class( t_species_gr ), intent(inout) :: this
  real( p_double ), intent(in) :: dt
  integer, intent(in) :: ptrcur, np

  integer :: i, pp, ix
  real( p_double ) :: r, p2, g, ppair_mod, pp_gth
  real( p_k_part ), dimension(p_p_dim) :: ppair

  integer :: par_idx, num_par_max_new
  integer, dimension(:,:), pointer :: par_ix
  real( p_k_part ), dimension(:,:), pointer :: par_x
  real( p_k_part ), dimension(:,:), pointer :: par_p
  real( p_k_part ), dimension(:), pointer :: par_q
  
  if (this%num_par == 0) return
  
  call alloc( par_ix, (/p_x_dim, np/) )
  call alloc( par_x, (/p_x_dim+3, np/) )
  call alloc( par_p, (/p_p_dim, np/) )
  call alloc( par_q, (/np/) )
  par_idx = 0

  ppair_mod = sqrt(0.25*this%pp_gpair**2 - 1.)

  !loop through all particles in the bunch
  pp = ptrcur
  do i=1,np

    ix = this%ix(1,pp)
    p2 = this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2
    g = sqrt(1.0d0+p2)
    r = this%geometry%r%f1(2,ix) + this%x(1,pp)*this%geometry%dr%f1(1,ix)

    pp_gth = this%pp_gth * min(1+.01/this%geometry%st%f1(2, this%ix(2,pp))**1.5_p_double,&
                               40.0_p_double)

    if ( r > this%pp_rmin .and. r < this%pp_rmax .and. g > pp_gth ) then
      ppair = ppair_mod * this%p(:,pp) / sqrt(p2)
      this%p(:,pp) = (1. - this%pp_gpair/g) * this%p(:,pp)
      
      ! save properties of pair to be created
      par_idx = par_idx + 1
      par_ix(:,par_idx) = this%ix(1:p_x_dim,pp)
      par_x(:,par_idx) = this%x(1:p_x_dim+3,pp)
      par_p(:,par_idx) = ppair
      par_q(par_idx) = this%q(pp)

      if ( this%track_nemit ) then
        this%nemit%f2(1, this%ix(1,pp), this%ix(2,pp)) = this%nemit%f2(1, this%ix(1,pp), this%ix(2,pp)) + 1
      endif
    endif
    pp = pp + 1
  end do

  if (par_idx > 0) then
    ! if ( this%num_par + par_idx > this%num_par_max ) then
    !   num_par_max_new = max(this%num_par_max / 4, this%num_par + par_idx)
    !   call this%grow_buffer( num_par_max_new )
    ! endif

    if ( this%grp_recip%num_par + par_idx > this%grp_recip%num_par_max ) then
      num_par_max_new = max(this%grp_recip%num_par_max / 4, this%grp_recip%num_par + par_idx)
      call this%grp_recip%grow_buffer( num_par_max_new )
    endif

    do i=1,par_idx
      call this%create_particle( &
        par_ix(:,i), & 
        par_x(:,i), &
        par_p(:,i), &
        par_q(i) &
      )

      call this%grp_recip%create_particle( &
        par_ix(:,i), & 
        par_x(:,i), &
        par_p(:,i), &
        -par_q(i) &
      )
    enddo
  endif

  call freemem( par_ix )
  call freemem( par_x )
  call freemem( par_p )
  call freemem( par_q )

end subroutine emit_qed_red
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine emit_qed_red_gr( this, dt, ptrcur, np )

  implicit none

  class( t_species_gr ), intent(inout) :: this
  real( p_double ), intent(in) :: dt
  integer, intent(in) :: ptrcur, np

  integer :: i, pp, ix
  real( p_double ) :: r, p2, g, ppair_mod, pp_gth
  real( p_k_part ), dimension(p_p_dim) :: ppair

  integer :: par_idx, num_par_max_new
  integer, dimension(:,:), pointer :: par_ix
  real( p_k_part ), dimension(:,:), pointer :: par_x
  real( p_k_part ), dimension(:,:), pointer :: par_p
  real( p_k_part ), dimension(:), pointer :: par_q
  
  if (this%num_par == 0) return
  
  call alloc( par_ix, (/p_x_dim, np/) )
  call alloc( par_x, (/p_x_dim+1, np/) )
  call alloc( par_p, (/p_p_dim, np/) )
  call alloc( par_q, (/np/) )
  par_idx = 0

  ppair_mod = sqrt(0.25*this%pp_gpair**2 - 1.)

  !loop through all particles in the bunch
  pp = ptrcur
  do i=1,np

    ix = this%ix(1,pp)
    p2 = this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2
    
    r = this%geometry%r%f1(2,ix) + this%x(1,pp) * this%geometry%dr%f1(1,ix)
    g = sqrt(1.0d0+p2)

    pp_gth = this%pp_gth * min(1+.01/this%geometry%st%f1(2, this%ix(2,pp))**1.5_p_double,&
                               40.0_p_double)

    if ( r > this%pp_rmin .and. r < this%pp_rmax .and. g > pp_gth ) then
      ppair = ppair_mod * this%p(:,pp) / sqrt(p2)
      this%p(:,pp) = (1.-this%pp_gpair/g) * this%p(:,pp)
      
      ! save properties of pair to be created
      par_idx = par_idx + 1
      par_ix(:,par_idx) = this%ix(1:p_x_dim,pp)
      par_x(:,par_idx) = this%x(1:p_x_dim+1,pp)
      par_p(:,par_idx) = ppair
      par_q(par_idx) = this%q(pp)

      if ( this%track_nemit ) then
        this%nemit%f2(1, this%ix(1,pp), this%ix(2,pp)) = this%nemit%f2(1, this%ix(1,pp), this%ix(2,pp)) + 1
      endif
    endif
    pp = pp + 1
  end do

  if (par_idx > 0) then
    ! if ( this%num_par + par_idx > this%num_par_max ) then
    !   num_par_max_new = max(this%num_par_max / 4, this%num_par + par_idx)
    !   call this%grow_buffer( num_par_max_new )
    ! endif

    if ( this%grp_recip%num_par + par_idx > this%grp_recip%num_par_max ) then
      num_par_max_new = max(this%grp_recip%num_par_max / 4, this%grp_recip%num_par + par_idx)
      call this%grp_recip%grow_buffer( num_par_max_new )
    endif

    do i=1,par_idx
      call this%create_particle( &
        par_ix(:,i), & 
        par_x(:,i), &
        par_p(:,i), &
        par_q(i) &
      )

      call this%grp_recip%create_particle( &
        par_ix(:,i), & 
        par_x(:,i), &
        par_p(:,i), &
        -par_q(i) &
      )
    enddo
  endif

  call freemem( par_ix )
  call freemem( par_x )
  call freemem( par_p )
  call freemem( par_q )

end subroutine emit_qed_red_gr

end module m_species_emit_gr

!-------------------------------------------------------------------------------
subroutine emit_species_gr( this, emf, dt, i0, i1 )
  
  use m_system
  use m_parameters

  use m_species_define_gr, only : t_species_gr
  use m_species_group_gr
  use m_emf_define, only : t_emf
  use m_species_emit_gr

  implicit none

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  real(p_double),                 intent(in) :: dt
  integer, intent(in) :: i0, i1
  
  ! cooling type
  select case ( this%cooling_type )
  case(p_cooling_qed_red)
    select case (this%push_type)
    case (p_std)
      call emit_qed_red( this, dt, i0, i1 )
    case default
      call emit_qed_red_gr( this, dt, i0, i1 )
    end select
  case(p_cooling_qed_red_soft)
    ! call emit_qed_red_soft( this, dt, i0, i1 )
  case(p_cooling_qed)
    ! call emit_qed( this, dt, i0, i1 )
  case default
    ! do nothinng
  end select

end subroutine emit_species_gr
!-------------------------------------------------------------------------------