# 1 "gr/os-species-emit-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-species-emit-gr.f03" 2
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
# 2 "gr/os-species-emit-gr.f03" 2
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
# 3 "gr/os-species-emit-gr.f03" 2

module m_species_emit_gr

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
# 7 "gr/os-species-emit-gr.f03" 2

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

  call alloc(par_ix, (/p_x_dim, np/),"gr/os-species-emit-gr.f03",55)
  call alloc(par_x, (/p_x_dim+3, np/),"gr/os-species-emit-gr.f03",56)
  call alloc(par_p, (/p_p_dim, np/),"gr/os-species-emit-gr.f03",57)
  call alloc(par_q, (/np/),"gr/os-species-emit-gr.f03",58)
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
    ! num_par_max_new = max(this%num_par_max / 4, this%num_par + par_idx)
    ! call this%grow_buffer( num_par_max_new )
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

  call freemem(par_ix,"gr/os-species-emit-gr.f03",121)
  call freemem(par_x,"gr/os-species-emit-gr.f03",122)
  call freemem(par_p,"gr/os-species-emit-gr.f03",123)
  call freemem(par_q,"gr/os-species-emit-gr.f03",124)

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

  call alloc(par_ix, (/p_x_dim, np/),"gr/os-species-emit-gr.f03",150)
  call alloc(par_x, (/p_x_dim+1, np/),"gr/os-species-emit-gr.f03",151)
  call alloc(par_p, (/p_p_dim, np/),"gr/os-species-emit-gr.f03",152)
  call alloc(par_q, (/np/),"gr/os-species-emit-gr.f03",153)
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
    ! num_par_max_new = max(this%num_par_max / 4, this%num_par + par_idx)
    ! call this%grow_buffer( num_par_max_new )
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

  call freemem(par_ix,"gr/os-species-emit-gr.f03",217)
  call freemem(par_x,"gr/os-species-emit-gr.f03",218)
  call freemem(par_p,"gr/os-species-emit-gr.f03",219)
  call freemem(par_q,"gr/os-species-emit-gr.f03",220)

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
  class( t_emf ), intent( in ) :: emf
  real(p_double), intent(in) :: dt
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
