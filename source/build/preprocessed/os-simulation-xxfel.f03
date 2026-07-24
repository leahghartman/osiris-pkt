# 1 "xxfel/os-simulation-xxfel.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "xxfel/os-simulation-xxfel.f03" 2

!
! This simulation module allowes more precise external field calculations
! which aid in tracking single-particle (or small number) motion particle motion
! like when studing XXFEL systems.
!
module m_simulation_xxfel

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
# 10 "xxfel/os-simulation-xxfel.f03" 2
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
# 11 "xxfel/os-simulation-xxfel.f03" 2

use m_parameters
use m_simulation, only : t_simulation
use m_species_define, only : t_species
use m_particles_define, only : t_particles
use m_space, only : t_space
use m_emf_define, only : t_emf, p_emfsmooth_local, p_emfsmooth_none, p_emfsmooth_stand, &
                          p_emfsmooth_nci, fsmoothev
use m_emf_ncifilter, only : filter_nci
use m_logprof, only : begin_event, end_event
use m_emf_gridcenter, only : grid_center
use m_vdf_smooth, only : smooth
use m_emf_interpolate, only : get_emf
use m_fparser, only : t_fparser, eval, p_k_fparse


implicit none

private

type, extends( t_simulation ) :: t_simulation_xxfel
  ! No member-specific data
  contains
  procedure :: allocate_objs => allocate_objs_sim_xxfel
end type t_simulation_xxfel


type, extends( t_particles ) :: t_particles_xxfel
  contains
  procedure :: allocate_objs => allocate_objs_part_xxfel
end type t_particles_xxfel


type, extends( t_species ) :: t_species_xxfel
  contains
  procedure :: set_dudt => set_dudt_spec_xxfel
  procedure :: get_emf => get_emf_spec_xxfel
end type t_species_xxfel


type, extends( t_emf ) :: t_emf_xxfel
  contains
  procedure :: update_particle_fld => update_particle_fld_emf_xxfel
end type t_emf_xxfel

public :: t_simulation_xxfel

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_sim_xxfel( sim )

  implicit none

  class( t_simulation_xxfel ), intent(inout) :: sim

  ! allocate custom emf and particles objects
  allocate( t_emf_xxfel :: sim%emf )
  allocate( t_particles_xxfel :: sim%part )

  ! call parent class
  call sim%t_simulation%allocate_objs()

end subroutine allocate_objs_sim_xxfel


!-----------------------------------------------------------------------------------------
! Allocate particle species, cathodes, etc.
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_xxfel( this )

  implicit none

  class( t_particles_xxfel ), intent(inout) :: this

  integer :: i
  class(t_species), pointer :: spec, tail

  ! call parent class
  call this % t_particles % allocate_objs ( )

  ! Kill the spcies objects made in the parent class
  ! A little messy but these types of simulations rarely
  ! have much associated species data.
  ! No need to call cleanup since read or init hasn't been called
  if ( this%num_species > 0 ) then
    spec => this % species
    tail => spec % next
    deallocate(spec)
    do i = 2, this%num_species
      spec => tail
      tail => spec % next
      deallocate(spec)
    enddo

    ! make the ones we want
    allocate( t_species_xxfel :: spec )
    this % species => spec
    do i = 2, this%num_species
      tail => spec
      allocate( t_species_xxfel :: spec )
      tail % next => spec
    enddo
  endif

end subroutine allocate_objs_part_xxfel


!-----------------------------------------------------------------------------------------
! Set the correct dudt
!-----------------------------------------------------------------------------------------
subroutine set_dudt_spec_xxfel( this )

  implicit none

  class( t_species_xxfel ), intent(inout) :: this

  call this % t_species % set_dudt()

  ! check that pusher type is valid for fel
  select case( this % push_type )
  case(p_std,p_vay,p_fullrot,p_euler,p_cond_vay,p_cary,p_exact,p_exact_rr)
    ! These pushers are valid. More can be added to this list once verified
  case default
    write(err_buf__,*) "Invalid pusher type for simulation mode 'xxfel'";call err__("xxfel/os-simulation-xxfel.f03",137)
    write(err_buf__,*) "Valid pushers are 'standard', 'vay', 'fullrot', 'euler',";call err__("xxfel/os-simulation-xxfel.f03",138)
    write(err_buf__,*) "'cond_vay', 'cary', 'exact'/'analytic' and 'exact-rr'/'analytic-rr'.";call err__("xxfel/os-simulation-xxfel.f03",139)
    call abort_program( p_err_invalid )
  end select

end subroutine set_dudt_spec_xxfel


!-----------------------------------------------------------------------------------------
! Interpolate fields at particle positions for cell based positions.
! in the FEL case, the external field is evaluated @ the center of mass by
! evaluating the parser directly
!-----------------------------------------------------------------------------------------
subroutine get_emf_spec_xxfel( this, emf, bp, ep, np, ptrcur, t )

  implicit none

  class( t_species_xxfel ), intent(in) :: this
  class( t_emf ), intent(in), target :: emf
  real(p_k_part), dimension(:,:), intent(out) :: bp, ep
  integer, intent(in) :: np, ptrcur
  real(p_double), intent(in), optional :: t

  integer :: ptrend, i, j
  real(p_double) , dimension(p_x_dim,p_cache_size) :: global_x
  real(p_k_fparse) :: eval_pos(p_max_dim+1)

  ! This temp_e/b stuff is here to get around the fact that t_emf
  ! is passed in as intent(in), but the 'eval' function needs
  ! intent(inout) to be passed.
  type( t_fparser ), pointer :: temp_e, temp_b

  ptrend = ptrcur + np - 1
  call this % get_position( ptrcur, ptrend, global_x )

  if ( present(t) ) then

    call get_emf( emf, bp, ep, this%ix(:,ptrcur:), this%x(:,ptrcur:), np, &
                  this%interpolation )

    ! add in exact evaluation of the external fields to E/B. For XXFEL
    do i = 1, np
      do j = 1, p_f_dim
        eval_pos(1:p_x_dim) = global_x(1:p_x_dim,i)
        eval_pos(p_x_dim+1) = t
        temp_e => emf%ext_emf%mfunc_e(j)
        temp_b => emf%ext_emf%mfunc_b(j)
        ep(j,i) = ep(j,i) + real( eval(temp_e , eval_pos) , p_k_fld)
        bp(j,i) = bp(j,i) + real( eval(temp_b , eval_pos) , p_k_fld)
        ! ep(j,i) = ep(j,i) + real( eval(emf%ext_emf%mfunc_e(j) , eval_pos) , p_k_fld)
        ! bp(j,i) = bp(j,i) + real( eval(emf%ext_emf%mfunc_b(j) , eval_pos) , p_k_fld)
      enddo
    enddo

  else
    write(err_buf__,*) 'This pusher is not correctly implemented for xxfel.';call err__("xxfel/os-simulation-xxfel.f03",193)
    write(err_buf__,*) 'Argument t is missing.';call err__("xxfel/os-simulation-xxfel.f03",194)
    call abort_program( p_err_invalid )
  endif

end subroutine get_emf_spec_xxfel


!-----------------------------------------------------------------------------------------
subroutine update_particle_fld_emf_xxfel( this, n, t, dt, space, nx_p_min )
!-----------------------------------------------------------------------------------------
! Update fields to be used by particles
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_xxfel ), intent( inout ) :: this
  integer, intent(in) :: n
  real( p_double ), intent(in) :: t, dt
  type( t_space ), intent(in) :: space
  integer, dimension(:), intent(in) :: nx_p_min

  ! Smoothing and external fields

  ! XXFEL: Don't call set_fld_values here, ext fields handled in dudt

  ! do local smooth of fields if necessary
  if ( this%smooth_type == p_emfsmooth_local ) then
    ! Stop smoothing after given iteration
    if (this%smooth_nmax > 0 .and. n >= this%smooth_nmax) then
      this%smooth_type = p_emfsmooth_none
    else if ( mod( n, this%smooth_niter ) == 0) then

      ! smooth e and b fields
      call begin_event( fsmoothev )

      call smooth( this%e, this%smooth )
      call smooth( this%b, this%smooth )

      call end_event( fsmoothev )
    endif
  endif

  ! Update fields for particle interpolation
  if ( this%part_fld_alloc ) then

    ! Copy current emf field values to e_part and b_part
    this%e_part = this%e
    this%b_part = this%b

    ! smooth e and b fields using standard smoothing
    select case ( this%smooth_type )
    case( p_emfsmooth_stand )
      call begin_event( fsmoothev )

      call smooth( this%e_part, this%smooth )
      call smooth( this%b_part, this%smooth )

      call end_event( fsmoothev )

    case( p_emfsmooth_nci )
      call begin_event( fsmoothev )
      call filter_nci( this%e_part, this%b_part, dt, this%part_grid_center )
      call end_event( fsmoothev )

    end select

    ! XXFEL: Don't add fields here, ext fields handled in dudt

    ! Grid center fields if necessary
    if ( this%part_grid_center ) then
      call grid_center( this )
    endif

  endif

end subroutine update_particle_fld_emf_xxfel
!-----------------------------------------------------------------------------------------


end module
