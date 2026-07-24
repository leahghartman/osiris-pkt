# 1 "cyl_modes/os-current-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-current-cyl-modes.f03" 2

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
# 3 "cyl_modes/os-current-cyl-modes.f03" 2
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
# 4 "cyl_modes/os-current-cyl-modes.f03" 2

module m_current_cyl_modes

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
# 8 "cyl_modes/os-current-cyl-modes.f03" 2

use m_system
use m_parameters
use m_cyl_modes
use m_current_define, only : t_current, jayboundev, smoothev
use m_current_diag, only : t_current_diag, diag_current_ev
use m_current_boundary, only : spec_bc_2d
use m_vdf_memory, only : alloc, freemem
use m_restart, only : t_restart_handle
use m_grid_define, only: t_grid
use m_node_conf, only: t_node_conf, n_threads
use m_grid_cyl_modes, only : t_grid_cyl_modes
use m_vdf, only : link_vdf
use m_logprof, only : begin_event, end_event, create_event
use m_vdf_smooth, only : smooth
use m_vdf_define, only : t_vdf, t_vdf_report
use m_vdf_comm, only : t_vdf_msg
use m_time_step, only : t_time_step
use m_vdf_report, only : if_report
use m_space, only : t_space

implicit none

private

!-----------------------------------------------------------------------------------------
! t_current_cyl_modes class definition
!-----------------------------------------------------------------------------------------
type, extends( t_current ) :: t_current_cyl_modes

  integer :: n_cyl_modes = -1

  ! This will point to jay_cyl_m_arr(1)
  type( t_cyl_modes ) :: jay_cyl_m

  ! We need an array for using multiple threads
  type( t_cyl_modes ), dimension(:), pointer :: jay_cyl_m_arr => null()

  contains
  procedure :: allocate_objs => allocate_objs_current_cyl_modes
  procedure :: init => init_current_cyl_modes
  procedure :: cleanup => cleanup_current_cyl_modes
  procedure :: smooth => smooth_current_cyl_modes
  procedure :: normalize_cyl => normalize_cyl_current_cyl_modes

  procedure :: update_boundary => update_boundary_current_cyl_modes

  ! diagnostic report
  procedure :: report => report_current_cyl_modes

  procedure :: reshape => reshape_current_cyl_modes

end type t_current_cyl_modes

type, extends( t_current_diag ) :: t_current_diag_cyl_modes

  integer :: n_cyl_modes

  contains
  procedure :: avail_report_quants => avail_report_quants_current_cyl_modes
  procedure :: init_report_quants => init_report_quants_current_cyl_modes
  procedure :: init => init_diag_current_cyl_modes

end type t_current_diag_cyl_modes

! diagnostic quantities specific to the quasi-3D geometry
character(len=11), dimension(4), parameter :: p_report_quants = &
   (/ 'j1_cyl_m   ', 'j2_cyl_m   ', 'j3_cyl_m   ', 'div_j_cyl_m'/)

integer, parameter :: p_j1_cyl_m = 1, p_j2_cyl_m = 2, p_j3_cyl_m = 3, p_div_j_cyl_m = 4

public :: t_current_cyl_modes
public :: t_current_diag_cyl_modes

!-----------------------------------------------------------------------------------------
contains

subroutine allocate_objs_current_cyl_modes( this )

  implicit none

  class( t_current_cyl_modes ), intent(inout) :: this

  ! Allocate default t_diag_emf object class
  if ( .not. associated( this%diag ) ) then
    allocate( t_current_diag_cyl_modes :: this%diag )
  endif

  ! call superclass's allocate
  call this % t_current % allocate_objs( )

end subroutine allocate_objs_current_cyl_modes

!-----------------------------------------------------------------------------------------
subroutine init_current_cyl_modes( this, grid, gc_min, dx, no_co, if_move, &
                          bcemf_type, interpolation, restart, restart_handle, sim_options )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_current_cyl_modes ), intent( inout ) :: this
  class( t_grid ), intent(in) :: grid
  integer, dimension(:,:), intent(in) :: gc_min
  real(p_double), dimension(:), intent(in) :: dx
  class( t_node_conf ), intent(in) :: no_co
  logical, dimension(:), intent(in) :: if_move
  integer, dimension(:,:), intent(in) :: bcemf_type
  integer, intent(in) :: interpolation

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type( t_options ), intent(in) :: sim_options
  integer :: nt, i

  ! setup superclass data
  call this%t_current%init(grid,gc_min,dx,no_co,if_move,bcemf_type, &
                            interpolation,restart,restart_handle,sim_options)

  ! get n_cyl_modes from grid
  select type( grid )
  class is ( t_grid_cyl_modes )
    this%n_cyl_modes = grid%n_cyl_modes
  end select

  nt = n_threads( no_co )

  ! allocate and set up array jay_cyl_m_arr, one for each thread
  allocate( this%jay_cyl_m_arr(nt) )
  do i = 1, n_threads(no_co)
    call setup( this%jay_cyl_m_arr(i), this%n_cyl_modes, this%pf(i) )
  enddo

  ! allocate and link jay_cyl_m to jay_cyl_m_arr(1) for use everywhere besides the pusher
  call alloc(this%jay_cyl_m%pf_re, (/ 0 /), (/ this%n_cyl_modes /),"cyl_modes/os-current-cyl-modes.f03",141)
  if (this%n_cyl_modes>0) then
    call alloc(this%jay_cyl_m%pf_im, (/ 1 /), (/ this%n_cyl_modes /),"cyl_modes/os-current-cyl-modes.f03",143)
  endif

  call link_vdf( this%jay_cyl_m%pf_re(0), this%jay_cyl_m_arr(1)%pf_re(0) )
  do i = 1, this%n_cyl_modes
    call link_vdf( this%jay_cyl_m%pf_re(i), this%jay_cyl_m_arr(1)%pf_re(i) )
    call link_vdf( this%jay_cyl_m%pf_im(i), this%jay_cyl_m_arr(1)%pf_im(i) )
  enddo

end subroutine init_current_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_current_cyl_modes(this)
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_current_cyl_modes ), intent( inout ) :: this

  integer :: i

  call this%t_current%cleanup()

  ! cleanup local data structures
  do i = 1, size( this%jay_cyl_m_arr )
    call cleanup( this%jay_cyl_m_arr(i) )
  enddo
  deallocate( this%jay_cyl_m_arr )

  ! we only need to freemem jay_cyl_m since it links to jay_cyl_m_arr(1)
  call freemem(this%jay_cyl_m%pf_re,"cyl_modes/os-current-cyl-modes.f03",174)
  call freemem(this%jay_cyl_m%pf_im,"cyl_modes/os-current-cyl-modes.f03",175)

end subroutine cleanup_current_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine smooth_current_cyl_modes( this )
!-----------------------------------------------------------------------------------------
! smoothing for electric current
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_current_cyl_modes ), intent(inout) :: this

  integer :: mode

  call begin_event(smoothev)

  do mode = 0, ubound(this%jay_cyl_m%pf_re,1)
    call smooth( this%jay_cyl_m%pf_re(mode) , this%curr_smooth)
    if (mode > 0) call smooth( this%jay_cyl_m%pf_im(mode) , this%curr_smooth)
  enddo

  call end_event(smoothev)

end subroutine smooth_current_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine norm_ring_grid_cyl_modes( j_re, r_dr, gshift_i2, mode, j_im )

  implicit none

  type( t_vdf ), intent(inout) :: j_re
  real(p_k_fld), intent(in) :: r_dr
  integer, intent(in) :: gshift_i2
  integer, intent(in) :: mode
  type( t_vdf ), intent(inout), optional :: j_im

  real(p_k_fld) :: vol_fac_cor, vol_fac_mid, mode_factor
  integer :: i1, i2

  if (mode > 0 .and. (.not.present(j_im))) then
    write(err_buf__,*) 'Must include j_im for mode > 0';call err__("cyl_modes/os-current-cyl-modes.f03",219)
    call abort_program(p_err_invalid)
  endif

  if ( mode == 0 ) then
    mode_factor = 1.0_p_k_fld
  else
    ! The higher Fourier modes actually have an extra factor of 2
    mode_factor = 2.0_p_k_fld
  endif

  do i2 = lbound(j_re%f2, 3), ubound(j_re%f2, 3)

    ! Get normalization factors 1/r at cell corner and middle
    vol_fac_cor = r_dr / ABS((i2 + gshift_i2 - 0.5_p_k_fld))
    if ( i2 + gshift_i2 /= 0 ) then
      vol_fac_mid = r_dr / ABS((i2 + gshift_i2))
    else
      ! Turns out we will re-define the j_r on-axis component again below, but we keep
      ! this here to show how we would properly normalize it
      if (mode == 1) then
        vol_fac_mid = real(8.0, p_k_fld)*r_dr ! this is 2 x 1/(dr/2)**2
      else
        vol_fac_mid = 0.0_p_k_fld
      endif
    endif

    vol_fac_cor = vol_fac_cor * mode_factor
    vol_fac_mid = vol_fac_mid * mode_factor

    if (mode == 0) then

      do i1 = lbound(j_re%f2, 2), ubound(j_re%f2, 2)
        j_re%f2(1,i1,i2) = j_re%f2(1,i1,i2) * vol_fac_cor
        j_re%f2(2,i1,i2) = j_re%f2(2,i1,i2) * vol_fac_mid
        j_re%f2(3,i1,i2) = j_re%f2(3,i1,i2) * vol_fac_cor
      enddo

    else

      do i1 = lbound(j_re%f2, 2), ubound(j_re%f2, 2)
        j_re%f2(1,i1,i2) = j_re%f2(1,i1,i2) * vol_fac_cor
        j_re%f2(2,i1,i2) = j_re%f2(2,i1,i2) * vol_fac_mid
        j_re%f2(3,i1,i2) = j_re%f2(3,i1,i2) * vol_fac_cor
        j_im%f2(1,i1,i2) = j_im%f2(1,i1,i2) * vol_fac_cor
        j_im%f2(2,i1,i2) = j_im%f2(2,i1,i2) * vol_fac_mid
        j_im%f2(3,i1,i2) = j_im%f2(3,i1,i2) * vol_fac_cor
      enddo

    endif
  enddo

  ! process axial boundary
  ! This could be moved into update_boundary_current to overlap this calculation with
  ! communication, but this is cleaner and makes this routine interchangeable with
  ! (update_boundary and smooth).

  if ( gshift_i2 < 0 ) then

    ! note that pf%f2(1,:,1) and pf%f2(3,:,1) are off axis, but that pf%f2(2,:,1) is on axis
    ! first take care of the "physical" cell on axis - this is not an actual guard cell
    ! The factors for folding in j_z, j_r and j_phi are all different, but verified with
    ! charge conservation tests.
    ! Note, the final values in the guard cells are really meaningless, so we just set
    ! them to the same value as for |r|.
    ! if quartic interpolation is ever implemented, also handle cell -1 in r

    if (mode==0) then

      do i1 = lbound( j_re%f2, 2 ), ubound( j_re%f2, 2 )
        j_re%f2(1,i1,2) = j_re%f2(1,i1,2) + j_re%f2(1,i1,1)
        j_re%f2(1,i1,3) = j_re%f2(1,i1,3) + j_re%f2(1,i1,0)

        ! The j%f2(2,i1,1) value was already set to 0 above
        j_re%f2(2,i1,2) = j_re%f2(2,i1,2) - j_re%f2(2,i1,0)

        j_re%f2(3,i1,2) = j_re%f2(3,i1,2) - j_re%f2(3,i1,1)
        j_re%f2(3,i1,3) = j_re%f2(3,i1,3) - j_re%f2(3,i1,0)

        ! Set values for guard cells, won't really affect anything
        j_re%f2(1,i1,1) = j_re%f2(1,i1,2)
        j_re%f2(1,i1,0) = j_re%f2(1,i1,3)

        j_re%f2(2,i1,0) = j_re%f2(2,i1,2)

        j_re%f2(3,i1,1) = j_re%f2(3,i1,2)
        j_re%f2(3,i1,0) = j_re%f2(3,i1,3)
      enddo

    elseif (mode==1) then

      do i1 = lbound( j_re%f2, 2 ), ubound( j_re%f2, 2 )
        j_re%f2(1,i1,2) = j_re%f2(1,i1,2) + j_re%f2(1,i1,1)
        j_re%f2(1,i1,3) = j_re%f2(1,i1,3) + j_re%f2(1,i1,0)
        j_im%f2(1,i1,2) = j_im%f2(1,i1,2) + j_im%f2(1,i1,1)
        j_im%f2(1,i1,3) = j_im%f2(1,i1,3) + j_im%f2(1,i1,0)

        ! Rather than using the standard value of j_r(r=0) from the normalization, we
        ! instead ensure that dj_r/dr(r=0) = 0. This quantity doesn't actually show up
        ! in the field solver, so it doesn't matter that much.
        j_re%f2(2,i1,2) = j_re%f2(2,i1,2) - j_re%f2(2,i1,0)
        j_re%f2(2,i1,1) = 4.0*j_re%f2( 2, i1, 2 )/3.0 - j_re%f2( 2, i1, 3 )/3.0
        j_im%f2(2,i1,2) = j_im%f2(2,i1,2) - j_im%f2(2,i1,0)
        j_im%f2(2,i1,1) = 4.0*j_im%f2( 2, i1, 2 )/3.0 - j_im%f2( 2, i1, 3 )/3.0

        j_re%f2(3,i1,2) = j_re%f2(3,i1,2) + j_re%f2(3,i1,1)
        j_re%f2(3,i1,3) = j_re%f2(3,i1,3) + j_re%f2(3,i1,0)
        j_im%f2(3,i1,2) = j_im%f2(3,i1,2) + j_im%f2(3,i1,1)
        j_im%f2(3,i1,3) = j_im%f2(3,i1,3) + j_im%f2(3,i1,0)

        ! Set values for guard cells, won't really affect anything
        j_re%f2(1,i1,1) = j_re%f2(1,i1,2)
        j_re%f2(1,i1,0) = j_re%f2(1,i1,3)
        j_im%f2(1,i1,1) = j_im%f2(1,i1,2)
        j_im%f2(1,i1,0) = j_im%f2(1,i1,3)

        j_re%f2(2,i1,0) = j_re%f2(2,i1,2)
        j_im%f2(2,i1,0) = j_im%f2(2,i1,2)

        j_re%f2(3,i1,1) = j_re%f2(3,i1,2)
        j_re%f2(3,i1,0) = j_re%f2(3,i1,3)
        j_im%f2(3,i1,1) = j_im%f2(3,i1,2)
        j_im%f2(3,i1,0) = j_im%f2(3,i1,3)
      enddo

    else ! mode >= 2

      do i1 = lbound( j_re%f2, 2 ), ubound( j_re%f2, 2 )
        j_re%f2(1,i1,2) = j_re%f2(1,i1,2) + j_re%f2(1,i1,1)
        j_re%f2(1,i1,3) = j_re%f2(1,i1,3) + j_re%f2(1,i1,0)
        j_im%f2(1,i1,2) = j_im%f2(1,i1,2) + j_im%f2(1,i1,1)
        j_im%f2(1,i1,3) = j_im%f2(1,i1,3) + j_im%f2(1,i1,0)

        ! The j%f2(2,i1,1) value was already set to 0 above
        j_re%f2(2,i1,2) = j_re%f2(2,i1,2) - j_re%f2(2,i1,0)
        j_im%f2(2,i1,2) = j_im%f2(2,i1,2) - j_im%f2(2,i1,0)

        j_re%f2(3,i1,2) = j_re%f2(3,i1,2) + j_re%f2(3,i1,1)
        j_re%f2(3,i1,3) = j_re%f2(3,i1,3) + j_re%f2(3,i1,0)
        j_im%f2(3,i1,2) = j_im%f2(3,i1,2) + j_im%f2(3,i1,1)
        j_im%f2(3,i1,3) = j_im%f2(3,i1,3) + j_im%f2(3,i1,0)

        ! Set values for guard cells, won't really affect anything
        j_re%f2(1,i1,1) = j_re%f2(1,i1,2)
        j_re%f2(1,i1,0) = j_re%f2(1,i1,3)
        j_im%f2(1,i1,1) = j_im%f2(1,i1,2)
        j_im%f2(1,i1,0) = j_im%f2(1,i1,3)

        j_re%f2(2,i1,0) = j_re%f2(2,i1,2)
        j_im%f2(2,i1,0) = j_im%f2(2,i1,2)

        j_re%f2(3,i1,1) = j_re%f2(3,i1,2)
        j_re%f2(3,i1,0) = j_re%f2(3,i1,3)
        j_im%f2(3,i1,1) = j_im%f2(3,i1,2)
        j_im%f2(3,i1,0) = j_im%f2(3,i1,3)
      enddo

    endif

  endif


end subroutine norm_ring_grid_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Normalize deposited current for ring charges and take care of axial boundary.
! Note that current beyond axial boundary is reversed since r is considered to be < 0.
! This algorithm assumes that B1 is defined on the cylindrical axis.
!-----------------------------------------------------------------------------------------
subroutine normalize_cyl_current_cyl_modes( this )

  implicit none

  class( t_current_cyl_modes ), intent( inout ), target :: this

  real(p_k_fld) :: r_dr
  integer :: gshift_i2, mode





  r_dr = real( 1.0_p_double / this%pf(1)%dx_( p_r_dim), p_k_fld )
  gshift_i2 = this%gix_pos(2) - 2

  ! Normalize current for ring charges
  call norm_ring_grid_cyl_modes( this%jay_cyl_m%pf_re(0), r_dr, gshift_i2, 0 )
  do mode = 1, this%n_cyl_modes
    call norm_ring_grid_cyl_modes( this%jay_cyl_m%pf_re(mode), r_dr, gshift_i2, mode, &
                                   this%jay_cyl_m%pf_im(mode) )
  enddo





end subroutine normalize_cyl_current_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_boundary_current_cyl_modes( this, no_co, send_msg, recv_msg )

  implicit none

  class( t_current_cyl_modes ), intent(inout) :: this
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, j





  ! call superclass boundary update, takes care of mode 0
  call this % t_current % update_boundary( no_co, send_msg, recv_msg )

  call begin_event( jayboundev )

  ! call boundary update for cylindrical modes
  call update_boundary( this%jay_cyl_m, p_vdf_add, no_co, send_msg, recv_msg )

  ! TODO: possibly correct processing of specular boundaries here below
  ! This affects "conducting" and "reflecting" field boundaries

  ! Process specular boundaries, only handle 2-dimensional case
  do i = 1, 2
    if ( this%bc_type( p_lower, i ) == p_bc_specular ) then
      do j = 1, this%n_cyl_modes
        call spec_bc_2d( this%jay_cyl_m%pf_re(j), i, p_lower, this%interpolation )
        call spec_bc_2d( this%jay_cyl_m%pf_im(j), i, p_lower, this%interpolation )
      enddo
    endif

    if ( this%bc_type( p_upper, i ) == p_bc_specular ) then
      do j = 1, this%n_cyl_modes
        call spec_bc_2d( this%jay_cyl_m%pf_re(j), i, p_upper, this%interpolation )
        call spec_bc_2d( this%jay_cyl_m%pf_im(j), i, p_upper, this%interpolation )
      enddo
    endif
  enddo

  call end_event( jayboundev )





end subroutine update_boundary_current_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! overridden functions related to the diagnostic
!-----------------------------------------------------------------------------------------
subroutine report_current_cyl_modes( this, space, grid, no_co, tstep, t )

  implicit none

  class( t_current_cyl_modes ), intent(inout) :: this
  type( t_space ), intent(in) :: space
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t

  type( t_cyl_modes ) :: cyl_m_a
  type( t_vdf_report ), pointer :: rep

  call begin_event(diag_current_ev)

  ! run cyl_modes specific diagnostics
  rep => this%diag%reports
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then

      select case ( rep%quant )
      case ( p_j1_cyl_m, p_j2_cyl_m, p_j3_cyl_m )
        call report_cyl_modes( rep, this%jay_cyl_m, rep%quant - p_j1_cyl_m + 1, space, &
                               grid, no_co, tstep, t )

      case ( p_div_j_cyl_m )
        ! calculate electric field divergence
        call new_copy( cyl_m_a, this%jay_cyl_m, f_dim=1, zero=.true. )
        call div( this%jay_cyl_m, cyl_m_a, this%gix_pos(2), .false. )

        ! report it
        call report_cyl_modes( rep, cyl_m_a, 1, space, grid, no_co, tstep, t )

        ! free temporary memory
        call cleanup( cyl_m_a, cleanup_0=.true. )

      end select

    endif

    rep => rep%next
  enddo

  call end_event(diag_current_ev)

end subroutine report_current_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! overridden function for reshaping
!-----------------------------------------------------------------------------------------
subroutine reshape_current_cyl_modes( this, old_grid, new_grid, no_co, send_msg, recv_msg)

  implicit none

  class( t_current_cyl_modes ), intent(inout) :: this
  class( t_grid ), intent(in) :: old_grid, new_grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_vdf_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, j

  call this % t_current % reshape( old_grid, new_grid, no_co, send_msg, recv_msg )

  ! reshape the vdf objects
  do i = 1, size( this%pf )
    call reshape_nocopy( this%jay_cyl_m_arr(i), new_grid, this%pf(i) )
    do j = 1, this%n_cyl_modes
      call this % jay_cyl_m_arr(i) % pf_re(j) % zero()
      call this % jay_cyl_m_arr(i) % pf_im(j) % zero()
    enddo
  enddo

  call link_vdf( this%jay_cyl_m%pf_re(0), this%jay_cyl_m_arr(1)%pf_re(0) )
  do i = 1, this%n_cyl_modes
    call link_vdf( this%jay_cyl_m%pf_re(i), this%jay_cyl_m_arr(1)%pf_re(i) )
    call link_vdf( this%jay_cyl_m%pf_im(i), this%jay_cyl_m_arr(1)%pf_im(i) )
  enddo

end subroutine reshape_current_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function avail_report_quants_current_cyl_modes( this )

  class( t_current_diag_cyl_modes ), intent(in) :: this
  integer :: avail_report_quants_current_cyl_modes

  avail_report_quants_current_cyl_modes = size( p_report_quants )

end function avail_report_quants_current_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_report_quants_current_cyl_modes( this )

  class( t_current_diag_cyl_modes ), intent(inout) :: this

  if ( .not. associated( this % report_quants ) ) then
    allocate( this % report_quants( this % avail_report_quants( ) ) )
  endif

  ! These reports, and only these reports
  this % report_quants( 1 : size(p_report_quants) ) = p_report_quants

end subroutine init_report_quants_current_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_diag_current_cyl_modes( this, interpolation )

  implicit none

  class( t_current_diag_cyl_modes ), intent(inout) :: this
  integer, intent(in) :: interpolation

  type( t_vdf_report ), pointer :: report
  integer, parameter :: izero = ichar('0')
  real(p_double) :: interp_offset

  ! Current quantities only move in the x1 direction based on interpolation
  select case( interpolation )
  case( p_linear, p_cubic )
    interp_offset = 0.0_p_double
  case( p_quadratic, p_quartic )
    interp_offset = 0.5_p_double
  case default
    interp_offset = 0.0_p_double
    write(err_buf__,*) 'Interpolation value not supported';call err__("cyl_modes/os-current-cyl-modes.f03",605)
    call abort_program(p_err_invalid)
  end select

  ! setup cyl_modes specific diagnostics
  report => this%reports
  do
    if ( .not. associated( report ) ) exit

    report%xname = (/'x1', 'x2', 'x3'/)
    report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    ! these are just dummy values for now
    report%time_units = '1 / \omega_p'
    report%dt = 1.0

    ! current is 1/2 time step behind
    report%offset_t = -0.5_p_double

    report%fileLabel = ''
    report%basePath = trim(path_mass) // 'FLD' // p_dir_sep

    select case ( report%quant )
    case ( p_j1_cyl_m, p_j2_cyl_m, p_j3_cyl_m )
      report%label = 'j_'//char(izero + 1 + report%quant - p_j1_cyl_m)
      report%units = 'e \omega_p^2 / c'
      report%offset_x = (/ 0.0_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)
      report%offset_x(report%quant-p_j1_cyl_m+1) = report%offset_x(report%quant-p_j1_cyl_m+1) + 0.5_p_double

    case ( p_div_j_cyl_m )
      report%label = '\bf{\nabla}\cdot\bf{j}'
      report%units = 'e \omega_p^2 / c^2'
      report%offset_x = (/ 0.0_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)

    end select

    report => report%next
  enddo

  if (diag_current_ev==0) diag_current_ev = create_event('electric current diagnostics')

end subroutine init_diag_current_cyl_modes
!-----------------------------------------------------------------------------------------


end module m_current_cyl_modes
