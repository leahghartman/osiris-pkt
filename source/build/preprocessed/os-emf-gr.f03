# 1 "gr/os-emf-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-emf-gr.f03" 2
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
# 2 "gr/os-emf-gr.f03" 2
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
# 3 "gr/os-emf-gr.f03" 2

module m_emf_gr

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
# 7 "gr/os-emf-gr.f03" 2

use m_system
use m_parameters

use m_geometry_gr
use m_emf_define_gr

contains

!-------------------------------------------------------------------------------
! Advance the B field
!-------------------------------------------------------------------------------
subroutine dbdt_gr( this, dt, step )

  implicit none

  class( t_emf_gr ), intent( inout ) :: this
  real(p_double), intent(in) :: dt
  integer, intent(in) :: step

  select case ( p_x_dim )
    case (2)
      call dbdt_2d_gr( this, dt, step )
  end select

end subroutine

!-------------------------------------------------------------------------------
! Advance the B field (2D)
!-------------------------------------------------------------------------------
subroutine dbdt_2d_gr( this, dt, step )

  implicit none
  integer, parameter :: rank = 2

  ! dummy variables
  class( t_emf_gr ), intent( inout ) :: this
  real(p_double), intent(in) :: dt
  integer, intent(in) :: step

  ! local variables
  integer :: i1, i2
  real(p_k_fld) :: curl_er, curl_et, curl_ep
  real(p_k_fld) :: beta0
  real( p_k_fld ), pointer :: e(:,:,:), b(:,:,:)
  real( p_k_fld ), pointer :: lr(:,:), lt(:,:), ap(:,:), alpha(:,:)
  real( p_k_fld ), pointer :: lp(:,:,:), ar(:,:,:), at(:,:,:), beta(:,:,:)

  ! executable statements
  ! setup pointers to improve readibility
  e => this%e%f2
  b => this%b%f2

  lr => this%geometry%lr%f1
  lt => this%geometry%lt%f1
  lp => this%geometry%lp%f2

  ar => this%geometry%ar%f2
  at => this%geometry%at%f2
  ap => this%geometry%ap%f1

  alpha => this%geometry%alpha%f1
  beta => this%geometry%beta%f2
  beta0 = this%geometry%beta0

  if ( step == 1 ) then
    !$omp parallel do private(i1)
    do i2 = 0, this%e%nx_(p_tc_dim) + 1
      do i1 = 0, this%e%nx_(p_rc_dim) + 1

        curl_er = alpha(1, i1) * ( &
            lp(1, i1, i2+1) * e(3, i1, i2+1) - &
            lp(1, i1, i2) * e(3, i1, i2) &
          ) / ar(2, i1, i2)

        curl_et = ( &
            lp(1, i1, i2) * alpha(1, i1) * e(3, i1, i2) - &
            lp(1, i1+1, i2) * alpha(1, i1+1) * e(3, i1+1, i2) &
          ) / at(1, i1, i2)

        curl_ep = ( &
            lr(2, i1) * ( &
              alpha(2, i1) * ( e(1, i1, i2) - e(1, i1, i2+1) ) + &
              beta0 * beta(1, i1, i2+1) * b(2, i1, i2+1) - &
              beta0 * beta(1, i1, i2) * b(2, i1, i2) &
            ) + &
            alpha(1, i1+1) * lt(1, i1+1) * e(2, i1+1, i2) - &
            alpha(1, i1) * lt(1, i1) * e(2, i1, i2) + &
            lt(1, i1+1) * beta0 * beta(2, i1+1, i2) * b(1, i1+1, i2) - &
            lt(1, i1) * beta0 * beta(2, i1, i2) * b(1, i1, i2) &
          ) / ap(2, i1)


        b(1, i1, i2) = b(1, i1, i2) - dt * curl_er ! * this%l_abs(i1)
        b(2, i1, i2) = b(2, i1, i2) - dt * curl_et ! * this%l_abs(i1)
        b(3, i1, i2) = b(3, i1, i2) - dt * curl_ep ! * this%l_abs(i1)

      enddo
    enddo
    !$omp end parallel do

  elseif ( step == 2 ) then
    ! split advance loops, since curl_ep requires advanced br, bt

    ! advance br, bt
    !$omp parallel do private(i1)
    do i2 = 0, this%e%nx_(p_tc_dim) + 1
      do i1 = 0, this%b%nx_(p_rc_dim) + 1

        curl_er = alpha(1, i1) * ( &
            lp(1, i1, i2+1) * e(3, i1, i2+1) - &
            lp(1, i1, i2) * e(3, i1, i2) &
          ) / ar(2, i1, i2)

        curl_et = ( &
            lp(1, i1, i2) * alpha(1, i1) * e(3, i1, i2) - &
            lp(1, i1+1, i2) * alpha(1, i1+1) * e(3, i1+1, i2) &
          ) / at(1, i1, i2)


        b(1, i1, i2) = b(1, i1, i2) - dt * curl_er ! * this%l_abs(i1)
        b(2, i1, i2) = b(2, i1, i2) - dt * curl_et ! * this%l_abs(i1)

      enddo
    enddo
    !$omp end parallel do

    ! advance Bp
    !$omp parallel do private(i1)
    do i2 = 0, this%e%nx_(p_tc_dim) + 1
      do i1 = 0, this%e%nx_(p_rc_dim) + 1

        curl_ep = ( &
            lr(2, i1) * ( &
              alpha(2, i1) * ( e(1, i1, i2) - e(1, i1, i2+1) ) + &
              beta0 * beta(1, i1, i2+1) * b(2, i1, i2+1) - &
              beta0 * beta(1, i1, i2) * b(2, i1, i2) &
            ) + &
            alpha(1, i1+1) * lt(1, i1+1) * e(2, i1+1, i2) - &
            alpha(1, i1) * lt(1, i1) * e(2, i1, i2) + &
            lt(1, i1+1) * beta0 * beta(2, i1+1, i2) * b(1, i1+1, i2) - &
            lt(1, i1) * beta0 * beta(2, i1, i2) * b(1, i1, i2) &
          ) / ap(2, i1)

        b(3, i1, i2) = b(3, i1, i2) - dt * curl_ep ! * this%l_abs(i1)

      enddo
    enddo
    !$omp end parallel do

  endif

  ! clear pointers
  nullify(e)
  nullify(b)

  nullify(lr)
  nullify(lt)
  nullify(lp)

  nullify(ar)
  nullify(at)
  nullify(ap)

  nullify(alpha)
  nullify(beta)

end subroutine

!-------------------------------------------------------------------------------
! Advance the E field
!-------------------------------------------------------------------------------
subroutine dedt_gr( this, jay, dt )

  use m_current_define, only : t_current

  implicit none

  class( t_emf_gr ), intent( inout ) :: this
  class( t_current ), intent(in) :: jay
  real(p_double), intent(in) :: dt

  select case ( p_x_dim )
    case (2)
      call dedt_2d_gr( this, jay, dt )
  end select

end subroutine

!-------------------------------------------------------------------------------
! Advance the E field (2D)
!-------------------------------------------------------------------------------
subroutine dedt_2d_gr( this, jay, dt )

  use m_current_define, only : t_current

  implicit none

  class( t_emf_gr ), intent( inout ) :: this
  class( t_current ), intent(in) :: jay
  real(p_double), intent(in) :: dt

  ! local variables

  integer :: i1, i2
  real(p_k_fld) :: curl_br, curl_bt, curl_bp
  real(p_k_fld) :: beta0, ersh, ers, etsh, ets
  real( p_k_fld ), pointer :: e(:,:,:), b(:,:,:), j(:,:,:)
  real( p_k_fld ), pointer :: lr(:,:), lt(:,:), ap(:,:), alpha(:,:)
  real( p_k_fld ), pointer :: lp(:,:,:), ar(:,:,:), at(:,:,:), beta(:,:,:)

  ! executable statements
  ! setup pointers to improve readibility
  e => this%e%f2
  b => this%b%f2
  j => jay%pf(1)%f2

  lr => this%geometry%lr%f1
  lt => this%geometry%lt%f1
  lp => this%geometry%lp%f2

  ar => this%geometry%ar%f2
  at => this%geometry%at%f2
  ap => this%geometry%ap%f1

  alpha => this%geometry%alpha%f1
  beta => this%geometry%beta%f2
  beta0 = this%geometry%beta0

  !$omp parallel do private(i1)
  do i2 = 1, this%e%nx_(p_tc_dim) + 1
    do i1 = 1, this%e%nx_(p_rc_dim) + 1

      ! save values of E_r and E_t at t = t_n
      ersh = e(1, i1-1, i2)
      ers = e(1, i1, i2)
      etsh = e(2, i1, i2-1)
      ets = e(2, i1, i2)

      curl_br = alpha(2, i1) * ( &
          lp(2, i1, i2) * b(3, i1, i2) - &
          lp(2, i1, i2-1) * b(3, i1, i2-1) &
        ) / ar(1, i1, i2)

      curl_bt = ( &
          lp(2, i1-1, i2) * alpha(2, i1-1) * b(3, i1-1, i2) - &
          lp(2, i1, i2) * alpha(2, i1) * b(3, i1, i2) &
        ) / at(2, i1, i2)

      e(1, i1, i2) = e(1, i1, i2) + dt * ( curl_br - j(1, i1, i2) ) ! * this%l_abs(i1)
      e(2, i1, i2) = e(2, i1, i2) + dt * ( curl_bt - j(2, i1, i2) ) ! * this%l_abs(i1)

      ! get E_r and E_t at t = t_{n+1/2}
      ersh = 0.5 * ( ersh + e(1, i1-1, i2) )
      ers = 0.5 * ( ers + e(1, i1, i2) )
      etsh = 0.5 * ( etsh + e(2, i1, i2-1) )
      ets = 0.5 * ( ets + e(2, i1, i2) )

      curl_bp = ( &
          lr(1, i1) * ( &
            alpha(1, i1) * ( b(1, i1, i2-1) - b(1, i1, i2) ) + &
            beta0 * beta(2, i1, i2-1) * etsh - &
            beta0 * beta(2, i1, i2) * ets &
          ) + &
          lt(2, i1) * alpha(2, i1) * b(2, i1, i2) - &
          lt(2, i1-1) * alpha(2, i1-1) * b(2, i1-1, i2) + &
          lt(2, i1-1) * beta0 * beta(1, i1-1, i2) * ersh - &
          lt(2, i1) * beta0 * beta(1, i1, i2) * ers &
        ) / ap(1, i1)

      e(3, i1, i2) = e(3, i1, i2) + dt * ( curl_bp - j(3, i1, i2) ) ! * this%l_abs(i1)
    enddo
  enddo
  !$omp end parallel do

  ! clear pointers
  nullify(e)
  nullify(b)
  nullify(j)

  nullify(lr)
  nullify(lt)
  nullify(lp)

  nullify(ar)
  nullify(at)
  nullify(ap)

  nullify(alpha)
  nullify(beta)

end subroutine

end module m_emf_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Allocate any objects contained within the EMF (t_emf) object
!-------------------------------------------------------------------------------
subroutine allocate_objs_gr( this )

    use m_emf_define_gr
    use m_emf_bound_gr, only : t_emf_bound_gr
    use m_emf_gridval_gr, only : t_emf_gridval_gr

    implicit none

    class( t_emf_gr ), intent(inout) :: this

    ! Allocate t_emf_bound_gr object instead of t_emf_bound
    if ( .not. associated( this%bnd_con ) ) then
      allocate( t_emf_bound_gr :: this%bnd_con )
    endif

    ! Allocate t_emf_gridval_gr objects instead of t_emf_gridval
    if ( .not. associated( this%init_emf ) ) then
      allocate( t_emf_gridval_gr :: this%init_emf )
    endif
    if ( .not. associated( this%ext_emf ) ) then
      allocate( t_emf_gridval_gr :: this%ext_emf )
    endif

    ! Allocate superclass objects
    call this%t_emf%allocate_objs()

end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read input file data for gr algorithm
! - this routine also calls the superclass (t_emf) read_nml routine to read
! that section first
!-------------------------------------------------------------------------------
subroutine read_input_gr( this, input_file, periodic, if_move, grid, &
                                   dx, dt, gamma )

  use m_emf_define_gr
  use m_emf_gr
  use m_input_file, only : t_input_file
  use m_grid_define, only : t_grid

  implicit none

  class( t_emf_gr ), intent( inout ) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  real(p_double), dimension(:), intent(in) :: dx
  real(p_double), intent(in) :: dt
  real(p_double), intent(in) :: gamma

  ! -------------- Allocate any sub-objects needed --------------------
  call this%allocate_objs()

  ! read the superclass (emf) section
  call this%t_emf%read_input(input_file, periodic, if_move, grid, dx, dt, gamma)

end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Initialize t_emf_gr object
! - this routine also calls the superclass (t_emf) setup routine first to
! setup the electric magnetic fields
! - ...
!-------------------------------------------------------------------------------
subroutine init_gr( this, part_grid_center, part_interpolation, &
                     g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                     no_co, send_msg, recv_msg, restart, restart_handle, sim_options )

  use m_emf_define_gr
  use m_emf_gr
  use m_space, only : t_space
  use m_grid_define, only : t_grid
  use m_node_conf, only : t_node_conf
  use m_restart, only : t_restart_handle
  use m_vdf_comm, only : update_boundary, t_vdf_msg
  use m_emf_bound_gr, only : t_emf_bound_gr
  use m_emf_gridval_gr, only : t_emf_gridval_gr
  use m_time_step, only : t_time_step

  implicit none

  ! input parameters
  class( t_emf_gr ), intent( inout ), target :: this
  logical, intent(in) :: part_grid_center
  integer, intent(in) :: part_interpolation
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  integer, dimension(:,:), intent(in) :: gc_min
  real( p_double ), dimension(:), intent(in) :: dx
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax
  class( t_node_conf ), intent(in), target :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type( t_options ), intent(in) :: sim_options


  ! set pointer to geometry in boundary conditions object
  select type( bnd_con => this%bnd_con )
    class is ( t_emf_bound_gr )
      bnd_con % geometry => this % geometry
    class default
      write(err_buf__,*) 'Pointer to geometry must be set with a t_emf_bound_gr object';call err__("gr/os-emf-gr.f03",412)
      call abort_program( p_err_invalid )
  end select

  ! set pointer to geometry in init_emf object
  select type( init_emf => this%init_emf )
    class is ( t_emf_gridval_gr )
    init_emf % geometry => this % geometry
    class default
      write(err_buf__,*) 'Pointer to geometry must be set with a t_emf_gridval_gr object';call err__("gr/os-emf-gr.f03",421)
      call abort_program( p_err_invalid )
  end select

  ! set pointer to geometry in ext_emf object
  select type( ext_emf => this%ext_emf )
    class is ( t_emf_gridval_gr )
    ext_emf % geometry => this % geometry
    class default
      write(err_buf__,*) 'Pointer to geometry must be set with a t_emf_gridval_gr object';call err__("gr/os-emf-gr.f03",430)
      call abort_program( p_err_invalid )
  end select

  ! read any t_emf_gr quantity before calling superclass
  if ( restart ) then
    call this % restart_read( restart_handle )
  endif

  ! setup superclass data - this will also setup the diagnostics
  call this % t_emf % init( part_grid_center, part_interpolation, &
                            g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                            no_co, send_msg, recv_msg, &
                            restart, restart_handle, sim_options )

end subroutine

!-------------------------------------------------------------------------------
! Advance the EMF (t_emf_gr) object
!-------------------------------------------------------------------------------
subroutine advance_gr( this, jay, dt, send_msg, recv_msg )

  use m_emf_define_gr
  use m_emf_gr
  use m_current_define, only : t_current
  use m_vdf_comm, only : t_vdf_msg
  use m_emf_define, only : fsolverev
  use m_logprof

  implicit none

  class( t_emf_gr ), intent( inout ) :: this
  class( t_current ), intent(inout) :: jay
  real(p_double), intent(in) :: dt
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  !local variables
  real(p_double) :: dt_b, dt_e

  ! Advance the emf fields
  call begin_event( fsolverev )

  dt_b = dt / 2.0_p_double
  dt_e = dt

  ! Advance B half time step
  call dbdt_gr( this, dt_b, step = 1 )
  call this%bnd_con%update_boundary_b( this%e, this%b, step = 1 )

  ! Advance E one full time step
  call dedt_gr( this, jay, dt_e)
  call this%bnd_con%update_boundary_e( this%e, this%b )

  ! Advance B another half time step
  call dbdt_gr( this, dt_b, step = 2 )
  call this%bnd_con%update_boundary_b( this%e, this%b, step = 2 )

  call end_event( fsolverev )

end subroutine

!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_gr( this, restart_handle )
!-----------------------------------------------------------------------------------------
! write object information into a restart file
!-----------------------------------------------------------------------------------------

  use m_emf_gr
  use m_restart
  use m_parameters

  implicit none

  class( t_emf_gr ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for t_emf_gr object.'
  integer :: ierr

  call restart_io_write("p_emf_gr_rst_id", p_emf_gr_rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"gr/os-emf-gr.f03",510)

  ! write superclass checkpoint data first
  call this % t_emf % write_checkpoint( restart_handle )

end subroutine write_checkpoint_gr
!-----------------------------------------------------------------------------------------

subroutine restart_read_gr( this, restart_handle )

  use m_emf_gr
  use m_restart
  use m_parameters

  implicit none

  class( t_emf_gr ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for emf_gr object.'
  character(len=len(p_emf_gr_rst_id)) :: rst_id
  integer :: ierr

  call restart_io_read(rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"gr/os-emf-gr.f03",534)

  !do nothing
  !restart read is done in the init function: init_gr


end subroutine restart_read_gr
