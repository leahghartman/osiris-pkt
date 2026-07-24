# 1 "gr/os-species-define-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-species-define-gr.f03" 2
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
# 2 "gr/os-species-define-gr.f03" 2
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
# 3 "gr/os-species-define-gr.f03" 2

module m_species_define_gr

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
# 7 "gr/os-species-define-gr.f03" 2

use m_system
use m_parameters
use m_input_file, only : t_input_file
use stringutil
use m_node_conf, only : t_node_conf
use m_vdf_comm, only : t_vdf_msg
use m_species_comm, only : t_spec_msg
use m_vdf_define, only : t_vdf
use m_grid_define, only : t_grid
use m_space, only : t_space
use m_time_step, only : t_time_step
use m_emf_define, only : t_emf
use m_current_define, only : t_current
use m_species_define, only : t_species, p_spec_rst_id, t_part_idx

use m_geometry_gr
use m_restart

implicit none

private

! string to id restart data
character(len=*), parameter, public :: p_spec_gr_rst_id = "species_gr rst data - 0x0006"

! type of pusher
integer, parameter :: p_rk4 = 44
integer, parameter :: p_rk6 = 66


type, extends( t_species ) :: t_species_gr

  class( t_geometry_gr ), pointer :: geometry

  ! cooling
  integer :: cooling_type

  ! pair production
  real(p_double) :: pp_rmin, pp_rmax, pp_gth, pp_gpair

  ! pointer to receive produced reciprocal particles
  class( t_species_gr ), pointer :: grp_recip

  ! pair/photon emission diagnostic
  logical :: track_nemit
  type( t_vdf ) :: nemit

  ! injection
  type( t_vdf ) :: charge

  contains

  procedure :: get_n_x_dims => get_n_x_dims_species_gr
  procedure :: init_buffer => init_buffer_species_gr
  procedure :: read_input => read_input_species_gr
  procedure :: init => init_species_gr

  procedure :: push => push_species_gr
  procedure :: emit => emit_species_gr

  procedure :: deposit_density => deposit_density_species_gr

  procedure :: position_single_4 => position_single_4_gr
  procedure :: position_single_8 => position_single_8_gr
  procedure :: position_idx_comp_4 => position_idx_comp_4_gr
  procedure :: position_idx_comp_8 => position_idx_comp_8_gr
  procedure :: position_ref_box_idx_comp_4 => position_ref_box_idx_comp_4_gr
  procedure :: position_ref_box_idx_comp_8 => position_ref_box_idx_comp_8_gr

  procedure :: restart_write => restart_write_gr
  procedure :: restart_read => restart_read_gr

end type t_species_gr

public t_species_gr
public p_rk4, p_rk6

interface
  subroutine read_input_species_gr( this, input_file, def_name, &
    periodic, if_move, grid, dt, read_prof, sim_options )
  import t_species_gr, t_input_file, t_options, p_double, t_grid
  class( t_species_gr ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  character(len = *), intent(in) :: def_name
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  real(p_double), intent(in) :: dt
  logical, intent(in) :: read_prof
  type(t_options), intent(in) :: sim_options
  end subroutine
end interface

interface
  subroutine init_species_gr( this, sp_id, interpolation, grid_center, &
    grid, g_space, emf, jay, no_co, send_vdf, recv_vdf, bnd_cross, node_cross, &
    send_spec, recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, &
    tmax )

  import t_species_gr, t_emf, t_current, t_space, t_grid, t_node_conf, t_restart_handle, &
         t_options, t_vdf_msg, t_part_idx, t_spec_msg, p_double, t_time_step

  class( t_species_gr ), intent(inout) :: this

  integer, intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in) :: grid
  type( t_space ), intent(in) :: g_space
  class( t_emf ), intent(inout) :: emf
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_vdf, recv_vdf
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_spec, recv_spec
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type(t_options), intent(in) :: sim_options
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax

  end subroutine
end interface

interface
  subroutine push_species_gr( this, emf, current, t, tstep, tid, n_threads, options )

  import t_species_gr, t_emf, t_vdf, p_double, t_time_step, t_options

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( inout ) :: emf
  type( t_vdf ), intent(inout) :: current

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  integer, intent(in) :: tid ! local thread id
  integer, intent(in) :: n_threads ! total number of threads
  type( t_options ), intent(in) :: options

  end subroutine
end interface

interface
subroutine emit_species_gr( this, emf, dt, ptrcur, np )

  import t_species_gr, t_emf, p_double

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  real(p_double), intent(in) :: dt
  integer, intent(in) :: ptrcur, np

end subroutine emit_species_gr
end interface

interface
  subroutine deposit_density_species_gr( this, charge, i1, i2, q )

  import t_species_gr, t_vdf, p_k_part

  class( t_species_gr ), intent(in) :: this
  type( t_vdf ), intent(inout) :: charge
  integer, intent(in) :: i1, i2
  real(p_k_part), dimension(:), intent(in) :: q

  end subroutine
end interface

contains

!-------------------------------------------------------------------------------
! Return the number of spatial dimensions on the x buffer
!-------------------------------------------------------------------------------
pure function get_n_x_dims_species_gr( this )

  implicit none

  class( t_species_gr ), intent(in) :: this
  integer :: get_n_x_dims_species_gr

  select case ( this%geometry%metric )
  case ( p_geometry_minkowski )
    get_n_x_dims_species_gr = p_x_dim+3
  case ( p_geometry_schwarzschild )
    get_n_x_dims_species_gr = p_x_dim+1
  case ( p_geometry_kerr_slow )
    get_n_x_dims_species_gr = p_x_dim+1
  case ( p_geometry_kerr )
    get_n_x_dims_species_gr = p_x_dim+1
  end select

end function get_n_x_dims_species_gr
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Allocate buffers for particle quantities
!-------------------------------------------------------------------------------
subroutine init_buffer_species_gr( this, num_par_req )

  implicit none

  class(t_species_gr), intent(inout) :: this
  integer, intent(in) :: num_par_req

  ! The buffer size must always be a multiple of vector width in size because of SIMD code
  this%num_par_max = (( num_par_req + p_vecwidth - 1 ) / p_vecwidth) * p_vecwidth

  ! Under some compilers / configurations (e.g. gfortran 4.9.1, OS X, single precision)
  ! if the total q buffer size is below 1Kb it may not be allocated to a 32bit boundary
  ! which is required by AVX code
  if ( this%num_par_max < 256 ) this%num_par_max = 256

  ! setup position buffer
  call freemem(this%x,"gr/os-species-define-gr.f03",223)
  select case ( this%geometry%metric )
  case ( p_geometry_minkowski )
    call alloc(this%x, (/ p_x_dim + 3, this%num_par_max /),"gr/os-species-define-gr.f03",226)
  case ( p_geometry_schwarzschild )
    call alloc(this%x, (/ p_x_dim + 1, this%num_par_max /),"gr/os-species-define-gr.f03",228)
  case ( p_geometry_kerr_slow )
    call alloc(this%x, (/ p_x_dim + 1, this%num_par_max /),"gr/os-species-define-gr.f03",230)
  case ( p_geometry_kerr )
    call alloc(this%x, (/ p_x_dim + 1, this%num_par_max /),"gr/os-species-define-gr.f03",232)
  case default
    write(err_buf__,*) 'Invalid geometry in init_buffer_species_gr';call err__("gr/os-species-define-gr.f03",234)
    call abort_program( p_err_invalid )
  end select

  ! initialize particle cell information
  call freemem(this%ix,"gr/os-species-define-gr.f03",239)
  call alloc(this%ix, (/ p_x_dim, this%num_par_max /),"gr/os-species-define-gr.f03",240)

  ! setup momenta buffer
  call freemem(this%p,"gr/os-species-define-gr.f03",243)
  call alloc(this%p, (/ p_p_dim, this%num_par_max /),"gr/os-species-define-gr.f03",244)

  ! setup particle charge buffer
  call freemem(this%q,"gr/os-species-define-gr.f03",247)
  call alloc(this%q, (/ this%num_par_max /),"gr/os-species-define-gr.f03",248)

  ! initialize tracking data if necessary
  if ( this%add_tag ) then
     call freemem(this%tag,"gr/os-species-define-gr.f03",252)
     call alloc(this%tag, (/ 2, this%num_par_max /),"gr/os-species-define-gr.f03",253)
  endif

end subroutine init_buffer_species_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Get particle position in single precision
!-------------------------------------------------------------------------------
subroutine position_idx_comp_4_gr( this, comp, idx, np, pos )

  implicit none

  class( t_species_gr ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_single ), dimension(:), intent(out) :: pos

  integer :: j, ixmin, ix
  real( p_single ) :: xmin, dx

  if (comp <= p_x_dim) then
    select case ( comp )
    case( p_rc_dim )
      do j = 1, np
        ix = this%ix( p_rc_dim, idx(j) )
        pos(j) = real( this%geometry%r%f1(2, ix) + &
          this%geometry%dr%f1(1, ix) * &
          this%x ( comp, idx(j) ) , p_single )

      enddo

    case( p_tc_dim )
      xmin = real( this%g_box( p_lower, comp ), p_single )
      ixmin = this%my_nx_p( 1, comp ) - 2
      dx = real( this%dx( comp ), p_single )

      do j = 1, np
        pos(j) = xmin + dx * real( (this%ix( comp, idx(j) ) + ixmin) + &
          this%x ( comp, idx(j) ) , p_single )
      enddo
    end select
  else
    ! extra coordinates do not have a corresponding ix
    do j = 1, np
      pos(j) = real(this%x(comp, idx(j)), p_single)
    enddo
  endif

end subroutine position_idx_comp_4_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Get particle position in double precision
!-------------------------------------------------------------------------------
subroutine position_idx_comp_8_gr( this, comp, idx, np, pos )
!-------------------------------------------------------------------------------

  implicit none

  class( t_species_gr ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_double ), dimension(:), intent(out) :: pos

  integer :: j, ixmin, ix
  real( p_double ) :: xmin, dx

  if (comp <= p_x_dim) then
    select case ( comp )
    case( p_rc_dim )
      do j = 1, np
        ix = this%ix( p_rc_dim, idx(j) )
        pos(j) = real( this%geometry%r%f1(2, ix) + &
          this%geometry%dr%f1(1, ix) * &
          this%x ( comp, idx(j) ) , p_double )

      enddo

    case( p_tc_dim )
      xmin = real( this%g_box( p_lower, comp ), p_double )
      ixmin = this%my_nx_p( 1, comp ) - 2
      dx = real( this%dx( comp ), p_double )

      do j = 1, np
        pos(j) = xmin + dx * real( (this%ix( comp, idx(j) ) + ixmin) + &
          this%x ( comp, idx(j) ) , p_double )
      enddo
    end select
  else
    ! extra coordinates do not have a corresponding ix
    do j = 1, np
      pos(j) = real(this%x(comp, idx(j)), p_double)
    enddo
  endif

end subroutine position_idx_comp_8_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine position_single_4_gr( this, idx, pos )

  implicit none

  class( t_species_gr ), intent(in) :: this
  integer, intent(in) :: idx
  real( p_single ), dimension(:), intent(out) :: pos

  integer :: j, ixmin, ix
  real( p_double ) :: xmin, dx

  do j = 1, size(this%x,1)
    if (j <= p_x_dim) then
      select case ( j )
      case( p_rc_dim )

        ix = this%ix( p_rc_dim, idx )
        pos(j) = real( this%geometry%r%f1(2, ix) + &
          this%geometry%dr%f1(1, ix) * &
          this%x ( j, idx ) , p_single )

      case( p_tc_dim )
        xmin = real( this%g_box( p_lower, j ), p_double )
        ixmin = this%my_nx_p( 1, j ) - 2
        dx = real( this%dx( j ), p_double )


        pos(j) = real( xmin + dx * ( (this%ix( j, idx ) + ixmin) + &
          this%x ( j, idx ) ) , p_single )

      end select
    else
      ! extra coordinates do not have a corresponding ix
      pos(j) = real(this%x(j, idx), p_single)
    endif
  enddo

end subroutine position_single_4_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine position_single_8_gr( this, idx, pos )

  implicit none

  class( t_species_gr ), intent(in) :: this
  integer, intent(in) :: idx
  real( p_double ), dimension(:), intent(out) :: pos

  integer :: j, ixmin, ix
  real( p_double ) :: xmin, dx

  do j = 1, size(this%x,1)
    if (j <= p_x_dim) then
      select case ( j )
      case( p_rc_dim )

        ix = this%ix( p_rc_dim, idx )
        pos(j) = real( this%geometry%r%f1(2, ix) + &
          this%geometry%dr%f1(1, ix) * &
          this%x ( j, idx ) , p_double )

      case( p_tc_dim )
        xmin = real( this%g_box( p_lower, j ), p_double )
        ixmin = this%my_nx_p( 1, j ) - 2
        dx = real( this%dx( j ), p_double )


        pos(j) = xmin + dx * real( (this%ix( j, idx ) + ixmin) + &
          this%x ( j, idx ) , p_double )

      end select
    else
      ! extra coordinates do not have a corresponding ix
      pos(j) = real(this%x(j, idx), p_double)
    endif
  enddo

end subroutine position_single_8_gr
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_ref_box_idx_comp_4_gr( this, comp, idx, np, pos, if_pos_ref_box )

  implicit none

  class( t_species_gr ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_single ), dimension(:), intent(out) :: pos
  logical, intent(in) :: if_pos_ref_box

  integer :: j, ixmin, ix
  real( p_double ) :: xmin, dx

  ! we ignore if_pos_ref_box here
  ! this function does the same as position_idx_comp_4_gr

  if (comp <= p_x_dim) then
    select case ( comp )
    case( p_rc_dim )
      do j = 1, np
        ix = this%ix( p_rc_dim, idx(j) )
        pos(j) = real( this%geometry%r%f1(2, ix) + &
          this%geometry%dr%f1(1, ix) * &
          this%x ( comp, idx(j) ) , p_single )

      enddo

    case( p_tc_dim )
      xmin = real( this%g_box( p_lower, comp ), p_double )
      ixmin = this%my_nx_p( 1, comp ) - 2
      dx = real( this%dx( comp ), p_double )

      do j = 1, np
        pos(j) = real( xmin + dx * ( (this%ix( comp, idx(j) ) + ixmin) + &
          this%x ( comp, idx(j) ) ) , p_single )
      enddo
    end select
  else
    ! extra coordinates do not have a corresponding ix
    do j = 1, np
      pos(j) = real(this%x(comp, idx(j)), p_single)
    enddo
  endif

end subroutine position_ref_box_idx_comp_4_gr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_ref_box_idx_comp_8_gr( this, comp, idx, np, pos, if_pos_ref_box )

  implicit none

  class( t_species_gr ), intent(in) :: this
  integer, intent(in) :: comp
  integer, intent(in), dimension(:) :: idx
  integer, intent(in) :: np
  real( p_double ), dimension(:), intent(out) :: pos
  logical, intent(in) :: if_pos_ref_box

  integer :: j, ixmin, ix
  real( p_double ) :: xmin, dx

  ! we ignore if_pos_ref_box here
  ! this function does the same as position_idx_comp_8_gr

  if (comp <= p_x_dim) then
    select case ( comp )
    case( p_rc_dim )
      do j = 1, np
        ix = this%ix( p_rc_dim, idx(j) )
        pos(j) = real( this%geometry%r%f1(2, ix) + &
          this%geometry%dr%f1(1, ix) * &
          this%x ( comp, idx(j) ) , p_double )

      enddo

    case( p_tc_dim )
      xmin = real( this%g_box( p_lower, comp ), p_double )
      ixmin = this%my_nx_p( 1, comp ) - 2
      dx = real( this%dx( comp ), p_double )

      do j = 1, np
        pos(j) = xmin + dx * real( (this%ix( comp, idx(j) ) + ixmin) + &
          this%x ( comp, idx(j) ) , p_double )
      enddo
    end select
  else
    ! extra coordinates do not have a corresponding ix
    do j = 1, np
      pos(j) = real(this%x(comp, idx(j)), p_double)
    enddo
  endif

end subroutine position_ref_box_idx_comp_8_gr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_write_gr( this, restart_handle )
!-----------------------------------------------------------------------------------------
! write object information into a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species_gr ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for species_gr object.'
  integer :: ierr

  ! write superclass restart data first
  call this % t_species % restart_write( restart_handle )

  ! write t_species_gr rst_id check
  call restart_io_write("p_spec_gr_rst_id", p_spec_gr_rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"gr/os-species-define-gr.f03",553)

  ! write t_species_gr charge restart data
  call this % charge % write_checkpoint( restart_handle )


end subroutine restart_write_gr
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine restart_read_gr( this, restart_handle )
!-----------------------------------------------------------------------------------------
! read object information from a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species_gr ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for species_gr object.'
  character(len=len(p_spec_gr_rst_id)) :: rst_id
  integer :: ierr

  ! write superclass restart data first
  call this % t_species % restart_read( restart_handle )

  ! new t_species_gr additions
  ! t_species_gr rst_id check
  call restart_io_read(rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"gr/os-species-define-gr.f03",584)

  ! check if restart file is compatible
  if ( rst_id /= p_spec_gr_rst_id ) then
    write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("gr/os-species-define-gr.f03",588)
    write(err_buf__,*) 'from incompatible binary (species_gr)';call err__("gr/os-species-define-gr.f03",589)
    call abort_program(p_err_rstrd)
  endif

  ! read t_species_gr charge restart data
  call this % charge % read_checkpoint( restart_handle )

end subroutine restart_read_gr
!-----------------------------------------------------------------------------------------


end module m_species_define_gr
