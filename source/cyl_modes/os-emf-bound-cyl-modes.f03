#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_bound_cyl_modes

use m_system
use m_parameters
use m_cyl_modes,   only: t_cyl_modes, p_max_cyl_modes, p_real
use m_emf_define,  only: t_emf_bound, t_lindman, t_vpml
use m_vdf_define,  only: t_vdf
use m_vdf_comm,    only: irecv, isend, irecv_wait, isend_wait, t_vdf_msg, update_periodic_vdf
use m_emf_bound,   only: emfboundev
use m_space,       only: t_space, x_dim, nx_move, if_move
use m_node_conf,   only: t_node_conf, neighbor, periodic
use m_grid_define, only: t_grid
use m_restart,     only: t_restart_handle, restart_io_write, restart_io_read
use m_vpml,        only: restart_write, restart_read, is_active, update_boundary, setup, cleanup
use m_vpml,        only: update_interface_x1_2d, update_e_x1_2d, update_b_x1_2d
use m_vpml,        only: update_interface_x2_2d, update_e_x2_2d, update_b_x2_2d
use m_vpml,        only: update_b_2d, update_e_2d, move_window, reshape_obj
use m_lindman,     only: restart_write, restart_read, is_active, update_boundary, setup, cleanup
use m_lindman,     only: update_b_2d, move_window, update_e_2d, update_emfbound, reshape_obj
use m_logprof,     only: create_event, begin_event, end_event
use m_wall_define, only: range
use m_emf_pec
use m_emf_pmc

implicit none

private

! string to id restart data
character(len=*), parameter :: p_bcemf_cyl_modes_rst_id = "bcemf_cyl_modes rst data - 0x0001"

! special boundary condition type for quasi-3D geometry
type, extends( t_emf_bound ), public :: t_emf_bound_cyl_modes

  type (t_lindman), dimension(2,p_max_cyl_modes) :: lindman_re
  type (t_lindman), dimension(2,p_max_cyl_modes) :: lindman_im

  type (t_vpml) :: vpml_all_re(p_max_cyl_modes)
  type (t_vpml) :: vpml_all_im(p_max_cyl_modes)

contains

  procedure :: init                 => init_emf_bound_cyl_modes
  procedure :: cleanup              => cleanup_emf_bound_cyl_modes
  procedure :: write_checkpoint     => write_checkpoint_emf_bound_cyl_modes
  procedure :: restart_read         => restart_read_emf_bound_cyl_modes
  procedure :: update_boundary_cm   => update_boundary_emf_bound_cyl_modes
  procedure :: update_boundary_b_cm => update_boundary_bfld_emf_bound_cyl_modes
  procedure :: update_boundary_e_cm => update_boundary_efld_emf_bound_cyl_modes
  procedure :: move_window          => move_window_emf_bound_cyl_modes
  procedure :: reshape              => reshape_emf_bound_cyl_modes

end type t_emf_bound_cyl_modes

!-----------------------------------------------------------------------------------------
contains

!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_emf_bound_cyl_modes( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  class (t_emf_bound_cyl_modes), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for emf_bound_cyl_modes object.'
  integer :: i, bnd, ierr, i_mode

  restart_io_wr( p_bcemf_cyl_modes_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%global_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! write lindman boundary data if needed
  do i = 1, p_x_dim
    do bnd = p_lower, p_upper
      if ( this%type( bnd, i ) == p_bc_lindman ) then
        call restart_write( this%lindman(bnd), restart_handle )
        do i_mode = 1, p_max_cyl_modes
          call restart_write( this%lindman_re(bnd,i_mode), restart_handle )
          call restart_write( this%lindman_im(bnd,i_mode), restart_handle )
        enddo
      endif
    enddo
  enddo

#ifdef __HAS_PML__

  ! write vpml boundary data if needed
  if(is_active(this%vpml_all)) then
    call restart_write( this%vpml_all, restart_handle )
    do i_mode = 1, p_max_cyl_modes
      call restart_write( this%vpml_all_re(i_mode), restart_handle )
      call restart_write( this%vpml_all_im(i_mode), restart_handle )
    enddo
  endif

#endif

end subroutine write_checkpoint_emf_bound_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_emf_bound_cyl_modes( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  class (t_emf_bound_cyl_modes), intent(inout)  ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for emf_bound_cyl_modes object.'
  character(len=len(p_bcemf_cyl_modes_rst_id)) :: rst_id
  integer :: ierr

  ! check if restart file is compatible
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  if ( rst_id /= p_bcemf_cyl_modes_rst_id ) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (emf_bound_cyl_modes)')
    ERROR('rst_id = "', rst_id, '"' )
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%global_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

end subroutine restart_read_emf_bound_cyl_modes
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine update_boundary_emf_bound_cyl_modes( this, b, e, g_space, no_co, i_mode, i_comp, send_msg, recv_msg )
!---------------------------------------------------
!       update boundary of electromagnetic field
!---------------------------------------------------

  implicit none

  class( t_emf_bound_cyl_modes ), intent(inout) :: this
  type( t_vdf ),       intent(inout) :: b, e
  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: i_mode, i_comp
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, my_node
  integer, dimension(p_max_dim) :: lnx_move
  integer :: lneighbor, uneighbor

  ! TODO: Make sure lindman and vpml are considered properly here

  call begin_event(emfboundev)

  lnx_move( 1 : x_dim( g_space ) ) = nx_move( g_space )

  ! update to the lindman boundaries must occur before updating em fields to
  ! get the corner value right

  if(i_comp .eq. p_real) then
    call update_boundary(this%lindman_re(p_lower,i_mode), no_co, lnx_move, send_msg, recv_msg )
    call update_boundary(this%lindman_re(p_upper,i_mode), no_co, lnx_move, send_msg, recv_msg )
  else
    call update_boundary(this%lindman_im(p_lower,i_mode), no_co, lnx_move, send_msg, recv_msg )
    call update_boundary(this%lindman_im(p_upper,i_mode), no_co, lnx_move, send_msg, recv_msg )
  end if

  ! Update VPML walls if necessary

  if (i_comp == p_real) then
    if (is_active(this%vpml_all_re(i_mode))) then
    call update_boundary(this%vpml_all_re(i_mode), no_co, lnx_move, send_msg, recv_msg )
    end if
  else
    if (is_active(this%vpml_all_im(i_mode))) then
    call update_boundary(this%vpml_all_im(i_mode), no_co, lnx_move, send_msg, recv_msg )
    end if
  end if

  ! Update other boundary types

  my_node = no_co%my_aid()

  do i = 1, e%x_dim_

    if ( my_node == neighbor(no_co,p_lower,i) ) then

      !single node periodic
      call update_periodic_vdf( e, i, p_vdf_replace )
      call update_periodic_vdf( b, i, p_vdf_replace )

    else

      ! communication with other nodes
      lneighbor = neighbor( no_co, p_lower, i )
      uneighbor = neighbor( no_co, p_upper, i )

      ! post receives
      if ( lneighbor > 0 ) call irecv( e, p_vdf_replace, i, p_lower, no_co, lnx_move(i), recv_msg, b )
      if ( uneighbor > 0 ) call irecv( e, p_vdf_replace, i, p_upper, no_co, lnx_move(i), recv_msg, b )

      ! post sends
      ! Note: sends MUST be posted in the opposite order of receives to account for a 2 node
      ! periodic partition, where 1 node sends 2 messages to the same node
      if ( uneighbor > 0 ) call isend( e, p_vdf_replace, i, p_upper, no_co, lnx_move(i), send_msg, b )
      if ( lneighbor > 0 ) call isend( e, p_vdf_replace, i, p_lower, no_co, lnx_move(i), send_msg, b )

      ! Process lindman boundaries if needed
      if (this%type( p_lower, i ) == p_bc_lindman) then
        if(i_comp .eq. p_real) then
          call update_emfbound( this%lindman_re( p_lower, i_mode ), b, e, i)
        else
          call update_emfbound( this%lindman_im( p_lower, i_mode ), b, e, i)
        end if
      endif
      if (this%type( p_upper, i ) == p_bc_lindman) then
        if(i_comp .eq. p_real) then
          call update_emfbound( this%lindman_re( p_upper, i_mode ), b, e, i)
        else
          call update_emfbound( this%lindman_im( p_upper, i_mode ), b, e, i)
        end if
      endif

      ! wait for receives and unpack data
      if ( lneighbor > 0 ) call irecv_wait( e, p_lower, recv_msg, b )
      if ( uneighbor > 0 ) call irecv_wait( e, p_upper, recv_msg, b )

      ! wait for sends to complete
      if ( uneighbor > 0 ) call isend_wait( e, p_upper, send_msg )
      if ( lneighbor > 0 ) call isend_wait( e, p_lower, send_msg )

    endif

  enddo

  call end_event(emfboundev)

end subroutine update_boundary_emf_bound_cyl_modes
!---------------------------------------------------

!---------------------------------------------------
subroutine init_emf_bound_cyl_modes( this, b, e, dt, part_interpolation, &
                            space, no_co, grid, restart, restart_handle )
!---------------------------------------------------
! sets up this data structure from the given information
!---------------------------------------------------

  implicit none

  class (t_emf_bound_cyl_modes), intent( inout )  ::  this
  type( t_vdf ), intent(in) :: b, e
  real( p_double ), intent(in) :: dt
  integer, intent(in) :: part_interpolation
  type( t_space ), intent(in) :: space
  class( t_node_conf ), intent(in) :: no_co
  class( t_grid ), intent(in) :: grid

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  logical, dimension(p_x_dim) :: ifpr_l, if_move_l
  logical, dimension(2,p_x_dim) :: if_vpml
  integer :: i_bnd, i_x_dim, previous, n_vpml, i_mode

  if ( restart ) then

    call this%restart_read( restart_handle )

  else

    ifpr_l = no_co%ifpr()
    if_move_l = if_move(space)

    ! overwrite b.c. flags if motion of simulation space
    do i_x_dim = 1, p_x_dim
      if ( if_move_l(i_x_dim) ) then
        this%type( p_lower , i_x_dim ) = p_bc_move_c
        this%type( p_upper , i_x_dim ) = p_bc_move_c
      endif
    enddo

    ! Store the global type of boundary condition used
    this%global_type = this%type
    do i_x_dim = 1, p_x_dim
      if ( periodic( no_co, i_x_dim) ) this%global_type( 1:2 , i_x_dim ) = p_bc_periodic
    enddo

    ! overwrite b.c.-flags if single node periodicity
    do i_x_dim = 1, p_x_dim
      if ( ifpr_l(i_x_dim) ) then
        this%type( : , i_x_dim ) = p_bc_periodic
      endif
    enddo

    ! overwrite b.c. flags for boundaries to other nodes
    do i_x_dim = 1, p_x_dim
      do i_bnd=1, 2
        if ( neighbor( no_co, i_bnd, i_x_dim ) > 0 ) then
          this%type( i_bnd, i_x_dim ) = p_bc_other_node
        endif
      enddo
    enddo

    ! do additional cylindrical coordinates setup
    if ( grid % coordinates == p_cylindrical_b ) then

      ! If node has the axial boundary set the bc accordingly
      if ( no_co%on_edge( p_r_dim, p_lower )) then
        this%type(p_lower, p_r_dim) = p_bc_cyl_axis
      endif

      ! check if we have periodic boundaries along radial direction
      if ((this%type(p_lower, p_r_dim) == p_bc_periodic) .or. &
          (this%type(p_upper, p_r_dim) == p_bc_periodic)) then
        ERROR('no periodic b.c. for radial direction')
        call abort_program( p_err_invalid )
      endif

    endif

    ! When using PEC or PMC boundaries adjust to the correct one according to the

    ! particle interpolation
    do i_x_dim = 1, p_x_dim
      do i_bnd=1, 2
        select case ( this%type(i_bnd,i_x_dim) )
        case ( p_bc_pec )
          if ( mod( part_interpolation, 2 ) == 1 ) then
            this%type(i_bnd,i_x_dim) = p_bc_pec_odd
          else
            this%type(i_bnd,i_x_dim) = p_bc_pec_even
          endif
        case ( p_bc_pmc )
          if ( mod( part_interpolation, 2 ) == 1 ) then
            this%type(i_bnd,i_x_dim) = p_bc_pmc_odd
          else
            this%type(i_bnd,i_x_dim) = p_bc_pmc_even
          endif
        case default
          continue
        end select
      enddo
    enddo

  endif

  ! setup lindman boundaries
  previous = -1

  do i_x_dim=1, p_x_dim
    do i_bnd = p_lower, p_upper
      if ( this%type( i_bnd, i_x_dim ) == p_bc_lindman ) then
        ! check if lindman boundaries were already defined in another direction
        ! (this is redundant because it was already checked when reading the input file)
        if ( previous > 0 .and. previous /= i_x_dim ) then
          ERROR("Lindman boundaries are only allowed in 1 direction")
          call abort_program( p_err_invalid )
        endif

        ! setup lindman boundary object
        call setup( this%lindman(i_bnd), i_bnd, i_x_dim, b, e, dt, restart, restart_handle )

        do i_mode=1,p_max_cyl_modes
          call setup( this%lindman_re(i_bnd,i_mode), i_bnd, i_x_dim, b, e, dt, restart, restart_handle )
          call setup( this%lindman_im(i_bnd,i_mode), i_bnd, i_x_dim, b, e, dt, restart, restart_handle )
        enddo

        previous = i_x_dim
      endif
    enddo

  enddo

  ! setup VPML boundaries
  n_vpml = 0
  if_vpml = .false.
  do i_x_dim=1, p_x_dim
    if ( this%type( p_lower, i_x_dim ) == p_bc_vpml ) then
      if_vpml(p_lower,i_x_dim) = .true.
      n_vpml = n_vpml + 1
    endif

    if ( this%type( p_upper, i_x_dim ) == p_bc_vpml ) then
      if_vpml(p_upper,i_x_dim) = .true.
      n_vpml = n_vpml + 1
    endif
  enddo

  ! setup if at least one VPML
  if (n_vpml > 0 ) then

    call setup( this%vpml_all, e, b, grid, n_vpml, this%vpml_bnd_size, this%vpml_diffuse, &
      if_vpml, if_move_l, dt, restart, restart_handle )

    do i_mode = 1,p_max_cyl_modes
      call setup( this%vpml_all_re(i_mode), e, b, grid, n_vpml, &
                  this%vpml_bnd_size, this%vpml_diffuse, &
                  if_vpml, if_move_l, dt, restart, restart_handle )

      call setup( this%vpml_all_im(i_mode), e, b, grid, n_vpml, &
                  this%vpml_bnd_size, this%vpml_diffuse, &
                  if_vpml, if_move_l, dt, restart, restart_handle )
    enddo

  endif

  ! setup events
  if (emfboundev==0) emfboundev = create_event('update emf boundary')

end subroutine init_emf_bound_cyl_modes
!---------------------------------------------------

subroutine cleanup_emf_bound_cyl_modes( this )
!---------------------------------------------------
!       release memory used by ptr-components of this variable
!---------------------------------------------------

  implicit none

  class(t_emf_bound_cyl_modes), intent( inout ) :: this

  integer i

  ! call superclass cleanup function
  call this % t_emf_bound % cleanup( )

  ! TODO CHECK IF LINDMAN initialized properly
  do i=1,p_max_cyl_modes
    call cleanup( this%lindman_re(p_lower,i))
    call cleanup( this%lindman_re(p_upper,i))
    call cleanup( this%lindman_im(p_lower,i))
    call cleanup( this%lindman_im(p_upper,i))
  enddo

  do i=1,p_max_cyl_modes
    call cleanup(this%vpml_all_re(i))
    call cleanup(this%vpml_all_im(i))
  enddo

end subroutine cleanup_emf_bound_cyl_modes


!-----------------------------------------------------------------------------------------
! Update physical emf boundaries after each B field advance. Used for:
!  - Open boundaries (Lindman or PML)
!  - Conducting (to be removed)
!  - Perfect electric conductor (PEC)
!-----------------------------------------------------------------------------------------
subroutine update_boundary_bfld_emf_bound_cyl_modes( this, e_cyl_m, b_cyl_m, step )

  implicit none

  class( t_emf_bound_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: e_cyl_m, b_cyl_m
  integer, intent(in) :: step

  integer :: bnd, i_dim, mode, n_modes

  n_modes = ubound( e_cyl_m%pf_re, 1 )

  do i_dim = 1, 2
    do bnd = p_lower, p_upper
      select case (this%type( bnd, i_dim ))
      case (p_bc_lindman)
        call update_b_2d(this%lindman(bnd), e_cyl_m%pf_re(0))
        do mode = 1, n_modes
          call update_b_2d(this%lindman_re(bnd,mode), e_cyl_m%pf_re(mode))
          call update_b_2d(this%lindman_im(bnd,mode), e_cyl_m%pf_im(mode))
        enddo

      ! I'll have to investigate how these relate to cylindrical modal geometries
      case ( p_bc_pec_odd )
        call pec_bc_b_2d_odd( b_cyl_m%pf_re(0), i_dim, bnd )
        do mode = 1, n_modes
          call pec_bc_b_2d_odd( b_cyl_m%pf_re(mode), i_dim, bnd )
          call pec_bc_b_2d_odd( b_cyl_m%pf_im(mode), i_dim, bnd )
        enddo

      case ( p_bc_pmc_odd )
        call pmc_bc_b_2d_odd( b_cyl_m%pf_re(0), i_dim, bnd )
        do mode = 1, n_modes
          call pmc_bc_b_2d_odd( b_cyl_m%pf_re(mode), i_dim, bnd )
          call pmc_bc_b_2d_odd( b_cyl_m%pf_im(mode), i_dim, bnd )
        enddo

      case ( p_bc_pec_even )
        call pec_bc_b_2d_even( b_cyl_m%pf_re(0), i_dim, bnd )
        do mode = 1, n_modes
          call pec_bc_b_2d_even( b_cyl_m%pf_re(mode), i_dim, bnd )
          call pec_bc_b_2d_even( b_cyl_m%pf_im(mode), i_dim, bnd )
        enddo

      case ( p_bc_pmc_even )
        call pmc_bc_b_2d_even( b_cyl_m%pf_re(0), i_dim, bnd )
        do mode = 1, n_modes
          call pmc_bc_b_2d_even( b_cyl_m%pf_re(mode), i_dim, bnd )
          call pmc_bc_b_2d_even( b_cyl_m%pf_im(mode), i_dim, bnd )
        enddo

      case default
      end select
    enddo
  enddo

  ! advance vpml E if necessary (cylindrical require all emf for radii correction)
  ! the update_e_2d function calls a special cylindrical version.  This works for m = 0
  ! but produces errors for m > 0.
  if ( (step == 1) .and. is_active(this%vpml_all) ) then
    call update_e_2d( this%vpml_all, b_cyl_m%pf_re(0) )
    do mode = 1, n_modes

      ! quasi-3D specific vpml routine
      call update_e_vpml_2d_cyl_modes( this%vpml_all_re(mode), b_cyl_m%pf_re(mode) )
      call update_e_vpml_2d_cyl_modes( this%vpml_all_im(mode), b_cyl_m%pf_im(mode) )

    enddo
  endif

end subroutine update_boundary_bfld_emf_bound_cyl_modes
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine update_e_vpml_2d_cyl_modes(this, bfld)
!---------------------------------------------------

  implicit none

  type (t_vpml), intent(inout) :: this
  type( t_vdf ), intent(inout) :: bfld

  integer, dimension(2) :: vpml_range, update_range
  integer :: i_vpml

  ! loop through vpml boundaries to update interfaces
  do i_vpml=1, this%n_vpml

    select case (this%dir_vpml(i_vpml))

    ! x1 direction
    case (1)

      call update_interface_x1_2d(this%wall_array_b, i_vpml, this%pos_corner(i_vpml,2,:),&
                                  bfld, this%loc_vpml(i_vpml) )

    ! x2 direction
    case (2)

      call update_interface_x2_2d(this%wall_array_b, i_vpml, &
                                  bfld, this%loc_vpml(i_vpml) )

    end select

  enddo

  ! loop through vpml boundaries to update fields
  do i_vpml=1, this%n_vpml

    ! E and B must have the same range
    vpml_range = range(this%wall_array_e(i_vpml))

    select case (this%dir_vpml(i_vpml))

    ! x1 direction
    case (1)

      ! Left wall
      if (this%loc_vpml(i_vpml) == p_lower) then
        update_range(1) = vpml_range(1)+1
        update_range(2) = vpml_range(2)-1

        ! TODO: We see significant errors when using the _cyl_ version commented below.
        ! Make an accurate cylindrical version.
        call update_e_x1_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
          this%dtdx(2), update_range, vpml_range(1), &
          this%coef_e_low, this%coef_e_up, &
          this%pos_corner(i_vpml,2,:), p_lower, &
          this%if_diffuse(p_lower, 1))

        ! call update_cyl_e_x1_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
        !                         this%dtdx(2), update_range, vpml_range(1), &
        !                         this%coef_e_low, this%coef_e_up, &
        !                         this%pos_corner(i_vpml,2,:), p_lower, this%gir_pos )

      ! Right wall
      else
        update_range(1) = vpml_range(1)+1
        update_range(2) = vpml_range(2)

        ! TODO: We see significant errors when using the _cyl_ version commented below.
        ! Make an accurate cylindrical version.
        call update_e_x1_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                            this%dtdx(2), update_range, vpml_range(1), &
                            this%coef_e_low, this%coef_e_up, &
                            this%pos_corner(i_vpml,2,:), p_upper, &
                            this%if_diffuse(p_upper, 1) )

        ! call update_cyl_e_x1_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
        !                         this%dtdx(2), update_range, vpml_range(1), &
        !                         this%coef_e_low, this%coef_e_up, &
        !                         this%pos_corner(i_vpml,2,:), p_upper, this%gir_pos )

      endif

    ! x2 direction
    ! for now, cylindrical use the standard solver: solver in cyl. coord differs by dr/2r factor
    ! which is always small, since only upper boundary pml is possible
    case (2)

      ! Bottom wall
      if (this%loc_vpml(i_vpml) == p_lower) then
        update_range(1) = vpml_range(1)+1
        update_range(2) = vpml_range(2)-1

        call update_e_x2_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                            this%coef_e_low(2,:,:), this%dtdx(1), &
                            update_range, vpml_range(1)-1, &
                            this%if_diffuse(p_lower, 2) )

      ! Top wall
      else
        update_range(1) = vpml_range(1)+1
        update_range(2) = vpml_range(2)

        call update_e_x2_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                            this%coef_e_up(2,:,:), this%dtdx(1), &
                            update_range, vpml_range(1), &
                            this%if_diffuse(p_upper, 2) )
      endif

    end select

  enddo

end subroutine update_e_vpml_2d_cyl_modes
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update physical emf boundaries after each E field advance. Used for:
!  - Open boundaries (Lindman or PML)
!  - Conducting (to be removed)
!  - Perfect electric conductor (PEC)
!-----------------------------------------------------------------------------------------
subroutine update_boundary_efld_emf_bound_cyl_modes( this, e_cyl_m, b_cyl_m )

  implicit none

  class( t_emf_bound_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: e_cyl_m, b_cyl_m

  integer :: bnd, i_dim, mode, n_modes

  n_modes = ubound( e_cyl_m%pf_re, 1 )

  do i_dim = 1, 2
    do bnd = p_lower, p_upper
      select case (this%type( bnd, i_dim ))
      case (p_bc_lindman)
        call update_e_2d(this%lindman(bnd), e_cyl_m%pf_re(0))
        do mode = 1, n_modes
          call update_e_2d(this%lindman_re(bnd,mode), e_cyl_m%pf_re(mode))
          call update_e_2d(this%lindman_im(bnd,mode), e_cyl_m%pf_im(mode))
        enddo

      case ( p_bc_pec_odd )
        call pec_bc_e_2d_odd( e_cyl_m%pf_re(0), i_dim, bnd )
        do mode = 1, n_modes
          call pec_bc_e_2d_odd( e_cyl_m%pf_re(mode), i_dim, bnd )
          call pec_bc_e_2d_odd( e_cyl_m%pf_im(mode), i_dim, bnd )
        enddo

      case ( p_bc_pmc_odd )
        call pmc_bc_e_2d_odd( e_cyl_m%pf_re(0), i_dim, bnd )
        do mode = 1, n_modes
          call pmc_bc_e_2d_odd( e_cyl_m%pf_re(mode), i_dim, bnd )
          call pmc_bc_e_2d_odd( e_cyl_m%pf_im(mode), i_dim, bnd )
        enddo

      case ( p_bc_pec_even )
        call pec_bc_e_2d_even( e_cyl_m%pf_re(0), i_dim, bnd )
        do mode = 1, n_modes
          call pec_bc_e_2d_even( e_cyl_m%pf_re(mode), i_dim, bnd )
          call pec_bc_e_2d_even( e_cyl_m%pf_im(mode), i_dim, bnd )
        enddo

      case ( p_bc_pmc_even)
        call pmc_bc_e_2d_even( e_cyl_m%pf_re(0), i_dim, bnd )
        do mode = 1, n_modes
          call pmc_bc_e_2d_even( e_cyl_m%pf_re(mode), i_dim, bnd )
          call pmc_bc_e_2d_even( e_cyl_m%pf_im(mode), i_dim, bnd )
        enddo

      case default
      end select
    enddo
  enddo

  ! advance vpml B if necessary (cylindrical require all emf for radii correction)
  ! the update_b_2d function calls a special cylindrical version.  This works for m = 0
  ! but produces errors for m > 0.
  if ( is_active(this%vpml_all) ) then
    call update_b_2d( this%vpml_all, e_cyl_m%pf_re(0) )
    do mode = 1, n_modes

      ! quasi-3D specific vpml routine
      call update_b_vpml_2d_cyl_modes( this%vpml_all_re(mode), e_cyl_m%pf_re(mode) )
      call update_b_vpml_2d_cyl_modes( this%vpml_all_im(mode), e_cyl_m%pf_im(mode) )

    enddo
  endif

end subroutine update_boundary_efld_emf_bound_cyl_modes
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine update_b_vpml_2d_cyl_modes(this, efld)
!---------------------------------------------------

  implicit none

  type (t_vpml), intent(inout) :: this
  type( t_vdf ), intent(inout) :: efld

  integer, dimension(2) :: vpml_range, update_range
  integer :: i_vpml

  ! loop through vpml boundaries and call corresponding push
  do i_vpml=1, this%n_vpml

    select case (this%dir_vpml(i_vpml))

    ! x1 direction
    case (1)

      call update_interface_x1_2d(this%wall_array_e, i_vpml, this%pos_corner(i_vpml,2,:), &
        efld, this%loc_vpml(i_vpml) )

    ! x2 direction
    case (2)

      call update_interface_x2_2d(this%wall_array_e, i_vpml, efld, this%loc_vpml(i_vpml) )

    end select

  enddo

  ! loop through vpml boundaries and call corresponding push
  do i_vpml=1, this%n_vpml

    ! E and B must have the same range
    vpml_range = range(this%wall_array_e(i_vpml))

    select case (this%dir_vpml(i_vpml))

    ! x1 direction
    case (1)

      ! Left wall
      if (this%loc_vpml(i_vpml) == p_lower) then
        update_range(1) = vpml_range(1)
        update_range(2) = vpml_range(2)-1

        ! TODO: We see significant errors when using the _cyl_ version commented below.
        ! Make an accurate cylindrical version.
        call update_b_x1_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
          this%dtdx(2), update_range, vpml_range(1), &
          this%coef_b_low, this%coef_b_up, &
          this%pos_corner(i_vpml,2,:), p_lower )

        ! call update_cyl_b_x1_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
        !                         this%dtdx(2), update_range, vpml_range(1), &
        !                         this%coef_b_low, this%coef_b_up, &
        !                         this%pos_corner(i_vpml,2,:), p_lower, this%gir_pos )

      ! Right wall
      else
        update_range(1) = vpml_range(1)+1
        update_range(2) = vpml_range(2)-1

        ! TODO: We see significant errors when using the _cyl_ version commented below.
        ! Make an accurate cylindrical version.
        call update_b_x1_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
          this%dtdx(2), update_range, vpml_range(1), &
          this%coef_b_low, this%coef_b_up, &
          this%pos_corner(i_vpml,2,:), p_upper )

        ! call update_cyl_b_x1_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
        !                         this%dtdx(2), update_range, vpml_range(1), &
        !                         this%coef_b_low, this%coef_b_up, &
        !                         this%pos_corner(i_vpml,2,:), p_upper, this%gir_pos )

      endif

    ! x2 direction
    case (2)

      ! Bottom wall
      if (this%loc_vpml(i_vpml) == p_lower) then
        update_range(1) = vpml_range(1)
        update_range(2) = vpml_range(2)-1

        call update_b_x2_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                            this%coef_b_low(2,:,:), this%dtdx(1), &
                            update_range, vpml_range(1)-2 )

      ! Top wall
      else
        update_range(1) = vpml_range(1)+1
        update_range(2) = vpml_range(2)-1

        call update_b_x2_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                            this%coef_b_up(2,:,:), this%dtdx(1), &
                            update_range, vpml_range(1) )
      endif

    end select

  enddo


end subroutine update_b_vpml_2d_cyl_modes
!---------------------------------------------------


!---------------------------------------------------
subroutine move_window_emf_bound_cyl_modes( this, space )
!---------------------------------------------------
!       move boundary of electromagnetic field
!---------------------------------------------------

  implicit none

  class (t_emf_bound_cyl_modes), intent(inout) :: this
  type( t_space ),  intent(in) :: space ! local or global, only nx_move is needed

  integer :: bnd
  integer :: i_mode

  ! move lindman boundary data
  do bnd = p_lower, p_upper
    if( is_active(this%lindman_re(bnd,1)) ) then
      do i_mode = 1,p_max_cyl_modes
        call move_window(this%lindman_re(bnd,i_mode), space)
        call move_window(this%lindman_im(bnd,i_mode), space)
      enddo
    end if

    if( is_active(this%lindman(bnd)) ) call move_window(this%lindman(bnd), space )
  enddo


  ! move vpml boundary data
  if ( is_active(this%vpml_all) ) then
    do i_mode = 1,p_max_cyl_modes
      call move_window(this%vpml_all_re(i_mode), space)
      call move_window(this%vpml_all_im(i_mode), space)
    enddo

    call move_window( this%vpml_all, space )
  endif


end subroutine move_window_emf_bound_cyl_modes
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reshape_emf_bound_cyl_modes( this, old_lb, new_lb, no_co, send_msg, recv_msg )

  implicit none

  class (t_emf_bound_cyl_modes), intent(inout) :: this
  class(t_grid), intent(in) :: old_lb, new_lb
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i_bnd, i_mode

  do i_bnd = p_lower, p_upper
    if ( is_active( this%lindman(i_bnd) ) ) then
      call reshape_obj( this%lindman(i_bnd), old_lb, new_lb, no_co, send_msg, recv_msg )
      do i_mode = 1, p_max_cyl_modes
        call reshape_obj( this%lindman_re(i_bnd,i_mode), old_lb, new_lb, no_co, send_msg,&
                          recv_msg )
        call reshape_obj( this%lindman_im(i_bnd,i_mode), old_lb, new_lb, no_co, send_msg,&
                          recv_msg )
      enddo
    endif
  enddo

  if ( is_active( this%vpml_all ) ) then
    call reshape_obj( this%vpml_all, old_lb, new_lb, no_co, send_msg, recv_msg )
    do i_mode = 1, p_max_cyl_modes
      call reshape_obj( this%vpml_all_re(i_mode), old_lb, new_lb, no_co, send_msg, &
                        recv_msg )
      call reshape_obj( this%vpml_all_im(i_mode), old_lb, new_lb, no_co, send_msg, &
                        recv_msg )
    enddo
  endif

end subroutine reshape_emf_bound_cyl_modes
!-----------------------------------------------------------------------------------------

end module m_emf_bound_cyl_modes
