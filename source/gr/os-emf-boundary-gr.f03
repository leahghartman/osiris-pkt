#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_bound_gr

#include "memory/memory.h"

use m_system
use m_parameters
use m_emf_define, only: t_emf_bound
use m_geometry_gr
use m_emf_mur_gr, only: t_mur_gr
use m_restart, only: t_restart_handle

implicit none

private

! type of supported boundary conditions
integer, parameter :: p_bc_pec_gr    = 1, &
                      p_bc_mur_gr    = 2, &
                      p_bc_star_gr   = 3, &
                      p_bc_axi_gr    = 4

! IBM XL compilers
! This has to be explicitly made public so that you can access the superclass
! of t_emf_bound_gr
public :: t_emf_bound_gr

! type :: t_star_gr
! end type t_star_gr

! string to id restart data
character(len=*), parameter :: p_bcemf_gr_rst_id = "bcemf_gr rst data - 0x0001"


!-------------------------------------------------------------------------------
! t_emf_bound_gr class definition
!-------------------------------------------------------------------------------
type, extends( t_emf_bound ) :: t_emf_bound_gr
  class( t_geometry_gr ), pointer :: geometry
  type( t_mur_gr ) :: mur

contains

  procedure :: read_input     => read_input_gr
  procedure :: init           => init_gr
  procedure :: list_algorithm => list_algorithm_gr
  procedure :: cleanup        => cleanup_gr

  procedure :: update_boundary_b  => update_boundary_b_gr
  procedure :: update_boundary_e  => update_boundary_e_gr

  procedure :: write_checkpoint     => write_checkpoint_gr
  procedure :: restart_read         => restart_read_gr

end type t_emf_bound_gr

contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine read_input_gr( this, input_file, periodic, if_move )

  use m_input_file
  use m_emf_bound, only: read_nml

  implicit none

  class (t_emf_bound_gr), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move

  ! local variables
  character(20), dimension(2, p_x_dim) :: type

  namelist /nl_emf_bound/ type

  integer :: ierr, i, j

  ! executable statements
  type = "---"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_emf_bound", ierr )

  if ( ierr /= 0 ) then
    if (ierr < 0) then
      print *, "Error reading emf_bound parameters"
    else
      print *, "Error: emf_bound parameters missing"
    endif
    print *, "aborting..."
    stop
  endif


  read (input_file%nml_text, nml = nl_emf_bound, iostat = ierr)
  if (ierr /= 0) then
    print *, "Error reading emf_bound parameters"
    print *, "aborting..."
    stop
  endif

  ! read boundary condition types
  do i = 1, p_x_dim
    do j = 1, 2
      select case ( trim( type(j,i) ) )
          case ( "pec", "conducting" )
            this%type(j,i) = p_bc_pec_gr
          case ( "star" )
            this%type(j,i) = p_bc_star_gr
          case ( "axial", "axisymmetric" )
            this%type(j,i) = p_bc_axi_gr
          case ( "open", "mur" )
            this%type(j,i) = p_bc_mur_gr
          case ( "---" )
            this%type(j,i) = p_bc_invalid
          case default
            if ( mpi_node() == 0 ) then
               print "(A,I0,A,I0)", "    Error reading emf_bound parameters, dir = ", i, ", wall = ", j
               print *, '   Unknown boundary condition type, "', trim(type(i,j)), &
                        '", valid values are "conducting", "star", "open", "axisymmetric"'
               print *, '   aborting...'
            endif
            stop
      end select
    enddo
  enddo

  ! verify boundary conditions
  do i = 1, p_x_dim
    ! periodic boundaries and moving window will override local settings
    if ((.not. periodic(i)) .and. (.not. if_move(i))) then
       do j = 1, 2
         if (this%type(j,i) == p_bc_invalid) then
           print *, ""
           print *, "   Error reading emf_bound parameters, dir = ", i, ", wall = ", j
           print *, "   Invalid or no boundary conditions specified."
           print *, "   aborting..."
           stop
         elseif (this%type(j,i) == p_bc_periodic) then
           print *, ""
           print *, "   Error reading emf_bound parameters, dir = ", i
           print *, "   Periodic boundaries cannot be specified unless global periodic"
           print *, "   boundaries were set in the node_conf section"
           print *, "   aborting..."
           stop
         endif
       enddo
    endif
  enddo

  ! check special spherical boundaries
  do i = 1, p_x_dim
     do j = 1, 2
       if ((this%type(j,i) == p_bc_star_gr) .and. &
           ((i /= p_rc_dim) .or. ( j /= p_lower ))) then
           print *,        ""
           print "(A,I0)", "   Error reading emf_bound parameters, dir = ", i
           print "(A)",    "   Star boundaries can only be specified for"
           print "(A,I0)", "   the lower boundary of dimension ",p_rc_dim
           print "(A)",    "   aborting..."
           stop
       endif

       if ((this%type(j,i) == p_bc_mur_gr) .and. &
           ((i /= p_rc_dim) .or. ( j /= p_upper ))) then
           print *,        ""
           print "(A,I0)", "   Error reading emf_bound parameters, dir = ", i
           print "(A)",    "   Open boundaries can only be specified for"
           print "(A,I0)", "   the upper boundary of dimension ",p_rc_dim
           print "(A)",    "   aborting..."
           stop
       endif

       if ((this%type(j,i) == p_bc_axi_gr) .and. &
           (i /= p_tc_dim)) then
           print *,        ""
           print "(A,I0)", "   Error reading emf_bound parameters, dir = ", i
           print "(A)",    "   Axisymmetric boundaries can only be specified for"
           print "(A,I0)", "   dimension ",p_tc_dim
           print "(A)",    "   aborting..."
           stop
       endif
     enddo
  enddo

end subroutine
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine init_gr( this, b, e, dt, part_interpolation, &
                            space, no_co, grid, restart, restart_handle )

  use m_emf_bound

  use m_restart
  use m_emf_define

  use m_logprof

  use m_space
  use m_node_conf
  use m_grid_define

  use m_parameters

  use m_vdf_define
  use m_vdf_comm
  use m_emf_mur_gr, only: setup

  implicit none

  ! dummy variables

  class (t_emf_bound_gr), intent( inout )  ::  this

  type( t_vdf ),       intent(in) :: b, e
  real( p_double ), intent(in) :: dt
  integer, intent(in) :: part_interpolation
  type( t_space ),     intent(in) :: space
  class( t_node_conf ), intent(in) :: no_co
  class( t_grid ), intent(in) :: grid

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables

  logical, dimension(p_x_dim) :: ifpr_l, if_move_l
  integer :: i_bnd, i_x_dim

  ! executable statements

  if ( restart ) then

    call this%restart_read( restart_handle )

  else

     ifpr_l = no_co%ifpr()
     if_move_l = if_move(space)

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

  endif

  ! setup Mur boundary
  do i_x_dim=1, p_x_dim
     if ( this%type( p_upper, i_x_dim ) == p_bc_mur_gr ) then
       call setup( this%mur, this%geometry, grid, dt, restart, restart_handle )
     endif
  enddo

  ! setup events
  if (emfboundev==0) emfboundev = create_event('update emf boundary')
end subroutine
!---------------------------------------------------

!---------------------------------------------------
subroutine list_algorithm_gr( this )

  implicit none

  class( t_emf_bound_gr ), intent(in) :: this

  integer :: i

  print *, ' '
  print *, '- EM Field Boundary condtions:'
  do i = 1, p_x_dim
    print '(A,I0,A,A,A,A,A)', '     x', i, ' : [', &
                trim(bnd_type_name_gr( this%global_type( p_lower, i ))), ', ', &
                trim(bnd_type_name_gr( this%global_type( p_upper, i ))), ']'
  enddo

end subroutine

!-------------------------------------------------------------------------------
! Get the name of a boundary type
!-------------------------------------------------------------------------------
function bnd_type_name_gr( type )

  implicit none

  integer, intent(in) :: type
  character(len = 39) :: bnd_type_name_gr

  select case ( type )
    case( p_bc_other_node )
      bnd_type_name_gr = "Other node"
    case( p_bc_periodic )
      bnd_type_name_gr = "Periodic"

    case( p_bc_pec_gr )
      bnd_type_name_gr = "Conducting (perfect electric conductor)"
    case( p_bc_star_gr )
      bnd_type_name_gr = "Star"
    case( p_bc_mur_gr )
      bnd_type_name_gr = "Open (Mur algorithm)"
    case( p_bc_axi_gr )
      bnd_type_name_gr = "Poloidal axis"
    case default
      bnd_type_name_gr = "unknown"
  end select

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine cleanup_gr( this )
!---------------------------------------------------
!release memory used by ptr-components of this variable
!---------------------------------------------------
  use m_emf_mur_gr, only: cleanup

  implicit none
  class (t_emf_bound_gr), intent( inout ) :: this

  ! cleanup Mur
  call cleanup( this%mur )

end subroutine

!-------------------------------------------------------------------------------
! Update B field boundary
!-------------------------------------------------------------------------------
subroutine update_boundary_b_gr( this, e, b, step )
  
  use m_vdf_define, only : t_vdf
  use m_emf_mur_gr, only : update_b_2d_mur
  use m_emf_star_gr, only : update_b_2d_star
  use m_emf_axi_gr, only : update_b_2d_axi
  use m_emf_pec_gr, only : update_b_2d_pec

  class( t_emf_bound_gr ), intent(inout) :: this
  class( t_vdf ), intent(inout) :: e, b
  integer, intent(in) :: step

  integer :: bnd, i_dim

  ! executable statements
  select case ( b%x_dim_ )
  case (2)
    do i_dim = 1, 2
      do bnd = p_lower, p_upper
        select case (this%type( bnd, i_dim ))
          case (p_bc_pec_gr)
            call update_b_2d_pec(this%geometry, b, bnd)
          case (p_bc_star_gr)
            call update_b_2d_star(this%geometry, b)
          case (p_bc_mur_gr)
            call update_b_2d_mur(this%mur, b, step)
          case (p_bc_axi_gr)
            call update_b_2d_axi(this%geometry, b, bnd)

          case default
        end select
      enddo
    enddo
  end select


end subroutine

!-------------------------------------------------------------------------------
! Update E field boundary
!-------------------------------------------------------------------------------
subroutine update_boundary_e_gr( this, e, b )
  
  use m_vdf_define, only : t_vdf
  use m_emf_mur_gr, only : update_e_2d_mur
  use m_emf_star_gr, only : update_e_2d_star
  use m_emf_axi_gr, only : update_e_2d_axi
  use m_emf_pec_gr, only : update_e_2d_pec

  class( t_emf_bound_gr ), intent(inout) :: this
  class( t_vdf ), intent(inout) :: e, b

  integer :: bnd, i_dim

  ! executable statements
  select case ( b%x_dim_ )
  case (2)
    do i_dim = 1, 2
      do bnd = p_lower, p_upper
        select case (this%type( bnd, i_dim ))
          case (p_bc_pec_gr)
            call update_e_2d_pec(this%geometry, e, bnd)
          case (p_bc_star_gr)
            call update_e_2d_star(this%geometry, e)
          case (p_bc_mur_gr)
            call update_e_2d_mur(this%mur, e)
          case (p_bc_axi_gr)
            call update_e_2d_axi(this%geometry, e, bnd)
          case default
        end select
      enddo
    enddo
  end select


end subroutine

!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_gr( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------
  use m_parameters
  use m_restart

  implicit none

  class (t_emf_bound_gr), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for emf_bound_gr object.'
  integer :: i_dim, bnd, ierr

  restart_io_wr( p_bcemf_gr_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  call this % t_emf_bound % write_checkpoint( restart_handle )

  ! write mur boundary data if needed
  do i_dim = 1, 2
    do bnd = p_lower, p_upper
      if (this%type( bnd, i_dim ) == p_bc_mur_gr) then
        call this % mur % fsave % write_checkpoint( restart_handle )
      endif
    enddo
  enddo
  
end subroutine write_checkpoint_gr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_gr( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file
!-----------------------------------------------------------------------------------------
  use m_parameters
  use m_restart

  implicit none

  class (t_emf_bound_gr), intent(inout)  ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for emf_bound_cyl_modes object.'
  character(len=len(p_bcemf_gr_rst_id)) :: rst_id
  integer :: i_dim, bnd, ierr

  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  call this % t_emf_bound % restart_read( restart_handle )
   
  do i_dim = 1, 2
    do bnd = p_lower, p_upper
      if (this%type( bnd, i_dim ) == p_bc_mur_gr) then
        call this % mur % fsave % read_checkpoint( restart_handle )
      endif
    enddo
  enddo
  
end subroutine restart_read_gr
!-----------------------------------------------------------------------------------------


end module m_emf_bound_gr
