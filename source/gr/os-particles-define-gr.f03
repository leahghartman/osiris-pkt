#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_define_gr

#include "memory/memory.h"

use m_system
use m_parameters

use m_input_file
use m_node_conf
use m_restart
use m_space
use m_grid_define
use m_emf_define
use m_bnd
use m_vdf_define
use m_vdf_comm
use m_current_define
use m_particles_define
use m_time_step
use m_node_conf

use m_geometry_gr
use m_species_define_gr
use m_species_group_gr
use m_zpulse, only: t_zpulse_list

implicit none

private

! string to id restart data
character(len=*), parameter, public :: p_part_gr_rst_id = "particles_gr rst data - 0x0006"

type, extends( t_particles ) :: t_particles_gr

  class( t_geometry_gr ), pointer :: geometry
  
  ! number of species groups
  integer :: num_groups

  type( t_species_group ), dimension(:), pointer :: spec_group  => null()

contains
  procedure :: read_input      => read_input_particles_gr
  procedure :: allocate_objs   => allocate_objs_particles_gr
  procedure :: list_algorithm  => list_algorithm_particles_gr
  procedure :: init            => init_particles_gr
  procedure :: cleanup         => cleanup_particles_gr
  procedure :: advance_deposit => advance_deposit_particles_gr
  procedure :: report          => report_particles_gr

  procedure :: write_checkpoint  => write_checkpoint_gr
  procedure :: restart_read      => restart_read_gr

end type t_particles_gr

interface
  subroutine read_input_particles_gr( this, input_file, periodic, &
    if_move, grid, dt, sim_options )

  import t_particles_gr, t_input_file, p_double, t_options, t_grid

  class( t_particles_gr ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ),               intent(in) :: grid
  real(p_double), intent(in) :: dt
  type( t_options ), intent(in) :: sim_options

  end subroutine
end interface

interface
  subroutine init_particles_gr( this, g_space, jay, emf, grid, no_co, bnd, zpulse_list, &
                                ndump_fac, restart, restart_handle, t, tstep, tmin, tmax,&
                                sim_options )

  import t_particles_gr, t_space, t_emf, t_vdf, t_grid, t_node_conf, &
    t_restart_handle, p_double, t_options, t_current, t_bnd, t_time_step, t_zpulse_list

  class( t_particles_gr ), intent(inout) :: this
  type( t_space ),     intent(in) :: g_space
  class( t_emf) ,intent(inout) :: emf
  class( t_current ) ,intent(inout) :: jay
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  class(t_bnd), intent(inout) :: bnd
  class( t_zpulse_list ), intent(in) :: zpulse_list
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  real(p_double), intent(in) :: t
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax
  type( t_options ), intent(in) :: sim_options

  end subroutine
end interface

interface
  subroutine advance_deposit_particles_gr( this, emf, jay, tstep, t, &
    no_co, options )

  import t_particles_gr, t_emf, t_current, p_double, t_time_step, &
    t_node_conf, t_options

  class( t_particles_gr ), intent(inout) :: this
  class( t_emf ), intent( inout )  ::  emf
  class( t_current ), intent(inout) :: jay

  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t
  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  end subroutine

end interface

interface
  subroutine report_particles_gr( this, emf, g_space, grid, no_co, &
    tstep, t, tmin, send_msg, recv_msg )

    import t_particles_gr, t_emf, t_space, t_grid, t_node_conf, &
      t_time_step, p_double, t_vdf_msg

    class( t_particles_gr ), intent(inout) :: this
    class( t_emf ),       intent(inout) :: emf

    type( t_space ),        intent(in) :: g_space
    class( t_grid ),         intent(in) :: grid
    class( t_node_conf ),    intent(in) :: no_co
    type( t_time_step ),    intent(in) :: tstep
    real(p_double),         intent(in) :: t, tmin
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

public t_particles_gr

contains

!-------------------------------------------------------------------------------
! Allocate particle species
!-------------------------------------------------------------------------------
subroutine allocate_objs_particles_gr( this )

  use m_species_define, only : t_species
  use m_species_define_gr, only : t_species_gr
  use m_photons_define_gr, only : t_photons_gr
  use m_species_group_gr, only : t_species_group

  implicit none

  class( t_particles_gr ), intent(inout) :: this

  integer :: i
  class(t_species_gr), pointer :: spec, tail
  class(t_photons_gr), pointer :: photons
  class(t_species), pointer :: spec_base, e, p, g
  logical :: check

  check = .false.
  
  ! allocate species list with proper types
  if ( this%num_species > 0 ) then

    allocate(t_species_gr :: spec)
    this % species => spec

    do i = 2, this%num_species
      if ( check ) then
        tail => photons
      else
        tail => spec
      endif

      if ( mod(this%num_species - 3*this%num_groups + i, 3) .ne. 0 ) then
        allocate(t_species_gr :: spec)
        tail % next => spec
        check = .false.
      else
        allocate(t_photons_gr :: photons)
        tail % next => photons
        check = .true.
      endif
    enddo
  endif


  if ( this%num_groups > 0 ) then
    ! allocate species groups
    allocate( t_species_group :: this%spec_group( this%num_groups ) )

    ! set group species
    ! skip standard species first
    spec_base => this % species
    do i = 1, this%num_species - 3*this%num_groups
      spec_base => spec_base % next
    enddo
    
    ! find and set pointers to group species
    do i = 1, this%num_groups
      e => spec_base
      spec_base => spec_base % next

      p => spec_base
      spec_base => spec_base % next

      g => spec_base
      spec_base => spec_base % next

      call this % spec_group(i) % set_species( e, p, g )
    enddo
  endif

end subroutine allocate_objs_particles_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Printout the algorithm used by each species pusher
!-------------------------------------------------------------------------------
subroutine list_algorithm_particles_gr( this )

  use m_species_define, only : t_species
  implicit none
  class( t_particles_gr ), intent(in) :: this
  class( t_species ), pointer :: species
  integer :: i

  print *, ' '
  print '(A)', 'Particles:'

  print '(A,I0)', '- Total number of species: ', this%num_species
  print '(A,I0)', '- Number of regular species: ', this%num_species - 3*this%num_groups

  if ( this%num_species > 0 ) then
  ! particle shape
  select case ( this%interpolation )
    case (p_linear)
    print '(A)', '- Linear (1st order) interpolation'
  end select

  endif

  if ( this%grid_center ) then
    print '(A)', '- Using spatially centered grid values for interpolation'
  else
    print '(A)', '- Using staggered Yee grid values for interpolation'
  endif

  ! Current deposition
  if ( this%low_jay_roundoff ) then
    print *, '- Depositing current using low roundoff algorithm'
  endif

  species => this % species
  do i=1, this%num_species - 3*this%num_groups
    if (.not. associated(species)) exit
    call species % list_algorithm()
    species => species % next
  enddo

  print '(A,I0)', '- Number of groups: ', this%num_groups
  if ( this%num_groups > 0 ) then

    do i=1, this%num_groups
      print '(A,I0)', ' - Species group: ', i
      call this % spec_group(i) % list_algorithm()
    enddo
  endif

end subroutine list_algorithm_particles_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Cleanup particles object
!-------------------------------------------------------------------------------
subroutine cleanup_particles_gr( this )
!-------------------------------------------------------------------------------

  use m_particles, only : cleanup

  implicit none

  class(t_particles_gr) , intent(inout) :: this

  integer :: i

  ! cleanup superclass
  call this % t_particles % cleanup()

  ! cleanup local data structures
  if ( this%num_groups > 0 ) then
    do i = 1, this%num_groups
      call this % spec_group(i) % cleanup()
    enddo
    deallocate( this%spec_group )
  endif
  

end subroutine cleanup_particles_gr
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! write checkpoint information
!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_gr( this, restart_handle )

use m_species_define, only : t_species
use m_particles
use m_particles_charge
use m_restart , only : t_restart_handle
use m_parameters

implicit none

  class ( t_particles_gr ), intent(in) ::  this
  type ( t_restart_handle ), intent(inout)    ::  restart_handle

  character(len=*), parameter :: err_msg = 'error writing checkpoint data for particles_gr object.'
  class( t_species ), pointer :: spec
  integer :: ierr
  
  restart_io_wr( p_part_gr_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  ! write local data to file
  restart_io_wr( this%interpolation, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  restart_io_wr( this%num_species, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  restart_io_wr( this%n_current, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  call restart_write( this%charge, restart_handle )
  
  spec => this % species
  do
    if (.not. associated(spec)) exit
    select type( spec )
      class is ( t_species_gr )
        call spec % restart_write( restart_handle )
      class default
        ERROR('write_checkpoint_gr must be called with t_species_gr objects')
        call abort_program( p_err_invalid )
    end select
    spec => spec % next
  enddo
  

end subroutine write_checkpoint_gr

!-----------------------------------------------------------------------------------------
! read object information from a restart file
!-----------------------------------------------------------------------------------------
subroutine restart_read_gr( this, restart_handle )
  
  use m_parameters
  use m_restart

  implicit none

  class( t_particles_gr ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables
  character(len=*), parameter :: err_msg = 'error reading restart data for particles_gr object.'
  integer :: interpolation, num_species
  integer :: ierr

  character(len=len(p_part_gr_rst_id)) :: rst_id

  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! check if restart file is compatible
  if ( rst_id /= p_part_gr_rst_id) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (part)')
    ERROR('rst_id = ', rst_id)
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( interpolation, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( num_species, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%n_current, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  if ( this%interpolation /= interpolation ) then
    ERROR('The interpolation level in the input deck is different')
    ERROR('from the interpolation level in the restart file')
    call abort_program(p_err_invalid)
  endif

  if (this%num_species /= num_species) then
    ERROR('The number of species specified in the input deck is different')
    ERROR('from the number of species in the restart file')
    call abort_program(p_err_invalid)
  endif

end subroutine restart_read_gr
!-----------------------------------------------------------------------------------------


end module m_particles_define_gr
