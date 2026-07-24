#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_cuda

#include "memory/memory.h"

use m_system
use m_parameters
use m_time_step,      only : t_time_step
use m_space,          only : t_space
use m_grid_define,    only : t_grid
use m_node_conf,      only : t_node_conf
use m_node_conf_cuda, only : t_node_conf_cuda
use m_restart,        only : t_restart_handle
use m_emf_define,     only : t_emf
use m_vdf_comm,       only : t_vdf_msg
use m_current_define, only : t_current
use m_species_define, only : t_species, t_part_idx, t_spec_msg

implicit none

private

type, extends( t_species ) :: t_species_cuda

  integer :: til_id
  integer :: til_id_min

contains

  procedure :: init                  => init_species_cuda
  procedure :: get_default_init_type => get_default_init_type_spec_cuda
  procedure :: list_algorithm        => list_algorithm_spec_cuda
  procedure :: grow_buffer           => grow_buffer_spec_cuda

end type

public :: t_species_cuda


interface
  subroutine sync_buffer_spec( til_id, spec_id, num_par_max, x, p, q, ix ) bind(C)
    use iso_c_binding
    use m_parameters
    integer(c_int), value :: til_id, spec_id, num_par_max
    type(c_ptr), value    :: x, p, q, ix
  end subroutine
end interface


contains


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_species_cuda( this, sp_id, interpolation, grid_center, grid, g_space, &
                              emf, jay, no_co, send_vdf, recv_vdf, bnd_cross, node_cross,&
                              send_spec, recv_spec, ndump_fac, restart, restart_handle, &
                              sim_options, tstep, tmin, tmax )

  implicit none

  class( t_species_cuda ), intent(inout) :: this

  integer, intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in)     :: grid
  type( t_space ),     intent(in) :: g_space
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

  integer :: g_til_aid_min, til_id

  select type( no_co )
  class is( t_node_conf_cuda )
    this % til_id     = no_co % ti_co % g_til_aid
    this % til_id_min = no_co % ti_co % g_til_aid_min_max( 1, no_co%local_rank()+1 )
  end select

  ! Call superclass init
  call this % t_species % init( sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
                     no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
                     recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, &
                     tmin, tmax )

  ! disable CPU sort when using CUDA 
  this%n_sort = 0

end subroutine init_species_cuda
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Returns the default init type for use in read_input
!-----------------------------------------------------------------------------------------
pure function get_default_init_type_spec_cuda( this )

  implicit none

  class( t_species_cuda ), intent(in) :: this

  character(len = 16) :: get_default_init_type_spec_cuda

  get_default_init_type_spec_cuda = "transposed"

end function get_default_init_type_spec_cuda
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_spec_cuda( this )

  implicit none
  class( t_species_cuda ), intent(in) :: this

  print *, ' '
  print *, trim(this%name),' :'
  print *, '- CUDA pusher'

  if ( this%free_stream ) then
    print *, '- Free streaming particles (no dudt)'
    write(0,*) '- (*warning*) ', trim(this%name), ' are free streaming!'
  endif

end subroutine list_algorithm_spec_cuda
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Grow the particle buffers
! There are two cuda-specific things which must happen here, hence the overload:
!    1. We need to modify grow_buffer_spec, see comments on grow_buffer_spec_cuda_, below.
!    2. We need to update the C pointers to the particle buffers
! The superclass grow_buffer_spec routine should never get called. If it does the code will crash.
!-----------------------------------------------------------------------------------------
subroutine grow_buffer_spec_cuda( this, num_par_req )

  implicit none

  class(t_species_cuda), intent(inout) :: this
  integer, intent(in) :: num_par_req

  ! call modified grow_buffer_spec routine
  call grow_buffer_spec_cuda_( this, num_par_req )

  ! Make sure C gets wind of where the new buffers are located
  call sync_buffer_spec( this%til_id-this%til_id_min+1, this%sp_id, this%num_par_max, &
                         c_loc(this%x), c_loc(this%p), c_loc(this%q), c_loc(this%ix) )

end subroutine grow_buffer_spec_cuda
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Grow the particle buffers
! The only difference between this version and the base class is that here we replace
! all instances of spec%num_par with spec%num_par_max. There are two problems which can 
! occur if num_par is used instead.
! * First, this routine uses num_par as a proxy for whether
!   particle buffers have already been allocated and just need to be resized, or whether they
!   need to be initialized. This is not a valid proxy in the cuda algorithm - num_par can be
!   greater than zero (indicating that there are particles on the GPU) before any CPU buffers 
!   are allocated. I'm not sure why the base class insists on using a proxy rather than actually
!   checking `if associated`, but anyway, num_par_max should be an accepatable proxy for our
!   use case.
! * Second, the base class routine assumes the old size of the buffer (num_par_old) is 
!   simply num_par. This is also not valid in the cuda algorithm. Sometimes num_par can
!   exceed num_par_max, which then results in a segfault when trying to copy the buffers.
!   Simply replacing num_par with num_par_max here has the small drawback of copying more 
!   data than necessary, but it should only be small.
! Final note: I suppose I could have proposed a change to the base class grow_buffer_spec()
! routine, but for now I decided to be conservative and just leave the base class alone and
! keep this change separate. This could always be changed later.
!-------------------------------------------------------------------------------
subroutine grow_buffer_spec_cuda_( this, num_par_req )

  implicit none

  class(t_species_cuda), intent(inout) :: this
  integer, intent(in) :: num_par_req

  real(p_k_part), dimension(:), pointer   :: temp1_r
  real(p_k_part), dimension(:,:), pointer :: temp2_r
  integer, dimension(:,:), pointer :: temp2_i

  integer :: num_par_old, num_par_new, n_x_dim

  ! The buffer size must always be a multiple of vector width in size because of SIMD code
  num_par_new = (( num_par_req + p_vecwidth - 1 ) / p_vecwidth) * p_vecwidth

  if ( this%num_par_max > 0 ) then

     num_par_old = this%num_par_max

     write(0,'(A,I0,A,A)') '[', mpi_node(), '] (* warning *) resizing particle buffers for this ', &
                            trim(this%name)
     write(0,'(A,I0,A,I0,A,I0)') '[', mpi_node(), '] (* warning *) Buffer size: ', this%num_par_max, ' -> ', num_par_new
     write(0,'(A,I0,A,I0)') '[', mpi_node(), &
             '] (* warning *) Number of particles currently in buffer: ', this%num_par_max

     if ( num_par_new <= num_par_old ) then
       ERROR('Invalid size for new buffer')
       call abort_program( p_err_invalid )
     endif

     ! particle positions (may not be p_x_dim)
     n_x_dim = this%get_n_x_dims()
     call alloc( temp2_r, (/ n_x_dim, num_par_new /))
     call memcpy( temp2_r, this%x, n_x_dim * num_par_old )
     call freemem( this%x )
     this%x => temp2_r

     ! particle cell index
     call alloc( temp2_i, (/ p_x_dim, num_par_new /))
     call memcpy( temp2_i, this%ix, p_x_dim * num_par_old )
     call freemem( this%ix )
     this%ix => temp2_i

     ! particle momenta
     call alloc( temp2_r, (/ p_p_dim, num_par_new /))
     call memcpy( temp2_r, this%p, p_p_dim * num_par_old )
     call freemem( this%p )
     this%p => temp2_r

     ! particle charge
     call alloc( temp1_r, (/ num_par_new /) )
     call memcpy( temp1_r, this%q, num_par_old )
     call freemem( this%q )
     this%q => temp1_r

#ifdef __HAS_SPIN__

    ! particle spin
   call alloc( temp2_r, (/ p_s_dim, num_par_new /))
   call memcpy( temp2_r, this%s, p_s_dim * num_par_old )
   call freemem( this%s )
   this%s => temp2_r

#endif

     ! Resize tracking data if necessary
     if ( this%add_tag ) then
        call alloc( temp2_i, (/ 2, num_par_new /) )
        call memcpy( temp2_i, this%tag, 2 * num_par_old )
        call freemem( this%tag )
        this%tag => temp2_i
     endif

     ! resize ionization data if necessary
     !(...)

     this%num_par_max = num_par_new
     write(0,'(A,I0,A)') '[', mpi_node(), '] (* warning *) resize successfull!'

  else

    ! no particles in buffer, simply reallocate the buffers
    call this % init_buffer( num_par_new )

  endif

end subroutine grow_buffer_spec_cuda_
!-------------------------------------------------------------------------------

end module m_species_cuda
