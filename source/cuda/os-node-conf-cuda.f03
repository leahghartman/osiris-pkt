!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     node configuration class for cuda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_node_conf_cuda

#include "memory/memory.h"

use m_parameters
use m_system
use m_restart
use m_node_conf, only : t_node_conf, p_topology_mpi
use m_node_conf_tiles, only : t_node_conf_tiles, p_default, p_heavy, p_light, &
                              p_spec_buf_block_t
use m_tile_conf
use m_input_file, only : t_input_file, get_namelist

implicit none

private

! string to id restart data
character(len=*), parameter :: p_nconf_cuda_rst_id = "node conf cuda rst data - 0x0001"

type, extends(t_node_conf_tiles) :: t_node_conf_cuda

  ! cuda configuration class
  integer :: chunk_pool_size, chunk_size, com_buf_size, ihole_size, nchunks_pad
  integer :: nchunks_max_init, nblocks_requested
  real(p_double) :: ihole_size_frac, nchunks_pad_frac, chunk_growth_factor

contains

  procedure :: read_input => read_input_node_conf_cuda
  procedure :: setup_tiles => setup_node_conf_tiles_cuda

end type t_node_conf_cuda

public :: t_node_conf_cuda

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_node_conf_cuda( this, input_file, x_dim )

  implicit none

  class( t_node_conf_cuda ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: x_dim

  integer, dimension(p_x_dim) ::  node_number
  logical, dimension(p_x_dim) ::  if_periodic
  integer :: n_threads, grow_spec_buff
  integer :: chunk_pool_size, chunk_size, com_buf_size, ihole_size, nchunks_pad
  integer :: nchunks_max_init, nblocks_requested

  real(p_double) :: ihole_size_frac, nchunks_pad_frac, chunk_growth_factor
  real(p_double) :: tiles_per_rank, num_par_per_tile

  character(len=16) :: topology, tile_topology, force_omp

  integer, dimension(p_x_dim) :: tile_number

  namelist /nl_node_conf/ node_number, if_periodic, n_threads, topology, tile_number, &
                          tile_topology, force_omp, grow_spec_buff, &
                          chunk_pool_size, chunk_size, com_buf_size, &
                          ihole_size_frac, ihole_size, &
                          nchunks_pad_frac, nchunks_pad, nchunks_max_init, chunk_growth_factor, &
                          nblocks_requested

  integer :: i, ierr

  this%x_dim = x_dim

  node_number = 1
  if_periodic = .false.
  n_threads = 1
  tile_number = 0
  grow_spec_buff = p_spec_buf_block_t

  ! cuda-specific parameters
  chunk_size = 512     ! size of a chunk in units of num particles
  chunk_pool_size = -1 ! per MPI rank, num chunks to allocate, in units of num chunks
  com_buf_size = 1024  ! per MPI rank - size of buffer for CPU-GPU communication, in units of MB, default 1GB

  ! per tile - size of buffer to store moving particles in sort
  ihole_size_frac = 0.05 ! specified as a fraction of num particles/tile
  ihole_size = -1        ! optionally can specify in units of num particles

  ! per tile - (1) how many empty chunks each tile should have to buffer particles moving
  ! locally during the sort, and (2) how many chunks each tile should keep to buffer particles
  ! moving externally (via MPI) during the sort
  nchunks_pad_frac = 0.05 ! specified as a fraction of num chunks/tile
  nchunks_pad = -1        ! optionally can specify in units of num chunks

  nchunks_max_init = 20
  chunk_growth_factor = 2.0_p_double
  nblocks_requested = 10000

  ! Default topology routines
  topology = "mpi"
  tile_topology = "snake"
  force_omp = "light"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_node_conf", ierr )

  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading node_conf parameters"
      else
        write(0,*) "Error: node_conf parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_node_conf, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "aborting..."
    endif
    stop
  endif

  this%no_num = 1
  this%ti_co%g_til_num = 1
  do i=1, x_dim
    this%nx(i)   = node_number(i)
    this%ifpr_(i) = if_periodic(i)
    this%no_num  = this%no_num * this%nx(i)
    this%ti_co%g_til(i) = tile_number(i)
    this%ti_co%g_til_num = this%ti_co%g_til_num * this%ti_co%g_til(i)
  enddo

  ! process node topology type
  select case ( trim( topology ))
  case ( "mpi" )
    this%topology_type = p_topology_mpi
  case ( "old" )
    if ( mpi_node() == 0 ) then
      write(0,*) "(*warning*) Defaulting to MPI topology with tiling"
    endif
    this%topology_type = p_topology_mpi

#ifdef __bgq__
  case ( "bgq" )
    if ( mpi_node() == 0 ) then
      write(0,*) "(*warning*) Defaulting to MPI topology with tiling"
    endif
    this%topology_type = p_topology_mpi
#endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) 'Error reading node_conf parameters'
      write(0,*) 'invalid topology type: "', trim( topology ), '"'
      write(0,*) 'With tiling, only valid value is "mpi".'
      write(0,*) 'aborting...'
    endif
    stop
  end select

  ! process tile topology type
  select case ( trim( tile_topology ))
  case ( "hilbert" )
    this%ti_co%topology_type = p_topology_hilbert

  case ( "snake" )
    this%ti_co%topology_type = p_topology_snake

  case ( "viktor" )
    this%ti_co%topology_type = p_topology_viktor

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) 'Error reading node_conf parameters'
      write(0,*) 'invalid tile_topology type: "', trim( tile_topology ), '"'
      write(0,*) 'Valid values are "hilbert", "snake" and "viktor".'
      write(0,*) 'aborting...'
    endif
    stop
  end select

  ! process forcing heavy/light tiles
  ! we force that all tiles be light for cuda, since there are only ever a very small
  ! number of particles on the host
  select case( trim(force_omp) )
  ! case( "default" )
  !   this%omp_patt%force_omp = p_default
  ! case( "heavy" )
  !   this%omp_patt%force_omp = p_heavy
  case( "light" )
    this%omp_patt%force_omp = p_light
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) 'Error reading node_conf parameters'
      write(0,*) 'invalid force_omp type: "', trim( force_omp ), '"'
      write(0,*) 'Only "light" is allowed when running with cuda.'
      write(0,*) 'aborting...'
    endif
    stop
  end select

  if ( this%no_num < 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) 'Number of nodes (node_number) must be > 0 in all directions'
      write(0,*) "aborting..."
    endif
    stop
  endif

  ! set how many particles to grow species buffer by in update boundary/dlb routines
  if ( grow_spec_buff > 0 ) then
    p_spec_buf_block_t = grow_spec_buff
  else
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "grow_spec_buff must be > 0"
      write(0,*) "aborting..."
    endif
    stop
 endif

  this%chunk_pool_size = chunk_pool_size
  this%chunk_size = chunk_size
  this%com_buf_size = com_buf_size
  this%ihole_size = ihole_size
  this%ihole_size_frac = ihole_size_frac
  this%nchunks_pad = nchunks_pad
  this%nchunks_pad_frac = nchunks_pad_frac
  this%nchunks_max_init = nchunks_max_init
  this%chunk_growth_factor = chunk_growth_factor
  this%nblocks_requested = nblocks_requested

  ! Multithread support
  this%n_threads = n_threads

  if ( n_threads < 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "Invalid number of threads per node (n_threads)."
      write(0,*) "aborting..."
    endif
    stop
  endif

#ifndef _OPENMP
  if ( n_threads /= 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "Multiple threads per node are only supported when using OpenMP."
      write(0,*) "aborting..."
    endif
    stop
  endif
#endif

  if ( this%ti_co%g_til_num < 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "Number of tiles (tile_number) must be > 0 in all directions."
      write(0,*) "aborting..."
    endif
    stop
  endif

  if ( this%ti_co%g_til_num < this%no_num ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "Number of tiles (tile_number) must be > number of nodes (node_number)."
      write(0,*) "aborting..."
    endif
    stop
  endif

end subroutine read_input_node_conf_cuda
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! sets up t_node_conf data structure for tiles simulation object (ie not for simulation group)
!-----------------------------------------------------------------------------------------
subroutine setup_node_conf_tiles_cuda( n_a, n_b )

  implicit none

  class( t_node_conf_cuda ), intent(inout) :: n_a
  class( t_node_conf ), intent(in) :: n_b

  call n_a % t_node_conf_tiles % setup_tiles( n_b )

  ! set n_threads to one for all tiles (to avoid extra current/species energy allocs)
  n_a % n_threads = 1

  select type(n_b)
  class is(t_node_conf_cuda)

    n_a%chunk_pool_size = n_b%chunk_pool_size
    n_a%chunk_size = n_b%chunk_size
    n_a%com_buf_size = n_b%com_buf_size
    n_a%ihole_size = n_b%ihole_size
    n_a%nchunks_pad = n_b%nchunks_pad
    n_a%nchunks_max_init = n_b%nchunks_max_init
    n_a%nblocks_requested = n_b%nblocks_requested
    n_a%chunk_growth_factor = n_b%chunk_growth_factor

  end select

end subroutine setup_node_conf_tiles_cuda
!-----------------------------------------------------------------------------------------

end module m_node_conf_cuda
