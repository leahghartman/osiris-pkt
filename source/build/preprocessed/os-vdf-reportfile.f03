# 1 "vdf/os-vdf-reportfile.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "vdf/os-vdf-reportfile.f03" 2
! if __TEMPLATE__ is defined then read the template definition at the end of the file


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
# 5 "vdf/os-vdf-reportfile.f03" 2
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
# 6 "vdf/os-vdf-reportfile.f03" 2

module m_vdf_reportfile

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
# 10 "vdf/os-vdf-reportfile.f03" 2

use m_vdf_define
use m_space
use m_system

private

interface init_diag_file_vdf
  module procedure init_diag_file_vdf
end interface

interface write_vdf_1d_r4
  module procedure write_vdf_1d_r4
end interface

interface write_vdf_2d_r4
  module procedure write_vdf_2d_r4
end interface

interface write_vdf_3d_r4
  module procedure write_vdf_3d_r4
end interface

interface write_vdf_1d_r8
  module procedure write_vdf_1d_r8
end interface

interface write_vdf_2d_r8
  module procedure write_vdf_2d_r8
end interface

interface write_vdf_3d_r8
  module procedure write_vdf_3d_r8
end interface

public :: init_diag_file_vdf
public :: write_vdf_1d_r4, write_vdf_2d_r4, write_vdf_3d_r4
public :: write_vdf_1d_r8, write_vdf_2d_r8, write_vdf_3d_r8

contains

!---------------------------------------------------------------------------------------------------
! Initializes t_diag_file file object with supplied information
!---------------------------------------------------------------------------------------------------
subroutine init_diag_file_vdf( this, report, g_space, grid, diagFile )

  use m_space
  use m_grid_define
  use m_diagnostic_utilities

  implicit none

  class( t_vdf ), intent(in) :: this
  type( t_vdf_report ), intent(in) :: report
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  class( t_diag_file ), allocatable :: diagFile

  integer :: i

  ! Initialize grid output file with default parameters
  call create_diag_file( diagFile )

  diagFile % ftype = p_diag_grid

  ! Fill in parameters from specific report options
  diagFile%name = report%name

  ! Iteration info
  diagFile%iter%n = report%n
  diagFile%iter%t = report%t
  diagFile%iter%time_units = report%time_units
  diagFile%grid%offset_t = report%offset_t

  ! Grid info
  diagFile%grid%ndims = this % x_dim_
  diagFile%grid%name = report % name
  diagFile%grid%label = report % label
  diagFile%grid%units = report % units

  do i = 1, this % x_dim_
    diagFile%grid%count(i) = grid % g_nx(i)
    diagFile%grid%axis(i)%type = diag_axis_linear
    diagFile%grid%axis(i)%min = xmin( g_space, i )
    diagFile%grid%axis(i)%max = xmax( g_space, i )
    diagFile%grid%axis(i)%name = report%xname(i)
    diagFile%grid%axis(i)%label = report%xlabel(i)
    diagFile%grid%axis(i)%units = report%xunits(i)
    diagFile%grid%offset_x(i) = report%offset_x(i)
  enddo

  ! Set file name and path
  diagFile%filename = trim(report%filename)
  diagFile%filepath = trim(report%path)

end subroutine init_diag_file_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Generate specific template functions for single and double precision datatypes.
! Note that the module interfaces are not generated automatically and must be explicity written
! in the module header above.
!
! - write_vdf_1d_r4
! - write_vdf_2d_r4
! - write_vdf_3d_r4
! - write_vdf_1d_r8
! - write_vdf_2d_r8
! - write_vdf_3d_r8
!---------------------------------------------------------------------------------------------------




! single precision real



# 1 "vdf/os-vdf-reportfile.f03" 1
! if is defined then read the template definition at the end of the file
# 208 "vdf/os-vdf-reportfile.f03"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Merge / write a 1D dataset
!---------------------------------------------------------------------------------------------------
subroutine write_vdf_1d_r4( vdf, report, fc, g_space, grid )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  use m_parameters
  use zdf
  use m_diagfile

  implicit none

  class( t_vdf ), intent(in) :: vdf
  type( t_vdf_report ), intent(in) :: report
  integer, intent(in) :: fc
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid

  ! Buffer for merging data - only group leaders will allocate this
  real( p_single ), dimension(:), pointer :: write_buffer

  ! Buffer for MPI communications - all nodes allocate this
  real( p_single ), dimension(:), pointer :: comm_buffer

  real( p_single ), dimension(:), pointer :: gather_buffer

  ! MPI datatype for receiving data on group root nodes
  integer :: recv_type

  ! Sizes and offsets for MPI datatype
  integer, dimension(1) :: sizes, starts

  ! Rank of node sending data
  integer :: rank

  ! Status variable for messages
  integer, dimension( MPI_STATUS_SIZE ) :: stat

  ! Ping variables
  integer :: ping, ping_handle

  ! Total data size being received
  integer :: data_size

  ! Sizes and displacements of data from each node
  integer, dimension(:), pointer :: recvcounts, displs

  ! Diagnostic file object
  class(t_diag_file), allocatable :: diagFile

  integer :: i, i1, k, ierr

  type( t_diag_dataset ) :: parDset
  type( t_diag_chunk ) :: parChunk

  ! Dummy receive buffer for MPI_GATHERV
  integer, dimension(1) :: tmp

  ! Merge data if necessary
  select case ( grid%io%merge_type )
  case (p_none)

    ! No merging, just copy local data to output buffer
    call alloc(write_buffer, [ vdf % nx_(1) ],"vdf/os-vdf-reportfile.f03",279)

    do i1 = 1, vdf % nx_(1)
      write_buffer( i1 ) = real( vdf % f1 ( fc, i1 ), p_single )
    enddo

  case (p_point2point)

    ! Group leaders receive data from other group members and store it in merge_buffer
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, 1:1 ),"vdf/os-vdf-reportfile.f03",291)

      ! Copy local data to merge buffer
      starts(1) = grid%io%tiles( 2, 1, 1 )
      do i1 = 1, vdf % nx_(1)
        write_buffer( starts(1)+i1 ) = real( vdf % f1 ( fc, i1 ), p_single )
      enddo

      ! Get data from other node
      do rank = 1, grid%io%group_size-1

        ! get size of grid tile on source node
        sizes(1) = grid%io%tiles( 1, 1, rank + 1)

        ! get start position of grid tile on merged data grid
        starts(1) = grid%io%tiles( 2, 1, rank + 1)

        ! Create dataype describing the tile to be received
        call mpi_type_create_subarray( 1, grid%io%group_tile( 1, : ), sizes, starts, &
              MPI_ORDER_FORTRAN, MPI_REAL, recv_type, ierr )
        call mpi_type_commit( recv_type, ierr )

        ! notify node that we are ready
        call mpi_isend( ping, 1, MPI_INTEGER, rank, 0, &
                        grid%io%group_comm, ping_handle, ierr )

        ! Receive tile data
        call mpi_recv( write_buffer, 1, recv_type, rank, 1, &
                       grid%io%group_comm, stat, ierr )

        ! Free the datatype
        call mpi_type_free( recv_type, ierr )

        ! Wait for ping message to complete (this is just for cleanup, the ping must
        ! have completed by now)
        call mpi_wait( ping_handle, stat, ierr )

      enddo

    else

      ! Allocate buffer for messaging
      call alloc(comm_buffer, [ vdf % nx_(1) ],"vdf/os-vdf-reportfile.f03",333)

      ! Pack data for sending
      do i1 = 1, vdf % nx_(1)
        comm_buffer( i1 ) = real( vdf % f1 ( fc, i1 ), p_single )
      enddo

      ! receive ping from group leader
      call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                     grid%io%group_comm, stat, ierr )

      ! send data to group leader
      call mpi_send( comm_buffer, size(comm_buffer), MPI_REAL, 0, &
                     1, grid%io%group_comm, ierr )

      ! free message buffer
      call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",349)
    endif

  case (p_gather)

    ! Find counts and displacements for gathering the data
    if ( grid%io%group_rank == 0 ) then
      call alloc(recvcounts, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",356)
      call alloc(displs, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",357)

      data_size = 0

      do i = 1, grid%io%group_size
        recvcounts(i) = grid%io%tiles( 1, 1, i )
        data_size = data_size + recvcounts(i)
      enddo

      displs(1) = 0
      do i = 2, grid%io%group_size
        displs(i) = displs(i-1) + recvcounts(i-1)
      enddo

      call alloc(gather_buffer, [ data_size ],"vdf/os-vdf-reportfile.f03",371)
    else
      recvcounts => null()
      displs => null()
      gather_buffer => null()
    endif

    ! Copy the data to communication buffer
    call alloc(comm_buffer, [ vdf % nx_(1) ],"vdf/os-vdf-reportfile.f03",379)
    do i1 = 1, vdf % nx_(1)
      comm_buffer( i1 ) = real( vdf % f1 ( fc, i1 ), p_single )
    enddo

    ! Gather the data on the group root
    if ( grid%io%group_rank == 0 ) then
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_REAL, &
                        gather_buffer, recvcounts, displs, MPI_REAL, &
                        0, grid%io%group_comm, ierr )
    else
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_REAL, &
                        tmp, tmp, tmp, MPI_REAL, &
                        0, grid%io%group_comm, ierr )
    endif

    ! free the communication buffer
    call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",396)

    ! Rearrange data into properly formed 2D arrays
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, : ),"vdf/os-vdf-reportfile.f03",402)

      ! Transpose data from merge_buffer into data_buffer
      k = 0
      do i = 1, grid%io%group_size
        do i1 = grid%io%tiles( 2, 1, i ) + 1, &
                grid%io%tiles( 2, 1, i ) + grid%io%tiles( 1, 1, i )
          k = k+1
          write_buffer( i1 ) = gather_buffer( k )
        enddo
      enddo

      ! Free allocated memory on group root
      call freemem(gather_buffer,"vdf/os-vdf-reportfile.f03",415)
      call freemem(displs,"vdf/os-vdf-reportfile.f03",416)
      call freemem(recvcounts,"vdf/os-vdf-reportfile.f03",417)

    endif

  end select

  ! Write merged data to disk

  ! Only group leaders participate in writing data to disk
  ! Note than when not using merging all nodes have group_rank = 0
  if ( grid%io%group_rank == 0 ) then

    ! create file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

    ! Open the file for (possible) parallel I/O
    call diagFile % open( p_diag_create, grid%io%comm )

    call diagFile % start_cdset( report%name, 1, grid % g_nx, freal_to_diagtype( p_single ), &
                    parDset, grid % io % chunk_size )

    parChunk % count(1) = grid % io % group_tile( 1 , 1 )
    parChunk % start(1) = grid % io % group_tile( 2 , 1 )
    parChunk % stride(1) = 1
    parChunk % data = c_loc( write_buffer )

    ! Write the data
    call diagFile % write_par_cdset( parDset, parChunk )

    ! Close dataset
    call diagFile % end_cdset( parDset )

    ! Close file
    call diagFile % close( )

    deallocate( diagFile )

    ! free merge data buffer
    call freemem(write_buffer,"vdf/os-vdf-reportfile.f03",455)

  endif

end subroutine write_vdf_1d_r4


!---------------------------------------------------------------------------------------------------
! Merge / write a 2D dataset
!---------------------------------------------------------------------------------------------------
subroutine write_vdf_2d_r4( vdf, report, fc, g_space, grid )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  use m_parameters
  use m_diagfile

  implicit none

  class( t_vdf ), intent(in) :: vdf ! vdf object
  type( t_vdf_report ), intent(in) :: report
  integer, intent(in) :: fc ! field component to report
  type( t_space ), intent(in) :: g_space ! spatial information
  class( t_grid ), intent(in) :: grid


  ! Buffer for merging data - only group leaders will allocate this
  real( p_single ), dimension(:,:), pointer :: write_buffer => null()

  ! Buffer for MPI communications - all nodes allocate this
  real( p_single ), dimension(:,:), pointer :: comm_buffer

  real( p_single ), dimension(:), pointer :: gather_buffer

  ! MPI datatype for receiving data on group root nodes
  integer :: recv_type

  ! Sizes and offsets for MPI datatype
  integer, dimension(2) :: sizes, starts

  ! Rank of node sending data
  integer :: rank

  ! Status variable for messages
  integer, dimension( MPI_STATUS_SIZE ) :: stat

  ! Ping variables
  integer :: ping, ping_handle

  ! Total data size being received
  integer :: data_size

  ! Sizes and displacements of data from each node
  integer, dimension(:), pointer :: recvcounts, displs

  ! Diagnostic file object
  class(t_diag_file), allocatable :: diagFile

  integer :: i, i1, i2, k, ierr

  type( t_diag_dataset ) :: parDset
  type( t_diag_chunk ) :: parChunk

  ! Dummy receive buffer for MPI_GATHERV
  integer, dimension(1) :: tmp

  ! Merge data if necessary
  select case ( grid%io%merge_type )
  case (p_none)

    ! No merging, just copy local data to output buffer
    call alloc(write_buffer, vdf % nx_(1:2),"vdf/os-vdf-reportfile.f03",528)

    do i2 = 1, vdf % nx_(2)
      do i1 = 1, vdf % nx_(1)
        write_buffer( i1, i2 ) = real( vdf % f2 ( fc, i1, i2 ), p_single )
      enddo
    enddo

  case (p_point2point)

    ! Group leaders receive data from other group members and store it in merge_buffer
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, 1:2 ),"vdf/os-vdf-reportfile.f03",542)

      ! Copy local data to merge buffer
      starts(1) = grid%io%tiles( 2, 1, 1 )
      starts(2) = grid%io%tiles( 2, 2, 1 )
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          write_buffer( starts(1)+i1, starts(2)+i2 ) = real( vdf % f2 ( fc, i1, i2 ), p_single )
        enddo
      enddo

      ! Get data from other node
      do rank = 1, grid%io%group_size-1

        ! get size of grid tile on source node
        sizes(1) = grid%io%tiles( 1, 1, rank + 1)
        sizes(2) = grid%io%tiles( 1, 2, rank + 1)

        ! get start position of grid tile on merged data grid
        starts(1) = grid%io%tiles( 2, 1, rank + 1)
        starts(2) = grid%io%tiles( 2, 2, rank + 1)

        ! Create dataype describing the tile to be received
        call mpi_type_create_subarray( 2, grid%io%group_tile( 1, 1:2 ), sizes, starts, &
              MPI_ORDER_FORTRAN, MPI_REAL, recv_type, ierr )
        call mpi_type_commit( recv_type, ierr )

        ! notify node that we are ready
        call mpi_isend( ping, 1, MPI_INTEGER, rank, 0, &
                        grid%io%group_comm, ping_handle, ierr )

        ! Receive tile data
        call mpi_recv( write_buffer, 1, recv_type, rank, 1, &
                       grid%io%group_comm, stat, ierr )

        ! Free the datatype
        call mpi_type_free( recv_type, ierr )

        ! Wait for ping message to complete (this is just for cleanup, the ping must
        ! have completed by now)
        call mpi_wait( ping_handle, stat, ierr )

      enddo

    else

      ! Allocate buffer for messaging
      call alloc(comm_buffer, vdf % nx_(1:2),"vdf/os-vdf-reportfile.f03",589)

      ! Pack data for sending
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          comm_buffer( i1, i2 ) = real( vdf % f2 ( fc, i1, i2 ), p_single )
        enddo
      enddo

      ! receive ping from group leader
      call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                     grid%io%group_comm, stat, ierr )

      ! send data to group leader
      call mpi_send( comm_buffer, size(comm_buffer), MPI_REAL, 0, &
                     1, grid%io%group_comm, ierr )

      ! free message buffer
      call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",607)
    endif

  case (p_gather)

    ! Find counts and displacements for gathering the data
    if ( grid%io%group_rank == 0 ) then
      call alloc(recvcounts, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",614)
      call alloc(displs, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",615)

      data_size = 0

      do i = 1, grid%io%group_size
        recvcounts(i) = grid%io%tiles( 1, 1, i ) * &
                        grid%io%tiles( 1, 2, i )
        data_size = data_size + recvcounts(i)
      enddo

      displs(1) = 0
      do i = 2, grid%io%group_size
        displs(i) = displs(i-1) + recvcounts(i-1)
      enddo

      call alloc(gather_buffer, [ data_size ],"vdf/os-vdf-reportfile.f03",630)
    else
      recvcounts => null()
      displs => null()
      gather_buffer => null()
    endif

    ! Copy the data to communication buffer
    call alloc(comm_buffer, vdf % nx_(1:2),"vdf/os-vdf-reportfile.f03",638)
    do i2 = 1, vdf % nx_(2)
      do i1 = 1, vdf % nx_(1)
        comm_buffer( i1, i2 ) = real( vdf % f2 ( fc, i1, i2 ), p_single )
      enddo
    enddo

    ! Gather the data on the group root
    if ( grid%io%group_rank == 0 ) then
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_REAL, &
                        gather_buffer, recvcounts, displs, MPI_REAL, &
                        0, grid%io%group_comm, ierr )
    else
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_REAL, &
                        tmp, tmp, tmp, MPI_REAL, &
                        0, grid%io%group_comm, ierr )
    endif

    ! free the communication buffer
    call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",657)

    ! Rearrange data into properly formed 2D arrays
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, : ),"vdf/os-vdf-reportfile.f03",663)

      ! Transpose data from merge_buffer into data_buffer
      k = 0
      do i = 1, grid%io%group_size
        do i2 = grid%io%tiles( 2, 2, i ) + 1, &
                grid%io%tiles( 2, 2, i ) + grid%io%tiles( 1, 2, i )
          do i1 = grid%io%tiles( 2, 1, i ) + 1, &
                  grid%io%tiles( 2, 1, i ) + grid%io%tiles( 1, 1, i )
            k = k+1
            write_buffer( i1, i2 ) = gather_buffer( k )
          enddo
        enddo
      enddo

      ! Free allocated memory on group root
      call freemem(gather_buffer,"vdf/os-vdf-reportfile.f03",679)
      call freemem(displs,"vdf/os-vdf-reportfile.f03",680)
      call freemem(recvcounts,"vdf/os-vdf-reportfile.f03",681)

    endif

  end select

  ! Write merged data to disk

  ! Only group leaders participate in writing data to disk
  ! Note than when not using merging all nodes have group_rank = 0
  if ( grid%io%group_rank == 0 ) then

    ! create file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

    ! Open the file for (possible) parallel I/O
    call diagFile % open( p_diag_create, grid%io%comm )

    call diagFile % start_cdset( report%name, 2, grid % g_nx, freal_to_diagtype( p_single ), &
                    parDset, grid % io % chunk_size )

    do i = 1, 2
      parChunk % count(i) = grid % io % group_tile( 1 , i )
      parChunk % start(i) = grid % io % group_tile( 2 , i )
      parChunk % stride(i) = 1
    enddo
    parChunk % data = c_loc( write_buffer )

    ! Write the data
    call diagFile % write_par_cdset( parDset, parChunk )

    ! Close dataset
    call diagFile % end_cdset( parDset )

    ! Close file
    call diagFile % close( )

    deallocate( diagFile )

    ! free merge data buffer
    call freemem(write_buffer,"vdf/os-vdf-reportfile.f03",721)

  endif

end subroutine write_vdf_2d_r4


!---------------------------------------------------------------------------------------------------
! Merge / write a 3D dataset
!---------------------------------------------------------------------------------------------------
subroutine write_vdf_3d_r4( vdf, report, fc, g_space, grid )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  use m_parameters
  use m_diagfile

  implicit none

  class( t_vdf ), intent(in) :: vdf ! vdf object
  type( t_vdf_report ), intent(in) :: report
  integer, intent(in) :: fc ! field component to report
  type( t_space ), intent(in) :: g_space ! spatial information
  class( t_grid ), intent(in) :: grid


  ! Buffer for merging data - only group leaders will allocate this
  real( p_single ), dimension(:,:,:), pointer :: write_buffer

  ! Buffer for MPI communications - all nodes allocate this
  real( p_single ), dimension(:,:,:), pointer :: comm_buffer

  real( p_single ), dimension(:), pointer :: gather_buffer

  ! MPI datatype for receiving data on group root nodes
  integer :: recv_type

  ! Sizes and offsets for MPI datatype
  integer, dimension(3) :: sizes, starts

  ! Rank of node sending data
  integer :: rank

  ! Status variable for messages
  integer, dimension( MPI_STATUS_SIZE ) :: stat

  ! Ping variables
  integer :: ping, ping_handle

  ! Total data size being received
  integer :: data_size

  ! Sizes and displacements of data from each node
  integer, dimension(:), pointer :: recvcounts, displs

  ! Diagnostic file object
  class(t_diag_file), allocatable :: diagFile

  integer :: i, i1, i2, i3, k, ierr

  type( t_diag_dataset ) :: parDset
  type( t_diag_chunk ) :: parChunk

  ! Dummy receive buffer for MPI_GATHERV
  integer, dimension(1) :: tmp

  ! Merge data if necessary
  select case ( grid%io%merge_type )
  case (p_none)

    ! No merging, just copy local data to output buffer
    call alloc(write_buffer, vdf % nx_(1:3),"vdf/os-vdf-reportfile.f03",794)

    do i3 = 1, vdf % nx_(3)
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          write_buffer( i1, i2, i3 ) = real( vdf % f3 ( fc, i1, i2, i3 ), p_single )
        enddo
      enddo
    enddo

  case (p_point2point)

    ! Group leaders receive data from other group members and store it in merge_buffer
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, 1:3 ),"vdf/os-vdf-reportfile.f03",810)

      ! Copy local data to merge buffer
      starts(1) = grid%io%tiles( 2, 1, 1 )
      starts(2) = grid%io%tiles( 2, 2, 1 )
      starts(3) = grid%io%tiles( 2, 3, 1 )

      do i3 = 1, vdf % nx_(3)
        do i2 = 1, vdf % nx_(2)
          do i1 = 1, vdf % nx_(1)
            write_buffer( starts(1)+i1, starts(2)+i2, starts(3)+i3 ) = &
              real( vdf % f3 ( fc, i1, i2, i3 ), p_single )
          enddo
        enddo
      enddo

      ! Get data from other node
      do rank = 1, grid%io%group_size-1

        ! get size of grid tile on source node
        sizes(1) = grid%io%tiles( 1, 1, rank + 1)
        sizes(2) = grid%io%tiles( 1, 2, rank + 1)
        sizes(3) = grid%io%tiles( 1, 3, rank + 1)

        ! get start position of grid tile on merged data grid
        starts(1) = grid%io%tiles( 2, 1, rank + 1)
        starts(2) = grid%io%tiles( 2, 2, rank + 1)
        starts(3) = grid%io%tiles( 2, 3, rank + 1)

        ! Create dataype describing the tile to be received
        call mpi_type_create_subarray( 3, grid%io%group_tile( 1, 1:3 ), sizes, starts, &
              MPI_ORDER_FORTRAN, MPI_REAL, recv_type, ierr )
        call mpi_type_commit( recv_type, ierr )

        ! notify node that we are ready
        call mpi_isend( ping, 1, MPI_INTEGER, rank, 0, &
                        grid%io%group_comm, ping_handle, ierr )

        ! Receive tile data
        call mpi_recv( write_buffer, 1, recv_type, rank, 1, &
                       grid%io%group_comm, stat, ierr )

        ! Free the datatype
        call mpi_type_free( recv_type, ierr )

        ! Wait for ping message to complete (this is just for cleanup, the ping must
        ! have completed by now)
        call mpi_wait( ping_handle, stat, ierr )

      enddo

    else

      ! Allocate buffer for messaging
      call alloc(comm_buffer, vdf % nx_(:),"vdf/os-vdf-reportfile.f03",864)

      ! Pack data for sending
      do i3 = 1, vdf % nx_(3)
        do i2 = 1, vdf % nx_(2)
          do i1 = 1, vdf % nx_(1)
            comm_buffer( i1, i2, i3 ) = real( vdf % f3 ( fc, i1, i2, i3 ), p_single )
          enddo
        enddo
      enddo

      ! receive ping from group leader
      call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                     grid%io%group_comm, stat, ierr )

      ! send data to group leader
      call mpi_send( comm_buffer, size(comm_buffer), MPI_REAL, 0, &
                     1, grid%io%group_comm, ierr )

      ! free message buffer
      call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",884)
    endif

  case (p_gather)

    ! Find counts and displacements for gathering the data
    if ( grid%io%group_rank == 0 ) then
      call alloc(recvcounts, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",891)
      call alloc(displs, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",892)

      data_size = 0

      do i = 1, grid%io%group_size
        recvcounts(i) = grid%io%tiles( 1, 1, i ) * &
                        grid%io%tiles( 1, 2, i ) * &
                        grid%io%tiles( 1, 3, i )
        data_size = data_size + recvcounts(i)
      enddo

      displs(1) = 0
      do i = 2, grid%io%group_size
        displs(i) = displs(i-1) + recvcounts(i-1)
      enddo

      call alloc(gather_buffer, [ data_size ],"vdf/os-vdf-reportfile.f03",908)
    else
      recvcounts => null()
      displs => null()
      gather_buffer => null()
    endif

    ! Copy the data to communication buffer
    call alloc(comm_buffer, vdf % nx_(1:3),"vdf/os-vdf-reportfile.f03",916)
    do i3 = 1, vdf % nx_(3)
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          comm_buffer( i1, i2, i3 ) = real( vdf % f3 ( fc, i1, i2, i3 ), p_single )
        enddo
      enddo
    enddo

    ! Gather the data on the group root
    if ( grid%io%group_rank == 0 ) then
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_REAL, &
                        gather_buffer, recvcounts, displs, MPI_REAL, &
                        0, grid%io%group_comm, ierr )
    else
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_REAL, &
                        tmp, tmp, tmp, MPI_REAL, &
                        0, grid%io%group_comm, ierr )
    endif

    ! free the communication buffer
    call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",937)

    ! Rearrange data into properly formed 2D arrays
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, : ),"vdf/os-vdf-reportfile.f03",943)

      ! Transpose data from merge_buffer into data_buffer
      k = 0
      do i = 1, grid%io%group_size
        do i3 = grid%io%tiles( 2, 3, i ) + 1, &
                grid%io%tiles( 2, 3, i ) + grid%io%tiles( 1, 3, i )
          do i2 = grid%io%tiles( 2, 2, i ) + 1, &
                  grid%io%tiles( 2, 2, i ) + grid%io%tiles( 1, 2, i )
            do i1 = grid%io%tiles( 2, 1, i ) + 1, &
                    grid%io%tiles( 2, 1, i ) + grid%io%tiles( 1, 1, i )
              k = k+1
              write_buffer( i1, i2, i3 ) = gather_buffer( k )
            enddo
          enddo
        enddo
      enddo

      ! Free allocated memory on group root
      call freemem(gather_buffer,"vdf/os-vdf-reportfile.f03",962)
      call freemem(displs,"vdf/os-vdf-reportfile.f03",963)
      call freemem(recvcounts,"vdf/os-vdf-reportfile.f03",964)

    endif

  end select

! print *, "[",mpi_node(),"] ",'Writing merged data...'

  ! Write merged data to disk

  ! Only group leaders participate in writing data to disk
  ! Note than when not using merging all nodes have group_rank = 0
  if ( grid%io%group_rank == 0 ) then

    ! create file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

! print *, "[",mpi_node(),"] ",'Opening file...'

    ! Open the file for (possible) parallel I/O
    call diagFile % open( p_diag_create, grid % io % comm )

! print *, "[",mpi_node(),"] ",'Starting CDSET...'

    call diagFile % start_cdset( report%name, 3, grid % g_nx, freal_to_diagtype( p_single ), &
                    parDset, grid % io % chunk_size )

    do i = 1, 3
      parChunk % count(i) = grid % io % group_tile( 1 , i )
      parChunk % start(i) = grid % io % group_tile( 2 , i )
      parChunk % stride(i) = 1
    enddo
    parChunk % data = c_loc( write_buffer )

! print *, "[",mpi_node(),"] ",'Writing CDSET...'

    ! Write the data
    call diagFile % write_par_cdset( parDset, parChunk )

    ! Close dataset
    call diagFile % end_cdset( parDset )

    ! Close file
    call diagFile % close( )

    deallocate( diagFile )

    ! free merge data buffer
    call freemem(write_buffer,"vdf/os-vdf-reportfile.f03",1012)

  endif

end subroutine write_vdf_3d_r4

!---------------------------------------------------------------------------------------------------
! End of template Functions, do not change below
!---------------------------------------------------------------------------------------------------
# 129 "vdf/os-vdf-reportfile.f03" 2

! double precision real



# 1 "vdf/os-vdf-reportfile.f03" 1
! if is defined then read the template definition at the end of the file
# 208 "vdf/os-vdf-reportfile.f03"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Merge / write a 1D dataset
!---------------------------------------------------------------------------------------------------
subroutine write_vdf_1d_r8( vdf, report, fc, g_space, grid )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  use m_parameters
  use zdf
  use m_diagfile

  implicit none

  class( t_vdf ), intent(in) :: vdf
  type( t_vdf_report ), intent(in) :: report
  integer, intent(in) :: fc
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid

  ! Buffer for merging data - only group leaders will allocate this
  real( p_double ), dimension(:), pointer :: write_buffer

  ! Buffer for MPI communications - all nodes allocate this
  real( p_double ), dimension(:), pointer :: comm_buffer

  real( p_double ), dimension(:), pointer :: gather_buffer

  ! MPI datatype for receiving data on group root nodes
  integer :: recv_type

  ! Sizes and offsets for MPI datatype
  integer, dimension(1) :: sizes, starts

  ! Rank of node sending data
  integer :: rank

  ! Status variable for messages
  integer, dimension( MPI_STATUS_SIZE ) :: stat

  ! Ping variables
  integer :: ping, ping_handle

  ! Total data size being received
  integer :: data_size

  ! Sizes and displacements of data from each node
  integer, dimension(:), pointer :: recvcounts, displs

  ! Diagnostic file object
  class(t_diag_file), allocatable :: diagFile

  integer :: i, i1, k, ierr

  type( t_diag_dataset ) :: parDset
  type( t_diag_chunk ) :: parChunk

  ! Dummy receive buffer for MPI_GATHERV
  integer, dimension(1) :: tmp

  ! Merge data if necessary
  select case ( grid%io%merge_type )
  case (p_none)

    ! No merging, just copy local data to output buffer
    call alloc(write_buffer, [ vdf % nx_(1) ],"vdf/os-vdf-reportfile.f03",279)

    do i1 = 1, vdf % nx_(1)
      write_buffer( i1 ) = real( vdf % f1 ( fc, i1 ), p_double )
    enddo

  case (p_point2point)

    ! Group leaders receive data from other group members and store it in merge_buffer
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, 1:1 ),"vdf/os-vdf-reportfile.f03",291)

      ! Copy local data to merge buffer
      starts(1) = grid%io%tiles( 2, 1, 1 )
      do i1 = 1, vdf % nx_(1)
        write_buffer( starts(1)+i1 ) = real( vdf % f1 ( fc, i1 ), p_double )
      enddo

      ! Get data from other node
      do rank = 1, grid%io%group_size-1

        ! get size of grid tile on source node
        sizes(1) = grid%io%tiles( 1, 1, rank + 1)

        ! get start position of grid tile on merged data grid
        starts(1) = grid%io%tiles( 2, 1, rank + 1)

        ! Create dataype describing the tile to be received
        call mpi_type_create_subarray( 1, grid%io%group_tile( 1, : ), sizes, starts, &
              MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_type, ierr )
        call mpi_type_commit( recv_type, ierr )

        ! notify node that we are ready
        call mpi_isend( ping, 1, MPI_INTEGER, rank, 0, &
                        grid%io%group_comm, ping_handle, ierr )

        ! Receive tile data
        call mpi_recv( write_buffer, 1, recv_type, rank, 1, &
                       grid%io%group_comm, stat, ierr )

        ! Free the datatype
        call mpi_type_free( recv_type, ierr )

        ! Wait for ping message to complete (this is just for cleanup, the ping must
        ! have completed by now)
        call mpi_wait( ping_handle, stat, ierr )

      enddo

    else

      ! Allocate buffer for messaging
      call alloc(comm_buffer, [ vdf % nx_(1) ],"vdf/os-vdf-reportfile.f03",333)

      ! Pack data for sending
      do i1 = 1, vdf % nx_(1)
        comm_buffer( i1 ) = real( vdf % f1 ( fc, i1 ), p_double )
      enddo

      ! receive ping from group leader
      call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                     grid%io%group_comm, stat, ierr )

      ! send data to group leader
      call mpi_send( comm_buffer, size(comm_buffer), MPI_DOUBLE_PRECISION, 0, &
                     1, grid%io%group_comm, ierr )

      ! free message buffer
      call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",349)
    endif

  case (p_gather)

    ! Find counts and displacements for gathering the data
    if ( grid%io%group_rank == 0 ) then
      call alloc(recvcounts, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",356)
      call alloc(displs, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",357)

      data_size = 0

      do i = 1, grid%io%group_size
        recvcounts(i) = grid%io%tiles( 1, 1, i )
        data_size = data_size + recvcounts(i)
      enddo

      displs(1) = 0
      do i = 2, grid%io%group_size
        displs(i) = displs(i-1) + recvcounts(i-1)
      enddo

      call alloc(gather_buffer, [ data_size ],"vdf/os-vdf-reportfile.f03",371)
    else
      recvcounts => null()
      displs => null()
      gather_buffer => null()
    endif

    ! Copy the data to communication buffer
    call alloc(comm_buffer, [ vdf % nx_(1) ],"vdf/os-vdf-reportfile.f03",379)
    do i1 = 1, vdf % nx_(1)
      comm_buffer( i1 ) = real( vdf % f1 ( fc, i1 ), p_double )
    enddo

    ! Gather the data on the group root
    if ( grid%io%group_rank == 0 ) then
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_DOUBLE_PRECISION, &
                        gather_buffer, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                        0, grid%io%group_comm, ierr )
    else
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_DOUBLE_PRECISION, &
                        tmp, tmp, tmp, MPI_DOUBLE_PRECISION, &
                        0, grid%io%group_comm, ierr )
    endif

    ! free the communication buffer
    call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",396)

    ! Rearrange data into properly formed 2D arrays
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, : ),"vdf/os-vdf-reportfile.f03",402)

      ! Transpose data from merge_buffer into data_buffer
      k = 0
      do i = 1, grid%io%group_size
        do i1 = grid%io%tiles( 2, 1, i ) + 1, &
                grid%io%tiles( 2, 1, i ) + grid%io%tiles( 1, 1, i )
          k = k+1
          write_buffer( i1 ) = gather_buffer( k )
        enddo
      enddo

      ! Free allocated memory on group root
      call freemem(gather_buffer,"vdf/os-vdf-reportfile.f03",415)
      call freemem(displs,"vdf/os-vdf-reportfile.f03",416)
      call freemem(recvcounts,"vdf/os-vdf-reportfile.f03",417)

    endif

  end select

  ! Write merged data to disk

  ! Only group leaders participate in writing data to disk
  ! Note than when not using merging all nodes have group_rank = 0
  if ( grid%io%group_rank == 0 ) then

    ! create file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

    ! Open the file for (possible) parallel I/O
    call diagFile % open( p_diag_create, grid%io%comm )

    call diagFile % start_cdset( report%name, 1, grid % g_nx, freal_to_diagtype( p_double ), &
                    parDset, grid % io % chunk_size )

    parChunk % count(1) = grid % io % group_tile( 1 , 1 )
    parChunk % start(1) = grid % io % group_tile( 2 , 1 )
    parChunk % stride(1) = 1
    parChunk % data = c_loc( write_buffer )

    ! Write the data
    call diagFile % write_par_cdset( parDset, parChunk )

    ! Close dataset
    call diagFile % end_cdset( parDset )

    ! Close file
    call diagFile % close( )

    deallocate( diagFile )

    ! free merge data buffer
    call freemem(write_buffer,"vdf/os-vdf-reportfile.f03",455)

  endif

end subroutine write_vdf_1d_r8


!---------------------------------------------------------------------------------------------------
! Merge / write a 2D dataset
!---------------------------------------------------------------------------------------------------
subroutine write_vdf_2d_r8( vdf, report, fc, g_space, grid )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  use m_parameters
  use m_diagfile

  implicit none

  class( t_vdf ), intent(in) :: vdf ! vdf object
  type( t_vdf_report ), intent(in) :: report
  integer, intent(in) :: fc ! field component to report
  type( t_space ), intent(in) :: g_space ! spatial information
  class( t_grid ), intent(in) :: grid


  ! Buffer for merging data - only group leaders will allocate this
  real( p_double ), dimension(:,:), pointer :: write_buffer => null()

  ! Buffer for MPI communications - all nodes allocate this
  real( p_double ), dimension(:,:), pointer :: comm_buffer

  real( p_double ), dimension(:), pointer :: gather_buffer

  ! MPI datatype for receiving data on group root nodes
  integer :: recv_type

  ! Sizes and offsets for MPI datatype
  integer, dimension(2) :: sizes, starts

  ! Rank of node sending data
  integer :: rank

  ! Status variable for messages
  integer, dimension( MPI_STATUS_SIZE ) :: stat

  ! Ping variables
  integer :: ping, ping_handle

  ! Total data size being received
  integer :: data_size

  ! Sizes and displacements of data from each node
  integer, dimension(:), pointer :: recvcounts, displs

  ! Diagnostic file object
  class(t_diag_file), allocatable :: diagFile

  integer :: i, i1, i2, k, ierr

  type( t_diag_dataset ) :: parDset
  type( t_diag_chunk ) :: parChunk

  ! Dummy receive buffer for MPI_GATHERV
  integer, dimension(1) :: tmp

  ! Merge data if necessary
  select case ( grid%io%merge_type )
  case (p_none)

    ! No merging, just copy local data to output buffer
    call alloc(write_buffer, vdf % nx_(1:2),"vdf/os-vdf-reportfile.f03",528)

    do i2 = 1, vdf % nx_(2)
      do i1 = 1, vdf % nx_(1)
        write_buffer( i1, i2 ) = real( vdf % f2 ( fc, i1, i2 ), p_double )
      enddo
    enddo

  case (p_point2point)

    ! Group leaders receive data from other group members and store it in merge_buffer
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, 1:2 ),"vdf/os-vdf-reportfile.f03",542)

      ! Copy local data to merge buffer
      starts(1) = grid%io%tiles( 2, 1, 1 )
      starts(2) = grid%io%tiles( 2, 2, 1 )
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          write_buffer( starts(1)+i1, starts(2)+i2 ) = real( vdf % f2 ( fc, i1, i2 ), p_double )
        enddo
      enddo

      ! Get data from other node
      do rank = 1, grid%io%group_size-1

        ! get size of grid tile on source node
        sizes(1) = grid%io%tiles( 1, 1, rank + 1)
        sizes(2) = grid%io%tiles( 1, 2, rank + 1)

        ! get start position of grid tile on merged data grid
        starts(1) = grid%io%tiles( 2, 1, rank + 1)
        starts(2) = grid%io%tiles( 2, 2, rank + 1)

        ! Create dataype describing the tile to be received
        call mpi_type_create_subarray( 2, grid%io%group_tile( 1, 1:2 ), sizes, starts, &
              MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_type, ierr )
        call mpi_type_commit( recv_type, ierr )

        ! notify node that we are ready
        call mpi_isend( ping, 1, MPI_INTEGER, rank, 0, &
                        grid%io%group_comm, ping_handle, ierr )

        ! Receive tile data
        call mpi_recv( write_buffer, 1, recv_type, rank, 1, &
                       grid%io%group_comm, stat, ierr )

        ! Free the datatype
        call mpi_type_free( recv_type, ierr )

        ! Wait for ping message to complete (this is just for cleanup, the ping must
        ! have completed by now)
        call mpi_wait( ping_handle, stat, ierr )

      enddo

    else

      ! Allocate buffer for messaging
      call alloc(comm_buffer, vdf % nx_(1:2),"vdf/os-vdf-reportfile.f03",589)

      ! Pack data for sending
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          comm_buffer( i1, i2 ) = real( vdf % f2 ( fc, i1, i2 ), p_double )
        enddo
      enddo

      ! receive ping from group leader
      call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                     grid%io%group_comm, stat, ierr )

      ! send data to group leader
      call mpi_send( comm_buffer, size(comm_buffer), MPI_DOUBLE_PRECISION, 0, &
                     1, grid%io%group_comm, ierr )

      ! free message buffer
      call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",607)
    endif

  case (p_gather)

    ! Find counts and displacements for gathering the data
    if ( grid%io%group_rank == 0 ) then
      call alloc(recvcounts, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",614)
      call alloc(displs, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",615)

      data_size = 0

      do i = 1, grid%io%group_size
        recvcounts(i) = grid%io%tiles( 1, 1, i ) * &
                        grid%io%tiles( 1, 2, i )
        data_size = data_size + recvcounts(i)
      enddo

      displs(1) = 0
      do i = 2, grid%io%group_size
        displs(i) = displs(i-1) + recvcounts(i-1)
      enddo

      call alloc(gather_buffer, [ data_size ],"vdf/os-vdf-reportfile.f03",630)
    else
      recvcounts => null()
      displs => null()
      gather_buffer => null()
    endif

    ! Copy the data to communication buffer
    call alloc(comm_buffer, vdf % nx_(1:2),"vdf/os-vdf-reportfile.f03",638)
    do i2 = 1, vdf % nx_(2)
      do i1 = 1, vdf % nx_(1)
        comm_buffer( i1, i2 ) = real( vdf % f2 ( fc, i1, i2 ), p_double )
      enddo
    enddo

    ! Gather the data on the group root
    if ( grid%io%group_rank == 0 ) then
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_DOUBLE_PRECISION, &
                        gather_buffer, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                        0, grid%io%group_comm, ierr )
    else
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_DOUBLE_PRECISION, &
                        tmp, tmp, tmp, MPI_DOUBLE_PRECISION, &
                        0, grid%io%group_comm, ierr )
    endif

    ! free the communication buffer
    call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",657)

    ! Rearrange data into properly formed 2D arrays
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, : ),"vdf/os-vdf-reportfile.f03",663)

      ! Transpose data from merge_buffer into data_buffer
      k = 0
      do i = 1, grid%io%group_size
        do i2 = grid%io%tiles( 2, 2, i ) + 1, &
                grid%io%tiles( 2, 2, i ) + grid%io%tiles( 1, 2, i )
          do i1 = grid%io%tiles( 2, 1, i ) + 1, &
                  grid%io%tiles( 2, 1, i ) + grid%io%tiles( 1, 1, i )
            k = k+1
            write_buffer( i1, i2 ) = gather_buffer( k )
          enddo
        enddo
      enddo

      ! Free allocated memory on group root
      call freemem(gather_buffer,"vdf/os-vdf-reportfile.f03",679)
      call freemem(displs,"vdf/os-vdf-reportfile.f03",680)
      call freemem(recvcounts,"vdf/os-vdf-reportfile.f03",681)

    endif

  end select

  ! Write merged data to disk

  ! Only group leaders participate in writing data to disk
  ! Note than when not using merging all nodes have group_rank = 0
  if ( grid%io%group_rank == 0 ) then

    ! create file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

    ! Open the file for (possible) parallel I/O
    call diagFile % open( p_diag_create, grid%io%comm )

    call diagFile % start_cdset( report%name, 2, grid % g_nx, freal_to_diagtype( p_double ), &
                    parDset, grid % io % chunk_size )

    do i = 1, 2
      parChunk % count(i) = grid % io % group_tile( 1 , i )
      parChunk % start(i) = grid % io % group_tile( 2 , i )
      parChunk % stride(i) = 1
    enddo
    parChunk % data = c_loc( write_buffer )

    ! Write the data
    call diagFile % write_par_cdset( parDset, parChunk )

    ! Close dataset
    call diagFile % end_cdset( parDset )

    ! Close file
    call diagFile % close( )

    deallocate( diagFile )

    ! free merge data buffer
    call freemem(write_buffer,"vdf/os-vdf-reportfile.f03",721)

  endif

end subroutine write_vdf_2d_r8


!---------------------------------------------------------------------------------------------------
! Merge / write a 3D dataset
!---------------------------------------------------------------------------------------------------
subroutine write_vdf_3d_r8( vdf, report, fc, g_space, grid )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  use m_parameters
  use m_diagfile

  implicit none

  class( t_vdf ), intent(in) :: vdf ! vdf object
  type( t_vdf_report ), intent(in) :: report
  integer, intent(in) :: fc ! field component to report
  type( t_space ), intent(in) :: g_space ! spatial information
  class( t_grid ), intent(in) :: grid


  ! Buffer for merging data - only group leaders will allocate this
  real( p_double ), dimension(:,:,:), pointer :: write_buffer

  ! Buffer for MPI communications - all nodes allocate this
  real( p_double ), dimension(:,:,:), pointer :: comm_buffer

  real( p_double ), dimension(:), pointer :: gather_buffer

  ! MPI datatype for receiving data on group root nodes
  integer :: recv_type

  ! Sizes and offsets for MPI datatype
  integer, dimension(3) :: sizes, starts

  ! Rank of node sending data
  integer :: rank

  ! Status variable for messages
  integer, dimension( MPI_STATUS_SIZE ) :: stat

  ! Ping variables
  integer :: ping, ping_handle

  ! Total data size being received
  integer :: data_size

  ! Sizes and displacements of data from each node
  integer, dimension(:), pointer :: recvcounts, displs

  ! Diagnostic file object
  class(t_diag_file), allocatable :: diagFile

  integer :: i, i1, i2, i3, k, ierr

  type( t_diag_dataset ) :: parDset
  type( t_diag_chunk ) :: parChunk

  ! Dummy receive buffer for MPI_GATHERV
  integer, dimension(1) :: tmp

  ! Merge data if necessary
  select case ( grid%io%merge_type )
  case (p_none)

    ! No merging, just copy local data to output buffer
    call alloc(write_buffer, vdf % nx_(1:3),"vdf/os-vdf-reportfile.f03",794)

    do i3 = 1, vdf % nx_(3)
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          write_buffer( i1, i2, i3 ) = real( vdf % f3 ( fc, i1, i2, i3 ), p_double )
        enddo
      enddo
    enddo

  case (p_point2point)

    ! Group leaders receive data from other group members and store it in merge_buffer
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, 1:3 ),"vdf/os-vdf-reportfile.f03",810)

      ! Copy local data to merge buffer
      starts(1) = grid%io%tiles( 2, 1, 1 )
      starts(2) = grid%io%tiles( 2, 2, 1 )
      starts(3) = grid%io%tiles( 2, 3, 1 )

      do i3 = 1, vdf % nx_(3)
        do i2 = 1, vdf % nx_(2)
          do i1 = 1, vdf % nx_(1)
            write_buffer( starts(1)+i1, starts(2)+i2, starts(3)+i3 ) = &
              real( vdf % f3 ( fc, i1, i2, i3 ), p_double )
          enddo
        enddo
      enddo

      ! Get data from other node
      do rank = 1, grid%io%group_size-1

        ! get size of grid tile on source node
        sizes(1) = grid%io%tiles( 1, 1, rank + 1)
        sizes(2) = grid%io%tiles( 1, 2, rank + 1)
        sizes(3) = grid%io%tiles( 1, 3, rank + 1)

        ! get start position of grid tile on merged data grid
        starts(1) = grid%io%tiles( 2, 1, rank + 1)
        starts(2) = grid%io%tiles( 2, 2, rank + 1)
        starts(3) = grid%io%tiles( 2, 3, rank + 1)

        ! Create dataype describing the tile to be received
        call mpi_type_create_subarray( 3, grid%io%group_tile( 1, 1:3 ), sizes, starts, &
              MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_type, ierr )
        call mpi_type_commit( recv_type, ierr )

        ! notify node that we are ready
        call mpi_isend( ping, 1, MPI_INTEGER, rank, 0, &
                        grid%io%group_comm, ping_handle, ierr )

        ! Receive tile data
        call mpi_recv( write_buffer, 1, recv_type, rank, 1, &
                       grid%io%group_comm, stat, ierr )

        ! Free the datatype
        call mpi_type_free( recv_type, ierr )

        ! Wait for ping message to complete (this is just for cleanup, the ping must
        ! have completed by now)
        call mpi_wait( ping_handle, stat, ierr )

      enddo

    else

      ! Allocate buffer for messaging
      call alloc(comm_buffer, vdf % nx_(:),"vdf/os-vdf-reportfile.f03",864)

      ! Pack data for sending
      do i3 = 1, vdf % nx_(3)
        do i2 = 1, vdf % nx_(2)
          do i1 = 1, vdf % nx_(1)
            comm_buffer( i1, i2, i3 ) = real( vdf % f3 ( fc, i1, i2, i3 ), p_double )
          enddo
        enddo
      enddo

      ! receive ping from group leader
      call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                     grid%io%group_comm, stat, ierr )

      ! send data to group leader
      call mpi_send( comm_buffer, size(comm_buffer), MPI_DOUBLE_PRECISION, 0, &
                     1, grid%io%group_comm, ierr )

      ! free message buffer
      call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",884)
    endif

  case (p_gather)

    ! Find counts and displacements for gathering the data
    if ( grid%io%group_rank == 0 ) then
      call alloc(recvcounts, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",891)
      call alloc(displs, [ grid % io % group_size ],"vdf/os-vdf-reportfile.f03",892)

      data_size = 0

      do i = 1, grid%io%group_size
        recvcounts(i) = grid%io%tiles( 1, 1, i ) * &
                        grid%io%tiles( 1, 2, i ) * &
                        grid%io%tiles( 1, 3, i )
        data_size = data_size + recvcounts(i)
      enddo

      displs(1) = 0
      do i = 2, grid%io%group_size
        displs(i) = displs(i-1) + recvcounts(i-1)
      enddo

      call alloc(gather_buffer, [ data_size ],"vdf/os-vdf-reportfile.f03",908)
    else
      recvcounts => null()
      displs => null()
      gather_buffer => null()
    endif

    ! Copy the data to communication buffer
    call alloc(comm_buffer, vdf % nx_(1:3),"vdf/os-vdf-reportfile.f03",916)
    do i3 = 1, vdf % nx_(3)
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          comm_buffer( i1, i2, i3 ) = real( vdf % f3 ( fc, i1, i2, i3 ), p_double )
        enddo
      enddo
    enddo

    ! Gather the data on the group root
    if ( grid%io%group_rank == 0 ) then
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_DOUBLE_PRECISION, &
                        gather_buffer, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                        0, grid%io%group_comm, ierr )
    else
      call mpi_gatherv( comm_buffer, size( comm_buffer ), MPI_DOUBLE_PRECISION, &
                        tmp, tmp, tmp, MPI_DOUBLE_PRECISION, &
                        0, grid%io%group_comm, ierr )
    endif

    ! free the communication buffer
    call freemem(comm_buffer,"vdf/os-vdf-reportfile.f03",937)

    ! Rearrange data into properly formed 2D arrays
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc(write_buffer, grid%io%group_tile( 1, : ),"vdf/os-vdf-reportfile.f03",943)

      ! Transpose data from merge_buffer into data_buffer
      k = 0
      do i = 1, grid%io%group_size
        do i3 = grid%io%tiles( 2, 3, i ) + 1, &
                grid%io%tiles( 2, 3, i ) + grid%io%tiles( 1, 3, i )
          do i2 = grid%io%tiles( 2, 2, i ) + 1, &
                  grid%io%tiles( 2, 2, i ) + grid%io%tiles( 1, 2, i )
            do i1 = grid%io%tiles( 2, 1, i ) + 1, &
                    grid%io%tiles( 2, 1, i ) + grid%io%tiles( 1, 1, i )
              k = k+1
              write_buffer( i1, i2, i3 ) = gather_buffer( k )
            enddo
          enddo
        enddo
      enddo

      ! Free allocated memory on group root
      call freemem(gather_buffer,"vdf/os-vdf-reportfile.f03",962)
      call freemem(displs,"vdf/os-vdf-reportfile.f03",963)
      call freemem(recvcounts,"vdf/os-vdf-reportfile.f03",964)

    endif

  end select

! print *, "[",mpi_node(),"] ",'Writing merged data...'

  ! Write merged data to disk

  ! Only group leaders participate in writing data to disk
  ! Note than when not using merging all nodes have group_rank = 0
  if ( grid%io%group_rank == 0 ) then

    ! create file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

! print *, "[",mpi_node(),"] ",'Opening file...'

    ! Open the file for (possible) parallel I/O
    call diagFile % open( p_diag_create, grid % io % comm )

! print *, "[",mpi_node(),"] ",'Starting CDSET...'

    call diagFile % start_cdset( report%name, 3, grid % g_nx, freal_to_diagtype( p_double ), &
                    parDset, grid % io % chunk_size )

    do i = 1, 3
      parChunk % count(i) = grid % io % group_tile( 1 , i )
      parChunk % start(i) = grid % io % group_tile( 2 , i )
      parChunk % stride(i) = 1
    enddo
    parChunk % data = c_loc( write_buffer )

! print *, "[",mpi_node(),"] ",'Writing CDSET...'

    ! Write the data
    call diagFile % write_par_cdset( parDset, parChunk )

    ! Close dataset
    call diagFile % end_cdset( parDset )

    ! Close file
    call diagFile % close( )

    deallocate( diagFile )

    ! free merge data buffer
    call freemem(write_buffer,"vdf/os-vdf-reportfile.f03",1012)

  endif

end subroutine write_vdf_3d_r8

!---------------------------------------------------------------------------------------------------
! End of template Functions, do not change below
!---------------------------------------------------------------------------------------------------
# 135 "vdf/os-vdf-reportfile.f03" 2

end module m_vdf_reportfile


!---------------------------------------------------------------------------------------------------
! Save vdf grid to file
!---------------------------------------------------------------------------------------------------
subroutine write_vdf( vdf, report, fc, g_space, grid )

  use m_vdf_define, only : t_vdf, t_vdf_report
  use m_space
  use m_grid_define
  use m_node_conf
  use m_system
  use m_parameters
  use m_vdf_reportfile

  implicit none

  class( t_vdf ), intent(in) :: vdf
  type( t_vdf_report ), intent(in) :: report
  integer, intent(in) :: fc
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid

  ! Precision to be used for diagnostics
  integer :: prec

  ! Don't do diagnostics with better precision than the one used for storing the values
  if ( report%prec > p_k_fld ) then
    prec = p_k_fld
  else
    prec = report%prec
  endif

  ! Call appropriate routine according to precision / dimensions
  select case( prec )
    case ( p_single )
      select case (vdf % x_dim_)
      case(1)
        call write_vdf_1d_r4( vdf, report, fc, g_space, grid )
      case(2)
        call write_vdf_2d_r4( vdf, report, fc, g_space, grid )
      case(3)
        call write_vdf_3d_r4( vdf, report, fc, g_space, grid )
      case default
        write(err_buf__,*) 'Invalid value for vdf % x_dim_';call err__("vdf/os-vdf-reportfile.f03",181)
        call abort_program( p_err_invalid )
      end select

    case ( p_double )
      select case (vdf % x_dim_)
      case(1)
        call write_vdf_1d_r8( vdf, report, fc, g_space, grid )
      case(2)
        call write_vdf_2d_r8( vdf, report, fc, g_space, grid )
      case(3)
        call write_vdf_3d_r8( vdf, report, fc, g_space, grid )
      case default
        write(err_buf__,*) 'Invalid value for vdf % x_dim_';call err__("vdf/os-vdf-reportfile.f03",194)
        call abort_program( p_err_invalid )
      end select

    case default
      write(err_buf__,*) 'Invalid precision selected for VDF diagnostics';call err__("vdf/os-vdf-reportfile.f03",199)
      call abort_program( p_err_invalid )
  end select

end subroutine write_vdf
