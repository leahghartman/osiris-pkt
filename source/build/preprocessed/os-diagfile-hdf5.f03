# 1 "dutil/os-diagfile-hdf5.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "dutil/os-diagfile-hdf5.f03" 2
# 13 "dutil/os-diagfile-hdf5.f03"
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
# 14 "dutil/os-diagfile-hdf5.f03" 2

module m_diagfile_hdf5

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
# 18 "dutil/os-diagfile-hdf5.f03" 2

use m_diagfile
use m_parameters
use hdf5
use mpi

implicit none

private

! TODO: use the ones that are in os-diagfile.f03 rather then cut-n-paste here
!integer, private, parameter :: p_single = kind(1.0e0)
!integer, private, parameter :: p_double = kind(1.0d0)
!integer, private, parameter :: p_byte = selected_int_kind(2)
!integer, private, parameter :: p_int64 = selected_int_kind(10)

type, extends( t_diag_file ) :: t_diag_file_hdf5

  ! File id
  integer(hid_t) :: id

  ! Parallel communicator
  integer :: comm

  ! Parallel I/O mode
  integer :: iomode

  ! flag to determine if local parallel node writes to the file
  logical :: write

contains

  procedure :: open_diag_file => open_file_hdf5
  procedure :: open_par_diag_file => open_par_file_hdf5

  procedure :: close => close_file_hdf5

  procedure :: add_dataset_1D_r4 => add_dataset_1D_r4_hdf5
  procedure :: add_dataset_2D_r4 => add_dataset_2D_r4_hdf5
  procedure :: add_dataset_3D_r4 => add_dataset_3D_r4_hdf5

  procedure :: add_dataset_1D_r8 => add_dataset_1D_r8_hdf5
  procedure :: add_dataset_2D_r8 => add_dataset_2D_r8_hdf5
  procedure :: add_dataset_3D_r8 => add_dataset_3D_r8_hdf5

  ! Gfortran throws the strangest error in 'os-paritcles-define.f90' if a enable this method.
  ! for disable for now since it is unused in the rest of the code
  !procedure :: add_dataset_loc_ => add_dataset_loc_hdf5_

  procedure :: start_cdset => start_cdset_hdf5
  procedure :: write_cdset => write_simple_cdset_hdf5
  procedure :: write_par_cdset => write_par_cdset_hdf5
  procedure :: end_cdset => end_cdset_hdf5

  procedure :: open_cdset => open_cdset_hdf5
  procedure :: extend_cdset => extend_cdset_hdf5
  procedure :: close_cdset => close_cdset_hdf5

  procedure :: file_extension => file_extension_hdf5

end type t_diag_file_hdf5

public :: t_diag_file_hdf5

contains

pure function file_extension_hdf5( this )

  implicit none
  class( t_diag_file_hdf5 ), intent(in) :: this
  character(len=3) :: file_extension_hdf5
  file_extension_hdf5 = 'h5'

end function

! --------------------------------------------------------------------------------------------------
! Close diagnostics file
! --------------------------------------------------------------------------------------------------
subroutine close_file_hdf5( this )

  implicit none
  class( t_diag_file_hdf5 ), intent(inout) :: this
  integer :: ierr

  if ( this % write ) then
    call h5fclose_f( this%id, ierr)
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",104 -1
  endif

end subroutine close_file_hdf5
! --------------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------------
! Adds general file attributes
! --------------------------------------------------------------------------------------------------
subroutine add_file_attr_hdf5( this )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this

  integer(hid_t) :: rootID
  integer :: ierr

  if ( this % write ) then

    ! Open root group to add metadata
    call h5gopen_f( this%id, '/', rootID, ierr )
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",128 -1

    ! Add name property
    call hdf5_add( rootID, 'NAME', this%name )

    ! Add general metadata
    select case ( this%ftype )
    case (p_diag_grid)
      call hdf5_add( rootID, "TYPE", "grid" )
      call hdf5_add( rootID, this % grid )
      call hdf5_add( rootID, this % iter )

    case (p_diag_particles)
      call hdf5_add( rootID, "TYPE", "particles");
      call hdf5_add( rootID, this % particles );
      call hdf5_add( rootID, this % iter );

    case (p_diag_tracks)
      call hdf5_add( rootID, "TYPE", "tracks-2");
      call hdf5_add( rootID, this % tracks );

    case default
      write(0,*) '(*error*) Unsupported file type, aborting.'
      call abort_program( p_err_invalid )
    end select

    ! Add simulation info
    call hdf5_add( rootID, "SIMULATION", sim_info )

    ! Close root group
    call h5gclose_f( rootID, ierr )
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",159 -1

  endif
end subroutine add_file_attr_hdf5

! --------------------------------------------------------------------------------------------------
! Open Diagnostic file
! --------------------------------------------------------------------------------------------------
subroutine open_file_hdf5( this, amode )

  use m_grid_define
  use m_system
  use m_parameters
  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  integer, intent(in) :: amode

  character(len=1024) :: path
  integer :: ierr
  integer(hid_t) :: plistID
  integer(hsize_t) :: threshold, alignment

  ! integer(size_t) :: rdcc_nslots, rdcc_nbytes
  ! real :: rdcc_w0

  ! Create directory
  call mkdir( this%filepath, ierr )

  this % write = .true.

  ! open the file
  path = trim(this%filepath)//trim(this%filename)//'.h5'

  ! set to true for all amode
  this % write = .true.

  if ( amode == p_diag_create ) then

    call h5pcreate_f(H5P_FILE_ACCESS_F, plistID, ierr)
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",201 -1

    ! Sets alignment properties
    threshold = 0
    alignment = 262144
    !call h5pset_alignment_f( plistID, threshold, alignment, ierr )
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",207 -1

    call h5fcreate_f( path, H5F_ACC_TRUNC_F, this % id, ierr, access_prp = plistID )
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",210 -1

    call add_file_attr_hdf5( this )

    ! Close the property list
    call h5pclose_f(plistID, ierr)
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",216 -1

  else
    call h5fopen_f( path, H5F_ACC_RDWR_F, this % id, ierr)
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",220 -1

  endif

  this % comm=MPI_COMM_NULL

end subroutine open_file_hdf5

subroutine open_par_file_hdf5( this, amode, comm, iomode )

  use m_grid_define
  use m_system
  use m_parameters
  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  integer, intent(in) :: amode
  integer, intent(in) :: comm
  integer, intent(in), optional :: iomode

  integer :: ierr, rank
# 251 "dutil/os-diagfile-hdf5.f03"
  ! Variables for setting the data chunk cache
  ! integer(hsize_t) :: rdcc_nslots, rdcc_nbytes
  ! real :: rdcc_w0

  if ( present(iomode)) then
    this % iomode = iomode
  else
    this % iomode = diag_file_default_iomode
  endif

  ! If invalid comm just use the serial file open
  ! Ricardo: Look here... does this fit with how you were thinking of things
  ! (I speak of the the ' comm == MPI_COMM_SELF' portion)
  if ( comm == MPI_COMM_NULL .or. comm == MPI_COMM_SELF) then
    call open_file_hdf5( this, amode )
    return
  endif

  ! Only root node creates the directory
  call mpi_comm_rank( comm, rank, ierr )
  if ( rank == 0 ) then
    call mkdir( this%filepath, ierr )
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",273 -1
  endif

  select case ( this % iomode )
  case ( DIAG_MPI )
    ! For MPI access mode only root node writes to the file
    if ( rank == 0 ) then
      this % write = .true.
      call open_file_hdf5( this, amode )
    else
      this % write = .false.
    endif

  case ( DIAG_MPIIO_INDEPENDENT, DIAG_MPIIO_COLLECTIVE )
# 329 "dutil/os-diagfile-hdf5.f03"
  if ( rank == 0 ) then
    write(0,'(A)') '(*error*) Unable to create HDF5 file, parallel i/o mode not supported.'
    write(0,'(A)') '(*error*) HDF5 library was not compiled with parallel i/o support please'
    write(0,'(A)') '(*error*) use "mpi" mode instead.'
  endif


  case default

    if ( rank == 0 ) then
      write(0,'(A)') '(*error*) Unable to create HDF5 file, "posix" parallel i/o mode not supported.'
    endif

  end select

  ! Store the communicator
  this % comm = comm

end subroutine open_par_file_hdf5

subroutine add_dataset_loc_hdf5_( this, name, buffer, ndims, count, data_type )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len=*), intent(in) :: name
  type(c_ptr), intent(in) :: buffer
  integer, intent(in) :: data_type, ndims
  integer, dimension(:), intent(in) :: count

  integer(p_byte), dimension(:), pointer :: tmp
  integer(hsize_t), dimension(ndims) :: dims
  integer(hid_t) :: dataspaceID, datasetID, datatypeID
  integer :: bsize, i, ierr

  bsize = diag_sizeof( data_type )
  do i = 1, ndims
    dims(i) = count(i)
    bsize = bsize * count(i)
  enddo
  call c_f_pointer( buffer, tmp, [bsize] )

  ! Convert data type to 1 equivalent
  datatypeID = diag_hdf5_type( data_type )

  ! create the dataset
  call h5screate_simple_f( ndims, dims, dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",378 -1

  call h5dcreate_f( this % id, name, datatypeID, dataspaceID, datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",381 -1

  ! write the data
  call h5dwrite_f( datasetID, datatypeID, tmp, dims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",385 -1

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",389 -1

  call h5dclose_f( datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",392 -1


end subroutine add_dataset_loc_hdf5_

! --------------------------------------------------------------------------------------------------
! Start chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine start_cdset_hdf5( this, name, ndims, count, data_type, dset, chunk_size )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len = *), intent(in) :: name
  integer, intent(in) :: ndims
  integer, dimension(:), intent(in) :: count
  integer, intent(in) :: data_type
  type( t_diag_dataset ), intent(inout) :: dset

  integer, dimension(:), intent(in), optional :: chunk_size

  integer(hsize_t), dimension(diag_max_dims) :: dims, maxDims, chunkDims
  integer(hid_t) :: filespaceID, dsetID, propertyID
  integer :: i, ierr

  dset % ndims = ndims
  dset % count(1:ndims) = count(1:ndims)
  dset % data_type = data_type

  if ( this % write ) then

    if ( present(chunk_size) ) then
      ! Create chunked dataset
      do i = 1, dset%ndims
        dims(i) = dset % count(i)
        maxDims(i) = dset % count(i)
        chunkDims(i) = chunk_size(i)
        if ( chunkDims(i) > maxDims(i) .and. maxDims(i) > 0 ) chunkDims(i) = maxDims(i)
      enddo
      maxDims( dset%ndims ) = H5S_UNLIMITED_F

      call h5screate_simple_f( dset % ndims, dims, filespaceID, ierr, &
                               maxdims = maxDims )
      if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",437 -1

      ! Using unlimited dataspaces requires chunking
      call h5pcreate_f( H5P_DATASET_CREATE_F, propertyID, ierr )
      if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",441 -1

      call h5pset_chunk_f( propertyID, dset%ndims, chunkDims, ierr )
      if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",444 -1

      call h5dcreate_f( this%id, trim(name), diag_hdf5_type( dset%data_type), filespaceID, dsetID, ierr, &
                        dcpl_id = propertyID )
      if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",448 -1

      call h5pclose_f( propertyID, ierr )
      if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",451 -1
    else
      ! Create contiguous dataset
      do i = 1, dset%ndims
        dims(i) = dset % count(i)
      enddo

      if ( dims(1) > 0 ) then
        call h5screate_simple_f( dset % ndims, dims, filespaceID, ierr )
        if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",460 -1
      else
        call h5screate_f( H5S_SCALAR_F, filespaceID, ierr )
        if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",463 -1
      endif

      call h5dcreate_f( this%id, trim(name), diag_hdf5_type( dset%data_type), filespaceID, dsetID, ierr )
      if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",467 -1
    endif

    call h5sclose_f(filespaceID, ierr)
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",471 -1

    ! Store dataset ID
    dset % id = dsetID


  endif

end subroutine start_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Write part of a chunked dataset (not parallel)
! --------------------------------------------------------------------------------------------------
subroutine write_simple_cdset_hdf5( this, dset, chunk )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  type( t_diag_chunk ), intent(inout) :: chunk

  integer(hid_t) :: filespaceID, memspaceID, datasetID
  integer(hsize_t), dimension(diag_max_dims) :: start, count
  integer(p_byte), dimension(:), pointer :: buffer
  integer :: bsize, i, ierr

  datasetID = dset % id

  bsize = diag_sizeof( dset % data_type )
  do i = 1, dset % ndims
    start(i) = chunk % start(i)
    count(i) = chunk % count(i)
    bsize = bsize * chunk % count(i)
  enddo
  call c_f_pointer( chunk % data, buffer, [bsize] )

  ! create memory dataspace
  call h5screate_simple_f( dset%ndims, count, memspaceID, ierr)
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",511 -1

  ! select hyperslab in the file
  call h5dget_space_f( datasetID, filespaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",515 -1
  call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",517 -1

  ! write data
  call h5dwrite_f(datasetID, diag_hdf5_type( dset % data_type ), buffer, count, ierr, &
      file_space_id = filespaceID, mem_space_id = memspaceID)

  ! close resources
  call h5sclose_f(filespaceID, ierr)
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",525 -1
  call h5sclose_f(memspaceID, ierr)
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",527 -1

end subroutine write_simple_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Write parallel chunked dataset using mpi/io
! --------------------------------------------------------------------------------------------------
# 611 "dutil/os-diagfile-hdf5.f03"
subroutine write_cdset_mpi( this, dset, chunk )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  type( t_diag_chunk ), intent(inout) :: chunk

  integer, parameter :: ping_tag = 1001, tile_tag = 1002, data_tag = 1003
  integer :: ping, data_size, comm_rank, comm_size, ping_req, data_req
  integer :: step, src, i, ierr, bsize
  integer, dimension( 2 * diag_max_dims ) :: cinfo0, cinfo1
  integer(p_byte), dimension(:), pointer :: buffer0, buffer1, bufferw
  integer :: bsize0, bsize1, bsizew
  integer(hsize_t), dimension(diag_max_dims) :: start, count
  integer(hid_t) :: datasetID, memspaceID, filespaceID

  datasetID = dset % id

  ping = 1234

  ! Get size of data element
  data_size = diag_sizeof( dset % data_type )

  ! Point buffer0 to local data
  bsize0 = 1
  do i = 1, dset % ndims
    cinfo0( i ) = chunk % count(i)
    cinfo0( dset % ndims + i ) = chunk % start(i)
    bsize0 = bsize0 * chunk % count(i)
  enddo
  call c_f_pointer( chunk % data, buffer0, [bsize0*data_size] )

  call MPI_COMM_RANK( this%comm, comm_rank, ierr )

  if ( comm_rank == 0 ) then

    ! Get filespace
    call h5dget_space_f(datasetID, filespaceID, ierr)
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",652 -1

    ! Loop over all nodes
    call MPI_COMM_SIZE( this%comm, comm_size, ierr )

    bsize0 = 0
    bsize1 = 0
    step = 0
    do src = 0, comm_size-1
      ! Start receiving data from next node
      if ( src < comm_size - 1 ) then
        ! send ping
        call MPI_ISEND( ping, 1, MPI_INTEGER, src+1, ping_tag, this % comm, ping_req, ierr )

        if ( step == 0 ) then
          ! Receive on buffer1
          call MPI_RECV( cinfo1, 2 * dset%ndims, MPI_INTEGER, src+1, tile_tag, this%comm, MPI_STATUS_IGNORE, ierr)
          bsize = cinfo1(1)
          do i = 2, dset % ndims
            bsize = bsize * cinfo1(i)
          enddo
          if ( bsize > bsize1 ) then
            if ( bsize1 > 0 ) call freemem(buffer1,"dutil/os-diagfile-hdf5.f03",674)
            call alloc(buffer1, [bsize * data_size],"dutil/os-diagfile-hdf5.f03",675)
            bsize1 = bsize
          endif

          if( bsize .ne. 0 ) call MPI_IRECV( buffer1, bsize, diag_mpi_type(dset%data_type), src+1, data_tag, this%comm, data_req, ierr)

        else
          ! Receive on buffer0
          call MPI_RECV( cinfo0, 2 * dset%ndims, MPI_INTEGER, src+1, tile_tag, this%comm, MPI_STATUS_IGNORE, ierr)
          bsize = cinfo0(1)
          do i = 2, dset % ndims
            bsize = bsize * cinfo0(i)
          enddo
          if ( bsize > bsize0 ) then
            if ( bsize0 > 0 ) call freemem(buffer0,"dutil/os-diagfile-hdf5.f03",689)
            call alloc(buffer0, [bsize * data_size],"dutil/os-diagfile-hdf5.f03",690)
            bsize0 = bsize
          endif

          if( bsize .ne. 0 ) call MPI_IRECV( buffer0, bsize, diag_mpi_type(dset%data_type), src+1, data_tag, this%comm, data_req, ierr)
        endif
      endif

      ! write data from current node
      if ( step == 0 ) then
        bufferw => buffer0
        do i = 1, dset % ndims
          count(i) = cinfo0( i )
          start(i) = cinfo0( dset % ndims + i )
        enddo
      else
        bufferw => buffer1
        do i = 1, dset % ndims
          count(i) = cinfo1( i )
          start(i) = cinfo1( dset % ndims + i )
        enddo
      endif

      bsizew = 1
      do i = 1, dset % ndims
        bsizew = bsizew * count(i)
      enddo

      if ( bsizew .ne. 0 ) then

        call h5screate_simple_f( dset % ndims, count, memspaceID, ierr)
        if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",721 -1

        call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )
        if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",724 -1

        call h5dwrite_f( datasetID, diag_hdf5_type( dset % data_type ), bufferw, count, ierr, &
          mem_space_id = memspaceID, file_space_id = filespaceID)
        if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",728 -1

        ! close resources
        call h5sclose_f(memspaceID, ierr)
        if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",732 -1

      endif

      ! Wait for messages to complete
      if ( src < comm_size - 1 ) then
        call MPI_WAIT( ping_req, MPI_STATUS_IGNORE, ierr )
        if ( bsize .ne. 0 ) call MPI_WAIT( data_req, MPI_STATUS_IGNORE, ierr )
      endif

      step = 1 - step
    enddo

    ! Free communication buffers
    if ( bsize0 > 0 ) call freemem(buffer0,"dutil/os-diagfile-hdf5.f03",746)
    if ( bsize1 > 0 ) call freemem(buffer1,"dutil/os-diagfile-hdf5.f03",747)

    ! close resources
    call h5sclose_f( filespaceID, ierr )

  else
    ! Wait for ping from root node
    call MPI_RECV( ping, 1, MPI_INTEGER, 0, ping_tag, this%comm, MPI_STATUS_IGNORE, ierr)

    ! send local chunk information to root node
    call MPI_SEND( cinfo0, 2 * dset%ndims, MPI_INTEGER, 0, tile_tag, this%comm, ierr)

    if( bsize0 > 0 ) then
      ! Send chunk data to root node
      call MPI_SEND( buffer0, bsize0, diag_mpi_type( dset % data_type ), 0, data_tag, this%comm, ierr )
    endif

  endif

end subroutine write_cdset_mpi

subroutine write_par_cdset_hdf5( this, dset, chunk, offset_ )

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  type( t_diag_chunk ), intent(inout) :: chunk

  ! Offset is not used by 1 routines
  integer, intent(in), optional :: offset_

  integer :: rank, ierr

  ! RICARDO LOOK HERE
  ! If comm is MPI_COMM_NULL switch to the serial write code,
  ! but the code should not get here so this is a bit of a hack.
  ! But not so much of a hack since the parallel open function
  ! switches back to serial... so it's kinda ok but maybe not..
  !
  ! Also, Ricardo: does the MPI_COMM_SELF fit with your thinking here?
  if ( this%comm == MPI_COMM_NULL .or. this%comm == MPI_COMM_SELF ) then
    call this%write_cdset( dset, chunk )
    return
  endif

  select case ( this % iomode )
  case(DIAG_MPI)
    call write_cdset_mpi( this, dset, chunk )






  case default
    call mpi_comm_rank( this % comm, rank, ierr )
    if ( rank == 0 ) then
      write(0,'(A)') '(*error*) Invalid parallel IO mode, unable to write dataset.'
    endif
  end select

end subroutine write_par_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Write end marker for a chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine end_cdset_hdf5( this, dset )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset

  integer :: ierr
  integer(hid_t) :: dsetID

  if ( this % write ) then
    dsetID = dset % id
    call h5dclose_f( dsetID, ierr )
    if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",829 -1
  endif

end subroutine end_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Open a chunked dataset for updating
! --------------------------------------------------------------------------------------------------
subroutine open_cdset_hdf5( this, dset )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset

  integer :: i, ierr
  integer(hsize_t), dimension(diag_max_dims) :: dims, maxdims
  integer(hid_t) :: datasetID, datatypeID, filespaceID

  if ( .not. this % write ) return

  call h5dopen_f( this%id, trim(dset%name), datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",853 -1

  ! Get datatype
  call h5dget_type_f(datasetID, datatypeID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",857 -1
  dset % data_type = hdf5_diag_type( datatypeID )
  call h5tclose_f( datatypeID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",860 -1

  ! Get dimensions
  call h5dget_space_f( datasetID, filespaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",864 -1

  call h5sget_simple_extent_ndims_f(filespaceID, dset % ndims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",867 -1

  call h5sget_simple_extent_dims_f( filespaceID, dims, maxdims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",870 -1

  do i = 1, dset % ndims
    dset % count(i) = dims(i)
  enddo

  ! Store dataset ID
  dset % id = datasetID

  ! The offset is not used for 1 files
  dset % offset = -1

end subroutine open_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Extend dimensions of chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine extend_cdset_hdf5( this, dset, dims )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  integer(p_int64), dimension(:), intent(in) :: dims

  integer :: i, ierr
  integer(hid_t) :: datasetID

  datasetID = dset % id
  call h5dset_extent_f( datasetID, dims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",902 -1

  do i = 1, dset % ndims
    dset % count(i) = dims(i)
  enddo

end subroutine extend_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Close previously open chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine close_cdset_hdf5( this, dset )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset

  integer :: ierr
  integer(hid_t) :: datasetID

  datasetID = dset % id
  call h5dclose_f( datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",927 -1

  dset % ndims = 0
  dset % id = -1

end subroutine close_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! add_dataset routines
! --------------------------------------------------------------------------------------------------
# 945 "dutil/os-diagfile-hdf5.f03"
# 1 "dutil/os-diagfile-hdf5.f03" 1
# 987 "dutil/os-diagfile-hdf5.f03"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------

subroutine add_dataset_1d_r4_hdf5( this, name, buffer )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(p_single), dimension(:), intent(inout), target :: buffer

  integer(hsize_t), dimension(1) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: i, ierr

  do i = 1, 1
    dims(i) = size(buffer, i)
  enddo

  ! create the dataset
  call h5screate_simple_f( 1, dims, dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1011 -1

  call h5dcreate_f( this % id, name, H5T_NATIVE_REAL, dataspaceID, datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1014 -1

  ! write the data
  call h5dwrite_f( datasetID, H5T_NATIVE_REAL, buffer, dims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1018 -1

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1022 -1

  call h5dclose_f( datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1025 -1


end subroutine add_dataset_1d_r4_hdf5
# 946 "dutil/os-diagfile-hdf5.f03" 2






# 1 "dutil/os-diagfile-hdf5.f03" 1
# 987 "dutil/os-diagfile-hdf5.f03"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------

subroutine add_dataset_2d_r4_hdf5( this, name, buffer )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(p_single), dimension(:,:), intent(inout), target :: buffer

  integer(hsize_t), dimension(2) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: i, ierr

  do i = 1, 2
    dims(i) = size(buffer, i)
  enddo

  ! create the dataset
  call h5screate_simple_f( 2, dims, dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1011 -1

  call h5dcreate_f( this % id, name, H5T_NATIVE_REAL, dataspaceID, datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1014 -1

  ! write the data
  call h5dwrite_f( datasetID, H5T_NATIVE_REAL, buffer, dims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1018 -1

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1022 -1

  call h5dclose_f( datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1025 -1


end subroutine add_dataset_2d_r4_hdf5
# 953 "dutil/os-diagfile-hdf5.f03" 2






# 1 "dutil/os-diagfile-hdf5.f03" 1
# 987 "dutil/os-diagfile-hdf5.f03"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------

subroutine add_dataset_3d_r4_hdf5( this, name, buffer )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(p_single), dimension(:,:,:), intent(inout), target :: buffer

  integer(hsize_t), dimension(3) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: i, ierr

  do i = 1, 3
    dims(i) = size(buffer, i)
  enddo

  ! create the dataset
  call h5screate_simple_f( 3, dims, dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1011 -1

  call h5dcreate_f( this % id, name, H5T_NATIVE_REAL, dataspaceID, datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1014 -1

  ! write the data
  call h5dwrite_f( datasetID, H5T_NATIVE_REAL, buffer, dims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1018 -1

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1022 -1

  call h5dclose_f( datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1025 -1


end subroutine add_dataset_3d_r4_hdf5
# 960 "dutil/os-diagfile-hdf5.f03" 2






# 1 "dutil/os-diagfile-hdf5.f03" 1
# 987 "dutil/os-diagfile-hdf5.f03"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------

subroutine add_dataset_1d_r8_hdf5( this, name, buffer )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(p_double), dimension(:), intent(inout), target :: buffer

  integer(hsize_t), dimension(1) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: i, ierr

  do i = 1, 1
    dims(i) = size(buffer, i)
  enddo

  ! create the dataset
  call h5screate_simple_f( 1, dims, dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1011 -1

  call h5dcreate_f( this % id, name, H5T_NATIVE_DOUBLE, dataspaceID, datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1014 -1

  ! write the data
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, buffer, dims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1018 -1

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1022 -1

  call h5dclose_f( datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1025 -1


end subroutine add_dataset_1d_r8_hdf5
# 967 "dutil/os-diagfile-hdf5.f03" 2






# 1 "dutil/os-diagfile-hdf5.f03" 1
# 987 "dutil/os-diagfile-hdf5.f03"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------

subroutine add_dataset_2d_r8_hdf5( this, name, buffer )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(p_double), dimension(:,:), intent(inout), target :: buffer

  integer(hsize_t), dimension(2) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: i, ierr

  do i = 1, 2
    dims(i) = size(buffer, i)
  enddo

  ! create the dataset
  call h5screate_simple_f( 2, dims, dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1011 -1

  call h5dcreate_f( this % id, name, H5T_NATIVE_DOUBLE, dataspaceID, datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1014 -1

  ! write the data
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, buffer, dims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1018 -1

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1022 -1

  call h5dclose_f( datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1025 -1


end subroutine add_dataset_2d_r8_hdf5
# 974 "dutil/os-diagfile-hdf5.f03" 2






# 1 "dutil/os-diagfile-hdf5.f03" 1
# 987 "dutil/os-diagfile-hdf5.f03"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------

subroutine add_dataset_3d_r8_hdf5( this, name, buffer )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(p_double), dimension(:,:,:), intent(inout), target :: buffer

  integer(hsize_t), dimension(3) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: i, ierr

  do i = 1, 3
    dims(i) = size(buffer, i)
  enddo

  ! create the dataset
  call h5screate_simple_f( 3, dims, dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1011 -1

  call h5dcreate_f( this % id, name, H5T_NATIVE_DOUBLE, dataspaceID, datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1014 -1

  ! write the data
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, buffer, dims, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1018 -1

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1022 -1

  call h5dclose_f( datasetID, ierr )
  if (ierr<0) write(0,"('(*io error*)',A,':',I0)") "dutil/os-diagfile-hdf5.f03",1025 -1


end subroutine add_dataset_3d_r8_hdf5
# 981 "dutil/os-diagfile-hdf5.f03" 2


end module m_diagfile_hdf5
