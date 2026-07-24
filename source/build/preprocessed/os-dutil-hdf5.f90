# 1 "dutil/os-dutil-hdf5.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "dutil/os-dutil-hdf5.f90" 2
!--------------------------------------------------------------------------------------------------
!
! Using parallel I/O:
!
! Chunking is a file property, so it is not possible to define different chunk sizes on different
! processes. This means that when the array size is different on some/any process the only options
! are to use normal (contiguous) storage, or (maybe) define a chunk size that is the maximum

! array size on every node.
!--------------------------------------------------------------------------------------------------

! if __TEMPLATE__ is defined then read the template definition at the end of the file


!#define DEBUG_IO 1
# 25 "dutil/os-dutil-hdf5.f90"
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
# 26 "dutil/os-dutil-hdf5.f90" 2

module hdf5_util

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
# 30 "dutil/os-dutil-hdf5.f90" 2

use hdf5
use mpi

implicit none

private

integer, private, parameter :: p_single = kind(1.0e0)
integer, private, parameter :: p_double = kind(1.0d0)
integer, private, parameter :: p_byte = selected_int_kind(2)
integer, private, parameter :: p_int64 = selected_int_kind(10)

integer, parameter :: p_lower = 1
integer, parameter :: p_upper = 2

integer, parameter :: p_hdf5_create = 0
integer, parameter :: p_hdf5_update = 1


type :: t_h5_tune

  ! Use gpfs hints / optimizations
  logical :: gpfs = .false.

  ! This was deprecated as of 1 v1.8.13
  ! use mpi+POSIX instead of MPI-IO for parallel I/O
  ! logical :: posix = .false.

  ! metadata block size
  integer(hsize_t) :: meta_block_size = -1

  ! cache size (I think this is only for read operations...)
  integer(size_t) :: cache = -1

  ! data alignment in file
  integer(hsize_t), dimension(2) :: alignment = -1

  ! sieve buffer size
  integer(size_t) :: sieve_buf_size = -1

  ! MPI-IO hints

  ! collective_buffering - Specifies whether the application may benefit from collective buffering.
  logical :: collective_buffering = .false.

  ! cb_block_size - Specifies the block size to be used for collective buffering. Data access happens
  ! in chunks of this size.
  integer :: cb_block_size = -1

  ! cb_buffer_size - Specifies the size of the buffer space that can be used on each target node.
  ! e.g. 16777216, 33554432
  integer :: cb_buffer_size = -1

  ! access_style - Specifies the manner in which the file will be accessed until
  ! it is closed or this info key is changed. Valid list elements are:
  ! read_once, write_once, read_mostly, write_mostly, sequential, reverse_sequential, random
  character( len = 1024 ) :: access_style = '-'

  ! romio_ds_write
  ! Specifies whether to use data sieving for write access.
  character( len = 1024 ) :: romio_ds_write = '-'

  ! romio_cb_write
  ! Specifies whether to use collective buffering for write access.
  ! enable, disable, automatic (default)
  character( len = 1024 ) :: romio_cb_write = '-'

  ! Striping parameters
  ! striping_factor - Specifies the number of I/O devices that the file should be striped across.
  ! Relevant only on file creation.
  integer :: striping_factor = -1

  ! striping_unit - Specifies the striping unit – the amount of consecutive data assigned to one
  ! I/O device – to be used for this file. Only relevant on file creation.
  integer :: striping_unit = -1


  ! dataset tune

  ! Internal buffer size
  integer(hsize_t) :: buffer_size = -1

  ! Use chunked dataset
  logical :: chunked = .true.
  logical :: no_chunk_last_dim = .false.

  ! use independent or collective transfers
  integer :: transferMode

end type t_h5_tune

type(t_h5_tune), save :: h5_tune

interface hdf5_add
  module procedure hdf5_add_str
  module procedure hdf5_add_str_v1
  module procedure hdf5_add_logical
  module procedure hdf5_add_logical_v1

  module procedure hdf5_add_single
  module procedure hdf5_add_v1_single

  module procedure hdf5_add_double
  module procedure hdf5_add_v1_double

  module procedure hdf5_add_int
  module procedure hdf5_add_v1_int

  ! Compound attributes
  module procedure add_grid_info
  module procedure add_part_info
  module procedure add_track_info
  module procedure add_iteration_info
  module procedure add_sim_info
end interface

interface open_hdf5
  module procedure open_hdf5
end interface

interface close_hdf5
  module procedure close_hdf5
end interface

interface diag_hdf5_type
  module procedure diag_hdf5_type
end interface

interface hdf5_diag_type
  module procedure hdf5_diag_type
end interface

public :: open_hdf5, close_hdf5, hdf5_add, diag_hdf5_type, hdf5_diag_type

contains

!---------------------------------------------------------------------------------------------------
! Get 1 type from diag type
!---------------------------------------------------------------------------------------------------
function diag_hdf5_type( data_type )

  use m_diagfile

  implicit none

  integer, intent(in) :: data_type
  integer(hid_t) :: diag_hdf5_type

  select case (data_type)
  case ( p_diag_float32 )
    diag_hdf5_type = H5T_NATIVE_REAL
  case ( p_diag_float64 )
    diag_hdf5_type = H5T_NATIVE_DOUBLE
  case ( p_diag_int32 )
    diag_hdf5_type = H5T_NATIVE_INTEGER
  case ( p_diag_int64 )
    diag_hdf5_type = h5kind_to_type( selected_int_kind(10), H5_INTEGER_KIND )
  case default
    diag_hdf5_type = -1
  end select

end function

!---------------------------------------------------------------------------------------------------
! Get diag type from 1 type
!---------------------------------------------------------------------------------------------------
function hdf5_diag_type( h5type )

  use m_diagfile

  implicit none

  integer :: hdf5_diag_type
  integer(hid_t), intent(in) :: h5type

  logical :: flag
  integer :: ierr

  call h5tequal_f( h5type, H5T_NATIVE_REAL, flag, ierr )

  if ( flag ) then
    hdf5_diag_type = p_diag_float32
    return
  endif

  call h5tequal_f( h5type, H5T_NATIVE_DOUBLE, flag, ierr )

  if ( flag ) then
    hdf5_diag_type = p_diag_float64
    return
  endif

  call h5tequal_f( h5type, H5T_NATIVE_INTEGER, flag, ierr )

  if ( flag ) then
    hdf5_diag_type = p_diag_int32
    return
  endif

  call h5tequal_f( h5type, h5kind_to_type( selected_int_kind(10), H5_INTEGER_KIND ), &
                   flag, ierr )

  if ( flag ) then
    hdf5_diag_type = p_diag_int64
    return
  endif

  write(0,'(A,I0)') '(*error*) Unable to determine diag type from H5 type: ', h5type
  hdf5_diag_type = p_diag_null

end


!---------------------------------------------------------------------------------------------------
! Start 1 interface
!---------------------------------------------------------------------------------------------------
subroutine open_hdf5( )

   implicit none

   integer :: ierr

   call h5open_f(ierr)



   ! this needs to be done after h5open_f othwerwise H5FD_MPIO_COLLECTIVE_F might not be
   ! initialized (strange implementation from the hdf guys)
# 267 "dutil/os-dutil-hdf5.f90"
end subroutine open_hdf5

!---------------------------------------------------------------------------------------------------
! Close 1 interface
!---------------------------------------------------------------------------------------------------
subroutine close_hdf5( ierr )

   implicit none

   integer, intent(out) :: ierr

   call h5close_f(ierr)


end subroutine close_hdf5


!---------------------------------------------------------------------------------------------------
! compound metadata tags
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine add_grid_info( objID, info )

  use m_diagfile

  implicit none
  integer(hid_t), intent(in) :: objID
  type( t_diag_grid_info ), intent(in) :: info
  real(p_double), dimension(2) :: axis_range
  integer(hid_t) :: axisGroupID, dataspaceID, datasetID
  integer(hsize_t), dimension(1) :: dims
  integer :: i, ierr
  integer, parameter :: izero = ichar('0')

  ! Add grid metadata
  call hdf5_add( objID, 'LABEL', info % label )
  call hdf5_add( objID, 'UNITS', info % units )

  ! Only add offset information if specified
  if ( info%offset_x(1) /= p_no_offset ) then
    call hdf5_add( objID, 'OFFSET_X', info % offset_x(1:info%ndims) )
  endif
  if ( info%offset_t /= p_no_offset ) then
    call hdf5_add( objID, 'OFFSET_T', info % offset_t )
  endif
  if ( info%offset_t_ax(1) /= p_no_offset ) then
    call hdf5_add( objID, 'OFFSET_T_AXES', info % offset_t_ax(1:info%ndims) )
  endif

  ! add axis information
  call h5gcreate_f( objID, 'AXIS', axisGroupID, ierr)


  dims(1) = 2
  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  do i = 1, info%ndims
    call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+i), H5T_NATIVE_DOUBLE, dataspaceID, &
      datasetID, ierr )


    call hdf5_add( datasetID, 'TYPE', 'linear' )
    call hdf5_add( datasetID, 'UNITS', info%axis(i)%units )
    call hdf5_add( datasetID, 'NAME', info%axis(i)%name )
    call hdf5_add( datasetID, 'LONG_NAME', info%axis(i)%label )

    axis_range(1) = info%axis(i)%min
    axis_range(2) = info%axis(i)%max

    call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, ierr )


    call h5dclose_f( datasetID, ierr )

  enddo

  call h5gclose_f( axisGroupID, ierr )


end subroutine add_grid_info

subroutine add_part_info( objID, info )

  use m_diagfile

  implicit none

  integer(hid_t), intent(in) :: objID
  type( t_diag_part_info ), intent(in) :: info

  integer :: nquants

  nquants = info % nquants

  call hdf5_add( objID, 'QUANTS', info % quants( 1:nquants))
  call hdf5_add( objID, 'LABELS', info % qlabels(1:nquants))
  call hdf5_add( objID, 'UNITS', info % qunits( 1:nquants))
  call hdf5_add( objID, 'OFFSET_T',info % offset_t(1:nquants))

  ! (*warning*) These are not currently supported
  ! add particle selection data
! call hdf5_add( rootID, 'SELECT GAMMA LIMIT', info%raw_gamma_limit )
! call hdf5_add( rootID, 'SELECT FRACTION', info%raw_fraction )
! call hdf5_add( rootID, 'SELECT MATH EXPR', info%raw_math_expr )

end subroutine add_part_info

subroutine add_track_info( objID, info )

  use m_diagfile

  implicit none

  integer(hid_t), intent(in) :: objID
  type( t_diag_track_info ), intent(in) :: info

  integer :: nquants

  nquants = info % nquants

  call hdf5_add( objID, 'NTRACKS', info % ntracks )
  call hdf5_add( objID, 'NDUMP', info % ndump )
  call hdf5_add( objID, 'NITER', info % niter )
  call hdf5_add( objID, 'QUANTS', info % quants(1:nquants))
  call hdf5_add( objID, 'LABELS', info % qlabels(1:nquants))
  call hdf5_add( objID, 'UNITS', info % qunits( 1:nquants))
  call hdf5_add( objID, 'OFFSET_T',info % offset_t(1:nquants))

end subroutine add_track_info

subroutine add_iteration_info( objID, info )

  use m_diagfile

  implicit none

  integer(hid_t), intent(in) :: objID
  type( t_diag_iteration ), intent(in) :: info

  call hdf5_add( objID, 'TIME', info%t )
  call hdf5_add( objID, 'ITER', info%n )
  call hdf5_add( objID, 'TIME UNITS', info%time_units )

end subroutine add_iteration_info

subroutine add_sim_info( objID, name, info )
  use m_system

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  type( t_sim_info ), intent(in) :: info

  integer :: ndims, ierr
  integer(hid_t) :: rootID

  call h5gcreate_f( objID, name, rootID, ierr)


  call hdf5_add( rootID, "GIT_VERSION", info % git_version )
  call hdf5_add( rootID, "COMPILE_TIME", info % compile_time )

  call hdf5_add( rootID, "INPUT_FILE", info % input_file )
  call hdf5_add( rootID, "INPUT_FILE_CRC32", info % input_file_crc32 )
  call hdf5_add( rootID, "TIMESTAMP", info % timestamp )

  call hdf5_add( rootID, "DT", info%dt )

  ndims = info % ndims
  call hdf5_add( rootID, "NDIMS", ndims )
  call hdf5_add( rootID, 'NX', info%box_nx ( 1:ndims) )
  call hdf5_add( rootID, 'XMIN', info%box_min ( 1:ndims) )
  call hdf5_add( rootID, 'XMAX', info%box_max ( 1:ndims) )
  call hdf5_add( rootID, 'PERIODIC', info%periodic( 1:ndims) )
  call hdf5_add( rootID, 'MOVE C', info%move_c ( 1:ndims) )

  ! add parallel partition information
  call hdf5_add( rootID, 'PAR_NODE_CONF', info%node_conf( 1:ndims) )
  call hdf5_add( rootID, 'PAR_NX_X1', info%nx_x1_node )
  if ( ndims > 1 ) then
    call hdf5_add( rootID, 'PAR_NX_X2', info%nx_x2_node )
    if ( ndims > 2 ) then
      call hdf5_add( rootID, 'PAR_NX_X3', info%nx_x3_node )
    endif
  endif

  call h5gclose_f( rootID, ierr )


end subroutine add_sim_info

!---------------------------------------------------------------------------------------------------
! Logical and string attributes need special treatment. Remaining interfaces are generated through
! template functions
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_str( objID, name, attribute )

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  character( len = * ), intent(in) :: attribute

  integer(hid_t) :: dataspaceID, typeID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(size_t) :: size

  integer :: ierr

  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr )

  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)


  size = len(attribute)
  call h5tset_size_f(typeID, size, ierr)

  call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, typeID, attribute, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5tclose_f( typeID, ierr )

  call h5sclose_f( dataspaceID, ierr )


end subroutine hdf5_add_str
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_str_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  character( len = * ), dimension(:), intent(in) :: attribute

  integer(hid_t) :: dataspaceID, typeID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(size_t) :: maxlen
  integer :: i, ierr

  dims(1) = size(attribute)
  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  maxlen = 0
  do i = 1, size(attribute)-1
    if (len(attribute(i)) > maxlen) maxlen = len(attribute(i))
  enddo

  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)

  call h5tset_size_f(typeID, maxlen, ierr)


  call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, typeID, attribute, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5tclose_f( typeID, ierr )

  call h5sclose_f( dataspaceID, ierr )


end subroutine hdf5_add_str_v1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_logical( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  logical, intent(in) :: attribute

  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: bool, ierr

  dims(1) = 1
  if ( attribute ) then
    bool = 1
  else
    bool = 0
  endif
  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5sclose_f( dataspaceID, ierr )


end subroutine hdf5_add_logical
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_logical_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  logical, dimension(:), intent(in) :: attribute

  integer(hid_t) :: dataspaceID, attrID
  integer :: i, ierr
  integer(hsize_t), dimension(1) :: dims
  integer, dimension(1) :: ldim
  integer, dimension(:), pointer :: bool

  ldim(1) = size(attribute)
  dims(1) = ldim(1)
  call alloc(bool, ldim,"dutil/os-dutil-hdf5.f90",602)
  do i = 1, int(dims(1))
    if ( attribute(i) ) then
      bool(i) = 1
    else
      bool(i) = 0
    endif
  enddo

  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5sclose_f( dataspaceID, ierr )


  call freemem(bool,"dutil/os-dutil-hdf5.f90",623)

end subroutine hdf5_add_logical_v1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Generate specific template functions for single, double and integer datatypes.
! Note that the module interfaces are not generated automatically and must be explicity written
! in the module header above.
!
! - hdf5_add
! - hdf5_add_v1
!---------------------------------------------------------------------------------------------------



! single precision real




# 1 "dutil/os-dutil-hdf5.f90" 1
!--------------------------------------------------------------------------------------------------
!
! Using parallel I/O:
!
! Chunking is a file property, so it is not possible to define different chunk sizes on different
! processes. This means that when the array size is different on some/any process the only options
! are to use normal (contiguous) storage, or (maybe) define a chunk size that is the maximum

! array size on every node.
!--------------------------------------------------------------------------------------------------

! if is defined then read the template definition at the end of the file
# 668 "dutil/os-dutil-hdf5.f90"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_single( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  real(p_single), intent(in) :: attribute

  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: ierr

  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  call h5acreate_f( objID, name, H5T_NATIVE_REAL, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, H5T_NATIVE_REAL, attribute, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5sclose_f( dataspaceID, ierr )


end subroutine hdf5_add_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_v1_single( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  real(p_single), dimension(:), intent(in) :: attribute

  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: ierr

  dims(1) = size(attribute)
  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  call h5acreate_f( objID, name, H5T_NATIVE_REAL, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, H5T_NATIVE_REAL, attribute, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5sclose_f( dataspaceID, ierr )


end subroutine hdf5_add_v1_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! End of template Functions, do not change below
!---------------------------------------------------------------------------------------------------
# 645 "dutil/os-dutil-hdf5.f90" 2

! double precision real




# 1 "dutil/os-dutil-hdf5.f90" 1
!--------------------------------------------------------------------------------------------------
!
! Using parallel I/O:
!
! Chunking is a file property, so it is not possible to define different chunk sizes on different
! processes. This means that when the array size is different on some/any process the only options
! are to use normal (contiguous) storage, or (maybe) define a chunk size that is the maximum

! array size on every node.
!--------------------------------------------------------------------------------------------------

! if is defined then read the template definition at the end of the file
# 668 "dutil/os-dutil-hdf5.f90"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_double( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  real(p_double), intent(in) :: attribute

  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: ierr

  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  call h5acreate_f( objID, name, H5T_NATIVE_DOUBLE, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, H5T_NATIVE_DOUBLE, attribute, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5sclose_f( dataspaceID, ierr )


end subroutine hdf5_add_double
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_v1_double( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  real(p_double), dimension(:), intent(in) :: attribute

  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: ierr

  dims(1) = size(attribute)
  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  call h5acreate_f( objID, name, H5T_NATIVE_DOUBLE, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, H5T_NATIVE_DOUBLE, attribute, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5sclose_f( dataspaceID, ierr )


end subroutine hdf5_add_v1_double
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! End of template Functions, do not change below
!---------------------------------------------------------------------------------------------------
# 652 "dutil/os-dutil-hdf5.f90" 2

! integer




# 1 "dutil/os-dutil-hdf5.f90" 1
!--------------------------------------------------------------------------------------------------
!
! Using parallel I/O:
!
! Chunking is a file property, so it is not possible to define different chunk sizes on different
! processes. This means that when the array size is different on some/any process the only options
! are to use normal (contiguous) storage, or (maybe) define a chunk size that is the maximum

! array size on every node.
!--------------------------------------------------------------------------------------------------

! if is defined then read the template definition at the end of the file
# 668 "dutil/os-dutil-hdf5.f90"
!---------------------------------------------------------------------------------------------------
! Template Function definitions
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_int( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  integer, intent(in) :: attribute

  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: ierr

  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, attribute, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5sclose_f( dataspaceID, ierr )


end subroutine hdf5_add_int
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine hdf5_add_v1_int( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  integer, dimension(:), intent(in) :: attribute

  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: ierr

  dims(1) = size(attribute)
  call h5screate_simple_f(1, dims, dataspaceID, ierr )


  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )

  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, attribute, dims, ierr)

  call h5aclose_f( attrID, ierr )

  call h5sclose_f( dataspaceID, ierr )


end subroutine hdf5_add_v1_int
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! End of template Functions, do not change below
!---------------------------------------------------------------------------------------------------
# 659 "dutil/os-dutil-hdf5.f90" 2

end module hdf5_util

!---------------------------------------------------------------------------------------------------
! end of hdf5_util module
!---------------------------------------------------------------------------------------------------
