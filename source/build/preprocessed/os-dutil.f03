# 1 "dutil/os-dutil.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "dutil/os-dutil.f03" 2

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
# 3 "dutil/os-dutil.f03" 2
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
# 4 "dutil/os-dutil.f03" 2

module m_diagnostic_utilities

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
# 8 "dutil/os-dutil.f03" 2

use m_parameters

use m_diagfile
use m_diagfile_zdf


use m_diagfile_hdf5


implicit none

private

! Precision for diagnostics
! some diagnostics (e.g. field grids) can be controlled in the input file using the prec
! parameter

integer, parameter :: p_diag_prec = p_single
!integer, parameter :: p_diag_prec = p_double

! Size of buffer for diagnostics (to be used by phasespaces & averaged grids)
# 48 "dutil/os-dutil.f03"
! (this sets it to 64 MB in single precision)
integer, parameter :: p_max_diag_buffer_size = 1024*1024*(64/4)



integer :: diag_buffer_size = -1
real(p_diag_prec), dimension(:), pointer :: diag_buffer


integer, parameter :: p_zdf_format = 0
integer, parameter :: p_hdf5_format = 1

! Default file format to use
integer, private :: file_format_default

interface get_filename
  module procedure get_filename
end interface

interface init_dutil
  module procedure init_diagutil
end interface

interface cleanup_dutil
  module procedure cleanup_diagutil
end interface

interface get_diag_buffer
  module procedure get_diag_buffer
end interface

interface diag_type_real
  module procedure diag_type_real
end interface

interface create_diag_file
  module procedure create_diag_file
end interface

public :: t_diag_file, create_diag_file
public :: p_diag_grid, p_diag_particles, p_diag_tracks
public :: get_filename

public :: init_dutil, cleanup_dutil, get_diag_buffer

public :: p_diag_prec, diag_type_real, p_max_diag_buffer_size

! pass though some parameters/symbols from os-diagfile for convenience
public :: diag_axis_linear, diag_axis_log10, diag_axis_log2
public :: t_diag_dataset, t_diag_chunk
public :: p_diag_create, p_diag_read, p_diag_update
public :: p_diag_float32, p_diag_float64, p_diag_int32, p_diag_int64, p_diag_null


public :: p_hdf5_format, p_zdf_format
public :: DIAG_MPI, DIAG_INDEPENDENT, DIAG_MPIIO_INDEPENDENT, DIAG_MPIIO_COLLECTIVE, DIAG_POSIX
public :: sim_info

contains

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
pure function diag_type_real( kind )

  implicit none

  integer :: diag_type_real
  integer, intent(in) :: kind

  select case (kind)
  case (p_double)
    diag_type_real = p_diag_float64
  case (p_single)
    diag_type_real = p_diag_float32
  case default
    diag_type_real = -1
  end select

end function diag_type_real

!--------------------------------------------------------------------------------------------------
! Initialize module variables
!--------------------------------------------------------------------------------------------------
subroutine init_diagutil( sim_options, buffer_size )


  use hdf5_util


  implicit none

  type( t_options ), intent(in) :: sim_options
  integer, intent(in) :: buffer_size


  ! Initialize hdf5 subsystem
  call open_hdf5( )


  ! Set default file format
  file_format_default = sim_options % file_format

  ! Set default parallel io algorithm
  call diag_file_set_def_iomode( sim_options % parallel_io )

  ! allocate diagnostics buffer
  if (buffer_size > p_max_diag_buffer_size) then
    diag_buffer_size = p_max_diag_buffer_size
  else
    diag_buffer_size = buffer_size
  endif

  if ( diag_buffer_size > 0 ) then
    call alloc(diag_buffer, (/diag_buffer_size/),"dutil/os-dutil.f03",161)
  endif

end subroutine init_diagutil
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
subroutine cleanup_diagutil( )
!--------------------------------------------------------------------------------------------------
! Cleanup module
!--------------------------------------------------------------------------------------------------
  implicit none

  if ( diag_buffer_size > 0 ) then
    call freemem(diag_buffer,"dutil/os-dutil.f03",175)
  endif

  call freemem(sim_info % nx_x1_node,"dutil/os-dutil.f03",178)
  call freemem(sim_info % nx_x2_node,"dutil/os-dutil.f03",179)
  call freemem(sim_info % nx_x3_node,"dutil/os-dutil.f03",180)

end subroutine cleanup_diagutil
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
subroutine get_diag_buffer( buffer, bsize )
!--------------------------------------------------------------------------------------------------
  implicit none

  real(p_diag_prec), dimension(:), pointer :: buffer
  integer, intent(out), optional :: bsize

  buffer => diag_buffer
  if ( present(bsize) ) then
    bsize = diag_buffer_size
  endif

end subroutine get_diag_buffer
!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
! Generate a file name from the information provided
!--------------------------------------------------------------------------------------------------
function get_filename( n, prefix, suffix, node )

  use stringutil

  implicit none

  ! arguments and return value

  character(len = 80) :: get_filename
  integer, intent(in) :: n
  character(len = *), intent(in) :: prefix

  character(len = *), intent(in), optional :: suffix
  integer, intent(in), optional :: node

  get_filename = trim(prefix)//'-'//trim(idx_string(n,p_time_length))

  if (present(node)) get_filename = trim(get_filename)// '.'//trim(idx_string(node,3))

  if (present(suffix)) get_filename = trim(get_filename)//trim(suffix)

end function get_filename
!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
! Return new diagfile
!--------------------------------------------------------------------------------------------------
subroutine create_diag_file(diagFile, type )

  implicit none

  class( t_diag_file ), allocatable :: diagFile
  integer, intent(in), optional :: type

  integer :: ltype

  if ( present(type) ) then
    ltype = type
  else
    ltype = file_format_default
  endif

  if ( allocated(diagFile) ) deallocate(diagFile)

  select case( ltype )
  case( p_zdf_format )
    allocate( t_diag_file_zdf :: diagFile )


  case( p_hdf5_format )
    allocate( t_diag_file_hdf5 :: diagFile )


  case default
    write(err_buf__,*) "Invalid diagnostics file format requested, aborting.";call err__("dutil/os-dutil.f03",260)
    call abort_program( p_err_invalid )
  end select

end subroutine
!--------------------------------------------------------------------------------------------------

end module m_diagnostic_utilities
