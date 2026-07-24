# 1 "tiles/os-spec-rawdiag-tiles.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-spec-rawdiag-tiles.f90" 2
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
# 2 "tiles/os-spec-rawdiag-tiles.f90" 2
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
# 3 "tiles/os-spec-rawdiag-tiles.f90" 2

module m_species_rawdiag_tiles

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
# 7 "tiles/os-spec-rawdiag-tiles.f90" 2

use m_parameters
use m_species_define, only : t_spec_arr
use m_species_rawdiag
use m_node_conf, only : t_node_conf
use m_diagnostic_utilities, only : t_diag_file, t_diag_dataset, t_diag_chunk, p_diag_prec

private

interface write_raw
  module procedure write_raw_tiles
end interface

public :: write_raw

contains

!-----------------------------------------------------------------------------------------
! write raw particle data for tiles
!-----------------------------------------------------------------------------------------
subroutine write_raw_tiles( spec, no_co, n, t, n_aux )

  implicit none

  type( t_spec_arr ), dimension(:), pointer :: spec
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: n
  real(p_double), intent(in) :: t
  integer, intent(in) :: n_aux

  class( t_diag_file ), allocatable :: diagFile
  type( t_diag_dataset ), dimension( spec(1)%s%get_n_x_dims() + p_p_dim + 3 ) :: datasets
  type( t_diag_chunk ) :: chunk, chunk_tag
  integer :: start, rank, num_par, total_par, color, comm, ierr
  integer :: ndims, i, j, k, l, ntil
  integer, dimension(:), pointer :: idx, num_par_til, np_til_scan
  real( p_diag_prec ), dimension(:), pointer :: write_buffer
  integer, dimension( :, : ), pointer :: write_tag_buffer
  real( p_k_part ) :: u2

  idx => null()

  ! total_par = number of particles totaled across all tiles and nodes
  ! num_par = number of particles on all tiles on this node
  ! num_par_til = array of number of particles on each tile on this node

  ndims = spec(1)%s%get_n_x_dims()
  ntil = size(spec)

  ! get total number of particles
  num_par = 0
  do i = 1, ntil
    num_par = num_par + spec(i)%s%num_par
  enddo

  ! Here num_par refers to total number of particles
  ! In this loop we set it to the number of raw particles to dump
  if (num_par > 0) then

    ! idx = array for all particles of all tiles on this node
    ! num_par_til = raw particle number for each tile
    ! np_til_scan = starting index in dix corresponding to each tile
    call alloc(idx, (/num_par/),"tiles/os-spec-rawdiag-tiles.f90",69)
    call alloc(num_par_til, (/ntil/),"tiles/os-spec-rawdiag-tiles.f90",70)
    call alloc(np_til_scan, (/ntil/),"tiles/os-spec-rawdiag-tiles.f90",71)

    ! get index of particles to output
    num_par = 0
    do i = 1, ntil
      np_til_scan(i) = num_par + 1
      if ( spec(i)%s%num_par > 0 ) then
        call write_raw_select_data( spec(i)%s, t, idx(num_par+1:), num_par_til(i) )
        num_par = num_par + num_par_til(i)
      else
        num_par_til(i) = 0
      endif
    enddo

  endif

  if ( no_co % comm_size() > 1 ) then
    call MPI_ALLREDUCE( num_par, total_par, 1, MPI_INTEGER, MPI_SUM, no_co % comm, ierr )

    if ( total_par > 0 ) then
      if ( num_par > 0 ) then
        color = 1
      else
        color = MPI_UNDEFINED
      endif
      call MPI_COMM_SPLIT( no_co % comm, color, 0, comm, ierr )

      if ( num_par > 0 ) then
        call MPI_COMM_RANK( comm, rank, ierr )
        call MPI_EXSCAN( num_par, start, 1, MPI_INTEGER, MPI_SUM, comm, ierr )
        if ( rank == 0 ) start = 0
      endif
    else
      comm = MPI_COMM_NULL
    endif
  else
    comm = MPI_COMM_SELF
    total_par = num_par
    start = 0
  endif

  ! Write data, if any
  if ( total_par > 0 ) then

    ! Only nodes with particles are involved in this section
    if ( num_par > 0 ) then

      ! Create the file
      call write_raw_create_file( spec(1)%s, total_par, n, t, n_aux, diagFile, comm )

      ! Create the datasets
      call create_datasets( diagFile, total_par, datasets )

      call alloc(write_buffer, [ num_par ],"tiles/os-spec-rawdiag-tiles.f90",124)
      chunk % start(1) = start
      chunk % count(1) = num_par
      chunk % stride(1) = 1
      chunk % data = c_loc( write_buffer )

      ! Positions
      k = 1
      do j = 1, ndims
        do l = 1, ntil
          if ( num_par_til(l) > 0 ) then
            call spec(l)%s % get_position( j, idx(np_til_scan(l):), num_par_til(l), &
                                           write_buffer(np_til_scan(l):), &
                                           spec(l)%s%diag%raw_ref_box(j) )
          endif
        enddo
        call diagFile % write_par_cdset( datasets(k), chunk, start )
        k = k+1
      enddo

      ! Momenta
      do j = 1, p_p_dim
        do l = 1, ntil
          do i = np_til_scan(l), np_til_scan(l) + num_par_til(l) - 1
            write_buffer(i) = real(spec(l)%s%p(j,idx(i)), p_diag_prec)
          enddo
        enddo
        call diagFile % write_par_cdset( datasets(k), chunk, start )
        k = k+1
      enddo

      ! Charge
      do l = 1, ntil
        do i = np_til_scan(l), np_til_scan(l) + num_par_til(l) - 1
          write_buffer(i) = real(spec(l)%s%q(idx(i)), p_diag_prec)
        enddo
      enddo
      call diagFile % write_par_cdset( datasets(k), chunk, start )
      k = k+1

      ! Energy
      do l = 1, ntil
        do i = np_til_scan(l), np_til_scan(l) + num_par_til(l) - 1
          u2 = spec(l)%s%p( 1, idx(i) )**2 + spec(l)%s%p( 2, idx(i) )**2 + spec(l)%s%p( 3, idx(i) )**2
          write_buffer(i) = real(u2 / ( sqrt(1 + u2) + 1) , p_diag_prec)
        enddo
      enddo
      call diagFile % write_par_cdset( datasets(k), chunk, start )

      call freemem(write_buffer,"tiles/os-spec-rawdiag-tiles.f90",173)

      ! Tags
      if ( spec(1)%s%add_tag ) then
        call alloc(write_tag_buffer, [2, num_par],"tiles/os-spec-rawdiag-tiles.f90",177)
        chunk_tag % count(1) = 2
        chunk_tag % count(2) = num_par
        chunk_tag % start(1) = 0
        chunk_tag % start(2) = start
        chunk_tag % stride(1) = 1
        chunk_tag % stride(2) = 1
        chunk_tag % data = c_loc(write_tag_buffer)

        do l = 1, ntil
          do i = np_til_scan(l), np_til_scan(l) + num_par_til(l) - 1
            write_tag_buffer(1,i) = spec(l)%s%tag(1,idx(i))
            write_tag_buffer(2,i) = spec(l)%s%tag(2,idx(i))
          enddo
        enddo

        k = k+1
        call diagFile % write_par_cdset( datasets(k), chunk_tag, 2*start )

        call freemem(write_tag_buffer,"tiles/os-spec-rawdiag-tiles.f90",196)

      endif

      call close_datasets( diagFile, datasets )

      ! close file
      call diagFile % close( )
      deallocate( diagFile )
    endif
  else
    ! No particles, root node writes an empty file
    if ( no_co % local_rank() == 0 ) then
      call write_raw_create_file( spec(1)%s, 0, n, t, n_aux, diagFile, MPI_COMM_NULL )
      call create_datasets( diagFile, num_par, datasets )
      call close_datasets( diagFile, datasets )
      call diagFile % close( )
      deallocate( diagFile )
  endif
  endif

  ! Free particle indexes
  if ( associated(idx) ) then
    call freemem(idx,"tiles/os-spec-rawdiag-tiles.f90",219)
    call freemem(num_par_til,"tiles/os-spec-rawdiag-tiles.f90",220)
    call freemem(np_til_scan,"tiles/os-spec-rawdiag-tiles.f90",221)
  endif

end subroutine write_raw_tiles
!-----------------------------------------------------------------------------------------

end module m_species_rawdiag_tiles
