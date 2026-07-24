# 1 "tiles/os-bnd-tiles.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-bnd-tiles.f03" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary class for tiles module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
# 6 "tiles/os-bnd-tiles.f03" 2
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
# 7 "tiles/os-bnd-tiles.f03" 2

module m_bnd_tiles

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
# 11 "tiles/os-bnd-tiles.f03" 2

use m_system
use m_parameters
use m_species_define, only : t_part_idx, t_spec_msg
use m_vdf_comm, only : t_vdf_msg, cleanup
use m_vdf_comm_tiles, only : t_bnd_msg
use m_species_comm, only : init_spec_comm, cleanup

implicit none

private

! In t_sim_group
type :: t_bnd_group

  ! Replaces module variables for vdf comm
  ! NOW A LINKED LIST of size max(unique_nodes over dimensions)
  type( t_vdf_msg ), pointer :: send_vdf => null(), recv_vdf => null()

  ! Replaces module variables for part comm
  ! NOW A LINKED LIST
  ! This 2nd particle buffers for this list will be purged at each dlb
  type( t_spec_msg ), pointer :: send_spec => null(), recv_spec => null()

  ! Size of linked lists
  integer :: n_msg = 0

  contains

  procedure :: init => init_bnd_group
  procedure :: cleanup => cleanup_bnd_group

end type t_bnd_group

! In t_sim_group
! To be deallocated and reallocated after dlb
! One for each dimension
type :: t_bnd_patt

  ! number of unique communication nodes in current dimension
  integer :: n_msg = 0

  ! description of tiles that are on a particular node's communication boundary
  ! The messages are in some random order, it doesn't really matter
  class(t_bnd_msg), dimension(:), pointer :: msg => null()

  ! description of tiles that are NOT on a node's communication boundary
  ! (i.e. tiles involved in a local update boundary)
  class(t_bnd_msg), pointer :: msg_loc => null()

  contains

  procedure :: map_node_to_idx
  procedure :: cleanup => cleanup_bnd_patt

end type t_bnd_patt

public :: t_bnd_group, t_bnd_patt

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_bnd_group( this, patt, buffer_size )

  implicit none

  class( t_bnd_group ), intent(inout) :: this
  class( t_bnd_patt ), dimension(:), pointer :: patt
  integer, intent(in), optional :: buffer_size

  integer :: n_msg, i
  type( t_vdf_msg ), pointer :: snd_vdf, rcv_vdf, snd_vdf_prev, rcv_vdf_prev
  type( t_spec_msg ), pointer :: snd_spc, rcv_spc, snd_spc_prev, rcv_spc_prev

  ! Find out the maximum number of messages in any dimension
  n_msg = 0
  do i = 1, p_x_dim
    n_msg = max(n_msg,patt(i)%n_msg)
  enddo

  ! Change size (grow or shrink) of linked list arrays if necessary

  ! Initialize if first time here
  if (this%n_msg == 0 .and. n_msg > 0) then
    ! Sanity check
    if (associated(this%send_vdf)) then
      write(err_buf__,*) 'Should not be associated at this point';call err__("tiles/os-bnd-tiles.f03",98)
    endif
    allocate(this%send_vdf)
    allocate(this%recv_vdf)
    allocate(this%send_spec)
    allocate(this%recv_spec)

    if (present(buffer_size)) then
      call init_spec_comm( this%send_spec, this%recv_spec, buffer_size )
    else
      call init_spec_comm( this%send_spec, this%recv_spec )
    endif

  endif

  ! Continue initializing the rest that are needed
  snd_vdf => this%send_vdf
  rcv_vdf => this%recv_vdf
  snd_spc => this%send_spec
  rcv_spc => this%recv_spec

  ! Clear 2nd particle buffer of current spec_msg to save memory
  if (associated(snd_spc)) then
    if (associated(snd_spc%buffer2)) call freemem(snd_spc%buffer2,"tiles/os-bnd-tiles.f03",121)
    if (associated(rcv_spc%buffer2)) call freemem(rcv_spc%buffer2,"tiles/os-bnd-tiles.f03",122)
  endif

  do i = 1, n_msg-1

    snd_vdf_prev => snd_vdf
    rcv_vdf_prev => rcv_vdf
    snd_spc_prev => snd_spc
    rcv_spc_prev => rcv_spc

    if (.not. associated(snd_vdf_prev%next)) then
      allocate(snd_vdf)
      allocate(rcv_vdf)
      allocate(snd_spc)
      allocate(rcv_spc)

      if (present(buffer_size)) then
        call init_spec_comm( snd_spc, rcv_spc, buffer_size )
      else
        call init_spec_comm( snd_spc, rcv_spc )
      endif

      snd_vdf_prev%next => snd_vdf
      rcv_vdf_prev%next => rcv_vdf
      snd_spc_prev%next => snd_spc
      rcv_spc_prev%next => rcv_spc

    else
      snd_vdf => snd_vdf_prev%next
      rcv_vdf => rcv_vdf_prev%next
      snd_spc => snd_spc_prev%next
      rcv_spc => rcv_spc_prev%next
    endif

    ! Clear 2nd particle buffer of current spec_msg to save memory
    if (associated(snd_spc%buffer2)) call freemem(snd_spc%buffer2,"tiles/os-bnd-tiles.f03",157)
    if (associated(rcv_spc%buffer2)) call freemem(rcv_spc%buffer2,"tiles/os-bnd-tiles.f03",158)

  enddo

  ! Deallocate extras to save memory
  if (this%n_msg > n_msg) then
    call cleanup( snd_vdf%next )
    deallocate( snd_vdf%next )
    snd_vdf%next => null()
    call cleanup( rcv_vdf%next )
    deallocate( rcv_vdf%next )
    rcv_vdf%next => null()

    call cleanup( snd_spc%next )
    deallocate( snd_spc%next )
    snd_spc%next => null()
    call cleanup( rcv_spc%next )
    deallocate( rcv_spc%next )
    rcv_spc%next => null()
  endif

  this%n_msg = n_msg

end subroutine init_bnd_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cleanup_bnd_group( this )

  implicit none

  class( t_bnd_group ), intent(inout) :: this

  if (this%n_msg>0) then

    call cleanup( this%send_vdf )
    call cleanup( this%recv_vdf )
    deallocate( this%send_vdf )
    deallocate( this%recv_vdf )

    call cleanup( this%send_spec )
    call cleanup( this%recv_spec )
    deallocate( this%send_spec )
    deallocate( this%recv_spec )

    this%n_msg = 0

  endif

end subroutine cleanup_bnd_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function map_node_to_idx( this, node )

  class(t_bnd_patt), intent(in) :: this
  integer, intent(in) :: node

  integer :: map_node_to_idx
  integer :: i

  map_node_to_idx = 0

  do i = 1, this%n_msg
    if (this%msg(i)%node == node) then
      map_node_to_idx = i
      exit
    endif
  enddo

  ! Sanity check
  if (map_node_to_idx == 0) then
    write(err_buf__,*) 'Node not found in map_node_to_idx';call err__("tiles/os-bnd-tiles.f03",232)
  endif

end function map_node_to_idx
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cleanup_bnd_patt( this )

  implicit none

  class(t_bnd_patt), intent(inout) :: this

  integer :: i

  ! cleanup this%msg
  if ( this%n_msg>0 ) then
    do i = 1, this%n_msg
      call this%msg(i)%cleanup()
    enddo
    deallocate(this%msg)
  endif

  ! cleanup this%msg_loc
  if ( associated(this%msg_loc) ) then
    call this%msg_loc%cleanup()
    deallocate(this%msg_loc)
  endif

end subroutine cleanup_bnd_patt
!-----------------------------------------------------------------------------------------

end module m_bnd_tiles
