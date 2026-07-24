# 1 "tiles/os-tconf.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-tconf.f03" 2
! ----
! tile configuration class
! ----

!#define DEBUG_FILE 1

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
# 8 "tiles/os-tconf.f03" 2
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
# 9 "tiles/os-tconf.f03" 2

module m_tile_conf

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
# 13 "tiles/os-tconf.f03" 2

use m_parameters
use m_system

implicit none

private

! topology calculation
integer, parameter :: p_topology_hilbert = 2
integer, parameter :: p_topology_snake = 3
integer, parameter :: p_topology_viktor = 4

type :: int_array_4d
  integer, dimension(:,:,:,:), pointer :: ii => null()
end type

type :: t_tile_conf

  ! total number of tiles in global simulation
  integer :: g_til_num

  ! number of tiles in each dimension over the entire simulation
  integer, dimension(p_x_dim) :: g_til

  ! global array identification number for this tile
  integer :: g_til_aid

  ! T_ile G_rid P_osition
  ! globally across the whole simulation
  integer, dimension(p_max_dim) :: g_tgp

  ! array identifying number of neighboring tiles
  ! Local and global numbering, respectively
  integer, dimension( 2, p_max_dim ) :: neighbor_g_til

  ! Type of topology to use for tiles ( hilbert, snake, viktor )
  integer :: topology_type

  ! Space-filling operator: g_til_aid maps to coordinates in the global grid of tiles
  ! dimensions are p_x_dim, g_til_num
  integer, dimension(:,:), pointer :: h => null()
  integer, dimension(:,:), pointer :: h_old => null()

  ! Inverse space-filling operator: tile's coordiantes map to g_til_aid
  integer, dimension(:), pointer :: hi1 => null()
  integer, dimension(:,:), pointer :: hi2 => null()
  integer, dimension(:,:,:), pointer :: hi3 => null()
  integer, dimension(:), pointer :: hi1_old => null()
  integer, dimension(:,:), pointer :: hi2_old => null()
  integer, dimension(:,:,:), pointer :: hi3_old => null()

  ! Min/max tile g_aid's for all nodes
  ! index 1 is dimension 3: lower/upper/num tiles
  ! index 2 is node aid
  integer, dimension(:,:), pointer :: g_til_aid_min_max => null()
  integer, dimension(:,:), pointer :: g_til_aid_min_max_old => null()

  ! Min/max tiles for each node in given simualtion direction
  type(int_array_4d), dimension(:), pointer :: cart_til_mm => null()
  type(int_array_4d), dimension(:), pointer :: cart_til_mm_old => null()

  ! describes rectangular region of tiles in inverse space filling operator (hi1/2/3)
  ! that node is responsible for for diagnostics
  ! - index 1 is lower/upper/#tiles
  ! - index 2 is spatial dimension
  integer, dimension(3, p_max_dim) :: my_nt_io

contains

  procedure :: cleanup => cleanup_tile_conf

end type t_tile_conf

interface vis_hid
  module procedure vis_hid
end interface

interface g_tgp
  module procedure g_tgp_tile
end interface

interface g_tgp_old
  module procedure g_tgp_old_tile
end interface

interface set_ti_co_old
  module procedure set_ti_co_old
end interface

public :: t_tile_conf, int_array_4d
public :: g_tgp, g_tgp_old
public :: vis_hid
public :: set_ti_co_old

public :: p_topology_hilbert
public :: p_topology_snake
public :: p_topology_viktor

contains

!-----------------------------------------------------------------------------------------
! print tile topology
!-----------------------------------------------------------------------------------------
subroutine vis_hid( this )

  implicit none

  type( t_tile_conf ), intent(in) :: this

  integer :: i, k

  select case (p_x_dim)
  case (1)
    print *, this%hi1(:)
  case (2)
    do i = this%g_til(2),1,-1
      print *, this%hi2(:,i)
    enddo
  case (3)
    do k = this%g_til(3),1,-1
      print *, 'k =', k
      do i = this%g_til(2),1,-1
        print *, this%hi3(:,i,k)
      enddo
    enddo
  end select

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! t_ile g_rid p_osition: result is coords of specified tile in the grid of tiles
!-----------------------------------------------------------------------------------------
function g_tgp_tile( this, til, dim )

  implicit none

  type( t_tile_conf ), intent(in) :: this
  integer, intent(in) :: til, dim
  integer :: g_tgp_tile

  g_tgp_tile=this%h(dim,til)

end function g_tgp_tile
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! T_ile G_rid P_osition: result is coords of specified tile in the old grid of tiles
!-----------------------------------------------------------------------------------------
function g_tgp_old_tile( this, til, dim )

  implicit none

  type( t_tile_conf ), intent( in ) :: this
  integer, intent(in) :: til, dim
  integer :: g_tgp_old_tile

  g_tgp_old_tile=this%h_old(dim,til)

end function g_tgp_old_tile
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_tile_conf( this )
!-----------------------------------------------------------------------------------------
! Object destructor, nullifies pointers (deallocated in cleanup_node_conf)
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_tile_conf ), intent( inout ) :: this
  integer :: i

  if ( associated(this%h) ) then
    call freemem(this%h,"tiles/os-tconf.f03",188)
  endif

  if ( associated(this%h_old) ) then
    call freemem(this%h_old,"tiles/os-tconf.f03",192)
  endif

  if ( associated(this%hi1) ) then
    call freemem(this%hi1,"tiles/os-tconf.f03",196)
  endif

  if ( associated(this%hi1_old) ) then
    call freemem(this%hi1_old,"tiles/os-tconf.f03",200)
  endif

  if ( associated(this%hi2) ) then
    call freemem(this%hi2,"tiles/os-tconf.f03",204)
  endif

  if ( associated(this%hi2_old) ) then
    call freemem(this%hi2_old,"tiles/os-tconf.f03",208)
  endif

  if ( associated(this%hi3) ) then
    call freemem(this%hi3,"tiles/os-tconf.f03",212)
  endif

  if ( associated(this%hi3_old) ) then
    call freemem(this%hi3_old,"tiles/os-tconf.f03",216)
  endif

  if ( associated(this%g_til_aid_min_max) ) then
    call freemem(this%g_til_aid_min_max,"tiles/os-tconf.f03",220)
  endif

  if ( associated(this%g_til_aid_min_max_old) ) then
    call freemem(this%g_til_aid_min_max_old,"tiles/os-tconf.f03",224)
  endif

  if ( associated(this%cart_til_mm) ) then
    do i = 1, p_x_dim
      deallocate( this%cart_til_mm(i)%ii )
    enddo
    deallocate( this%cart_til_mm )
    this%cart_til_mm => null()
  endif

  if ( associated(this%cart_til_mm_old) ) then
    do i = 1, p_x_dim
      deallocate( this%cart_til_mm_old(i)%ii )
    enddo
    deallocate( this%cart_til_mm_old )
    this%cart_til_mm_old => null()
  endif

end subroutine cleanup_tile_conf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Store old tile configuration data
!-----------------------------------------------------------------------------------------
subroutine set_ti_co_old( this )

  implicit none

  type( t_tile_conf ), intent(inout) :: this
  integer :: i

  this%g_til_aid_min_max_old = this%g_til_aid_min_max

  if ( this%topology_type==p_topology_viktor ) then
    this%h_old = this%h
    do i = 1, p_x_dim
      this%cart_til_mm_old(i)%ii = this%cart_til_mm(i)%ii
    enddo
    select case( p_x_dim )
      case(1)
        this%hi1_old = this%hi1
      case(2)
        this%hi2_old = this%hi2
      case(3)
        this%hi3_old = this%hi3
    end select
  endif

end subroutine set_ti_co_old
!-----------------------------------------------------------------------------------------

end module m_tile_conf
