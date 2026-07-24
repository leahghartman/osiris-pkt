! functions supporting parallel grid information

#include "os-preprocess.fpp"

module m_grid_tiles

#include "memory/memory.h"

use m_system
use m_parameters
use m_grid_define, only : t_grid, t_msg_patt, t_msg, dynlb_sum_up_load_ev, &
                          p_dynamic_load_balance, p_expr_load_balance, p_no_load_balance,&
                          p_static_load_balance, p_test_load_balance, p_none, &
                          p_point2point, p_gather, p_global, p_node, p_grid, alloc
use m_grid_parallel, only : even_dist, load_dist, add_expression_load
use m_fparser
use m_input_file
use m_tile_conf, only : p_topology_hilbert, p_topology_snake, p_topology_viktor, g_tgp
use m_tile_conf, only : int_array_4d
use m_logprof, only : create_event
use m_node_conf, only : t_node_conf, my_ngp, no_num, nx, p_max, p_sum, reduce_array, n_threads, root
use m_node_conf_tiles, only : t_node_conf_tiles, set_ti_co_viktor, p_light
use m_restart


implicit none

private

integer, parameter :: p_send = 1
integer, parameter :: p_recv = -1

integer, parameter :: p_tile_node = 4
integer, parameter :: p_tile_load = 5
integer, parameter :: p_hl_tils   = 6

type, extends(t_grid) :: t_grid_tiles

  ! frequency to dump tile node info (node for each tile)
  integer :: ndump_tile_node   = 0
  integer :: ndump_tile_load   = 0
  integer :: ndump_heavy_light = 0

contains

  procedure :: read_input      => read_input_grid_tiles
  procedure :: init            => init_grid_tiles
  procedure :: cleanup         => cleanup_grid_tiles
  procedure :: if_report       => if_report_grid_tiles
  procedure :: ndump           => ndump_report_grid_tiles
  procedure :: test_partition  => test_partition_grid_tiles
  procedure :: get_node_limits => get_node_limits_tiles
  procedure :: parallel_partition_tiles

end type t_grid_tiles

type, extends(t_grid_tiles) :: t_grid_group

  ! number of cells for each tile for all directions
  ! partition( bound, part_idx )
  ! - bound is 1:3 (lower/upper/#cells)
  ! - part_idx is organized as [ 1 .. nt1, 1 .. nt2, 1 .. nt3 ]
  !   where nt1 (2/3) is num tiles in dimension 1 (2/3)
  integer, pointer, dimension(:,:) :: nx_til => null()

  ! arrays holding the integral of the load
  real(p_double), pointer, dimension(:,:)   :: til_load_2 => null()
  real(p_double), pointer, dimension(:,:,:) :: til_load_3 => null()

contains

  procedure :: init => init_grid_group
  procedure :: cleanup => cleanup_grid_group
  procedure :: gather_load => gather_load_grid_group
  procedure :: parallel_partition_group
  procedure :: parallel_partition_io

end type t_grid_group

interface new_dlb_patt
  module procedure new_dlb_patt_
  module procedure new_dlb_patt_viktor
end interface

interface init_tile_load
  module procedure init_tile_load
end interface

interface add_load
  module procedure add_grid_load_tiles
end interface

interface conv_load_2_int_load
  module procedure conv_load_2_int_load_tiles
  module procedure conv_load_2_int_load_tiles_1
  module procedure conv_load_2_int_load_tiles_2
end interface

interface project
  module procedure project_2
  module procedure project_3
end interface

interface load_dist
  module procedure load_dist_viktor
end interface

interface even_dist
  module procedure even_dist_viktor
end interface

interface unique_count_1d
  module procedure unique_count_1d
end interface

interface unique_count_2d
  module procedure unique_count_2d
end interface

public :: new_dlb_patt, init_tile_load, add_load, conv_load_2_int_load
public :: load_dist, even_dist, unique_count_1d, unique_count_2d

public :: t_grid_tiles, t_grid_group, p_tile_node, p_tile_load, p_hl_tils

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_grid_tiles( this, input_file, partition )

  implicit none

  class( t_grid_tiles ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in), dimension(:) :: partition

  ! local variables

  integer, dimension(p_x_dim) :: nx_p
  character(len=16)               :: coordinates

  logical, dimension(p_x_dim) :: load_balance

  logical :: any_lb

  character(len=16) :: lb_type, lb_gather
  integer :: n_dynamic, start_load_balance
  logical :: balance_on_start
  real( p_single ) :: max_imbalance, cell_weight

  integer :: ndump_global_load, ndump_node_load, ndump_grid_load
  integer :: ndump_tile_node, ndump_tile_load, ndump_heavy_light
  character (len=p_max_expr_len) :: spatial_loaddensity

  integer, dimension(p_x_dim) :: io_nmerge
  character(len=16) :: io_merge_type


  namelist /nl_grid/ nx_p, coordinates, io_nmerge, io_merge_type, &
                     load_balance, lb_type, lb_gather, n_dynamic, start_load_balance, &
                     balance_on_start, max_imbalance, cell_weight, &
                     ndump_global_load, ndump_node_load, ndump_grid_load, &
                     ndump_tile_node, ndump_tile_load, ndump_heavy_light, &
                     spatial_loaddensity


  integer :: i, ierr

  ! executable statements

  this%x_dim = p_x_dim


  nx_p = 0
  coordinates = "cartesian"

  io_merge_type = "point2point"
  io_nmerge = 1

  load_balance = .false.
  lb_type   = "none"
  lb_gather = "sum"
  n_dynamic = -1
  spatial_loaddensity = "NO FUNCTION DEFINED"

  ! Default is to start load balance at timestep 0 (initialization)
  start_load_balance = -1
  balance_on_start = .false.

  max_imbalance = 0.0
  cell_weight   = 0.0

  ndump_global_load = 0
  ndump_node_load   = 0
  ndump_grid_load   = 0
  ndump_tile_node   = 0
  ndump_tile_load   = 0
  ndump_heavy_light = 0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_grid", ierr )

  if ( ierr /= 0 ) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading grid parameters"
      else
        write(0,*) "Error - grid parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_grid, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading grid parameters"
      write(0,*) "aborting..."
    endif
    stop
  endif

  ! validate global grid parameters
  do i = 1, this%x_dim
    if (nx_p(i) <= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "Invalid grid parameters, nx_p must be >= 0 for all dimensions."
        write(0,*) "nx_p(",p_x_dim,") = ", nx_p
        write(0,*) "(also check the number of dimensions of the binary)"
        write(0,*) "aborting..."
      endif
      stop
    endif
  enddo

  ! store global grid size
  this%g_nx(1:this%x_dim) = nx_p(1:this%x_dim)

  select case ( trim( coordinates ) )
  case ( "cylindrical" )
    if (disp_out(input_file)) then
      SCR_ROOT('cylindrical coordinates - B1 on axis')
    endif
    this%coordinates = p_cylindrical_b
    if ( p_x_dim /= 2 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "Invalid coordinates, cylindrical coordinates "
        write(0,*) "are only available in 2D."
      endif
      stop
    endif

  case ( "cartesian" )
    if (disp_out(input_file)) then
      SCR_ROOT('cartesian coordinates')
    endif
    this%coordinates = p_cartesian

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error in grid parameters"
      write(0,*) "Invalid coordinates, coordinates must be:"
      write(0,*) "'cartesian' (1D, 2D, 3D) or 'cylindrical' (2D)"
    endif
    stop

  end select

  this%start_load_balance = start_load_balance
  this%balance_on_start = balance_on_start

  select case (trim(lb_type))
  case("none")
    this%lb_type = p_no_load_balance
    this%load_balance = .false.

  case("test")
    this%lb_type = p_test_load_balance
    this%load_balance = .true.
    this%n_dynamic = n_dynamic

  case("static")
    this%lb_type = p_static_load_balance
    this%cell_weight = cell_weight

    ! The start load balance parameter is ignored, since this occurs a timestep 0
    this%start_load_balance = -1

  case("dynamic")
    this%lb_type = p_dynamic_load_balance
    if (n_dynamic < 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "n_dynamic must be >= 1 when using dynamic load balancing"
        write(0,*) "aborting..."
      endif
      stop
    endif
    this%n_dynamic = n_dynamic

    this%max_imbalance = max_imbalance
    this%cell_weight = cell_weight

    select case (trim( lb_gather ))
    case ("sum")
      this%lb_gather_max = .false.
    case ("max")
      if ( mpi_node() == 0 ) then
        write(0,*) "Warning, the lb_gather grid parameter is ignored when using tiles."
        write(0,*) "'sum' is the default since load is calculated for each tile."
        ! Note, using sum or max actually yeilds equivalent results for tiling
      endif
      this%lb_gather_max = .false.
    case default
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "lb_gather must be 'sum' or 'max' when using dynamic load balancing"
        write(0,*) "aborting..."
      endif
      stop
    end select

  case("expression")
    this%lb_type = p_expr_load_balance
    this%cell_weight = cell_weight

    select case (p_x_dim)
    case (1)
      call setup(this%spatial_loaddensity, trim(spatial_loaddensity), (/'x1'/), ierr)
    case (2)
      call setup(this%spatial_loaddensity, trim(spatial_loaddensity), (/'x1','x2'/), ierr)
    case (3)
      call setup(this%spatial_loaddensity, trim(spatial_loaddensity), (/'x1','x2','x3'/), ierr)
    end select

    ! check if function compiled ok
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error compiling spatial_loaddensity"
        write(0,*) "aborting..."
      endif
      stop
    endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error in grid parameters"
      write(0,*) "Invalid lb_type selected."
      write(0,*) "aborting..."
    endif
    stop
  end select

  if ( this%lb_type /= p_no_load_balance ) then

    ! Store load balance directions
    any_lb = .false.
    do i = 1, this%x_dim
      this%load_balance(i) = load_balance(i)
      if ( load_balance(i) ) any_lb = .true.
    enddo

    if ( .not. any_lb ) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "Error in grid parameters"
        write(0,*) "Load balance type '", trim(lb_type), "' selected, but no load balance"
        write(0,*) "direction specified"
        write(0,*) "aborting..."
      endif
      stop
    endif

    ! store cell weight
    this%cell_weight = cell_weight
    if ( cell_weight < 0.0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "Error in grid parameters"
        write(0,*) "cell_weight must be >= 0."
        write(0,*) "aborting..."
      endif
      stop
    endif
  endif


  ! grid reports
  this%ndump_global_load   = ndump_global_load
  this%ndump_node_load     = ndump_node_load
  this%ndump_grid_load     = ndump_grid_load
  this%ndump_tile_node     = ndump_tile_node
  this%ndump_tile_load     = ndump_tile_load
  this%ndump_heavy_light   = ndump_heavy_light

  ! Data merging for file output
  select case (trim(io_merge_type))
  case("none")
    this % io % merge_type = p_none
  case("point2point")
    this % io % merge_type = p_point2point
  case("gather")
    this % io % merge_type = p_gather
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error in grid parameters"
      write(0,*) "io_merge_type must be one of 'none', 'point2point', or 'gather'"
      write(0,*) "aborting..."
    endif
    stop
  end select

  do i = 1, this%x_dim

    if ( io_nmerge(i) < 1 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "io_nmerge(:) must be >= 1 in all directions."
        write(0,*) "aborting..."
      endif
      stop
    endif

    if ( mod(partition(i), io_nmerge(i) ) /= 0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "io_nmerge(:) must divide number of nodes in parallel partition exactly."
        write(0,*) "aborting..."
      endif
      stop
    endif

    this % io % n_merge(i) = io_nmerge(i)
  enddo

  if ( this % io % merge_type == p_none .and. product(io_nmerge(1:this%x_dim)) > 1 ) then
    if (mpi_node() == 0) then
      write(0,'(A,I0,A)') " (*warning*) io_nmerge(*) is > 1 but io_merge_type was set to 'none'"
      write(0,'(A)')      " (*warning*) no data merging for I/O will be performed"
    endif
  endif


end subroutine read_input_grid_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Setup grid object for t_simulation_tiles
!-----------------------------------------------------------------------------------------
subroutine init_grid_tiles( this, no_co, gc_num, g_box, restart, restart_handle )

  implicit none

  ! dummy variables
  class( t_grid_tiles ), intent(inout) :: this
  class( t_node_conf ), intent(in)  :: no_co    ! node configuration

  integer, dimension(:, :), intent(in) :: gc_num
  real(p_double), dimension(:,:), intent(in) :: g_box

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle


  ! local variables
  integer :: i, x_dim_

  if ( .not. restart ) then

    x_dim_ = this%x_dim
    select type(no_co); class is(t_node_conf_tiles)
      this%nnodes(1:x_dim_) = no_co%ti_co%g_til
    class default
      ERROR('Should only have t_node_conf_tiles type here')
      call abort_program()
    end select

    do i = 1, x_dim_
      this%g_box( :, i ) = g_box( :, i )

      ! Minimum number of cells per node
      this%min_nx(i) = max( 1, gc_num(p_lower,i)+gc_num(p_upper,i))
    enddo

  else

    call this % read_checkpoint( restart_handle )

  endif

  ! create event for timing
  if (dynlb_sum_up_load_ev < 0) dynlb_sum_up_load_ev = create_event('dynlb sum up load')

end subroutine init_grid_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! clear object data
!-----------------------------------------------------------------------------------------
subroutine cleanup_grid_tiles(this)

  implicit none

  class( t_grid_tiles ), intent(inout) :: this

  ! Parent class stuff, slightly modified
  ! call freemem( this%nx_node )
  call freemem( this%int_load )

  ! TODO: remove this?
  call this % io % cleanup()

end subroutine cleanup_grid_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get dump frequency of required diagnostic
!-----------------------------------------------------------------------------------------
function ndump_report_grid_tiles( this, rep_type )
  
  implicit none
  
  class( t_grid_tiles ), intent(in) :: this
  integer, intent(in) :: rep_type
  integer :: ndump_report_grid_tiles
  
  select case( rep_type )
    case ( p_global )
      ndump_report_grid_tiles = this%ndump_global_load
    case ( p_node )
      ndump_report_grid_tiles = this%ndump_node_load
    case ( p_grid )
      ndump_report_grid_tiles = this%ndump_grid_load
    case ( p_tile_node )
      ndump_report_grid_tiles = this%ndump_tile_node
    case ( p_tile_load )
      ndump_report_grid_tiles = this%ndump_tile_load
    case ( p_hl_tils )
      ndump_report_grid_tiles = this%ndump_heavy_light
    case default
      ndump_report_grid_tiles = 0
  end select
    
end function ndump_report_grid_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Check if load diagnostics are required at this timestep
!-----------------------------------------------------------------------------------------
function if_report_grid_tiles( this, n, rep_type )
  
  implicit none
  
  class( t_grid_tiles ), intent(in) :: this
  integer, intent(in) :: n, rep_type
  logical :: if_report_grid_tiles
  
  select case( rep_type )
    case ( p_global )
      if_report_grid_tiles = test_if_report( n, this%ndump_global_load )
    case ( p_node )
      if_report_grid_tiles = test_if_report( n, this%ndump_node_load )
    case ( p_grid )
      if_report_grid_tiles = test_if_report( n, this%ndump_grid_load )
    case ( p_tile_node )
      if_report_grid_tiles = test_if_report( n, this%ndump_tile_node )
    case ( p_tile_load )
      if_report_grid_tiles = test_if_report( n, this%ndump_tile_load )
    case ( p_hl_tils )
      if_report_grid_tiles = test_if_report( n, this%ndump_heavy_light )
    case default
      if_report_grid_tiles = .false.
  end select
  
contains

  function test_if_report( n, ndump )
    
    implicit none
    
    integer, intent(in) :: n, ndump
    logical :: test_if_report
    
    if ( ndump > 0 ) then
      test_if_report = ( mod( n, ndump ) == 0 )
    else
      test_if_report = .false.
    endif
  
  end function test_if_report
  
end function if_report_grid_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Test parallel partition
!  - The minimum number of grid cells per node is gc_num( p_lower ) + gc_num( p_upper )
!-----------------------------------------------------------------------------------------
subroutine test_partition_grid_tiles( this, no_co, gc_num )

  implicit none

  class( t_grid_tiles ), intent(inout) :: this
  class( t_node_conf ), intent(in)  :: no_co
  integer, dimension(:,:), intent(in) :: gc_num

  integer :: i, min_nx

  do i = 1, p_x_dim

    ! Minimum number of cells along given direction
    min_nx = max( 1, gc_num( p_lower, i ) + gc_num( p_upper, i ) )

    select type(no_co); class is(t_node_conf_tiles)
      if ( min_nx * no_co%ti_co%g_til(i) > this%g_nx(i) ) then
        write(0,'(A,I1)') '(*error*) Too many partitions along direction ', i
        write(0,       *) '(*error*) nx = ', this%g_nx(i), ' n_tiles = ', no_co%ti_co%g_til(i), &
                          ' mininum grid cells/tile = ', min_nx
        stop
      endif
    end select

  enddo

end subroutine test_partition_grid_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  Returns the grid limits for all nodes
!  The nx_p_nodes needs to be already allocated and has the following structure:
!  nx_p_nodes(bound, dim, pe)
!  where bound is 1:2 (lower/upper boundary), dim is 1:p_x_dim (dimension),
!  and pe is 1:g_til_num(no_co)+2 (global tile number plus two extra slots)
!-----------------------------------------------------------------------------------------
subroutine get_node_limits_tiles( this, no_co, nx_p_nodes  )

  implicit none

  class( t_grid_tiles ), intent(in)     :: this
  class( t_node_conf ), intent(in)        :: no_co
  integer, dimension(:,:,:), intent(out) :: nx_p_nodes

  ! local variables

  integer :: til, i, idx
  integer, dimension(3) :: g_til_aid_mm

  select type(no_co); class is(t_node_conf_tiles)

  do til = 1, product(no_co%ti_co%g_til)
    idx = 0
    do i = 1, this%x_dim
      nx_p_nodes(1,i,til) = this%nx_node( 1, idx + g_tgp( no_co%ti_co, til, i ) )
      nx_p_nodes(2,i,til) = this%nx_node( 2, idx + g_tgp( no_co%ti_co, til, i ) )
      idx = idx + this%nnodes(i)
    enddo
  enddo

  ! Add two extra rows with important tile information for diagnostic
  g_til_aid_mm = no_co%ti_co%g_til_aid_min_max( :, no_co%my_aid() )
  ! Min g_til_aid for current node
  nx_p_nodes(1,this%x_dim,til) = g_til_aid_mm(1)
  ! Max g_til_aid for current node
  nx_p_nodes(2,this%x_dim,til) = g_til_aid_mm(2)
  ! Number of tiles for current node
  nx_p_nodes(1,this%x_dim,til+1) = g_til_aid_mm(3)
  ! g_til_aid for this tile
  nx_p_nodes(2,this%x_dim,til+1) = no_co%ti_co%g_til_aid

  end select

end subroutine get_node_limits_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Sets number of cells per tile for all directions
!-----------------------------------------------------------------------------------------
subroutine parallel_partition_tiles( this, no_co, nx_til )

  implicit none

  class( t_grid_tiles ), intent(inout) :: this
  class( t_node_conf ), intent(in)     :: no_co
  integer, pointer, dimension(:,:)     :: nx_til

  integer :: i, j, idx_node, idx_nx

  idx_nx   = 0
  idx_node = 0

  ! We want each tile to have the same number of cells
  ! so we default to the even distribution

  ! May be allocated already due to restart
  if (associated(this%nx_node)) then
    call freemem(this%nx_node)
  endif

  ! get num cells for each tile in each direction
  this%nx_node => nx_til

  do i = 1, this%x_dim

    ! get maximum number of cells in partition (used for diagnostics output)
    this%max_nx(i) = maxval( this%nx_node( 3, idx_node + 1 : idx_node + this%nnodes(i) ) )

    ! store local tile info
    select type(no_co); class is(t_node_conf_tiles)
      this%my_nx(:,i) = this%nx_node( :, idx_node + no_co%ti_co%g_tgp(i) )
    end select

    idx_nx   = idx_nx   + this%g_nx(i)
    idx_node = idx_node + this%nnodes(i)

  enddo

  ! (* debug *) Validate the partition
  idx_node = 0
  do i = 1, this%x_dim

    ! Check partition sizes
    do j = 1, this%nnodes(i)
      if ( this%nx_node(3,j+idx_node) < this%min_nx(i) ) then
        ERROR('Partition too small')
        call abort_program()
      endif
    enddo

    ! Check overlap
    do j = 2, this%nnodes(i)
      if ( this%nx_node(1,j+idx_node) /= this%nx_node(2,j+idx_node-1) + 1 ) then
        ERROR('Partition overlap')
        call abort_program()
      endif
    enddo

    idx_node = idx_node + this%nnodes(i)
  enddo

end subroutine parallel_partition_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Setup grid object for t_simulation_tiles
!-----------------------------------------------------------------------------------------
subroutine init_grid_group( this, no_co, gc_num, g_box, restart, restart_handle )

  implicit none

  ! dummy variables
  class( t_grid_group ), intent(inout) :: this
  class( t_node_conf ), intent(in)  :: no_co    ! node configuration
  integer, dimension(:, :), intent(in) :: gc_num
  real(p_double), dimension(:,:), intent(in) :: g_box
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  ! TODO: check to see if there is anything extraneous in here
  call this % t_grid % init( no_co, gc_num, g_box, restart, restart_handle )

end subroutine init_grid_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! clear object data
!-----------------------------------------------------------------------------------------
subroutine cleanup_grid_group(this)

  implicit none

  class( t_grid_group ), intent(inout) :: this

  call freemem( this%nx_til )
  call freemem( this%til_load_2 )
  call freemem( this%til_load_3 )

  call this % t_grid % cleanup()

end subroutine cleanup_grid_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Calculate total load from load in each process
!-----------------------------------------------------------------------------------------
subroutine gather_load_grid_group( this, no_co )

  implicit none

  class( t_grid_group ), intent(inout) :: this
  class( t_node_conf ), intent(in)  :: no_co    ! node configuration

  integer :: topology_type

  select type(no_co); class is(t_node_conf_tiles)
    topology_type = no_co%ti_co%topology_type
  end select

  ! Note, lb_gather_max is ignored here (always false for tiles)
  select case(topology_type)
  case(p_topology_viktor)

    select case(p_x_dim)
    case(2)
      call reduce_array( no_co, this%til_load_2, operation = p_sum, all =.true. )
    case(3)
      call reduce_array( no_co, this%til_load_3, operation = p_sum, all =.true. )
    end select

  case default
    call reduce_array( no_co, this%int_load, operation = p_sum, all =.true. )
  end select

end subroutine gather_load_grid_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Sets number of tiles per node
!-----------------------------------------------------------------------------------------
subroutine parallel_partition_group( this, no_co, n )

  implicit none

  class( t_grid_group ), intent(inout) :: this
  class( t_node_conf ), intent(inout)  :: no_co
  integer, intent(in)                  :: n ! iteration

  ! local variables
  integer :: j, min_til

  integer, dimension(:,:), pointer :: g_til_aid_mm

  select type(no_co); class is(t_node_conf_tiles)

    ! Each node needs at least 1 tile, or n_threads tiles if forcing light tiles
    select case( no_co%omp_patt%force_omp )
    case( p_light )
      min_til = n_threads( no_co )
    case default
      min_til = 1
    end select

    select case( no_co%ti_co%topology_type )
    case( p_topology_hilbert, p_topology_snake )

      call alloc( g_til_aid_mm, (/ 3, no_num(no_co) /) )

      if ( any(this%load_balance) .and. (n >= this%start_load_balance )) then
        call load_dist( this%lb_type, g_til_aid_mm, &
          no_num(no_co), product(no_co%ti_co%g_til), &
          this%int_load, min_til )
      else
        call even_dist( g_til_aid_mm, no_num(no_co), product(no_co%ti_co%g_til) )
      endif

      ! (* debug *) Validate the partition
      ! Check partition sizes
      do j = 1, no_num(no_co)
        if ( g_til_aid_mm(3,j) < min_til ) then
          ERROR('Partition too small')
          call abort_program()
        endif
      enddo

      ! Check overlap
      do j = 2, no_num(no_co)
        if ( g_til_aid_mm(1,j) /= g_til_aid_mm(2,j-1) + 1 ) then
          ERROR('Partition overlap')
          call abort_program()
        endif
      enddo

      no_co%ti_co%g_til_aid_min_max =  g_til_aid_mm

      call freemem(g_til_aid_mm)

    case( p_topology_viktor )

      ! For viktor topology, we need to set space filling and inverse space
      ! filling operators here because in this topology they depend on how
      ! tiles are distributed accross nodes (see set_ti_co_viktor)

      if ( this%x_dim == 1 ) then
        ERROR('For 1d simulations, use topology="snake" or topology="hilbert"')
        call abort_program()
      endif

      ! Find partition boundaries in each dim (result is filling cart_til_mm)
      if ( any(this%load_balance) .and. (n >= this%start_load_balance )) then
        call load_dist( this, no_co, min_til, no_co%ti_co%cart_til_mm )
      else
        call even_dist( no_co, no_co%ti_co%cart_til_mm )
      endif

      ! Set other tile configuration data
      call set_ti_co_viktor( no_co, min_til )

    case default

      ERROR('Load balance not implemented for tile topology ',no_co%ti_co%topology_type,',')
      ERROR('Try tile_topology = hilbert, snake or viktor.')
      call abort_program(p_err_notimplemented)

    end select

    ! ! (*debug*)
    ! if (root(no_co)) then
    !   print *, '    --g_til_aid_min_max'
    !   print '    (3I6)', no_co%ti_co%g_til_aid_min_max
    ! endif

  end select

end subroutine parallel_partition_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize member data of group grid needed for diagnostics
! We use an even distribution for dumping diagnostics
!-----------------------------------------------------------------------------------------
subroutine parallel_partition_io(this, no_co)

  implicit none

  class( t_grid_group ), intent(inout) :: this
  class( t_node_conf ), intent(inout) :: no_co

  ! min/max/num elements of inverse space filling operator indexed node
  ! (set up analogous to nx_node/nx_til)
  integer, dimension(:,:), pointer :: nt_node

  integer :: idx_node, idx_til
  integer :: i, j

  select type(no_co); class is(t_node_conf_tiles)

    ! allocate group grid nx_til
    call alloc( this%nx_til, (/3, sum(no_co%ti_co%g_til)/) )

    ! allocate local pointer
    call alloc( nt_node, (/3, sum(this%nnodes(1:p_x_dim))/) )
    
    idx_til  = 0
    idx_node = 0

    this%nx_node = 0
    this%my_nx   = 0

    do i = 1,this%x_dim
      ! call even_dist to figure out how many tiles per node in each dimension
      call even_dist( nt_node( :, idx_node+1:idx_node+this%nnodes(i) ), &
        this%nnodes(i), no_co%ti_co%g_til(i) )

      ! store local node info - num tiles in this dimension
      no_co%ti_co%my_nt_io(:,i) = nt_node( :, idx_node + my_ngp(no_co,i) )

      ! call even_dist to figure out how many cells per tile in each dimension
      call even_dist( this%nx_til( :, idx_til+1:idx_til+no_co%ti_co%g_til(i) ), &
        no_co%ti_co%g_til(i), this%g_nx(i) )

      ! combine last two steps to find out how many cells per node in each dimension
      do j = idx_node+1, idx_node+this%nnodes(i)
        ! lower cell bnd of jth node is lower bnd of jth node's first tile
        this%nx_node(1,j) = this%nx_til( 1, idx_til+nt_node(1,j) )
        ! upper cell bnd of jth node is upper bnd of jth node's last tile
        this%nx_node(2,j) = this%nx_til( 2, idx_til+nt_node(2,j) )
        ! total number of cells is the difference (plus 1 obvs) 
        this%nx_node(3,j) = this%nx_node(2,j) - this%nx_node(1,j) + 1
      enddo

      ! get maximum number of cells on a node in this dimension (used for diagnostics)
      this%max_nx(i) = maxval( this%nx_node( 3, idx_node+1 : idx_node+this%nnodes(i) ) )

      ! store local node info - num cells in this dimension
      this%my_nx(:,i) = this%nx_node( :, idx_node + my_ngp(no_co,i) )

      ! increment idx_til, idx_node
      idx_til  = idx_til  + no_co%ti_co%g_til(i)
      idx_node = idx_node + this%nnodes(i)
    enddo

    call freemem( nt_node )

  end select

end subroutine parallel_partition_io
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Get dlb message pattern: outgoing/incoming tiles and their destinations
!-------------------------------------------------------------------------------
subroutine new_dlb_patt_( patt, no_co, g_til_aid_mm_old, g_til_aid_mm_new, diff )

  implicit none

  type( t_msg_patt ), intent(inout)   :: patt         ! who's patt?
  class( t_node_conf ), intent(in)     :: no_co
  integer, dimension(:,:), intent(in) :: g_til_aid_mm_old
  integer, dimension(:,:), intent(in) :: g_til_aid_mm_new
  integer, dimension(:)               :: diff

  integer, dimension(:), pointer :: msg_nodes
  integer :: my_node, num_nodes
  integer :: lb, ub


  ! Prepare local variables for t_msg_patt calculations
  my_node   = no_co%my_aid()
  num_nodes = no_num(no_co)

  ! worst case scenario we exchange messages with all other nodes
  call alloc( msg_nodes, (/ num_nodes-1 /) )

  ! Get send message pattern
  ! g_til_aid lower/upper bound for tils that were under my control before dlb
  lb = g_til_aid_mm_old(1,my_node)
  ub = g_til_aid_mm_old(2,my_node)
  call get_dlb_msg( patt%send, patt%n_send, diff, msg_nodes, my_node, lb, ub, 1 )

  ! Get receive message pattern
  ! g_til_aid lower/upper bound for tils that will be under my control after dlb
  lb = g_til_aid_mm_new(1,my_node)
  ub = g_til_aid_mm_new(2,my_node)
  call get_dlb_msg( patt%recv, patt%n_recv, diff, msg_nodes, my_node, lb, ub, -1 )

  ! (* debug *)
  ! if ( patt%n_send>0 )
  !   do i = 1, patt%n_send
  !       DEBUG( '--send msg num    :',i )
  !       DEBUG( '--send target node:',patt%send(i)%node )
  !       DEBUG( '--send tiles      :',patt%send(i)%tils(:) )
  !   enddo
  ! endif
  ! if ( patt%n_recv>0 ) then
  !   do i = 1, patt%n_recv
  !       DEBUG( '--recv msg num    :',i )
  !       DEBUG( '--recv target node:',patt%recv(i)%node )
  !       DEBUG( '--recv tiles      :',patt%recv(i)%tils(:) )
  !   enddo
  ! endif

  ! Cleanup
  call freemem( msg_nodes )


end subroutine new_dlb_patt_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Fill in recv's and send's of t_msg_patt
!-------------------------------------------------------------------------------
subroutine get_dlb_msg( msg, n_msg, diff, msg_nodes, my_node, lb, ub, const )

  implicit none

  type( t_msg ), dimension(:), pointer, intent(inout) :: msg
  integer, intent(inout)                              :: n_msg
  integer, dimension(:), intent(in)                   :: diff
  integer, dimension(:), pointer, intent(inout)       :: msg_nodes
  integer, intent(in)                                 :: my_node
  integer, intent(in)                                 :: lb, ub
  integer, intent(in)                                 :: const

  integer :: i, j, k, n

  j = 0
  n_msg = 0
  do i = lb, ub
    if ( diff(i) /= j .and. diff(i) /= 0 ) then
      ! increment number of messages
      n_msg = n_msg + 1
      ! find target or source node, const=+1 for sends, cont=-1 for recvs
      msg_nodes( n_msg ) = my_node + const*diff(i)
      ! store last unchanged value in diff
      j = diff(i)
    endif
  enddo

  if ( n_msg>0 ) then
    call alloc( msg, (/ n_msg /) )

    ! j takes on the last unchanged value of diff
    ! k is placeholder, marking g_til_aid lb of current message
    ! n counts which message we're on
    j = 0
    k = 0
    n = 0
    do i = lb, ub
      ! Handle the end of the array differently
      if (i==ub .and. diff(i)/=0) then
        if (diff(i)==j) then
          ! send only the current message
          n=n+1
          call new_msg( msg(n), msg_nodes(n), i, k )
        else
          if (j==0) then
            ! send only the last tile as a new message
            k = i
            n = n+1
            call new_msg( msg(n), msg_nodes(n), i, k )
          else
            ! send current message
            n = n+1
            call new_msg( msg(n), msg_nodes(n), i-1, k )

            ! then send the last tile
            k = i
            n = n+1
            call new_msg( msg(n), msg_nodes(n), i, k )
          endif
        endif

      elseif ( diff(i)/=j ) then
        ! Don't send message when entering the array or after 0's
        if ( i>lb .and. j/=0 ) then
          n = n+1
          call new_msg( msg(n), msg_nodes(n), i-1, k )
        endif

        k = i
        j = diff(i)
      endif
    enddo
  endif

end subroutine get_dlb_msg
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Fill in member data of nth msg
!-------------------------------------------------------------------------------
subroutine new_msg( msg_n, node, til_ub, til_lb )

  implicit none

  type( t_msg ), intent(inout) :: msg_n
  integer, intent(in)          :: node
  integer, intent(in)          :: til_ub, til_lb

  integer :: i

  ! set tile range
  call alloc( msg_n%tils, (/ til_ub - til_lb + 1 /) )
  do i = til_lb, til_ub
    msg_n%tils(i - til_lb + 1) = i
  enddo

  ! set target node
  msg_n%node = node

end subroutine new_msg
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Get message pattern: outgoing/incoming tiles and their destinations
!-------------------------------------------------------------------------------
subroutine new_dlb_patt_viktor( no_co, patt, diff, difff )

  implicit none

  class( t_node_conf ), intent(in)     :: no_co
  type( t_msg_patt ), intent(inout)   :: patt         ! who's patt?
  integer, dimension(:,:), optional   :: diff
  integer, dimension(:,:,:), optional :: difff

  integer, dimension(p_x_dim) :: lb, ub
  integer :: my_node

  my_node = no_co%my_aid()

  select type(no_co); class is(t_node_conf_tiles)

  ! Get send message pattern
  select case(p_x_dim)
    case(2)
      ! dimension 1 lb and ub
      lb(1) = no_co%ti_co%cart_til_mm_old(1)%ii( 1, my_ngp(no_co,1), 1, 1 )
      ub(1) = no_co%ti_co%cart_til_mm_old(1)%ii( 2, my_ngp(no_co,1), 1, 1 )
      ! dimension 2 lb and ub
      lb(2) = no_co%ti_co%cart_til_mm_old(2)%ii( 1, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
      ub(2) = no_co%ti_co%cart_til_mm_old(2)%ii( 2, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
      ! get msg pattern
      call get_dlb_msg_2d_viktor( no_co, patt%send, patt%n_send, diff, lb, ub, my_node, p_send )
    case(3)
      ! dimension 1 lb and ub
      lb(1) = no_co%ti_co%cart_til_mm_old(1)%ii( 1, my_ngp(no_co,1), 1, 1 )
      ub(1) = no_co%ti_co%cart_til_mm_old(1)%ii( 2, my_ngp(no_co,1), 1, 1 )
      ! dimension 2 lb and ub
      lb(2) = no_co%ti_co%cart_til_mm_old(2)%ii( 1, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
      ub(2) = no_co%ti_co%cart_til_mm_old(2)%ii( 2, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
      ! dimension 3 lb and ub
      lb(3) = no_co%ti_co%cart_til_mm_old(3)%ii( 1, my_ngp(no_co,1), my_ngp(no_co,2), my_ngp(no_co,3) )
      ub(3) = no_co%ti_co%cart_til_mm_old(3)%ii( 2, my_ngp(no_co,1), my_ngp(no_co,2), my_ngp(no_co,3) )
      ! get msg pattern
      call get_dlb_msg_3d_viktor( no_co, patt%send, patt%n_send, difff, lb, ub, my_node, p_send )
  end select

  ! Get recv message pattern
  select case(p_x_dim)
    case(2)
      ! dimension 1 lb and ub
      lb(1) = no_co%ti_co%cart_til_mm(1)%ii( 1, my_ngp(no_co,1), 1, 1 )
      ub(1) = no_co%ti_co%cart_til_mm(1)%ii( 2, my_ngp(no_co,1), 1, 1 )
      ! dimension 2 lb and ub
      lb(2) = no_co%ti_co%cart_til_mm(2)%ii( 1, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
      ub(2) = no_co%ti_co%cart_til_mm(2)%ii( 2, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
      ! get msg pattern
      call get_dlb_msg_2d_viktor( no_co, patt%recv, patt%n_recv, diff, lb, ub, my_node, p_recv )
    case(3)
      ! dimension 1 lb and ub
      lb(1) = no_co%ti_co%cart_til_mm(1)%ii( 1, my_ngp(no_co,1), 1, 1 )
      ub(1) = no_co%ti_co%cart_til_mm(1)%ii( 2, my_ngp(no_co,1), 1, 1 )
      ! dimension 2 lb and ub
      lb(2) = no_co%ti_co%cart_til_mm(2)%ii( 1, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
      ub(2) = no_co%ti_co%cart_til_mm(2)%ii( 2, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
      ! dimension 3 lb and ub
      lb(3) = no_co%ti_co%cart_til_mm(3)%ii( 1, my_ngp(no_co,1), my_ngp(no_co,2), my_ngp(no_co,3) )
      ub(3) = no_co%ti_co%cart_til_mm(3)%ii( 2, my_ngp(no_co,1), my_ngp(no_co,2), my_ngp(no_co,3) )
      ! get msg pattern
      call get_dlb_msg_3d_viktor( no_co, patt%recv, patt%n_recv, difff, lb, ub, my_node, p_recv )
  end select

  end select

  ! (* debug *)
  ! if ( patt%n_send>0 ) then
  !   do i = 1, patt%n_send
  !       DEBUG( '--send msg num    :',i )
  !       DEBUG( '--send target node:',patt%send(i)%node )
  !       DEBUG( '--send tiles      :',patt%send(i)%tils(:) )
  !   enddo
  ! endif
  ! if ( patt%n_recv>0 ) then
  !   do i = 1, patt%n_recv
  !       DEBUG( '--recv msg num    :',i )
  !       DEBUG( '--recv target node:',patt%recv(i)%node )
  !       DEBUG( '--recv tiles      :',patt%recv(i)%tils(:) )
  !   enddo
  ! endif

end subroutine new_dlb_patt_viktor
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Fill in recv's and send's of t_msg_patt
!-----------------------------------------------------------------------------------------
subroutine get_dlb_msg_2d_viktor( no_co, msg, n_msg, diff, lb, ub, my_node, const )

  implicit none

  class( t_node_conf ), intent(in)                     :: no_co
  type( t_msg ), dimension(:), pointer, intent(inout) :: msg
  integer, intent(inout)                              :: n_msg
  integer, intent(in)                                 :: my_node
  integer, dimension(:), intent(in)                   :: lb, ub
  integer, dimension(:,:), intent(in)                 :: diff
  integer, intent(in)                                 :: const

  integer, dimension(:,:), pointer :: unique
  integer :: val, count, n, t
  integer :: i, j

  if ( all(diff(lb(1):ub(1),lb(2):ub(2))==0) ) then
    n_msg = 0
    return
  endif

  unique => null()
  call unique_count_2d( diff(lb(1):ub(1),lb(2):ub(2)), unique )

  n_msg = size(unique,1)
  call alloc( msg, (/ n_msg /) )

  select type(no_co); class is(t_node_conf_tiles)

  do n = 1, n_msg
    ! value to look for in diff
    val = unique(n,1)
    ! pack node id into the msg
    msg(n)%node = my_node + const*val
    ! number of times that value will appear (number of tiles going to that node)
    call alloc( msg(n)%tils, (/unique(n,2)/) )
    ! find g_til_aid's and pack them into the msg
    count = 0
    do i = lb(1), ub(1)
      do j = lb(2), ub(2)
        if ( diff(i,j)==val ) then
          count = count+1
          t = no_co%ti_co%hi2( i, j)
          msg(n)%tils(count) = t
        endif
      enddo
    enddo
  enddo

  end select

  call freemem(unique)

end subroutine get_dlb_msg_2d_viktor
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Fill in recv's and send's of t_msg_patt
!-----------------------------------------------------------------------------------------
subroutine get_dlb_msg_3d_viktor( no_co, msg, n_msg, diff, lb, ub, my_node, const )

  implicit none

  class( t_node_conf ), intent(in)                     :: no_co
  type( t_msg ), dimension(:), pointer, intent(inout) :: msg
  integer, intent(inout)                              :: n_msg
  integer, intent(in)                                 :: my_node
  integer, dimension(:), intent(in)                   :: lb, ub
  integer, dimension(:,:,:), intent(in)               :: diff
  integer, intent(in)                                 :: const

  integer, dimension(:,:), pointer :: unique
  integer :: val, count, n, t
  integer :: i, j, k

  if ( all( diff(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) )==0) ) then
    n_msg = 0
    return
  endif

  unique => null()
  call unique_count_3d( diff(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), unique )

  n_msg = size(unique,1)
  call alloc( msg, (/ n_msg /) )

  select type(no_co); class is(t_node_conf_tiles)

  do n = 1, n_msg
    ! value to look for in diff
    val = unique(n,1)
    ! pack node id into the msg
    msg(n)%node = my_node + const*val
    ! number of times that value will appear (number of tiles going to that node)
    call alloc( msg(n)%tils, (/unique(n,2)/) )
    ! find g_til_aid's and pack them into the msg
    count = 0
    do i = lb(1), ub(1)
      do j = lb(2), ub(2)
        do k = lb(3), ub(3)
          if ( count==size(msg(n)%tils) ) exit
          if ( diff(i,j,k)==val ) then
            count = count+1
            t = no_co%ti_co%hi3( i, j, k)
            msg(n)%tils(count) = t
          endif
        enddo
      enddo
    enddo
  enddo

  end select

  call freemem(unique)

end subroutine get_dlb_msg_3d_viktor
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize tile load array
!-----------------------------------------------------------------------------------------
subroutine init_tile_load( this, no_co )

  implicit none

  class( t_grid ), intent(inout)   :: this
  class( t_node_conf ), intent(in) :: no_co

  integer :: size
  integer :: i

  select type(this); class is(t_grid_group)
  select type(no_co); class is(t_node_conf_tiles)

  select case( no_co%ti_co%topology_type )
  case( p_topology_viktor )

    select case( p_x_dim )
    case(2)
      call freemem( this%til_load_2 )
      call alloc( this%til_load_2, no_co%ti_co%g_til )
      this%til_load_2 = 0
    case(3)
      call freemem( this%til_load_3 )
      call alloc( this%til_load_3, no_co%ti_co%g_til )
      this%til_load_3 = 0
    end select

  case default

    call freemem( this%int_load )

    size = product(no_co%ti_co%g_til)
    call alloc( this%int_load, (/ size /))
    do i = 1, size
      this%int_load(i) = 0
    enddo

  end select

  end select
  end select

end subroutine init_tile_load
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Adds the tile calculation load (for local cells) to the load array and/or fixed expression
! This routine must be called before gather_load
!-----------------------------------------------------------------------------------------
subroutine add_grid_load_tiles( this, no_co, n, int_load, til_load_2, til_load_3 )

  implicit none

  class( t_grid ), intent(inout) :: this
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: n
  real( p_double ), dimension(:), pointer, optional :: int_load
  ! used for viktor's alg:
  real(p_double), dimension(:,:), pointer, optional   :: til_load_2
  real(p_double), dimension(:,:,:), pointer, optional :: til_load_3

  integer :: dim, i, j, t_lower, t_upper, til_num, nnodes, my_aid, l_nx_til, mod_res_til
  real(p_double) :: l_vol
  integer, dimension(p_x_dim) :: tgp

  ! Global tile number and cell number in each direction
  integer, dimension( p_x_dim ) :: gtil, g_nx

  ! Grids/node (l_nx) and remainder of the same (mod_res), and l_nx1 = l_nx + 1
  ! Patterned after even_dist in os-grid-parallel
  integer, dimension(p_max_dim) :: l_nx, l_nx1, mod_res

  ! Add cell calculation load
  if ( ( ( this % lb_type == p_static_load_balance ) .or. &
         ( this % lb_type == p_dynamic_load_balance ) ) .and. &
       ( this%cell_weight > 0 ) ) then

    select type(no_co); class is(t_node_conf_tiles)

    if ( n > 0 ) then

      ! the node partition has already been defined so use the local grid cells only

      l_vol = this%my_nx(3, 1)
      do dim = 2, this%x_dim
        l_vol = l_vol * this%my_nx(3, dim)
      enddo

      select case( no_co%ti_co%topology_type )
      case(p_topology_viktor)

        select case(p_x_dim)
        case(2)
          ! til_load_2 will be present
          til_load_2(no_co%ti_co%g_tgp(1),no_co%ti_co%g_tgp(2)) = &
            til_load_2(no_co%ti_co%g_tgp(1),no_co%ti_co%g_tgp(2)) + this%cell_weight * l_vol
        case(3)
          ! til_load_3 will be present
          til_load_3(no_co%ti_co%g_tgp(1),no_co%ti_co%g_tgp(2),no_co%ti_co%g_tgp(3)) = &
            til_load_3(no_co%ti_co%g_tgp(1),no_co%ti_co%g_tgp(2),no_co%ti_co%g_tgp(3)) + &
            this%cell_weight * l_vol
        end select

      case default

        ! int_load will be present for n>0
        int_load(no_co%ti_co%g_til_aid) = int_load(no_co%ti_co%g_til_aid) + this%cell_weight * l_vol

      end select

    else

      ! No node partition has been defined yet so use global volume

      ! get the sizes of each tile
      gtil = no_co%ti_co%g_til(1:p_x_dim)
      do i = 1, this%x_dim
        g_nx(i) = this%g_nx(i)
        l_nx(i) = g_nx(i)/gtil(i)
        mod_res(i) = mod(g_nx(i),gtil(i))
      enddo
      l_nx1 = l_nx + 1

      ! Find which tiles this node is responsible for. This is essentially even_dist_
      til_num = product(gtil); nnodes = no_num(no_co)
      l_nx_til = til_num / nnodes
      mod_res_til = mod(til_num, nnodes)
      my_aid = no_co % my_aid()

      if (my_aid > mod_res_til) then
        t_lower = 1 + (my_aid - 1) * l_nx_til + mod_res_til
        t_upper = t_lower + l_nx_til - 1
      else
        t_lower = 1 + (my_aid - 1) * (l_nx_til + 1)
        t_upper = t_lower + l_nx_til
      endif

      ! Loop over tiles and calculate grid load for each
      do i = t_lower, t_upper
        l_vol = 1
        do j = 1, this%x_dim
          tgp(j) = g_tgp( no_co%ti_co, i, j )
          ! Multiply local volume based on tile size
          if ( tgp(j) > mod_res(j) ) then
            l_vol = l_vol * l_nx(j)
          else
            l_vol = l_vol * l_nx1(j)
          endif
        enddo

        select case( no_co%ti_co%topology_type )
        case(p_topology_viktor)

          select type(this); class is(t_grid_group)

          select case(p_x_dim)
          case(2)
            ! til_load_2 will not be present
            this%til_load_2(tgp(1),tgp(2)) = &
              this%til_load_2(tgp(1),tgp(2)) + this%cell_weight * l_vol
          case(3)
            ! til_load_3 will not be present
            this%til_load_3(tgp(1),tgp(2),tgp(3)) = &
              this%til_load_3(tgp(1),tgp(2),tgp(3)) + this%cell_weight * l_vol
          end select

          end select

        case default

          ! int_load will not be present for n==0
          this%int_load(i) = this%int_load(i) + this%cell_weight * l_vol

        end select

      enddo

    endif

    end select

  endif

  ! Add expression load
  if ( this % lb_type == p_expr_load_balance ) then
    ! (*debug*) This can only be done at initialization (n=0)
    if ( n > 0 ) then

    endif

    ERROR('Expression load not yet implemented for tiles')
    call abort_program(p_err_notimplemented)

    call add_expression_load( this )
  endif


end subroutine add_grid_load_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Converts a load array into a cumulative load array for tiles
!-----------------------------------------------------------------------------------------
subroutine conv_load_2_int_load_tiles( int_load )

  implicit none

  real(p_double), pointer, dimension(:), intent(inout) :: int_load

  integer :: i

  ! if ( mpi_node() == 0 ) then
  ! print *, 'load'
  ! print *, int_load
  ! endif

  ! sanity check - checks for bad initial int_load
  do i = 1, size(int_load)
    if ( int_load(i) < 0 ) then
      SCR_MPINODE('Invalid value for int_load (A), tile = ', i)
      call abort_program( p_err_invalid )
    endif
  enddo

  do i = 2, size(int_load)
    int_load(i) = int_load(i) + int_load(i-1)
  enddo

  ! sanity check - checks for overflows
  do i = 1, size(int_load)
    if ( int_load(i) < 0 ) then
      SCR_MPINODE('Invalid value for int_load (B), tile = ', i)
      call abort_program( p_err_invalid )
    endif
  enddo

  ! if ( mpi_node() == 0 ) then
  ! print *, 'int_load'
  ! print *, int_load
  ! endif

end subroutine conv_load_2_int_load_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Converts a load array into a cumulative load array for tiles for local segment
!-----------------------------------------------------------------------------------------
subroutine conv_load_2_int_load_tiles_1( int_load, lb, ub )

  implicit none

  real(p_double), pointer, dimension(:), intent(inout) :: int_load
  integer, intent(in) :: lb, ub

  integer :: i

  ! if ( mpi_node() == 0 ) then
  ! print *, 'load'
  ! print *, int_load
  ! endif

  ! sanity check - checks for bad initial int_load
  do i = lb, ub
    if ( int_load(i) < 0 ) then
      SCR_MPINODE('Invalid value for int_load (A), tile = ', i)
      call abort_program( p_err_invalid )
    endif
  enddo

  do i = lb+1, ub
    int_load(i) = int_load(i) + int_load(i-1)
  enddo

  ! sanity check - checks for overflows
  do i = lb, ub
    if ( int_load(i) < 0 ) then
      SCR_MPINODE('Invalid value for int_load (B), tile = ', i)
      call abort_program( p_err_invalid )
    endif
  enddo

  ! if ( mpi_node() == 0 ) then
  ! print *, 'int_load'
  ! print *, int_load
  ! endif

end subroutine conv_load_2_int_load_tiles_1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Converts a load array into a cumulative load array for tiles after mpi reduce
!-----------------------------------------------------------------------------------------
subroutine conv_load_2_int_load_tiles_2( int_load, g_til_aid_mm )

  implicit none

  real(p_double), pointer, dimension(:), intent(inout) :: int_load
  integer, dimension(:,:), intent(in) :: g_til_aid_mm

  integer :: i, tot, lb, ub

  ! if ( mpi_node() == 0 ) then
  ! print *, 'load'
  ! print *, int_load
  ! endif

  ! sanity check - checks for bad initial int_load
  ! shouldn't need to do this since each tile checked
  ! do i = 1, size(int_load)
  !   if ( int_load(i) < 0 ) then
  !     SCR_MPINODE('Invalid value for int_load (A), tile = ', i)
  !     call abort_program( p_err_invalid )
  !   endif
  ! enddo

  ! int_load has been integrated by each process,
  ! so we just need to add the last point of each process to the next values
  do i = 2, size(g_til_aid_mm,2)
    tot = int( int_load( g_til_aid_mm(2,i-1) ) )
    lb = g_til_aid_mm(1,i)
    ub = g_til_aid_mm(2,i)
    int_load(lb:ub) = int_load(lb:ub) + tot
  enddo

  ! sanity check - checks for overflows
  do i = 1, size(int_load)
    if ( int_load(i) < 0 ) then
      SCR_MPINODE('Invalid value for int_load (B), tile = ', i)
      call abort_program( p_err_invalid )
    endif
  enddo

  ! if ( mpi_node() == 0 ) then
  ! print *, 'int_load'
  ! print *, int_load
  ! endif

end subroutine conv_load_2_int_load_tiles_2
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! projects a 2d load matrix onto dimension specified by dim
!-------------------------------------------------------------------------------
subroutine project_2( load, int_load, dim )

  implicit none

  real(p_double), dimension(:,:), intent(in) :: load
  real(p_double), dimension(:), pointer      :: int_load
  integer, intent(in)                        :: dim

  integer :: i

  do i = 1, size(load,dim)
    if (dim==1) int_load(i) = sum(load(i,:))
    if (dim==2) int_load(i) = sum(load(:,i))
  enddo


end subroutine project_2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! projects a 3d load matrix onto dimension specified by dim
!-------------------------------------------------------------------------------
subroutine project_3( load, int_load, dim )

  implicit none

  real(p_double), dimension(:,:,:), intent(in) :: load
  real(p_double), dimension(:), pointer        :: int_load
  integer, intent(in)                          :: dim

  integer :: i

  do i = 1, size(load,dim)
    if (dim==1) int_load(i) = sum(load(i,:,:))
    if (dim==2) int_load(i) = sum(load(:,i,:))
    if (dim==3) int_load(i) = sum(load(:,:,i))
  enddo


end subroutine project_3
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Sets the number of cells per node as specified by lb_type
!-----------------------------------------------------------------------------------------
subroutine load_dist_viktor( this, no_co, min_til, cart_til_mm )

  implicit none

  class( t_grid ), intent(in)               :: this
  class( t_node_conf ), intent(in)          :: no_co
  integer, intent(in)                      :: min_til
  type(int_array_4d), dimension(:), pointer :: cart_til_mm

  real(p_double), dimension(:), pointer :: int_load

  integer :: i, k

  select type(this); class is(t_grid_group)
  select type(no_co); class is(t_node_conf_tiles)

  select case( this%x_dim )
  case(2)
    ! Fill in cart_til_mm(1)%ii
    call alloc( int_load, (/no_co%ti_co%g_til(1)/) )
    ! project to first dimension
    call project( this%til_load_2, int_load, 1 )
    ! conv load array to integrated load array
    call conv_load_2_int_load( int_load )
    ! find distribution for 1st dimension (fill cart_til_mm(1)%ii
    call load_dist( this%lb_type, cart_til_mm(1)%ii(:,:,1,1), nx(no_co,1), &
      no_co%ti_co%g_til(1), int_load, min_til )
    ! cleanup objs
    call freemem( int_load )

    ! Fill in cart_til_mm(2)%ii
    do i = 1, nx(no_co,1)
      call alloc( int_load, (/no_co%ti_co%g_til(2)/) )
      ! project node's range of tiles in dim one to second dim
      call project( this%til_load_2(cart_til_mm(1)%ii(1,i,1,1): &
                                    cart_til_mm(1)%ii(2,i,1,1),:), &
                                    int_load, 2 )
      ! conv load array to integrated load array
      call conv_load_2_int_load( int_load )
      ! find distribution for 2nd dimension (fill cart_til_mm(2)%ii)
      call load_dist( this%lb_type, cart_til_mm(2)%ii(:,i,:,1), nx(no_co,2), &
        no_co%ti_co%g_til(2), int_load, min_til )
      ! cleanup objs
      call freemem( int_load )
    enddo
  case(3)
    ! Fill in cart_til_mm(1)%ii
    call alloc( int_load, (/no_co%ti_co%g_til(1)/) )
    ! project to first dimension
    call project( this%til_load_3, int_load, 1 )
    ! conv load array to integrated load array
    call conv_load_2_int_load( int_load )
    ! find distribution for 1st dimension (fill cart_til_mm(1)%ii)
    call load_dist( this%lb_type, cart_til_mm(1)%ii(:,:,1,1), nx(no_co,1), &
      no_co%ti_co%g_til(1), int_load, min_til )
    ! cleanup objs
    call freemem( int_load )

    ! Fill in cart_til_mm(2)%ii
    do i = 1, nx(no_co,1)
      call alloc( int_load, (/no_co%ti_co%g_til(2)/) )
      ! project node's range of tiles in dim one to second dim
      call project( this%til_load_3(cart_til_mm(1)%ii(1,i,1,1):&
                                    cart_til_mm(1)%ii(2,i,1,1),:,:), &
                                    int_load, 2 )
      ! conv load array to integrated load array
      call conv_load_2_int_load( int_load )
      ! find distribution for 1st dimension (fill cart_til_mm(2)%ii)
      call load_dist( this%lb_type, cart_til_mm(2)%ii(:,i,:,1), nx(no_co,2), &
        no_co%ti_co%g_til(2), int_load, min_til )
      ! cleanup objs
      call freemem( int_load )

      ! Fill in cart_til_mm(3)%ii
      do k = 1, nx(no_co,2)
        call alloc( int_load, (/no_co%ti_co%g_til(3)/) )
        ! project node's range of tiles in dim two to third dim
        call project( this%til_load_3(cart_til_mm(1)%ii(1,i,1,1):&
                      cart_til_mm(1)%ii(2,i,1,1),&
                      cart_til_mm(2)%ii(1,i,k,1):cart_til_mm(2)%ii(2,i,k,1),:),&
                      int_load, 3 )
        ! conv load array to integrated load array
        call conv_load_2_int_load( int_load )
        ! find distribution for 1st dimension (fill cart_til_mm(:,:,3))
        call load_dist( this%lb_type, cart_til_mm(3)%ii(:,i,k,:), nx(no_co,3), &
          no_co%ti_co%g_til(3), int_load, min_til )
        ! cleanup objs
        call freemem( int_load )
      enddo
    enddo
  end select

  end select
  end select

end subroutine load_dist_viktor
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Distributes the tiles as evenly as possible between all nodes, with viktor topology
!-----------------------------------------------------------------------------------------
subroutine even_dist_viktor( no_co, cart_til_mm )

  implicit none

  class( t_node_conf ), intent(in)          :: no_co
  type(int_array_4d), dimension(:), pointer :: cart_til_mm

  integer :: i, j

  select type(no_co); class is(t_node_conf_tiles)

  select case(p_x_dim)
  case (1)
    ! find even dist for dim 1
    call even_dist( cart_til_mm(1)%ii(:,:,1,1), nx(no_co,1), no_co%ti_co%g_til(1) )
  case (2)
    ! find even dist for dim 1
    call even_dist( cart_til_mm(1)%ii(:,:,1,1), nx(no_co,1), no_co%ti_co%g_til(1) )
    ! find even dist for dim 2
    do i = 1, nx(no_co,1)
      call even_dist( cart_til_mm(2)%ii(:,i,:,1), nx(no_co,2), no_co%ti_co%g_til(2) )
    enddo
  case (3)
    ! find even dist for dim 1
    call even_dist( cart_til_mm(1)%ii(:,:,1,1), nx(no_co,1), no_co%ti_co%g_til(1) )
    ! find even dist for dim 2
    do i = 1, nx(no_co,1)
      call even_dist( cart_til_mm(2)%ii(:,i,:,1), nx(no_co,2), no_co%ti_co%g_til(2) )
    enddo
    ! find even dist for dim 3
    do i = 1, nx(no_co,1)
      do j = 1, nx(no_co,2)
        call even_dist( cart_til_mm(3)%ii(:,i,j,:), nx(no_co,3), no_co%ti_co%g_til(3) )
      enddo
    enddo
  end select

  end select

end subroutine even_dist_viktor
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! input: unsorted integer array
! output: sorted integer array of unique values in first row (excluding 0),
! and count of each value in second row
!-----------------------------------------------------------------------------------------
subroutine unique_count_1d( list, unique )

  implicit none

  integer, dimension(:), intent(in) :: list
  integer, dimension(:,:), pointer :: unique

  integer :: min_val, max_val, i
  integer, dimension(:,:), pointer :: tmp_uniq

  call alloc( tmp_uniq, (/size(list),2/) )

  min_val = minval(list)-1
  max_val = maxval(list)

  i = 0
  do while (min_val<max_val)
    min_val = minval(list, mask=list>min_val)
    if (min_val /= 0) then
      i = i+1
      tmp_uniq(i,1) = min_val
      tmp_uniq(i,2) = count(list==min_val)
    endif
  enddo

  if ( associated(unique) ) then
    ERROR('pointer "unique" should not be initialized')
    call abort_program(p_err_invalid)
  endif

  if (i>0) then
    call alloc( unique, (/i,2/) )
    unique = tmp_uniq(1:i,:)
  endif

  call freemem(tmp_uniq)

end subroutine unique_count_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! input: unsorted integer array
! output: sorted integer array of unique values in first row (excluding 0),
! and count of each value in second row
!-----------------------------------------------------------------------------------------
subroutine unique_count_2d( list, unique )

  implicit none

  integer, dimension(:,:), intent(in) :: list
  integer, dimension(:,:), pointer :: unique

  integer :: min_val, max_val, i
  integer, dimension(:,:), pointer :: tmp_uniq

  call alloc( tmp_uniq, (/size(list),2/) )

  min_val = minval(list)-1
  max_val = maxval(list)

  i = 0
  do while (min_val<max_val)
    min_val = minval(list, mask=list>min_val)
    if (min_val /= 0) then
      i = i+1
      tmp_uniq(i,1) = min_val
      tmp_uniq(i,2) = count(list==min_val)
    endif
  enddo

  if ( associated(unique) ) then
    ERROR('pointer "unique" should not be initialized')
    call abort_program(p_err_invalid)
  endif

  if (i>0) then
    call alloc( unique, (/i,2/) )
    unique = tmp_uniq(1:i,:)
  endif

  call freemem(tmp_uniq)

end subroutine unique_count_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! input: unsorted integer array
! output: sorted integer array of unique values in first row (excluding 0),
! and count of each value in second row
!-----------------------------------------------------------------------------------------
subroutine unique_count_3d( list, unique )

  implicit none

  integer, dimension(:,:,:), intent(in) :: list
  integer, dimension(:,:), pointer :: unique

  integer :: min_val, max_val, i
  integer, dimension(:,:), pointer :: tmp_uniq

  call alloc( tmp_uniq, (/size(list),2/) )

  min_val = minval(list)-1
  max_val = maxval(list)

  i = 0
  do while (min_val<max_val)
    min_val = minval(list, mask=list>min_val)
    if (min_val /= 0) then
      i = i+1
      tmp_uniq(i,1) = min_val
      tmp_uniq(i,2) = count(list==min_val)
    endif
  enddo

  if ( associated(unique) ) then
    ERROR('pointer "unique" should not be initialized')
    call abort_program(p_err_invalid)
  endif

  if ( i>0 ) then
    call alloc( unique, (/i,2/) )
    unique = tmp_uniq(1:i,:)
  endif

  call freemem(tmp_uniq)

end subroutine unique_count_3d
!-----------------------------------------------------------------------------------------

end module m_grid_tiles
