# 1 "overdense/os-spec-splitting.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "overdense/os-spec-splitting.f03" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Particle splitting and recombination module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------------------------
! To-do/ to-fix (or at least not done/not fixed)
!-----------------------------------------------------------------------------------------
!
! - Recombination does not work with moving window (deeply so)
! - Linear momentum is conserved in recombinations, but not angular momentum, I think.
! Although this is clearly unphysical, I don't see how it matters.
! Let's say a neutral was involved somehow, thus making it physical again.
! - Store and output energy non-conservation in splitting and recombinations
!
!-----------------------------------------------------------------------------------------

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
# 18 "overdense/os-spec-splitting.f03" 2
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
# 19 "overdense/os-spec-splitting.f03" 2

module m_spec_splitting

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
# 23 "overdense/os-spec-splitting.f03" 2

use m_system
use m_parameters
use m_space, only : t_space, xmin, xmax, if_move_any
use m_grid_define, only : t_grid
use m_input_file, only : t_input_file, get_namelist

implicit none

private

! Maximum number of recombination regions (one on each cube face)
integer, parameter :: p_max_recombine = 6

! Maximum number of gammas for splitting
! necessary because of the limitations in the use of namelists
! namelist can not include variable-size data structures
integer, parameter :: p_num_split_max = 100

!-----------------------------------------------------------------------------------------
! recombination class
!-----------------------------------------------------------------------------------------
type :: t_recombine

  ! recombination related variables
  !-------------------
  ! region (in simulation units) over which recombination is happening
  real(p_k_part), dimension(p_x_dim,2) :: recombination_box
  ! recombination "cell" size (in simulation units), usually set to dx
  ! determines grid for binning particles to look to recombine
  real(p_k_part), dimension(p_x_dim) :: recombination_cell

  ! calculated variables
  real(p_k_part), dimension(p_x_dim,2) :: lr_box ! LocalRecombination_box (but with globally referenced positions)
  logical :: if_recombine_local = .true.
  integer, dimension(p_x_dim) :: l_grid_size ! number of recombination cells in each direction
  integer, dimension(:), pointer :: cells => null() ! Indecies of particles to recombine, one per cell
  integer :: cell_count

  type(t_recombine), pointer :: next => null()

end type t_recombine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! splitting class
!-----------------------------------------------------------------------------------------
type :: t_split

  ! splitting related variables
  !-------------------
  ! array of energies where splitting occurs and how many times it happens
  ! e.g., splitting_gammas(1:2) = 1.5,2.0,
  ! not split below 1.5, split once between 1.5-2.0, split twice above 2.0
  real(p_k_part), dimension(:), pointer :: splitting_gammas => null()
  ! flag for giving a push to each split particle to separate them slightly
  logical :: if_push
  ! fraction of transverse momentum to add to the particles. Should be small
  real(p_k_part) :: boost
  ! flag for multiplying this small boost by a random number between 0.5 and 1.5
  logical :: random_boost
  ! run splitting every n_split timesteps
  integer :: n_split
  ! flag to only do splitting past a certain x1 value
  logical :: use_min_x
  ! minimum x1 value for performing splitting
  real(p_k_part) :: min_x
  ! split in front of (.true.) or behind (.false.) min_x
  logical :: forward_of_min_x
  ! direction to check min_x
  integer :: dir_min_x

  ! calculated variables
  logical :: if_split = .false.
  logical :: if_split_local = .false.
  real(p_k_part) :: base_particle_q

  ! recombination related variables
  !-------------------
  ! flag to recombine particles
  logical :: if_recombine
  ! run recombination every n_recombine timesteps
  integer :: n_recombine
  real(p_double) :: recombine_start_time
  type(t_recombine), pointer :: recombine => null()

end type t_split
!-----------------------------------------------------------------------------------------

interface init
  module procedure init_split
end interface init

interface read_nml
  module procedure read_nml_splitting
end interface read_nml

interface local_boundaries
  module procedure local_boundaries_split
end interface local_boundaries

interface cleanup
  module procedure cleanup_split
end interface cleanup

public :: init, read_nml, local_boundaries, cleanup, if_split
public :: if_recombine_local, floor_log_2

public :: t_recombine, t_split

contains

!-----------------------------------------------------------------------------------------
subroutine read_nml_splitting( this, input_file, if_collide, rqm, ppc, density )
!-----------------------------------------------------------------------------------------

  type( t_split ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, intent(in) :: if_collide
  real(p_k_part), intent(in) :: rqm
  integer, intent(in) :: ppc
  real(p_k_part), intent(in) :: density

  integer :: num_split
  ! the charge of an unsplit particle
  real(p_k_part) :: unsplit_charge_density
  ! to get unsplit charge exactly right, specify unsplit_charge_density (numerator) as the
  ! density and unsplit_charge_ppc as the particles per cell
  ! divides unsplit_charge_density
  integer :: unsplit_charge_ppc
  real(p_k_part), dimension(p_num_split_max) :: splitting_gammas
  real(p_k_part) :: boost
  logical :: random_boost, if_push
  integer :: n_split
  logical :: use_min_x
  real(p_k_part) :: min_x
  logical :: forward_of_min_x
  integer :: dir_min_x

  real(p_k_part), dimension(p_x_dim,2,p_max_recombine) :: recomb_box_bounds
  real(p_k_part), dimension(p_x_dim,p_max_recombine) :: recomb_cell_dx
  logical :: if_recombine
  integer :: n_recombine
  integer :: n_recombine_regions
  real(p_double) :: recombine_start_time
  type(t_recombine), pointer :: recombine, tail

  namelist /nl_splitting/ unsplit_charge_density, unsplit_charge_ppc, splitting_gammas, &
                          boost, random_boost, n_split, if_push, recomb_box_bounds, &
                          recomb_cell_dx, if_recombine, n_recombine, use_min_x, &
                          min_x, forward_of_min_x, dir_min_x, recombine_start_time, &
                          n_recombine_regions

  integer :: i, j, ierr

  unsplit_charge_density = density
  unsplit_charge_ppc = ppc
  splitting_gammas = -1.0
  boost = 0.01
  random_boost = .true.
  n_split = 0
  if_push = .true.

  use_min_x = .false.
  min_x = - huge(min_x)
  forward_of_min_x = .true.
  dir_min_x = 1

  if_recombine = .false.
  recomb_box_bounds = 0.0
  recomb_cell_dx = 0.0
  n_recombine = 25
  n_recombine_regions = 1
  recombine_start_time = -1.0_p_double

  call get_namelist( input_file, "nl_splitting", ierr )

  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_splitting, iostat = ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "   Error reading splitting parameters"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
    this%if_split = .false.
    return
  endif

  num_split = 0
  this%if_split = .false.
  if (splitting_gammas(1) .ge. 1.0 .and. n_split .gt. 0) then

    num_split = 1
    this%if_split = .true.

    do_loop: do i = 2, p_num_split_max
      if (splitting_gammas(i) .lt. 1.0) exit do_loop

      if (splitting_gammas(i-1) .le. splitting_gammas(i)) then
        num_split = num_split + 1
      else
        if (mpi_node()==0) print *,'Warning: Splitting energies must be non-decreasing.'
        if (mpi_node()==0) print *,"I'm only using the first ", num_split, "values. The rest I'm ignoring."
        exit do_loop
      endif

    enddo do_loop

    call alloc(this%splitting_gammas, (/num_split/),"overdense/os-spec-splitting.f03",235)
    do i = 1, num_split
      this%splitting_gammas( i ) = splitting_gammas( i )
    enddo

    this%base_particle_q = sign(unsplit_charge_density, rqm)
    this%base_particle_q = this%base_particle_q / unsplit_charge_ppc
    this%boost = boost
    this%random_boost = random_boost
    this%n_split = n_split
    this%if_push = if_push

    this%use_min_x = use_min_x
    this%min_x = min_x
    this%forward_of_min_x = forward_of_min_x
    this%dir_min_x = dir_min_x

    if ( this%n_split .lt. 1 ) then
      write(err_buf__,*) "n_split must be a positive integer";call err__("overdense/os-spec-splitting.f03",253)
      call abort_program( p_err_invalid )
    endif

    if ( unsplit_charge_density .eq. 0. ) then
      write(err_buf__,*) "unsplit_charge must be non-zero";call err__("overdense/os-spec-splitting.f03",258)
      call abort_program( p_err_invalid )
    endif

    if ( unsplit_charge_ppc .le. 0 ) then
      write(err_buf__,*) "unsplit_charge_ppc must be a positive integer";call err__("overdense/os-spec-splitting.f03",263)
      call abort_program( p_err_invalid )
    endif

    if ( this%if_push .and. ( this%boost .eq. 0.0) ) then
      write(err_buf__,*) "If using pushing splitting, boost must be non-zero";call err__("overdense/os-spec-splitting.f03",268)
      write(err_buf__,*) "(There's no good reason for it to be less than zero, but it shouldn't matter)";call err__("overdense/os-spec-splitting.f03",269)
      call abort_program( p_err_invalid )
    endif

    if ( this%if_push .and. if_collide ) then
      write(err_buf__,*) "Not necessary to push split particles apart when collisions are on";call err__("overdense/os-spec-splitting.f03",274)
      write(err_buf__,*) "Please turn off one or the other and run again";call err__("overdense/os-spec-splitting.f03",275)
      call abort_program( p_err_invalid )
    endif

    if ( (.not. this%if_push) .and. (.not. if_collide) ) then
      write(err_buf__,*) "Particle splitting is on, but there is nothing I know of to divert the particle's phase trajectory";call err__("overdense/os-spec-splitting.f03",280)
      write(err_buf__,*) "If there's something I haven't thought of to do that you'll have to edit out this message and recompile";call err__("overdense/os-spec-splitting.f03",281)
      call abort_program( p_err_invalid )
    endif

    if ( this%if_split .and. (p_p_dim .ne. 3) ) then
      write(err_buf__,*) "Particle splitting is hardwired for three momentum dimensions";call err__("overdense/os-spec-splitting.f03",286)
      write(err_buf__,*) "There's no good reason for that, I apologize";call err__("overdense/os-spec-splitting.f03",287)
      call abort_program( p_err_invalid )
    endif

    if ( this%use_min_x .and. (dir_min_x < 1 .or. dir_min_x > p_x_dim) ) then
      write(err_buf__,*) "dir_min_x must be between 1 and p_x_dim";call err__("overdense/os-spec-splitting.f03",292)
      call abort_program( p_err_invalid )
    endif

    this%if_recombine = if_recombine
    this%n_recombine = n_recombine
    this%recombine_start_time = recombine_start_time

    if ( this%if_recombine .and. n_recombine_regions > p_max_recombine ) then
      write(err_buf__,*) "n_recombine_regions cannot exceed ",p_max_recombine;call err__("overdense/os-spec-splitting.f03",301)
      call abort_program( p_err_invalid )
    endif

    recombine => null()
    tail => null()
    if (this%if_recombine) then

      do j = 1, n_recombine_regions

        allocate(recombine)
        recombine%recombination_box = recomb_box_bounds(:,:,j)
        recombine%recombination_cell = recomb_cell_dx(:,j)

        do i = 1, p_x_dim
          if (recombine%recombination_box(i,p_upper) - recombine%recombination_box(i,p_lower) .le. 0.) then
            write(err_buf__,*) "Bad recombination box";call err__("overdense/os-spec-splitting.f03",317)
            call abort_program( p_err_invalid )
          endif
          if (recombine%recombination_cell(i) .le. 0.) then
            write(err_buf__,*) "Bad recombination cell";call err__("overdense/os-spec-splitting.f03",321)
            call abort_program( p_err_invalid )
          endif
        enddo

        if (associated(tail)) then
          tail%next => recombine
        else
          this%recombine => recombine
        endif
        tail => recombine

      enddo

      if (this%n_recombine .lt. 1) then
        write(err_buf__,*) "n_recombine must be a positive integer";call err__("overdense/os-spec-splitting.f03",336)
        call abort_program( p_err_invalid )
      endif

    endif

  elseif ( n_split .gt. 0) then

    write(err_buf__,*) "n_split appears to be set in the input deck, but something isn't right";call err__("overdense/os-spec-splitting.f03",344)
    write(err_buf__,*) "Please check the input deck and try launching again,";call err__("overdense/os-spec-splitting.f03",345)
    write(err_buf__,*) "Or else set n_split to 0";call err__("overdense/os-spec-splitting.f03",346)
    call abort_program( -1)
  endif

  if ( if_recombine .and. (.not. this%if_split) ) then
    write(err_buf__,*) "Recombination was requested without splitting";call err__("overdense/os-spec-splitting.f03",351)
    call abort_program( p_err_invalid )
  endif

end subroutine read_nml_splitting
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_split( this, g_space, grid, dx )
!-----------------------------------------------------------------------------------------
  type( t_split ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  real(p_double), dimension(p_x_dim), intent(in) :: dx

  integer :: cell_count, i, dir
  real(p_k_part), dimension(2,p_x_dim) :: ls_box
  type(t_recombine), pointer :: recombine

  if ( .not. this%if_split ) return

  if ( this%use_min_x ) then

    dir = this%dir_min_x

    if ( this%min_x < xmin(g_space,dir) ) then
      write(err_buf__,*) "Using min_x, but it is set to a value left of the simulation boundary";call err__("overdense/os-spec-splitting.f03",377)
      call abort_program( p_err_invalid )
    elseif ( this%min_x > xmax(g_space,dir) ) then
      write(err_buf__,*) "Using min_x, but it is set to a value right of the simulation boundary";call err__("overdense/os-spec-splitting.f03",380)
      call abort_program( p_err_invalid )
    endif

    if (this%forward_of_min_x) then
      if ( (xmin(g_space,dir) + dx(dir)*grid%my_nx(p_upper,dir)) .le. this%min_x ) then
        this%if_split_local = .false.
      else
        this%if_split_local = .true.
      endif
    else
      if ( (xmin(g_space,dir) + dx(dir)*(grid%my_nx(p_lower,dir)-1)) .ge. this%min_x ) then
        this%if_split_local = .false.
      else
        this%if_split_local = .true.
      endif
    endif

  else

    this%if_split_local = .true.

  endif

  if ( this%if_recombine ) then

    if ( if_move_any( g_space ) ) then
      write(err_buf__,*) "Recombinations don't work with moving window.";call err__("overdense/os-spec-splitting.f03",407)
      write(err_buf__,*) "Aborting.";call err__("overdense/os-spec-splitting.f03",408)
      call abort_program( p_err_invalid )
    endif

    do i = 1, p_x_dim
      ls_box(p_lower, i) = real(xmin(g_space, i) + dx(i)*grid%my_nx(p_lower, i), p_k_part)
      ls_box(p_upper, i) = real(xmin(g_space, i) + dx(i)*grid%my_nx(p_upper, i), p_k_part)
    enddo

    recombine => this%recombine
    do
      if (.not. associated(recombine)) exit
      call local_boundaries(recombine, ls_box)

      if ( .not. recombine%if_recombine_local ) return

      cell_count = 1
      do i = 1, p_x_dim
        recombine%l_grid_size(i) = int(&
          (recombine%lr_box(i,p_upper) - recombine%lr_box(i,p_lower))/recombine%recombination_cell(i) + 2 )
        ! The extra 2 added in shouldn't be needed if there are no errors
        ! but there were errors before and I haven't confirmed there aren't some still
        cell_count = cell_count * recombine%l_grid_size(i)
      enddo

      call alloc(recombine%cells, (/ cell_count /),"overdense/os-spec-splitting.f03",433)

      recombine%cells = 0

      recombine%cell_count = cell_count

      recombine => recombine%next

    enddo

  endif

end subroutine init_split
!-----------------------------------------------------------------------------------------

subroutine local_boundaries_split( this, ls_box )

  type( t_recombine ), intent(inout) :: this
  real(p_k_part), dimension(:,:), intent(in) :: ls_box !LocalSpace_box

  integer :: i

  ! LocalRecombination_box, a subset of ls_box
  ! this%lr_box(p_x_dim, 2)
  this%lr_box = 0.

  do i = 1, p_x_dim

    if ( ( ls_box(p_upper, i) .le. this%recombination_box(i, p_lower ) ) .or. &
         ( ls_box(p_lower, i) .ge. this%recombination_box(i, p_upper ) ) ) then
      this%if_recombine_local = .false.
      return
    endif

    if ( ls_box(p_upper, i) .le. this%recombination_box(i, p_upper) ) then
      this%lr_box(i, p_upper) = ls_box(p_upper, i)
    else
      this%lr_box(i, p_upper) = this%recombination_box(i, p_upper)
    endif

    if ( ls_box(p_lower, i) .ge. this%recombination_box(i, p_lower) ) then
      this%lr_box(i, p_lower) = ls_box(p_lower, i)
    else
      this%lr_box(i, p_lower) = this%recombination_box(i, p_lower)
    endif

  enddo

  this%if_recombine_local = .true.

end subroutine local_boundaries_split


!-----------------------------------------------------------------------------------------
subroutine cleanup_split( this )

  type(t_split) :: this
  type(t_recombine), pointer :: recombine, next

  call freemem(this%splitting_gammas,"overdense/os-spec-splitting.f03",492)

  recombine => this%recombine
  do
    if (.not. associated(recombine)) exit
    call freemem(recombine%cells,"overdense/os-spec-splitting.f03",497)
    next => recombine%next
    deallocate(recombine)
    recombine => next
  enddo

end subroutine cleanup_split
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function if_split( this, n )

  type( t_split ), intent(in) :: this
  integer, intent(in) :: n
  logical :: if_split

  if_split = .false.

  if (.not. this%if_split) return

  if (.not. this%if_split_local) return

  if ( this%n_split > 0 ) then
    if ( mod( n, this%n_split ) == 0 ) if_split = .true.
  endif

end function if_split
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function if_recombine_local( this, recombine, n )

  type( t_split ), intent(in) :: this
  type( t_recombine ), intent(in) :: recombine
  integer, intent(in) :: n
  logical :: if_recombine_local

  if_recombine_local = .false.

  if (.not. this%if_recombine) return

  if (.not. recombine%if_recombine_local ) return

  if ( (this%n_recombine > 0) .and. (n .gt. 0) ) then
    if ( mod( n, this%n_recombine ) == 0 ) if_recombine_local = .true.
  endif

end function if_recombine_local
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function floor_log_2(x)
!-----------------------------------------------------------------------------------------
! Returns the floor of the log base 2 of a real number
! Returns 0 for anything less than 1 (including negative numbers)
!-----------------------------------------------------------------------------------------
  real(p_k_part) :: x

  integer :: int_x

  floor_log_2 = 0
  int_x = int(x)
  do
    if ( int_x .lt. 1) exit
    floor_log_2 = floor_log_2 + 1
    int_x = ISHFT(int_x, -1)
  enddo

end function floor_log_2
!-----------------------------------------------------------------------------------------

end module m_spec_splitting
