# 1 "cyl_modes/os-particles-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-particles-cyl-modes.f03" 2
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
# 2 "cyl_modes/os-particles-cyl-modes.f03" 2
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
# 3 "cyl_modes/os-particles-cyl-modes.f03" 2

module m_particles_cyl_modes

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
# 7 "cyl_modes/os-particles-cyl-modes.f03" 2

use m_system
use m_parameters
use m_space, only: t_space, x_dim, nx_move, if_move
use m_grid_define, only: t_grid
use m_node_conf, only: t_node_conf, n_threads
use m_restart, only: t_restart_handle, restart_io_write, restart_io_read
use m_time_step, only: t_time_step, dt, n
use m_input_file, only: t_input_file, get_namelist
use m_particles_define, only: t_particles, pushev, cathodeev, reduce_current_ev
use m_particles, only: validate
use m_species_cyl_modes_define, only: t_species_cyl_modes, t_spec_cyl_arr
use m_logprof, only: begin_event, end_event
use m_neutral_cyl_modes, only: t_neutral_cyl_modes
use m_emf_define, only: t_emf
use m_current_define, only: t_current
use m_species_define, only: t_species
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_current_cyl_modes, only: t_current_cyl_modes
use m_cyl_modes, only: reduce, check_nan, t_cyl_modes
use m_cathode_cyl_modes, only: t_cathode_cyl_modes

implicit none

private

type, extends( t_particles ) :: t_particles_cyl_modes

  ! all species must be of type species_cyl_modes
  ! all cathodes must be of type cathode_cyl_modes
  ! all neutrals must be of type neutral_cyl_modes
  type(t_cyl_modes) :: jay_tmp_cyl_m

contains

  procedure :: allocate_objs => allocate_objs_part_cyl_modes
  procedure :: read_input => read_input_part_cyl_modes
  procedure :: advance_deposit => advance_deposit_part_cyl_modes
  procedure :: inject_cathode => inject_cathode_part_cyl_modes

end type t_particles_cyl_modes

! Add interface to superclass read_input function to avoid type casting
interface
  subroutine read_input_particles( this, input_file, periodic, if_move, grid, dt, &
                               sim_options )

  import t_particles, t_input_file, p_double, t_options, t_grid

  class( t_particles ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  real(p_double), intent(in) :: dt
  type( t_options ), intent(in) :: sim_options

  end subroutine
end interface

public :: t_particles_cyl_modes

contains

!-----------------------------------------------------------------------------------------
! Overrides allocate_objs method and allocates t_species_cyl_modes and t_neutral_cyl_modes objects
! instead
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_cyl_modes( this )

  implicit none

  class( t_particles_cyl_modes ), intent(inout) :: this

  integer :: i
  class(t_species_cyl_modes), pointer :: spec, tail

  if ( this%num_species > 0 ) then

    allocate( t_species_cyl_modes :: spec )
    this % species => spec
    do i = 2, this%num_species
      tail => spec
      allocate( t_species_cyl_modes :: spec )
      tail % next => spec
    enddo

    if ( this%num_cathode > 0 ) then
      allocate(t_cathode_cyl_modes :: this%cathode( this%num_cathode ))
    endif


    if ( this%num_neutral > 0 ) then
    allocate(t_neutral_cyl_modes :: this%neutral( this%num_neutral ))
    endif


  endif

end subroutine allocate_objs_part_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Don't allow quartic interpolation
!-----------------------------------------------------------------------------------------
subroutine read_input_part_cyl_modes( this, input_file, periodic, if_move, grid, dt, &
                                      sim_options )

  implicit none

  class( t_particles_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  real(p_double), intent(in) :: dt
  type( t_options ), intent(in) :: sim_options

  call read_input_particles( this, input_file, periodic, if_move, grid, dt, sim_options )

  if ( this%interpolation == p_quartic ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error reading particles parameters"
      write(0,*) "   Allowed interpolation orders are 'linear', 'quadratic', and 'cubic'"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  if ( this%grid_center ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error reading particles parameters"
      write(0,*) "   grid_center must be set to false when using cyl_modes"
      write(0,*) "   aborting..."
    endif
    stop
  endif

end subroutine read_input_part_cyl_modes
!-----------------------------------------------------------------------------------------

subroutine advance_deposit_part_cyl_modes( this, emf, jay, tstep, t, no_co, options )

  implicit none

  class( t_particles_cyl_modes ), intent(inout) :: this
  class( t_emf ), intent( inout ) :: emf
  class( t_current ), intent(inout) :: jay

  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t

  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  integer :: n_modes, nt, tid, i
  class( t_species ), pointer :: species
  type( t_spec_cyl_arr ), dimension(1:this%num_species) :: spec_arr
# 181 "cyl_modes/os-particles-cyl-modes.f03"
  select type( jay )
  class is ( t_current_cyl_modes )

  select type( emf )
  class is ( t_emf_cyl_modes )

    n_modes = jay%n_cyl_modes

    nt = n_threads( no_co )
    if ( nt > 1 ) then

      call begin_event( pushev )

      ! Multi-threaded pusher

      ! Get child-class pointers to all species outside of OpenMP loop
      species => this % species
      do i = 1, this%num_species
        if ( .not. associated(species)) exit
        select type( species )
        class is( t_species_cyl_modes )
          spec_arr(i)%s => species
        end select
        species => species % next
      enddo

      !$omp parallel do private(i)
      do tid = 1, nt
        ! Zero the current grid for this thread and all modes
        call jay%jay_cyl_m_arr(tid)%pf_re(0)%zero()
        do i = 1, n_modes
          call jay%jay_cyl_m_arr(tid)%pf_re(i)%zero()
          call jay%jay_cyl_m_arr(tid)%pf_im(i)%zero()
        enddo

        do i = 1, this%num_species
          if ( .not. associated(spec_arr(i) % s)) exit
          call spec_arr(i) % s % push_cyl( emf, jay % jay_cyl_m_arr(tid), t, tstep, &
                                            tid-1, nt, options )
        enddo
      enddo
      !$omp end parallel do

      call end_event( pushev )

      ! Add current from all threads for all modes
      ! This is also done with OpenMP parallelism
      call begin_event( reduce_current_ev )
      call reduce( jay % jay_cyl_m_arr )
      call end_event( reduce_current_ev )

    else

      call begin_event( pushev )

      ! Single-threaded pusher

      species => this % species

      call jay%jay_cyl_m%pf_re(0)%zero()
      do i = 1, n_modes
        call jay%jay_cyl_m%pf_re(i)%zero()
        call jay%jay_cyl_m%pf_im(i)%zero()
      enddo

      if ( this%low_jay_roundoff ) then
        if (mpi_node()==0) print *,"(*warning*) Low jay roundoff not implemented in cylindrical mode, using default deposit."
      endif

      do
        if ( .not. associated(species)) exit
        select type( species )
        class is ( t_species_cyl_modes )
          call species % push_cyl( emf, jay%jay_cyl_m, t, tstep, 0, 1, options )
        end select
        species => species % next
      enddo

      call end_event( pushev )

    endif

    ! Advance iteration information
    this%n_current = this%n_current + 1





  end select

  end select






end subroutine advance_deposit_part_cyl_modes

subroutine inject_cathode_part_cyl_modes( this, jay, t, dt, no_co, coordinates )

  implicit none

  class( t_particles_cyl_modes ), intent(inout) :: this
  class( t_current ), intent(inout), target :: jay
  real( p_double ), intent(in) :: t
  real( p_double ), intent(in) :: dt
  integer, intent(in) :: coordinates
  class( t_node_conf ), intent(in) :: no_co

  ! local variables

  integer :: i
  type( t_cyl_modes ), pointer :: jay_cyl_m

  ! executable statements

  call begin_event(cathodeev)

  ! Get pointer to jay_cyl_m
  select type(obj => jay)
  class is (t_current_cyl_modes)
    jay_cyl_m => obj % jay_cyl_m
  class default
    jay_cyl_m => null() ! Should never happen, just to suppress compiler warning
  end select

  if ( this%low_jay_roundoff ) then

    ! this version has a lower roundoff error because the current from each cathode
    ! is deposited onto a clean grid, and the grids are then added
    do i = 1, this%num_cathode
      select type (obj => this%cathode(i))
      class is (t_cathode_cyl_modes)
        call obj%inject_cyl_lowroundoff(jay_cyl_m, this%jay_tmp_cyl_m, t, dt, no_co, coordinates)
      end select
    enddo

  else

    do i = 1, this%num_cathode
      select type (obj => this%cathode(i))
      class is (t_cathode_cyl_modes)
      call obj%inject_cyl(jay_cyl_m, t, dt, no_co, coordinates)
    end select
  enddo

  endif

  call end_event(cathodeev)

end subroutine inject_cathode_part_cyl_modes

end module m_particles_cyl_modes
