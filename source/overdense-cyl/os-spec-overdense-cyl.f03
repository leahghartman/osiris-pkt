
#include "os-config.h"
#include "os-preprocess.fpp"

module m_spec_overdense_cyl

#include "memory/memory.h"

use m_system
use m_parameters
use m_species_cyl_modes_define
use m_spec_overdense
use m_species_define, only : t_part_idx, t_spec_msg, p_cell_low, p_cell_near
use m_random, only : rng
use m_input_file, only : t_input_file
use m_grid_define, only : t_grid, t_msg_patt
use m_spec_fluid_moments
use m_node_conf, only : t_node_conf, reduce
use m_space, only : t_space
use m_emf_define, only : t_emf
use m_current_define, only : t_current
use m_vdf_define, only : t_vdf
use m_vdf_comm, only : t_vdf_msg, reshape_nocopy
use m_restart, only : t_restart_handle
use m_time_step, only : t_time_step, n, dt
use m_spec_damp
use m_spec_splitting
use m_math, only : pi, pi_2
use m_fparser, only : eval, p_k_fparse
use m_psource_std, only : t_psource_std
use m_emf_cyl_modes, only : t_emf_cyl_modes
use m_cyl_modes, only : t_cyl_modes

#ifdef __HAS_TRACKS__
use m_species_tracks, only : missing_particles, inbound_particle
#endif

implicit none

private

! species type, extends the cyl_modes type and add arguments/methods from overdense
type, extends( t_species_cyl_modes ) :: t_species_overdense_cyl

  type(t_split) :: splitter
  type(t_damper), pointer :: damper => null()
  type(t_fluid_moments) :: fluid_moments

contains

  procedure :: read_input => read_input_spec_overdense_cyl
  procedure :: init => init_spec_overdense_cyl
  procedure :: push_cyl => push_species_overdense_cyl
  procedure :: cleanup => cleanup_spec_overdense_cyl
  procedure :: list_algorithm => list_algorithm_spec_overdense_cyl
  procedure :: reshape => reshape_spec_overdense_cyl

  ! damper routines
  procedure :: prep_damper
  procedure :: damp => damp_species

  ! splitting/recombining routines
  procedure :: push_particles_apart
  procedure :: split => split_species
  procedure :: recombine => recombine_species

end type t_species_overdense_cyl

public :: t_species_overdense_cyl

contains

!-----------------------------------------------------------------------------------------
! Read information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_spec_overdense_cyl( this, input_file, def_name, periodic, if_move, &
                                          grid, dt, read_prof, sim_options )

  implicit none

  class( t_species_overdense_cyl ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  character(len = *), intent(in) :: def_name
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  type(t_options), intent(in) :: sim_options
  real(p_double), intent(in) :: dt
  logical, intent(in) :: read_prof

  type(t_damper), pointer :: item, tail
  real(p_k_part) :: density

  call this % t_species_cyl_modes % read_input( input_file, def_name, periodic, if_move, &
                                                grid, dt, read_prof, sim_options )

  ! Read possibly multiple damper namelists
  tail => null()
  do
    select case (trim(input_file % get_section_name()))
    case('nl_damper')
      allocate(item)
    case default
      ! Not a damper section, exit
      exit
    end select

    ! Read input for damper
    call read_nml( item, input_file )
    if (item%if_cold_perp) then
      if (mpi_node()==0) then
        write(0,*) '(*warning*) The damper if_cold_perp parameter is ignored in cyl geometry.'
        write(0,*) '(*warning*) Setting if_cold_perp = .false.'
      endif
      item%if_cold_perp = .false.
    endif

    if (associated(tail)) then
      tail % next => item
    else
      this % damper => item
    endif

    tail => item

  enddo

  ! Read additional namelist for splitting
  ! First, supply density from default value if source is or inherits from t_psource_std.
  ! Note, this density value may not be desired by the user.
  select type( source => this%source )
  class is( t_psource_std )
    density = source % density
  class default
    density = 0.0_p_k_part
  end select

  call read_nml( this%splitter, input_file, this%if_collide, this%rqm, &
                 product(this%num_par_x(1:p_x_dim)), density )

end subroutine read_input_spec_overdense_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       sets up this data structure from the given information
!-----------------------------------------------------------------------------------------
subroutine init_spec_overdense_cyl( this, sp_id, interpolation, grid_center, grid, &
                                    g_space, emf, jay, no_co, send_vdf, recv_vdf, &
                                    bnd_cross, node_cross, send_spec, recv_spec, &
                                    ndump_fac, restart, restart_handle, sim_options, &
                                    tstep, tmin, tmax )

  implicit none

  class( t_species_overdense_cyl ), intent(inout) :: this

  integer, intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in) :: grid
  type( t_space ), intent(in) :: g_space
  class( t_emf ), intent(inout) :: emf
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_vdf, recv_vdf
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_spec, recv_spec
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type(t_options), intent(in) :: sim_options
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax

  type(t_damper), pointer :: damper

  call this % t_species_cyl_modes % init( sp_id, interpolation, grid_center, grid, &
                                          g_space, emf, jay, no_co, send_vdf, recv_vdf, &
                                          bnd_cross, node_cross, send_spec, recv_spec, &
                                          ndump_fac, restart, restart_handle, &
                                          sim_options, tstep, tmin, tmax )

  damper => this%damper
  do
    if (.not. associated(damper)) exit
    call init( damper, this%fluid_moments, grid, no_co, this%dx, g_space )
    damper => damper%next
  enddo

  call init( this%splitter, g_space, grid, this%dx )

end subroutine init_spec_overdense_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Particle pusher for cyl_modes simulations. Note that when using these simulations the
! push_type parameter is ignored
!-----------------------------------------------------------------------------------------
subroutine push_species_overdense_cyl( this, emf, jay_cyl_m, t, tstep, tid, n_threads, &
                                       options )

  implicit none

  class( t_species_overdense_cyl ), intent(inout) :: this
  class( t_emf_cyl_modes ), intent( inout )  ::  emf
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  type( t_options ), intent(in) :: options
  integer, intent(in) :: tid        ! local thread id
  integer, intent(in) :: n_threads  ! total number of threads

  integer :: chunk, i0, i1
  type(t_damper), pointer :: damper

  ! Skip if before push start time
  if ( t >= this%push_start_time ) then

    ! Extra routine for damping hot particles
    damper => this%damper
    do
      if (.not. associated(damper)) exit
      if ( if_damp( damper, n(tstep), t ) ) then

        ! range of particles for each thread
        chunk = ( this%num_par + n_threads - 1 ) / n_threads
        i0    = tid * chunk + 1
        i1    = min( (tid+1) * chunk, this%num_par )

        call this % damp( damper, dt(tstep), i0, i1 )

      endif
      damper => damper%next
    enddo

  endif

  ! Call superclass push
  call this % t_species_cyl_modes % push_cyl( emf, jay_cyl_m, t, tstep, tid, n_threads, &
                                              options )

end subroutine push_species_overdense_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Call superclass cleanup, then cleanup extra members
!-----------------------------------------------------------------------------------------
subroutine cleanup_spec_overdense_cyl( this )

  implicit none

  class( t_species_overdense_cyl ), intent(inout) :: this

  type(t_damper), pointer :: damper, next

  ! Call superclass cleanup
  call this % t_species_cyl_modes % cleanup()

  ! Cleanup extra damping member data
  call cleanup( this%fluid_moments )
  damper => this%damper
  do
    if (.not. associated(damper)) exit
    next => damper%next
    deallocate(damper)
    damper => next
  enddo

  ! Cleanup extra splitting member data
  call cleanup( this%splitter )

end subroutine cleanup_spec_overdense_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! List damping algorithm
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_spec_overdense_cyl( this )

  implicit none

  class( t_species_overdense_cyl ), intent(in) :: this

  ! Call superclass list_algorithm
  call this % t_species_cyl_modes % list_algorithm()

  ! List damping algorithm
  if (associated(this%damper)) then
    if ( this%damper%if_damp ) then
      print *, '- Particle damping enabled'
    endif
  endif

  ! List splitting algorithm
  if ( this%splitter%if_split .and. this%splitter%n_split > 0 ) then
    if ( this%splitter%if_recombine .and. this%splitter%n_recombine > 0 ) then
      print *, '- Particle splitting and recombination enabled'
    else
      print *, '- Particle splitting enabled'
    endif
  endif

end subroutine list_algorithm_spec_overdense_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! redistribute the species particle data through all nodes when node grids change and
! store new grid boundaries
!-----------------------------------------------------------------------------------------
subroutine reshape_spec_overdense_cyl( this, old_grid, new_grid, msg_patt, no_co, &
                                       send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species_overdense_cyl ), intent(inout) :: this
  type( t_msg_patt ), intent(in) :: msg_patt
  class( t_grid ), intent(in) :: old_grid, new_grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, dir, cell_count
  real(p_k_part), dimension(2,p_x_dim) :: ls_box
  type(t_damper), pointer :: damper
  type(t_recombine), pointer :: recombine

  call this % t_species_cyl_modes % reshape( old_grid, new_grid, msg_patt, no_co, &
                                             send_msg, recv_msg )

  ! Fill in local xmax using g_box instead of g_space
  do i = 1, p_x_dim
    ls_box(p_lower,i) = real( this%g_box(p_lower,i) + &
                        this%dx(i)*(new_grid%my_nx(p_lower,i)-0.5_p_double), p_k_part )
    ls_box(p_upper,i) = real( this%g_box(p_lower,i) + &
                        this%dx(i)*(new_grid%my_nx(p_upper,i)-0.5_p_double), p_k_part )
  enddo
  if ( this%coordinates == p_cylindrical_b .and. this%pos_type == p_cell_near ) then
    ls_box(p_lower,p_r_dim) = real( this%g_box(p_lower,p_r_dim) + &
                              this%dx(p_r_dim)*new_grid%my_nx(p_lower,p_r_dim), p_k_part )
    ls_box(p_upper,p_r_dim) = real( this%g_box(p_lower,p_r_dim) + &
                              this%dx(p_r_dim)*new_grid%my_nx(p_upper,p_r_dim), p_k_part )
  endif

  ! Take care of damper logic
  damper => this%damper
  do
    if (.not. associated(damper)) exit
    if ( damper%if_damp ) then

      if( damper%reemit_with_local_temperature .or. damper%local_vth_cutoff) then
        if ( damper%fourth_root_temperature ) then
          call reshape_nocopy( this%fluid_moments%density, new_grid )
          call reshape_nocopy( this%fluid_moments%fourth_root_abs_momentum, new_grid )
        else
          call reshape_nocopy( this%fluid_moments%density, new_grid )
          call reshape_nocopy( this%fluid_moments%abs_momentum, new_grid )
        endif
      endif

      dir = damper%damper_direction

      ! Recalculate if_damp_local since the grid has changed
      if ( (ls_box(p_upper,dir) .le. damper%x_start .and. &
            damper%orientation==p_forward) .or. &
           (ls_box(p_lower,dir) .ge. damper%x_start .and. &
            damper%orientation==p_backward) ) then
        damper%if_damp_local = .false.
      else
        damper%if_damp_local = .true.
      endif

    endif

    damper => damper%next

  enddo

  ! Take care of splitting logic
  if ( this%splitter%if_split ) then

    if ( this%splitter%use_min_x ) then

      dir = this%splitter%dir_min_x

      if ( this%splitter%forward_of_min_x ) then
        if ( ls_box(p_upper,dir) .le. this%splitter%min_x ) then
          this%splitter%if_split_local = .false.
        else
          this%splitter%if_split_local = .true.
        endif
      else
        if ( ls_box(p_lower,dir) .ge. this%splitter%min_x ) then
          this%splitter%if_split_local = .false.
        else
          this%splitter%if_split_local = .true.
        endif
      endif

    endif

  endif

  ! Take care of recombination logic
  if ( this%splitter%if_recombine ) then

    recombine => this%splitter%recombine
    do
      if (.not. associated(recombine)) exit

      if ( recombine%if_recombine_local ) then
        call freemem( recombine%cells )
      endif

      call local_boundaries(recombine, ls_box)

      if ( recombine%if_recombine_local ) then

        cell_count = 1
        do i = 1, p_x_dim
          recombine%l_grid_size(i) = int((recombine%lr_box(i,p_upper) - &
            recombine%lr_box(i,p_lower))/recombine%recombination_cell(i) + 2 )
          cell_count = cell_count * recombine%l_grid_size(i)
        enddo

        call alloc( recombine%cells, (/ cell_count /) )

        recombine%cells = 0

        recombine%cell_count = cell_count

      endif

      recombine => recombine%next

    enddo

  endif

end subroutine reshape_spec_overdense_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Hazard function for electron "survival" based on a probability density function of
! pdf(t) = sin(pi*t)**2, where t = x / l_damp, where l_damp is the length of the damping
! region.  The mean free path will be l_damp/2.
! The hazard function is pdf(t) / (1 - cdf(t)), for the cumulative distribution function
! cdf(t).
!-----------------------------------------------------------------------------------------
function hazard_sin2( x, l_damp )

  implicit none

  real(p_k_part) :: hazard_sin2
  real(p_k_part), intent(in) :: x, l_damp

  real(p_k_part) :: t, pi_

  t = x / l_damp
  pi_ = real(pi,p_k_part)
  hazard_sin2 = 4 * pi_ * sin(pi_*t)**2 / ( 2*pi_*(1-t) + sin(2*pi_*t) )

end function hazard_sin2
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Calculate the necessary fluid moments if needed for damper
!-----------------------------------------------------------------------------------------
subroutine prep_damper( this, n, t )

  class( t_species_overdense_cyl ), intent(inout) :: this
  integer, intent(in) :: n
  real(p_double), intent(in) :: t

  type(t_damper), pointer :: damper

  damper => this%damper

  do
    if (.not. associated(damper)) exit

    if ( if_damp( damper, n, t ) ) then

      if (damper%reemit_with_local_temperature .or. damper%local_vth_cutoff) then
        if( damper%fourth_root_temperature ) then
          call calculate( this%fluid_moments, this, 'fourth_root_abs_proper_velocity', n)
        else
          call calculate( this%fluid_moments, this, 'abs_proper_velocity', n)
        endif
      endif

    endif

    damper => damper%next

  enddo

end subroutine prep_damper
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Particle energy absorption subroutine
!-----------------------------------------------------------------------------------------
subroutine damp_species( this, damper, dt, i0, i1 )

  class( t_species_overdense_cyl ), intent(inout) :: this
  type( t_damper ), intent(inout) :: damper
  real(p_double), intent(in) :: dt
  integer, intent(in) :: i0, i1

  integer :: i, j
  real(p_k_part), dimension(4) :: x_vector
  real(p_k_part) :: x ! position of the particle
  real(p_k_part) :: psq
  real(p_k_part) :: gamma ! gamma of the particle

  real(p_k_part) :: velocity ! beta in the direction of the absober
  ! could be changed to just beta, which would be *slightly* more physical
  real(p_k_part) :: pos_interpol ! position interpolation factor, to find the local
                                  ! stopping distance etc.
  real(p_k_part) :: stopping_distance ! stopping distance at particle position
  real(p_k_part), dimension(p_p_dim) :: reemit_uth ! reemision thermal velocity at
                                                    ! particle position
  real(p_k_part) :: prob ! stopping probabilty
  real(p_k_part) :: rand ! named for the senator, of course
  integer :: dir

  ! average abs value of proper velocity (fourth-root proper velocity) in a given cell
  real(p_k_part), dimension(p_p_dim) :: local_avg_u, local_avg_fourth_root_u
  ! average proper thermal velocity, calculated from the above
  real(p_k_part), dimension(p_p_dim) :: local_avg_uth
  logical :: skip_part, use_local_vth

  dir = damper%damper_direction ! just a nickname

  do i = i0, i1

    call this % get_position( i, x_vector )

    x = x_vector(dir)

    if (damper%orientation == p_forward) then

      if ( x .lt. damper%x_start ) cycle

      select case( damper%stop_dir )
      case(p_all)
        ! do nothing
      case( p_forward )
        if ( this%p(dir,i) < 0.0_p_k_part ) cycle
      case( p_backward )
        if ( this%p(dir,i) > 0.0_p_k_part ) cycle
      end select

      if ( x .lt. damper%x_full ) then
        pos_interpol = (damper%x_full - x)/(damper%x_range)
      else
        pos_interpol = 0.0_p_k_part
      endif

    else ! backward

      if ( x .gt. damper%x_start ) cycle

      select case( damper%stop_dir )
      case(p_all)
        ! do nothing
      case( p_forward )
        if ( this%p(dir,i) > 0.0_p_k_part ) cycle
      case( p_backward )
        if ( this%p(dir,i) < 0.0_p_k_part ) cycle
      end select

      if ( x .gt. damper%x_full ) then
        pos_interpol = (damper%x_full - x)/(damper%x_range)
      else
        pos_interpol = 0.0_p_k_part
      endif

    endif

    psq = 0.0_p_k_part
    do j = 1, p_p_dim
      psq = psq + this%p(j,i)**2
    enddo
    gamma = sqrt( 1.0_p_k_part + psq )

    if( damper%local_vth_cutoff ) then

      ! Check if we actually want the local vth based on the expression
      use_local_vth = .true.

      if ( damper%use_local_vth_func ) then
        if ( eval( damper%local_vth_func, real(x_vector,p_k_fparse) ) > 0.0_p_k_fparse ) then
          use_local_vth = .true.
        else
          use_local_vth = .false.
        endif
      endif

      if ( use_local_vth ) then

        if ( damper%fourth_root_temperature ) then
          call read_fluid_moment( this%fluid_moments, this, i, &
            'fourth_root_abs_proper_velocity', local_avg_fourth_root_u )

          ! This comes from solving for uth after integrating |u|^1/4*exp(-u^2/(2uth^2))
          local_avg_uth = 1.64801_p_k_part * local_avg_fourth_root_u**4

        else
          call read_fluid_moment( this%fluid_moments, this, i, 'abs_proper_velocity', &
                                  local_avg_u )

          ! This comes from solving for uth after integrating |u|*exp(-u^2/(2uth^2))
          local_avg_uth = real( sqrt(pi_2), p_k_part ) * local_avg_u

        endif

        ! Check if momentum in any direction is greater than the cutoff
        skip_part = .true.
        do j = 1, p_p_dim
          if( abs(this%p(j,i)) .ge. (local_avg_uth(j)*(damper%vth_mult_full + &
                          pos_interpol * damper%vth_mult_range)) ) skip_part = .false.
        enddo
        ! If p is not high in any direction, then cycle
        if (skip_part) cycle

      else

        ! use gamma cutoff
        if( gamma .lt. (damper%gamma_full + &
                        pos_interpol * damper%gamma_range) ) cycle

      endif

    else
      ! Otherwise we use a gamma cutoff
      if( gamma .lt. (damper%gamma_full + &
                      pos_interpol * damper%gamma_range) ) cycle
    endif

    ! Compute absolute value of velocity
    velocity = sqrt( this%p(1,i)**2 + this%p(2,i)**2 + this%p(3,i)**2 ) / gamma

    stopping_distance = damper%stopping_dist_full + &
                        pos_interpol * damper%stopping_dist_range

    if ( damper%linear ) then
      prob = damper%n_damp * real(dt,p_k_part) * velocity / stopping_distance
    else
      prob = damper%n_damp * real(dt,p_k_part) * velocity / stopping_distance * &
             hazard_sin2( abs(x-damper%x_start), stopping_distance )
      if ( prob < 0.0_p_k_part ) then
        prob = 1.0_p_k_part
      endif
    endif
    rand = real( rng % genrand_real3(), p_k_part )

    if( rand .gt. prob ) cycle

    if (damper%reemit_with_local_temperature) then

      ! Only need to evaluate use_local_vth if it wasn't already defined
      if ( .not. damper%local_vth_cutoff ) then
        use_local_vth = .true.

        if ( damper%use_local_vth_func ) then
          if ( eval( damper%local_vth_func, real(x_vector,p_k_fparse) ) > 0.0_p_k_fparse ) then
            use_local_vth = .true.
          else
            use_local_vth = .false.
          endif
        endif
      endif

      ! Check whether to use local vth or not
      if ( use_local_vth ) then

        ! This lookup is expensive, so put in some logic for it
        ! The lookup may have been done already
        if ( .not. damper%local_vth_cutoff ) then

          if (damper%fourth_root_temperature) then
            call read_fluid_moment( this%fluid_moments, this, i, &
              'fourth_root_abs_proper_velocity', local_avg_fourth_root_u )

            local_avg_uth = 1.64801_p_k_part * local_avg_fourth_root_u**4

          else
            call read_fluid_moment( this%fluid_moments, this, i, 'abs_proper_velocity', &
                                    local_avg_u )

            local_avg_uth = real(sqrt(pi_2),p_k_part) * local_avg_u

          endif

        endif

        reemit_uth = local_avg_uth

        ! if_cold_perp does not apply in quasi-3D geometry
        ! if ( damper%if_cold_perp ) then
        !   do j = p_x_dim + 1, p_p_dim
        !     reemit_uth(j) = 0.0_p_k_part
        !   enddo
        ! endif

      else

        reemit_uth = damper%reemit_temp_full + &
                     pos_interpol * damper%reemit_temp_range

      endif

    else
      reemit_uth = damper%reemit_temp_full + &
                   pos_interpol * damper%reemit_temp_range
    endif

    do j = 1, p_p_dim
      this%p(j,i) = reemit_uth(j) * real( rng % genrand_gaussian(), p_k_part )
    enddo

    ! Reemit in correct direction if specified
    if (damper%orientation == p_forward) then
      select case(damper%reemit_dir)
      case(p_all)
        ! do nothing
      case(p_forward)
        this%p(dir,i) = abs(this%p(dir,i))
      case(p_backward)
        this%p(dir,i) = -abs(this%p(dir,i))
      end select

    else ! backward
      select case(damper%reemit_dir)
      case(p_all)
        ! do nothing
      case(p_forward)
        this%p(dir,i) = -abs(this%p(dir,i))
      case(p_backward)
        this%p(dir,i) = abs(this%p(dir,i))
      end select

    endif

  enddo

end subroutine damp_species
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine push_particles_apart( this, i, j )

  class( t_species_overdense_cyl ), intent(inout) :: this
  integer, intent(in) :: i, j

  real( p_k_part ) :: mag_i !magnitude of i (and also j, I hope)
  real( p_k_part ), dimension(3) :: u, v
  real( p_k_part ) :: a
  real( p_k_part ) :: phase
  real( p_k_part ), parameter :: zero = 0.0, one = 1.0, small = 0.00000001

  a = this%splitter%boost
  u = 0 !why does it cast these correctly? oh well
  v = 0

  mag_i = magnitude(this%p(:,i))
  ! this seems like a bad idea, right? this will be linear around zero, that's not good
  ! if ( this%splitter%random_boost ) a = a*(1 + genrand_gaussian() )
  ! how about this, it multiplies a evenly by between 0.5 and 1.5
  if (this%splitter%random_boost) a = a*(1.0_p_k_part + &
                                    0.5_p_k_part * real( rng%genrand_real1(), p_k_part ) )

  ! Only use 3-D geometry
  if (this%p(2,i) .ne. 0.0) then
    v = (/ one, zero, zero /)
    u = crossproduct( this%p(:,i), v )
  else if ( (this%p(1,i) .ne. 0.0) .or. (this%p(3,i) .ne. 0.0) ) then
    v = (/ zero, one, zero /)
    u = crossproduct( this%p(:,i), v )
  else
    u = (/ small, small, small /)
    this%p(:,i) = this%p(:,i) - u
    this%p(:,j) = this%p(:,j) + u
    return
  endif

  u = u * a * mag_i/magnitude(u)
  v = crossproduct( this%p(:,i), u )
  v = v * a * mag_i/magnitude(v)

  phase = 2.0_p_k_part * real( pi * rng%genrand_real2(), p_k_part )
  u = u*sin(phase) + v*cos(phase)

  if ( magnitude(u) .lt. small )  then
    u = small
  endif

  this%p(:,i) = this%p(:,i) - u
  this%p(:,j) = this%p(:,j) + u

  contains

  function crossproduct(u,v)
    real( p_k_part ), dimension(3) :: u,v
    real( p_k_part ), dimension(3) :: crossproduct

    crossproduct(1) = u(2)*v(3) - u(3)*v(2)
    crossproduct(2) = -u(1)*v(3) + u(3)*v(2)
    crossproduct(3) = u(1)*v(2) - u(2)*v(1)
  end function crossproduct

  function magnitude(a)
    real( p_k_part ), dimension(3) :: a
    real( p_k_part ) :: magnitude

    magnitude = sqrt( a(1)**2 + a(2)**2 + a(3)**2 )
  end function magnitude

  function dot3( a, b )
    implicit none

    real(p_double), dimension(3), intent(in) :: a, b
    real(p_double) :: dot3

    dot3 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

  end function dot3

end subroutine push_particles_apart
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine split_species( this, n )

  class( t_species_overdense_cyl ), intent(inout), target :: this
  integer, intent(in) :: n

  type( t_split ), pointer :: split

  integer :: original_num_par, num_split, splitting_state_i, split_state_i_power2
  integer :: i
  real(p_k_part) :: gamma_i, q_copy
  real(p_k_part) :: real_splitting_state_i
  real(p_k_part), dimension(4) :: position, x_copy
  integer, dimension(p_x_dim) :: ix_copy
  real(p_k_part), dimension(p_p_dim) :: p_copy

  if ( .not. if_split( this%splitter, n ) ) return

  split => this%splitter

  original_num_par = this%num_par
  num_split = size( split%splitting_gammas )

  do i = 1, original_num_par

    if ( split%use_min_x ) then
      call this % get_position( i, position )
      if ( split%forward_of_min_x ) then
        if (position(split%dir_min_x) .lt. split%min_x) cycle
      else
        if (position(split%dir_min_x) .gt. split%min_x) cycle
      endif
    endif

    gamma_i = sqrt( 1.0_p_k_part + this%p(1,i)**2 + this%p(2,i)**2 + this%p(3,i)**2 )
    real_splitting_state_i = split%base_particle_q/this%q(i) + 0.5
    split_state_i_power2 = int(real_splitting_state_i) ! base q = 1, split once = 2, twice = 4, thrice = 8
    splitting_state_i = 1
    do while (split_state_i_power2 .gt. 1)
      splitting_state_i = splitting_state_i + 1
      split_state_i_power2 = ISHFT(split_state_i_power2, -1)
    end do
    ! split_state_i_power2 = 1 if unsplit, 2 if split once, 3 if twice, 4 if thrice...

    if ( splitting_state_i .le. num_split ) then !I'm pretty sure that should be .le. ....
      if (gamma_i .gt. split%splitting_gammas(splitting_state_i) ) then
        this%q(i) = this%q(i)/2.0_p_k_part

        q_copy = this%q(i)
        p_copy(:) = this%p(:,i)
        x_copy(:) = this%x(:,i)
        ix_copy(:) = this%ix(:,i)

        ! Create the new particle
        call this % create_particle( ix_copy, x_copy, p_copy, q_copy )

        ! Push the particle (otherwise, hopefully collisions or something else is being used)
        if ( split%if_push ) call this % push_particles_apart( i, this%num_par )

      endif ! gamma_i .gt. ...
    endif ! splitting_state_i ...

  enddo

end subroutine split_species
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine recombine_species( species, recombine, n, t )

  class( t_species_overdense_cyl ), intent(inout), target :: species
  type( t_recombine ), intent(inout) :: recombine
  integer, intent(in) :: n
  real(p_double), intent(in) :: t

  type( t_split ), pointer :: split

  integer :: i, fac_grid_size
  integer :: idx ! particle InDeX
  real( p_k_part ), dimension(4) :: x ! Position vector for particle #idx
  real( p_k_part ) :: gamma ! gamma of particle #idx
  integer :: s_s_l2 ! Splitting_State_Log_2 of particle #idx
  integer :: or_idx ! Old_Resident_IDX (or ORiginal_IDX, whatever)
  integer :: or_s_s_l2 ! Old_Resident_Splitting_State_Log_2
  integer :: num_split ! Total number of splitting gammas
  ! real( p_k_part ) :: or_gamma ! Old_Resident_gamma
  ! real( p_k_part ) :: E_i, E_f ! Energy before, after recombination
  real( p_k_part ) :: index_r ! recombination cell index as a real number
  integer :: rc_cell_index ! recombination cell index
  integer :: dim ! x, p index iterator
  real( p_k_part ) :: q_sum ! The sum of the 2 q's (mostly unnessisary)

  if ( t < species%splitter%recombine_start_time ) return

  split => species%splitter

  if ( .not. if_recombine_local( split, recombine, n ) ) return

  num_split = ubound(split%splitting_gammas,1)
  recombine%cells = 0
  idx = 0
  particle_loop: do
    idx = idx + 1 ! This is safer here than at the cycle statements (as I stupidly had been doing)
    if ( idx .gt. species%num_par) exit particle_loop

    ! Find particle position, and cycle if the particle isn't in the recombination box
    call species % get_position( idx, x )
    do dim = 1, p_x_dim
      ! exclude points on the boundary to avoid rounding error then say just under zero below
      ! probably not possible, equality is equality, but
      if ( (x(dim) .le. recombine%lr_box(dim, p_lower)) .or. &
           (x(dim) .ge. recombine%lr_box(dim, p_upper)) ) then
        cycle particle_loop
      endif
    enddo

    ! find the integer log_2 of the inverse particle size
    ! (the 'splitting state' )
    ! ( `-1' because '1' is the base state in split_species(), this should be refactored)
    ! FIX THIS
    s_s_l2 = floor_log_2( split%base_particle_q/species%q(idx) + 0.5 ) - 1

    if (s_s_l2>num_split) then
      s_s_l2 = num_split
    endif

    ! Check if particle is even split
    if ( s_s_l2 .le. 0 ) cycle particle_loop

    ! check if particle is cold enough to recombine
    gamma = sqrt( 1.0_p_k_part + species%p(1,idx)**2 + species%p(2,idx)**2 + species%p(3,idx)**2)
    ! This hardcodes a factor of 4 hysteresis loop with splitting
    ! The factor should probably be moved to the input deck
    if ( (1.+(gamma-1.)*4.) .gt. split%splitting_gammas(s_s_l2) ) cycle particle_loop

    ! Find the index for the recombination cell the particle is in
    index_r = x(1) - recombine%lr_box(1, p_lower)
    index_r = index_r/recombine%recombination_cell(1)
    rc_cell_index = min(ceiling(index_r), recombine%l_grid_size(1))
    ! Using ceiling here avoids adding 1 for fortan indexing overall
    ! min() is then in case it's on the edge and there's a rounding issue (I believe that's excluded from earlier but)

    fac_grid_size = 1
    do i = 2, p_x_dim

      fac_grid_size = fac_grid_size * recombine%l_grid_size(i-1)

      index_r = ( x(i) - recombine%lr_box(i,p_lower) ) / recombine%recombination_cell(i)

      rc_cell_index = rc_cell_index + &
        min( floor(index_r), recombine%l_grid_size(i)-1 ) * fac_grid_size

    enddo

    if (rc_cell_index .gt. recombine%cell_count) then
      write(0,"(A,I0,A,I0,A,I0)") "[", mpi_node(), "] cell count ", rc_cell_index, &
                                  "is larger than maximum, ", recombine%cell_count
      write(0,"(A,I0,A,I0)") "[", mpi_node(), "] particle index: ", idx
      write(0,*) "[", mpi_node(), "] x(:) = ", x
      write(0,*) "[", mpi_node(), "] lr_box(:,:) = ", recombine%lr_box
      call abort_program( p_err_invalid )
    endif

    if (rc_cell_index .lt. 1) then
      write(0,"(A,I0,A,I0,A)") "[", mpi_node(), "] cell count ", rc_cell_index, &
                               " is less than minimum, 1"
      write(0,"(A,I0,A,I0)") "[", mpi_node(), "] particle index: ", idx
      write(0,*) "[", mpi_node(), "] x(:) = ", x
      write(0,*) "[", mpi_node(), "] lr_box(:,:) = ", recombine%lr_box
      call abort_program( p_err_invalid )
    endif

    or_idx = recombine%cells(rc_cell_index) !Just a nickname


    ! Do various things depending on who (if anyone) is already in the cell
    if ( or_idx .eq. 0 ) then
      ! If cell is currently empty, store particle array index and cycle
      recombine%cells(rc_cell_index) = idx
      cycle particle_loop
    else
      ! else find charge state of stored particle, then decide what to do
      ! again here, -1 because of poor definitions here and in splitting
      or_s_s_l2 = floor_log_2( split%base_particle_q/species%q( or_idx ) + 0.5 ) - 1
      if (or_s_s_l2>num_split) then
        or_s_s_l2 = num_split
      endif
      if ( or_s_s_l2 .eq.  s_s_l2 ) then
        ! Only combine particles in the same octant in momentum space (in the lab frame)
        ! This is a very poor-man's way of sectioning phase space (both the method and the hard-cording)
        ! but the right way to do this is to re-write this all as a sort rather than a rejection method,
        ! so until(?) I do that I'm going to just leave this as is
        if ( sign(1.0_p_k_part, species%p(1,idx)) .ne. sign(1.0_p_k_part, species%p(1,or_idx)) .or. &
             sign(1.0_p_k_part, species%p(2,idx)) .ne. sign(1.0_p_k_part, species%p(2,or_idx)) .or. &
             sign(1.0_p_k_part, species%p(3,idx)) .ne. sign(1.0_p_k_part, species%p(3,or_idx)) ) then
          cycle particle_loop
        endif

        ! Particles match, recombine!

        ! calculate original energy in lab frame (this was for optional diagnostic)
        !or_gamma = sqrt( 1.0_p_double + &
        !               species%p(1,or_idx)**2 + species%p(2,or_idx)**2 + species%p(3,or_idx)**2)
        !E_i = ( gamma + or_gamma - 2.0 ) * species%q(idx) * species%rqm

        !Update momentum
        q_sum = species%q(or_idx) + species%q(idx)
        do dim = 1, p_p_dim
          ! TODO
          species%p(dim, or_idx) = species%p(dim, or_idx)*( species%q(or_idx) / q_sum ) + &
                                   species%p(dim, idx)   *( species%q(idx)    / q_sum )
        enddo

        ! Enbiggen charge (and hence mass)
        species%q(or_idx) = species%q(or_idx) + species%q(idx)

        ! Calculate new lab frame energy (this was for optional diagnostic)
        !gamma = sqrt( 1.0_p_double + &
        !              species%p(1,or_idx)**2 + species%p(2,or_idx)**2 + species%p(3,or_idx)**2)
        !E_f = ( gamma - 1.0 ) * species%q(or_idx) * species%rqm

        ! Move particle from the end of the particle array to the 'hole' currently occupied by the 'ghost' electron
        if (idx .lt. species%num_par) then
          do dim = 1, species%get_n_x_dims()
            species%x(dim, idx) = species%x(dim, species%num_par)
          enddo
          do dim = 1, p_x_dim
            species%ix(dim, idx) = species%ix(dim, species%num_par)
          enddo
          do dim = 1, p_p_dim
            species%p(dim, idx) = species%p(dim, species%num_par)
          enddo
          species%q(idx) = species%q( species%num_par )

#ifdef __HAS_SPIN__

          do dim = 1, p_s_dim
            species%s(dim, idx) = species%s(dim, species%num_par)
          enddo

#endif

          ! Zap tracks for annihilated particle and moved particle, then recreate tracks for moved particle
#ifdef __HAS_TRACKS__

          if ( species%diag%ndump_fac_tracks > 0 ) then
            call missing_particles( species%diag%tracks, (/idx,species%num_par/), 2, species%tag )
          endif

#endif

          if ( species%add_tag ) then

            species%tag(:,idx) = species%tag( :, species%num_par )
#ifdef __HAS_TRACKS__
            ! tell tracks that the particle index has moved
            ! tracks doesn't do this test itself, I have to do it by hand at every entry point?
            ! That seems dumb, I should probably fix that.
            if ( (species%diag%ndump_fac_tracks > 0) .and. (species%tag(1,idx) < 0) )  then
              call inbound_particle( species%diag%tracks, species%tag(:,idx), idx )
            endif

#endif

          endif

        else

#ifdef __HAS_TRACKS__
          if ( species%diag%ndump_fac_tracks > 0 ) then
            call missing_particles( species%diag%tracks, (/species%num_par/), 1, species%tag )
          endif
#endif

        endif

        !An electron dies not with a bang, but with a variable increment
        species%num_par = species%num_par - 1

        ! clear the recombination cell
        recombine%cells(rc_cell_index) = 0

        ! counter the idx iterate because we want to check the new particle as well
        idx = idx - 1
        cycle particle_loop

      elseif (  s_s_l2  .gt. or_s_s_l2 ) then
        ! Swap in smaller particle
        recombine%cells(rc_cell_index) = idx
        cycle particle_loop

      else
        ! leave particle in cell, drop current particle, and cycle
        cycle particle_loop

      endif ! or_s_s_l2 .eq.  s_s_l2

    endif ! split%cells(rc_cell_index) .eq. 0

    ! We should never get here, so I'll enforce that
    ERROR("Something went wrong in particle recombinations.")
    ERROR("Hopefully Josh was still debugging them.")
    call abort_program( p_err_invalid )

  enddo particle_loop

end subroutine recombine_species
!-----------------------------------------------------------------------------------------

end module m_spec_overdense_cyl
