#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_rad

#include "memory/memory.h"

use m_system
use m_parameters
use m_species_radiat_define
use m_particles_define, only : t_particles, p_report_quants
use m_input_file, only : t_input_file, get_namelist
use m_vdf_define, only : p_max_reports_len, p_max_reports, p_n_report_type, &
                         p_full, p_savg, p_senv, p_line, p_slice
use m_species_define, only : t_species, p_cell_near, p_cell_low, p_max_spname_len
use m_diagnostic_utilities,  only : p_diag_prec
use m_vdf_report, only : new
use m_species_collisions, only : read_nml
use m_grid_define, only : t_grid
use m_logprof, only : create_event, begin_event, end_event
use m_space, only : t_space
use m_emf_define, only : t_emf
use m_current_define, only : t_current
use m_node_conf, only : t_node_conf
use m_bnd, only : t_bnd
use m_restart, only : t_restart_handle
use m_time_step, only : t_time_step
use m_zpulse, only: t_zpulse_list

implicit none

private

type, extends( t_particles ) :: t_particles_rad

  ! radiat species groups
  integer :: num_radiat = 0
  integer :: num_neutral_rad = 0
  integer :: num_neutral_mov_ions_rad = 0
  integer :: num_neutral_mov_ions = 0

contains

  ! Allocate species_radiat objects instead
  procedure :: allocate_objs => allocate_objs_part_rad
  procedure :: read_input => read_input_part_rad
  procedure :: init => init_part_rad
  procedure :: advance_deposit => advance_deposit_part_rad

end type t_particles_rad

integer :: rad_calc_ev

public :: t_particles_rad, rad_calc_ev

contains

!-----------------------------------------------------------------------------------------
subroutine read_input_part_rad( this, input_file, periodic, if_move, grid, dt, sim_options )

  implicit none

  class( t_particles_rad ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ),             intent(in) :: grid
  real(p_double), intent(in) :: dt
  type( t_options ), intent(in) :: sim_options

  integer  :: num_species, i
  class( t_species ), pointer :: species
  integer  :: num_cathode
  integer  :: num_neutral, num_neutral_mov_ions
  logical  :: low_jay_roundoff

  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports
  integer  :: ndump_fac, ndump_fac_ave, ndump_fac_lineout, ndump_fac_ene
  integer, dimension(p_x_dim) :: n_ave
  integer                     :: prec, n_tavg

  character(len=20) :: interpolation
  logical :: grid_center

  integer :: num_radiat, num_neutral_rad, num_neutral_mov_ions_rad

  namelist /nl_particles/ num_species, num_cathode, num_neutral, &
  num_neutral_mov_ions, low_jay_roundoff, &
  ndump_fac, ndump_fac_ave, ndump_fac_lineout, n_ave, prec, &
  reports, n_tavg, interpolation, grid_center, ndump_fac_ene, &
  num_radiat, num_neutral_rad, num_neutral_mov_ions_rad

  integer, dimension( p_n_report_type ) :: ndump_fac_all
  integer :: ierr

  character(len = p_max_spname_len) :: spname

  num_species = 0
  num_cathode = 0
  num_neutral = 0
  num_neutral_mov_ions = 0
  num_radiat=0
  num_neutral_rad = 0
  num_neutral_mov_ions_rad = 0

  low_jay_roundoff = .false.

  ! Diagnostics
  reports = "-"
  ndump_fac = 0
  ndump_fac_ave = 0
  ndump_fac_lineout = 0
  ndump_fac_ene = 0

  n_ave = -1
  n_tavg = -1
  prec = p_diag_prec

  interpolation = "quadratic"
  grid_center = .false.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_particles", ierr )

  if ( ierr /= 0 ) then
    if (ierr < 0) then
      print *, "Error reading particles parameters"
    else
      print *, "Error: particles parameters missing"
    endif
    print *, "aborting..."
    stop
  endif

  read (input_file%nml_text, nml = nl_particles, iostat = ierr)
  if (ierr /= 0) then
    print *, "Error reading particles parameters"
    print *, "aborting..."
    stop
  endif

  ! process interpolation scheme
  select case ( trim( interpolation ))
    ! case ( "ngp" ) ! not implemented
    !  this%interpolation = p_ngp
  case ( "linear" )
    this%interpolation = p_linear
  case ( "quadratic" )
    this%interpolation = p_quadratic
  case ( "cubic" )
    this%interpolation = p_cubic
  case ( "quartic" )
    this%interpolation = p_quartic

  case default
    print *, '   Error reading species parameters'
    print *, '   invalid interpolation: "', trim( interpolation ), '"'
    print *, '   Valid values are "linear", "quadratic", "cubic" and "quartic".'
    print *, '   aborting...'
    stop
  end select

  this%grid_center = grid_center

  ! process position type
  select case ( this%interpolation )
  case ( p_linear, p_cubic )
    this%pos_type = p_cell_low
  case ( p_quadratic, p_quartic )
    this%pos_type = p_cell_near
  end select

  this%low_jay_roundoff = low_jay_roundoff


  ! Particle kinetic energy diagnostics
  this%ndump_fac_ene = ndump_fac_ene

  ! global charge diagnostics
  ndump_fac_all(p_full)  = ndump_fac
  ndump_fac_all(p_savg)  = ndump_fac_ave
  ndump_fac_all(p_senv)  = ndump_fac_ave
  ndump_fac_all(p_line)  = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  ! process normal reports
  call new( this%reports, reports, p_report_quants, &
    ndump_fac_all, n_ave, n_tavg, prec, &
    p_x_dim, ierr )
  if ( ierr /= 0 ) then
    print *, "(*error*) Invalid report"
    print *, "(*error*) aborting..."
    stop
  endif


  ! Allocate sub-objects; species are ordered as listed here
  this%num_species = num_radiat + num_species + num_cathode + num_neutral_rad + &
                     num_neutral + 2*num_neutral_mov_ions_rad + 2*num_neutral_mov_ions
  this%num_cathode = num_cathode
  this%num_neutral = num_neutral_rad + num_neutral + num_neutral_mov_ions_rad + &
                     num_neutral_mov_ions

  this%num_radiat = num_radiat
  this%num_neutral_rad = num_neutral_rad
  this%num_neutral_mov_ions_rad = num_neutral_mov_ions_rad
  this%num_neutral_mov_ions = num_neutral_mov_ions_rad + num_neutral_mov_ions

  call this % allocate_objs()

  species => this % species

  ! Read the input files for standard species
  if ( num_species + num_radiat > 0 ) then

    do i=1, num_species + num_radiat
      write(spname, '(A,I0)') 'species ',i
      if ( mpi_node() == 0 ) then
        print '(A,I0,A)', " - species (",i,") configuration..."
      endif
      call species % read_input( input_file, spname, periodic, if_move, grid, &
        dt, .true.,  sim_options )

      species => species % next
    enddo
  endif


  ! Take care of cathodes
  if ( num_cathode > 0 ) then

    do i=1, num_cathode
      write(spname, '(A,I0)') 'cathode ',i
      if ( mpi_node() == 0 ) then
        print '(A,I0,A)', " - cathode (",i,") configuration..."
      endif

      ! read the cathode configuration, and the associated species
      ! configuration
      call this % cathode(i) % read_input( input_file, species, spname,  periodic, &
                                           if_move, grid, dt, sim_options )

      species => species % next

    end do

  endif

#ifdef __HAS_IONIZATION__

  ! Take care of neutrals

  if ( num_neutral + num_neutral_rad > 0 ) then

    if (sim_options%omega_p0 <= 0.0) then
      print *, "(*error*) When using ionization (num_neutral > 0) omega_p0 or n0 must be"
      print *,"(*error*) set in the simulation section at the beggining of the input file"
      print *, "(*error*) bailing out..."
      call abort_program()
    endif

    do i=1, num_neutral + num_neutral_rad
      write(spname, '(A,I0)') 'neutral ',i

      if ( mpi_node() == 0 ) then
        print '(A,I0,A)', " - neutral (",i,") configuration..."
      endif

      call this%neutral( i ) % read_input( input_file, species, spname,  .false., &
                                           periodic, if_move, grid, dt, &
                                           sim_options )

      species => species % next
    end do


  endif

  ! Take care of neutrals with moving ions

  if ( num_neutral_mov_ions + num_neutral_mov_ions_rad > 0 ) then

    if (sim_options%omega_p0 <= 0.0) then
      print *, "(*error*) When using ionization (num_neutral_mov_ions > 0) omega_p0 or n0"
      print *, "(*error*) must be set in the simulation section at the beggining of the"
      print *, "(*error*) input file bailing out..."
      call abort_program()
    endif

    do i=1, num_neutral_mov_ions + num_neutral_mov_ions_rad

      write(spname, '(A,I0)') "neutral_mov_ions ", i

      SCR_ROOT(" - neutral with moving ions (",i,") configuration...")

      call this%neutral( num_neutral + i )%read_input( input_file, species, spname, &
                                                       .true., periodic, if_move, grid, &
                                                       dt, sim_options )

      species => species % next

    end do

  endif

#else

  if (( num_neutral + num_neutral_rad > 0 ) .or. &
      ( num_neutral_mov_ions + num_neutral_mov_ions_rad > 0 )) then
    print *, "(*error*) Ionization is not supported in this version"
    stop
  endif

#endif

  ! validate species names
  call this % validate_names()

#ifdef __HAS_COLLISIONS__

  ! Read collision data
  call read_nml( this%coll, input_file, this%species, this%num_species, sim_options )

#endif

end subroutine read_input_part_rad
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_part_rad( this, g_space, jay, emf, grid, no_co, bnd, zpulse_list, &
                          ndump_fac, restart, restart_handle, t, tstep, tmin, tmax, &
                          sim_options)

  implicit none

  class( t_particles_rad ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  class( t_emf) ,intent(inout) :: emf
  class( t_current ) ,intent(inout) :: jay
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  class(t_bnd), intent(inout) :: bnd
  class( t_zpulse_list ), intent(in) :: zpulse_list
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  real(p_double), intent(in) :: t
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax
  type( t_options ), intent(in) :: sim_options

  call this % t_particles % init( g_space, jay, emf, grid, no_co, bnd, zpulse_list, &
                                  ndump_fac, restart, restart_handle, t, tstep, tmin, &
                                  tmax, sim_options )

  ! setup events in this file
  if ( rad_calc_ev == 0 ) then
    rad_calc_ev = create_event('RaDiO calculation')
  endif

end subroutine init_part_rad
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Overrides allocate_objs method and allocates t_species_qed objects instead
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_rad( this )

  implicit none

  class( t_particles_rad ), intent(inout) :: this
  integer :: i, n_spec, i0_rad, i1_rad, i0_neutral_rad, i1_neutral_rad
  integer :: i0_neutral_mov_ions_rad, i1_neutral_mov_ions_rad
  class(t_species), pointer :: spec, tail

  ! number of vanilla species, including num_radiat and regular species
  n_spec = this%num_species - (this%num_cathode + this%num_neutral + &
                               this%num_neutral_mov_ions)
  i0_rad = 1
  i1_rad = this%num_radiat
  i0_neutral_rad = n_spec + this%num_cathode + 1
  i1_neutral_rad = n_spec + this%num_cathode + this%num_neutral_rad
  i0_neutral_mov_ions_rad = n_spec + this%num_cathode + this%num_neutral &
                            - this%num_neutral_mov_ions + 1
  i1_neutral_mov_ions_rad = n_spec + this%num_cathode + this%num_neutral &
                            - this%num_neutral_mov_ions + this%num_neutral_mov_ions_rad

  if ( this%num_species > 0 ) then

    ! Check if first species should be t_species_radiat
    if ( (i0_rad==1 .and. i1_rad>=i0_rad) .or. &
         (i0_neutral_rad==1 .and. i1_neutral_rad>=i0_neutral_rad) .or. &
         (i0_neutral_mov_ions_rad==1 .and. &
          i1_neutral_mov_ions_rad>=i0_neutral_mov_ions_rad) ) then
      allocate(t_species_radiat :: spec)
    else
      allocate(spec)
    endif

    this % species => spec
    do i = 2, this%num_species
      tail => spec
      if( (i >= i0_rad .and. i <= i1_rad) .or. &
          (i >= i0_neutral_rad .and. i <= i1_neutral_rad) .or. &
          (i >= i0_neutral_mov_ions_rad .and. i <= i1_neutral_mov_ions_rad) ) then
        allocate( t_species_radiat :: spec )
      else
        allocate(spec)
      endif
      tail % next => spec
    enddo

    if ( this%num_cathode > 0 ) then
      allocate( this%cathode( this%num_cathode ) )
    endif

#ifdef __HAS_IONIZATION__
    if ( this%num_neutral > 0 ) then
      allocate( this%neutral( this%num_neutral ) )
    endif
#endif

  endif

end subroutine allocate_objs_part_rad
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Advance/deposit particles using multiple tasks per node (OpenMP)
!  - Low jay roundoff is currently not supported
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_part_rad( this, emf, jay, tstep, t, no_co, options )

  implicit none

  class( t_particles_rad ), intent(inout) :: this
  class( t_emf ), intent( inout )  ::  emf
  class( t_current ), intent(inout) :: jay
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  class(t_species), pointer :: species

  call this % t_particles % advance_deposit( emf, jay, tstep, t, no_co, options )

  call begin_event( rad_calc_ev )
  species => this % species
  do
    if (.not. associated(species)) exit
    select type(species)
    class is( t_species_radiat )
      call species % rad_calc( t, tstep, no_co )
    end select
    species => species % next
  enddo
  call end_event( rad_calc_ev )

end subroutine advance_deposit_part_rad
!-----------------------------------------------------------------------------------------

end module m_particles_rad
