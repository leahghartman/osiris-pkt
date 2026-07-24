#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_gr

#include "memory/memory.h"

use m_system
use m_parameters

use m_geometry_gr
use m_emf_define_gr

interface sort
  module procedure sort_particles_gr
end interface

contains

!-----------------------------------------------------------------------------------------
! sorts all particles
!-----------------------------------------------------------------------------------------
subroutine sort_particles_gr( this, n, t )
!-----------------------------------------------------------------------------------------
  
  use m_system
  use m_particles_define, only : t_particles
  use m_species_define, only : t_species
  use m_species, only : sort

  implicit none

  type( t_particles ), intent(inout) :: this
  integer, intent(in) :: n
  real(p_double), intent(in) :: t

  class( t_species ), pointer :: species

  ! executable statements
  species => this % species
  do
      if (.not. associated(species)) exit
      call sort( species, n, t )
      species => species % next
  enddo


end subroutine sort_particles_gr

end module m_particles_gr

!-------------------------------------------------------------------------------
subroutine read_input_particles_gr( this, input_file, periodic, &
    if_move, grid, dt, sim_options )

  use m_system
  use m_parameters

  use m_input_file, only : get_namelist, disp_out
  use m_particles_define, only : p_report_quants
  use m_species_define, only : t_species, p_max_spname_len, p_cell_low, p_cell_near
  use m_input_file, only : t_input_file
  use m_grid_define, only : t_grid
  use m_vdf_define, only : p_max_reports_len, p_max_reports, p_n_report_type, &
    p_full, p_savg, p_senv, p_line, p_slice
  use m_vdf_report, only : new
  use m_diagnostic_utilities, only : p_diag_prec
  use stringutil
  use m_species_collisions, only : read_nml

  use m_particles_define_gr, only : t_particles_gr

  implicit none

  ! dummy variables
  class( t_particles_gr ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  real(p_double), intent(in) :: dt
  type( t_options ), intent(in) :: sim_options

  class( t_species ), pointer :: species
  integer :: num_species, i
  integer :: num_groups
  integer :: e_id, p_id
  logical :: low_jay_roundoff

  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports
  integer :: ndump_fac, ndump_fac_ave, ndump_fac_lineout, ndump_fac_ene
  integer, dimension(p_x_dim) :: n_ave
  integer :: prec, n_tavg

  character(len=20) :: interpolation
  logical :: grid_center

  namelist /nl_particles/ num_species, num_groups, low_jay_roundoff, &
  ndump_fac, ndump_fac_ave, ndump_fac_lineout, n_ave, prec, &
  reports, n_tavg, interpolation, grid_center, ndump_fac_ene

  integer, dimension( p_n_report_type ) :: ndump_fac_all
  integer :: ierr

  character(len = p_max_spname_len) :: spname

  ! executable statements
  num_species = 0
  num_groups = 0

  ! Diagnostics
  reports = "-"
  ndump_fac = 0
  ndump_fac_ave = 0
  ndump_fac_lineout = 0
  ndump_fac_ene = 0

  n_ave = -1
  n_tavg = -1
  prec = p_diag_prec

  low_jay_roundoff = .false.
  interpolation = "linear"
  grid_center = .false.

  ! get namelist text from input file
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

  case ( "linear" )
    this%interpolation = p_linear

  case default
    print *, '   Error reading species parameters'
    print *, '   invalid interpolation: "', trim( interpolation ), '"'
    print *, '   Only "linear" interpolation is available.'
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

  this%num_species = num_species + 3*num_groups
  this%num_groups  = num_groups

  ! allocate sub-objects
  call this % allocate_objs()
  
  ! read the input file for standard species
  species => this % species

  if ( num_species>0 ) then
    do i=1, num_species
    write(spname, '(A,I0)') 'species ',i
    if ( mpi_node() == 0 .and. disp_out(input_file) ) then
      print '(A,I0,A)', " - species (",i,") configuration..."
    endif
    call species % read_input( input_file, spname, periodic, if_move, grid, &
      dt, .true., sim_options )
    
    species => species % next
    enddo
  endif    
! 
  ! read the input file for species groups
  if ( this%num_groups > 0 ) then
    SCR_ROOT("Species groups:")
    SCR_ROOT("- Number of groups: ", trim(tostring(this%num_groups)))

    do i=1, this%num_groups
      SCR_ROOT(" - Species group (",trim(tostring(i)),"):")
      call this % spec_group(i) % read_input( input_file, grid )

      SCR_ROOT("  - electrons configuration...")
      e_id = 3*num_species + 3*(i-1) + 1
      write(spname, '(A,I0)') 'species ', e_id
      call species % read_input( input_file, spname, periodic, if_move, grid, &
        dt, .true., sim_options )
      species => species % next

      SCR_ROOT("  - positrons/ions configuration...")
      p_id = e_id + 1
      write(spname, '(A,I0)') 'species ', p_id
      call species % read_input( input_file , spname, periodic, if_move, grid, &
        dt, .true., sim_options )
      species => species % next

      SCR_ROOT("  - photons configuration...")
      call species % read_input( input_file, spname, periodic, if_move, grid, &
        dt, .false., sim_options )
      species => species % next

    enddo
  endif


  ! validate species names
  call this % validate_names()

end subroutine read_input_particles_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine init_particles_gr( this, g_space, jay, emf, grid, no_co, bnd, zpulse_list, &
                              ndump_fac, restart, restart_handle, t, tstep, tmin, tmax, &
                              sim_options )

  use m_system
  use m_parameters
  use m_particles_define_gr, only : t_particles_gr
  use m_space, only : t_space
  use m_emf_define, only : t_emf
  use m_current_define, only : t_current
  use m_grid_define, only : t_grid
  use m_node_conf, only : t_node_conf
  use m_bnd, only : t_bnd
  use m_restart, only : t_restart_handle
  use m_vdf_define, only : t_vdf_report
  use m_species_define, only : t_species, p_cell_low, p_cell_near
  use m_particles_define
  use m_particles, only : restart_read
  use m_particles_charge, only : setup
  use m_species
  use m_logprof, only : create_event
  use m_time_step, only : t_time_step
  use m_zpulse, only : t_zpulse_list

  use m_species_define_gr, only : t_species_gr

  implicit none

  ! dummy variables
  class( t_particles_gr ), intent(inout) :: this
  type( t_space ),     intent(in) :: g_space
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

  ! local variables
  type( t_vdf_report ), pointer :: report
  integer, parameter :: izero = iachar('0')
  logical :: keep_previous_charge

  class(t_species), pointer :: species
  integer :: i

  ! check that restart values match input deck
  if ( restart ) then
    call this % restart_read( restart_handle )
  else
    this % n_current = 0
  endif

  ! setup diagnostics
  keep_previous_charge = .false.

  report => this%reports
  do
    if ( .not. associated( report ) ) exit

    report%xname  = (/'x1', 'x2', 'x3'/)
    report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    ! this are just dummy values for now
    report%time_units = '1 / \omega_p'
    report%dt         = 1.0

    report%fileLabel = ''
    report%basePath  = trim(path_mass) // 'FLD' // p_dir_sep

    select case ( report%quant )
    case ( p_charge )
      report%label = '\rho'
      if ( p_x_dim == 1 ) then
        report%units  = 'e \omega_p / c'
      else
        report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
        '/ c^'//char(izero+p_x_dim)
      endif
      
    case ( p_charge_htc )
      report%label = '\rho_{n-1/2}'
      if ( p_x_dim == 1 ) then
        report%units  = 'e \omega_p / c'
      else
        report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
        '/ c^'//char(izero+p_x_dim)
      endif
      keep_previous_charge = .true.
      
    case ( p_dcharge_dt )
      report%label = 'd\rho/dt'
      if ( p_x_dim == 1 ) then
        report%units  = 'e \omega_p / c'
      else
        report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
        '/ c^'//char(izero+p_x_dim)
      endif
      keep_previous_charge = .true.
      
    end select

    report => report%next
  enddo

  ! Setup charge
  call setup( this%charge, keep_previous_charge, restart, restart_handle )
  ! setup species
  if ( this%num_species > 0 ) then
    
    ! Initialize buffers for communications and current deposition
    call init_buffers_spec( )
    
    species => this % species; i = 1
    do
      if (.not. associated(species)) exit

      select type( s => species )
      class is ( t_species_gr )
        s%geometry => this%geometry
      class default
        ERROR('pointer to geometry must be set with a t_species_geometry object')
        call abort_program( p_err_invalid )
      end select

      call species%init( i, this%interpolation, this%grid_center, grid, &
      g_space, emf, jay, no_co, bnd%send_vdf, bnd%recv_vdf, bnd%bnd_cross, bnd%node_cross, &
      bnd%send_spec, bnd%recv_spec, ndump_fac, &
      restart, restart_handle, sim_options, tstep, tmin, tmax )
      species => species % next; i = i+1
    enddo
    
    ! setup additional buffers for low roundoff current deposition
    if (( this%num_species > 1 ) .and. (this%low_jay_roundoff)) then
      call this % jay_tmp % new( jay % pf(1) )
    else
      this%low_jay_roundoff = .false.
    endif
  endif

  if (this%num_groups > 0) then
    do i = 1, this%num_groups
      call this % spec_group(i) % init( this%geometry, grid )
    enddo
  endif

  call this % copy_zpulse_to_spec( zpulse_list )

  ! setup events in this file
  if (pushev==0) then
    pushev            = create_event('advance deposit')
    reduce_current_ev = create_event('reduce current')
    
    partboundev       = create_event('update particle boundary')
    diag_part_ev      = create_event('particle diagnostics')
    
    ! setup events in the species class
    sortev            = create_event('particle sort (total)')
    sort_genidx_ev    = create_event('particle sort, gen. idx')
    sort_rearrange_ev = create_event('particle sort, rearrange particles')
  endif

end subroutine init_particles_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine advance_deposit_particles_gr( this, emf, jay, tstep, t, &
  no_co, options )

  use m_system
  use m_parameters

  use m_particles_define_gr, only : t_particles_gr
  use m_emf_define, only : t_emf
  use m_current_define, only : t_current
  use m_node_conf, only : t_node_conf
  use m_time_step, only : t_time_step, dt

  class( t_particles_gr ), intent(inout) :: this
  class( t_emf ), intent( inout )  ::  emf
  class( t_current ), intent(inout) :: jay

  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t
  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  integer :: i

  call this%t_particles%advance_deposit(emf, jay, tstep, t, &
    no_co, options )

  ! create pairs and advance photons
  do i = 1, this%num_groups
    ! call this % spec_group(i) % photons %  push ( emf, jay % pf(1), t, tstep, 0, 1, options )
      
    ! inject particles
    call this % spec_group(i) % inject ( emf%b, emf%e, dt( tstep ), no_co )
  enddo

end subroutine advance_deposit_particles_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine report_particles_gr( this, emf, g_space, grid, no_co, &
  tstep, t, tmin, send_msg, recv_msg )
    
  use m_system
  use m_parameters
  use m_species_define
  use m_particles_define_gr
  use m_particles_charge
  use m_particles
  use m_emf_define
  use m_space
  use m_grid_define
  use m_node_conf
  use m_time_step
  use m_vdf_define
  use m_vdf_report
  use m_vdf_comm, only : t_vdf_msg, update_boundary
  use m_particles_define, only: p_charge
  use m_logprof
  use m_neutral
  use m_species_group_gr
  
  implicit none
  
  class( t_particles_gr ), intent(inout) :: this
  class( t_emf ),       intent(inout) :: emf
  
  type( t_space ),      intent(in) :: g_space
  class( t_grid ),      intent(in) :: grid
  class( t_node_conf ),  intent(in) :: no_co
  type( t_time_step ),  intent(in) :: tstep
  real(p_double),       intent(in) :: t, tmin
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, i1, i2
  integer, parameter :: p_part_block = 4096
  real(p_k_part), dimension(p_part_block) :: q

  call this % t_particles % report( emf, g_space, grid, no_co, tstep, t, tmin, &
    send_msg, recv_msg )

  do i = 1, this%num_groups
    call this % spec_group(i) % report( g_space, grid, no_co, tstep, t )

    if ( this % spec_group(i) % inj_limit == p_inj_limit_fracnGJpole .or. &
         this % spec_group(i) % inj_limit == p_inj_limit_fracnGJpolealt .or. &
         this % spec_group(i) % inj_limit == p_inj_limit_fracnGJlocal .or. &
         this % spec_group(i) % inj_limit == p_inj_limit_fracnGJlocalalt .or. &
         this % spec_group(i) % inj_limit == p_inj_limit_fracEBteo .or. &
         this % spec_group(i) % inj_limit == p_inj_limit_fracEBteolim .or. &
         this % spec_group(i) % inj_limit == p_inj_limit_sigma ) then
          
      ! deposit
      do i1 = 1, this % spec_group(i) % electrons % num_par, p_part_block
        i2 = i1 + p_part_block - 1
        if ( i2 > this % spec_group(i) % electrons % num_par ) i2 = this % spec_group(i) % electrons % num_par
        call this % spec_group(i) % electrons % get_quant( i1, i2, p_charge, q )
        call this % spec_group(i) % electrons % deposit_density( this % spec_group(i) % electrons % charge, i1, i2, q  )
      enddo
      do i1 = 1, this % spec_group(i) % positrons % num_par, p_part_block
        i2 = i1 + p_part_block - 1
        if ( i2 > this % spec_group(i) % positrons%num_par ) i2 = this % spec_group(i) % positrons % num_par
        call this % spec_group(i) % positrons % get_quant( i1, i2, p_charge, q )
        call this % spec_group(i) % positrons % deposit_density( this % spec_group(i) % positrons % charge, i1, i2, q  )
      enddo

      ! update boundaries
      call update_boundary(this % spec_group(i) % electrons % charge, p_vdf_add, no_co, send_msg, recv_msg )
      call update_boundary(this % spec_group(i) % positrons % charge, p_vdf_add, no_co, send_msg, recv_msg )
    endif

  enddo

end subroutine report_particles_gr
!-------------------------------------------------------------------------------