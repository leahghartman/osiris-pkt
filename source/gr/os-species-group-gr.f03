#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_group_gr

#include "memory/memory.h"

use m_system
use m_parameters

use m_geometry_gr
use m_species_define_gr
use m_photons_define_gr
use m_math, only : pi, pi_2

implicit none

private

integer, parameter, public :: p_cooling_none         = 1, &
                              p_cooling_rr           = 2, &
                              p_cooling_qed_red      = 3, &
                              p_cooling_qed_red_soft = 4, &
                              p_cooling_qed          = 5

! type of supported boundary conditions
integer, parameter, public :: p_inj_local_none   = 1, &
                              p_inj_local_surf   = 2, &
                              p_inj_local_vol    = 3

integer, parameter, public :: p_inj_type_none        = 0, &
                              p_inj_type_const       = 1, &
                              p_inj_type_fracnGJ     = 2, &
                              p_inj_type_fracnGJpole = 3, &
                              p_inj_type_fracEB      = 4, &
                              p_inj_type_fracEBteo   = 5

integer, parameter, public :: p_inj_limit_none            = 0, &
                              p_inj_limit_fracnGJpole     = 1, &
                              p_inj_limit_fracnGJpolealt  = 2, &
                              p_inj_limit_fracnGJlocal    = 3, &
                              p_inj_limit_fracnGJlocalalt = 4, &
                              p_inj_limit_fracEB          = 5, &
                              p_inj_limit_fracEBteo       = 6, &
                              p_inj_limit_fracEBteolim    = 7, &
                              p_inj_limit_sigma           = 8

type :: t_diag_group

  integer :: ndump_fac_nemit

contains

  procedure :: read_input => read_input_diag_group
  procedure :: init       => init_diag_group
  
end type t_diag_group

public t_species_group

type :: t_species_group

  class( t_geometry_gr ), pointer :: geometry

  integer :: cooling_type

  ! particle injection variables
  integer           :: inj_local, inj_type
  real(p_k_fld)     :: ks
  real(p_k_fld)     :: vs

  integer           :: inj_limit
  real(p_k_fld)     :: k_limit
  real(p_k_fld)     :: sig_limit

  ! auxiliary flags
  logical :: star_bnd

  ! pointers to group species
  class( t_species_gr ), pointer :: electrons
  class( t_species_gr ), pointer :: positrons
  class( t_photons_gr ), pointer :: photons

  class( t_diag_group ), pointer  :: diag  => null()

contains
  
  procedure :: allocate_objs   => allocate_objs_group
  procedure :: read_input      => read_input_group
  procedure :: init            => init_group
  procedure :: set_species     => set_species_group
  procedure :: list_algorithm  => list_algorithm_group
  procedure :: inject          => inject_group
  procedure :: report          => report_group
  procedure :: cleanup         => cleanup_group

end type t_species_group

contains

!-------------------------------------------------------------------------------
! Read information from input file
!-------------------------------------------------------------------------------
subroutine read_input_diag_group( this, input_file, group, grid )

  use m_system
  use m_parameters
  use m_input_file
  use m_grid_define, only : t_grid

  implicit none

  class( t_diag_group ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  class( t_species_group ), intent(inout) :: group
  class( t_grid ), intent(in) :: grid
  
  ! local variables
  integer :: ndump_fac_nemit

  integer :: ierr

  namelist /nl_group_diag/ ndump_fac_nemit

  ndump_fac_nemit = 0
  
  ! Get namelist text from input file
  call get_namelist( input_file, "nl_group_diag", ierr )
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading species group diagnostic parameters"
      else
        write(0,*) "Error: species group diagnostic parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif


  read (input_file%nml_text, nml = nl_group_diag, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error reading species group diagnostic parameters"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  if ( .not. ANY( group % cooling_type == (/p_cooling_none, p_cooling_rr/) ) ) then
    this%ndump_fac_nemit = ndump_fac_nemit
  else
    this%ndump_fac_nemit = 0
  endif
  
end subroutine read_input_diag_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine init_diag_group( this, group, grid )

  use m_system
  use m_parameters
  use m_grid_define, only : t_grid

  implicit none

  class( t_diag_group ), intent(inout) :: this
  class( t_species_group ), intent(inout) :: group
  class( t_grid ), intent(in) :: grid
  
  ! local variables

  ! do nothing for now

end subroutine init_diag_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine allocate_objs_group( this )

  use m_system
  use m_parameters
  use m_input_file

  implicit none

  class( t_species_group ), intent(inout) :: this

  ! allocate diagnostic object class
  if ( .not. associated( this%diag ) ) then
    allocate( t_diag_group :: this%diag )
  endif

end subroutine allocate_objs_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read information from input file
!-------------------------------------------------------------------------------
subroutine read_input_group( this, input_file, grid )

  use m_system
  use m_parameters
  use m_input_file
  use m_grid_define, only : t_grid

  implicit none

  class( t_species_group ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  class( t_grid ), intent(in) :: grid
  
  ! local variables
  ! cooling
  character(len = 10 ) :: cooling_type
  real(p_double) :: pp_rmin, pp_rmax, pp_gth, pp_gpair
  ! real(p_double)    :: qed_g_cutoff, qed_p_cutoff, qed_zeta

  ! particle injection
  character(len=20) :: inj_local, inj_type, inj_limit
  real(p_k_fld) :: ks, vs, k_limit, sig_limit

  integer :: ierr

  namelist /nl_group/ inj_local, inj_type, &
    ks, vs, inj_limit, k_limit, sig_limit, &
    cooling_type, pp_gth, pp_gpair, pp_rmin, pp_rmax

    ! qed_g_cutoff, qed_p_cutoff, qed_zeta
  
  ! cooling
  cooling_type = "none"
  pp_rmin    = 1.d0
  pp_rmax    = 10.d0
  pp_gth    = 30.d0
  pp_gpair  = 18.d0

  ! ! QED
  ! qed_g_cutoff = 10.0d0
  ! qed_p_cutoff = 2.0d0
  ! qed_zeta     = 1.0d0

  ! particle injection
  inj_local     = "none"
  inj_type      = "none"
  ks            = 1.0d0
  vs            = 0.0d0
  inj_limit     = "none"
  k_limit       = 1.0d0
  sig_limit     = 1.0d3
  
  ! Get namelist text from input file
  call get_namelist( input_file, "nl_group", ierr )
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading species group parameters"
      else
        write(0,*) "Error: species group parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif


  read (input_file%nml_text, nml = nl_group, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error reading species group parameters"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  ! allocate objects
  call this%allocate_objs()

  ! save particle injection variables
  select case ( trim( inj_local ) )
  case ( "surface" )
    this%inj_local = p_inj_local_surf
  case ( "volume" )
    this%inj_local = p_inj_local_vol
  case default
    this%inj_local = p_inj_local_none
  end select

  select case ( trim( inj_type ) )
  case ( "fracnGJ" )
    this%inj_type = p_inj_type_fracnGJ
  case ( "fracnGJpole" )
    this%inj_type = p_inj_type_fracnGJpole
  case ( "fracEB" )
    this%inj_type = p_inj_type_fracEB
  case ( "fracEBteo" )
    this%inj_type = p_inj_type_fracEBteo
  case ( "const" )
    this%inj_type = p_inj_type_const
  case default
    this%inj_type = p_inj_type_none
  end select
  
  this%ks = ks
  this%vs = vs

  select case ( trim( inj_limit ) )
  case ( "nGJpole" )
    this%inj_limit = p_inj_limit_fracnGJpole
  case ( "altnGJpole" )
    this%inj_limit = p_inj_limit_fracnGJpolealt
  case ( "nGJlocal" )
    this%inj_limit = p_inj_limit_fracnGJlocal
  case ( "altnGJlocal" )
    this%inj_limit = p_inj_limit_fracnGJlocalalt
  case ( "fracEB" )
    this%inj_limit = p_inj_limit_fracEB
  case ( "fracEBteo" )
    this%inj_limit = p_inj_limit_fracEBteo
  case ( "fracEBteolim" )
    this%inj_limit = p_inj_limit_fracEBteolim
  case ( "sigma" )
    this%inj_limit = p_inj_limit_sigma
  case default
    this%inj_limit = p_inj_limit_none
  end select

  this%k_limit   = k_limit
  this%sig_limit = sig_limit

  ! cooling type
  select case ( trim( cooling_type ) )
  case('none')
    this%cooling_type = p_cooling_none
  case('rr')
    this%cooling_type = p_cooling_rr
  case('qed_red')
    this%cooling_type = p_cooling_qed_red
  case('qed_red_soft')
    this%cooling_type = p_cooling_qed_red_soft
  case('qed')
    this%cooling_type = p_cooling_qed
  case default
    this%cooling_type = p_cooling_none
  end select

  this%electrons%cooling_type  = this%cooling_type
  this%positrons%cooling_type  = this%cooling_type

  ! pair production variables
  this%electrons%pp_rmin  = pp_rmin
  this%electrons%pp_rmax  = pp_rmax
  this%electrons%pp_gth   = pp_gth
  this%electrons%pp_gpair = pp_gpair
  
  this%positrons%pp_rmin  = pp_rmin
  this%positrons%pp_rmax  = pp_rmax
  this%positrons%pp_gth   = pp_gth
  this%positrons%pp_gpair = pp_gpair

  ! save QED variables
  ! this%electrons%qed_g_cutoff = qed_g_cutoff
  ! this%electrons%qed_p_cutoff = qed_p_cutoff
  ! this%electrons%qed_zeta = qed_zeta

  ! this%positrons%qed_g_cutoff = qed_g_cutoff
  ! this%positrons%qed_p_cutoff = qed_p_cutoff
  ! this%positrons%qed_zeta = qed_zeta

  ! read group diagnostics input
  call this % diag % read_input(input_file, this, grid)

end subroutine read_input_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine init_group( this, geometry, grid )

  use m_system
  use m_parameters
  use m_grid_define, only : t_grid

  implicit none

  class( t_species_group ), intent(inout) :: this
  class( t_geometry_gr ), intent(inout), pointer :: geometry
  class( t_grid ), intent(in) :: grid

  integer, dimension (2, p_x_dim) :: gc_num
  
  this%geometry => geometry
  call this % diag % init(this, grid)

  if (this % diag % ndump_fac_nemit > 0) then
    this % electrons % track_nemit = .true.
    this % positrons % track_nemit = .true.

    gc_num = 0
    call this % electrons % nemit % new( p_x_dim, 1, grid%my_nx(3,:), gc_num, this % electrons % dx(1:2), .true. )
    call this % positrons % nemit % new( p_x_dim, 1, grid%my_nx(3,:), gc_num, this % positrons % dx(1:2), .true. )    
    
    if (this%cooling_type == p_cooling_qed) then
      this % photons % track_nemit = .true.
      call this % photons % nemit % new( p_x_dim, 1, grid%my_nx(3,:), gc_num, this % electrons % dx(1:2), .true. )
    endif
  endif

  if ( this % inj_limit == p_inj_limit_fracnGJpole .or. &
       this % inj_limit == p_inj_limit_fracnGJpolealt .or. &
       this % inj_limit == p_inj_limit_fracnGJlocal .or. &
       this % inj_limit == p_inj_limit_fracnGJlocalalt .or. &
       this % inj_limit == p_inj_limit_fracEBteo .or. &
       this % inj_limit == p_inj_limit_fracEBteolim .or. &
       this % inj_limit == p_inj_limit_sigma ) then
    gc_num = 1
    call this % electrons % charge % new( p_x_dim, 1, grid%my_nx(3,:), gc_num, this % electrons % dx(1:2), .true. )
    call this % positrons % charge % new( p_x_dim, 1, grid%my_nx(3,:), gc_num, this % positrons % dx(1:2), .true. )
  endif
  
end subroutine init_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set the pointers to species in the group
!-------------------------------------------------------------------------------
subroutine set_species_group( this, electrons, positrons, photons )

  use m_system
  use m_parameters
  use m_species_define, only : t_species

  implicit none

  class( t_species_group ), intent(inout) :: this
  class( t_species ), intent(inout), target :: electrons, positrons, photons

  select type( electrons )
    class is ( t_species_gr )
      ! set pointer to group photons
      ! electrons%grp_photons => this%photons

      this%electrons => electrons
    class default
      ERROR('set_species_group must be called with t_species_gr objects')
      call abort_program( p_err_invalid )
  end select

  select type( positrons )
    class is ( t_species_gr )
      ! set pointer to group photons
      ! positrons%grp_photons => this%photons

      this%positrons => positrons
    class default
      ERROR('set_species_group must be called with t_species_gr objects')
      call abort_program( p_err_invalid )
  end select

  select type( photons )
    class is ( t_photons_gr )
      this%photons => photons
    class default
      ERROR('set_species_group must be called with t_species_gr objects')
      call abort_program( p_err_invalid )
  end select

  ! set pointers to group reciprocal species
  this%electrons%grp_recip => this%positrons
  this%positrons%grp_recip => this%electrons

end subroutine set_species_group

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Printout the algorithm used by each group
!-------------------------------------------------------------------------------
subroutine list_algorithm_group( this )

  use m_species_define, only : t_species
  implicit none
  class( t_species_group ), intent(in) :: this
  
  call this % electrons % list_algorithm()
  call this % positrons % list_algorithm()
  call this % photons % list_algorithm()

end subroutine list_algorithm_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine inject_group( this, b, e, dt, no_co )

  use m_vdf_define, only : t_vdf
  use m_node_conf, only : t_node_conf
  use m_species_define_gr
  implicit none
  class( t_species_group ), intent(inout) :: this
  class( t_vdf ), intent(in) :: b, e
  real( p_double ) :: dt
  class( t_node_conf ), intent(in) :: no_co
  
  select case ( this%electrons%push_type )
  case ( p_std )
    select case ( this%inj_local )
    case ( p_inj_local_surf )
      call inject_surface_group( this, b, e, dt, no_co )
    case ( p_inj_local_vol )
      call inject_volume_group( this, b, e, dt, no_co )
    case default
      ! do nothing
    end select
  case ( p_rk4 )
    select case ( this%inj_local )
    case ( p_inj_local_surf )
      call inject_surface_group_gr( this, b, e, dt, no_co )
    case ( p_inj_local_vol )
      call inject_volume_group_gr( this, b, e, dt, no_co )
    case default
      ! do nothing
    end select
  case ( p_rk6 )
    select case ( this%inj_local )
    case ( p_inj_local_surf )
      call inject_surface_group_gr( this, b, e, dt, no_co )
    case ( p_inj_local_vol )
      call inject_volume_group_gr( this, b, e, dt, no_co )
    case default
      ! do nothing
    end select
  end select

  if ( this % inj_limit == p_inj_limit_fracnGJpole .or. &
       this % inj_limit == p_inj_limit_fracnGJpolealt .or. &
       this % inj_limit == p_inj_limit_fracnGJlocal .or. &
       this % inj_limit == p_inj_limit_fracnGJlocalalt .or. &
       this % inj_limit == p_inj_limit_fracEBteo .or. &
       this % inj_limit == p_inj_limit_fracEBteolim .or. &
       this % inj_limit == p_inj_limit_sigma ) then
    call this % electrons % charge % zero()
    call this % positrons % charge % zero()
  endif

end subroutine inject_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine inject_surface_group( this, b, e, dt, no_co )

  use m_vdf_define, only : t_vdf
  use m_node_conf, only : t_node_conf
  implicit none

  ! dummy variables
  class( t_species_group ), intent(inout) :: this
  class( t_vdf ), intent(in) :: b, e
  real( p_double ) :: dt
  class( t_node_conf ), intent(in) :: no_co

  ! local variables
  real(p_k_part) :: r, t, ct, st, ph, cph, sph
  real(p_k_part) :: pr, pt, pph

  integer :: i, i1, i2, ppcell

  integer, parameter :: p_max_x_dim = 3
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2
  real(p_k_part), dimension(:,:), pointer :: ppos_cell
  real(p_k_part)  :: x1

  ! attributes for the new particles
  integer          , dimension(p_x_dim  ) :: ix
  real(p_k_part)   , dimension(p_x_dim+3) :: x
  real(p_k_part)   , dimension(p_p_dim  ) :: p
  real(p_k_part) :: q_norm, q_aux, q_inj, n_aux
  real(p_k_part) :: q, gamma, aveposden, aveeleden, epar
  real(p_k_part) :: Erave, Epave, Btave, Bpave
  real(p_k_part) :: nGJpole, nGJlocal, EdotB, B2
  real(p_k_part) :: tmin, tmax
  logical        :: inject
  
  class( t_geometry_gr ), pointer :: g

  ! executable statements
  g => this % geometry

  if ( no_co%my_ngp(1) .ne. 1 )  return

  do i = 1, p_x_dim
    ! get half distance between particles
    dxp_2(i) = 0.5_p_k_part/this%electrons%num_par_x(i)
  enddo  

  ! find total number of particles per cell
  ppcell = this%electrons%num_par_x(1)
  do i = 2, p_x_dim
    ppcell = ppcell * this%electrons%num_par_x(i)
  enddo

  ! initialize temp buffers
  call alloc( ppos_cell, (/p_x_dim, ppcell/) )

  ! get particle positions inside the cell depending on the number of particles per cell
  i = 0
  do i1 = 1, this%electrons%num_par_x(1)
    x1 = (2*i1 - 1 - this%electrons%num_par_x(1)) * dxp_2(1)
    do i2 = 1, this%electrons%num_par_x(2)
      i = i + 1
      ppos_cell(1,i) = x1
      ppos_cell(2,i) = (2*i2 - 1 - this%electrons%num_par_x(2)) * dxp_2(2)
    enddo
  enddo

  ! find normalization factor (normalized by the CFL condition)
  q_norm = real( ( (1.+1./g%A)**3 - (1.-1./g%A)**3 ) * 4.188790204786391 * dt / g%dtcfl, p_k_part ) ! 4/3 * pi

  i1 = 1
  do i2 = 1, b%nx_(2)

    Btave = (b%f2(2, i1, i2) + b%f2(2, i1-1, i2) + b%f2(2, i1, i2+1) +&
             b%f2(2, i1-1, i2+1))/4._p_double
    Bpave = (b%f2(3, i1, i2) + b%f2(3, i1-1, i2))/2._p_double

    ! total magnetic field squared (centered)
    B2 = b%f2(1, i1, i2)**2 + Btave**2 + Bpave**2 + 1.

    ! local GJ density
    select case( g%bprofile )
    case( p_bprofile_monopole )
      nGJlocal = 0._p_double
      epar     = 0._p_double
      nGJpole  = 0._p_double
      
    case( p_bprofile_dipole )
      !nGJlocal is valid only on the stellar surface
      nGJlocal = g%bs * abs( g%omega ) * abs( 1._p_double - 3._p_double * g%ct%f1(2, i2)**2 )
      !parallel E field is valid only on the stellar surface (exact)
      epar     = g%bs * abs( g%omega ) * abs( 2._p_double * g%ct%f1(2, i2)**3 ) /&
                          sqrt( 3._p_double * g%ct%f1(2, i2)**2 + 1._p_double)
      !nGJpole is valid only on the polar surface (exact)
      nGJpole  =  2._p_double * g%bs * abs( g%omega )

    case default 
      nGJlocal = 0._p_double
      epar     = 0._p_double
      nGJpole  = 0._p_double

    end select

    ! determine if particles should be created according to injection criterion
    select case( this%inj_limit )
    case( p_inj_limit_none )
      inject = .true.

    case( p_inj_limit_fracnGJpole )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden + aveeleden ) < (this%k_limit * nGJpole )
      
    case( p_inj_limit_fracnGJpolealt )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden - aveeleden ) < (this%k_limit * nGJpole )

    case( p_inj_limit_fracnGJlocal )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden + aveeleden ) < (this%k_limit * nGJlocal )
      
    case( p_inj_limit_fracnGJlocalalt )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden - aveeleden ) < (this%k_limit * nGJlocal )
      
    case( p_inj_limit_fracEB )
      Erave = (e%f2(1, i1, i2) + e%f2(1, i1-1, i2) + e%f2(1, i1, i2+1) + e%f2(1, i1-1, i2+1))/4.
      Epave = (e%f2(3, i1, i2) + e%f2(3, i1, i2+1))/2.
      
      EdotB = Erave * b%f2(1, i1, i2) + &
              e%f2(2, i1, i2) * Btave + &
              Epave * Bpave

      inject = abs( EdotB ) / sqrt(B2) > this%k_limit

    case( p_inj_limit_fracEBteo )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden - aveeleden ) < (this%k_limit * 2._p_double * epar )

    case( p_inj_limit_sigma )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = real( B2 / ( abs( aveposden - aveeleden ) + 1._p_double ) , p_k_part) > this%sig_limit
      
    case default
      inject = .false.

    end select


    if ( .not. inject ) cycle

    ! determine charge density to inject depending on chosen injection mechanism
    select case( this%inj_type )
    case( p_inj_type_const )
      n_aux = real( this%ks, p_k_part )

    case( p_inj_type_fracnGJ )
      n_aux = real( this%ks * nGJlocal, p_k_part )

    case( p_inj_type_fracnGJpole )
      n_aux = real( this%ks * nGJpole, p_k_part )

    case( p_inj_type_fracEB )
      Erave = (e%f2(1, i1, i2) + e%f2(1, i1-1, i2) + e%f2(1, i1, i2+1) +&
               e%f2(1, i1-1, i2+1))/4._p_double
      Epave = (e%f2(3, i1, i2) + e%f2(3, i1, i2+1))/2._p_double

      EdotB = Erave * b%f2(1, i1, i2) + &
              e%f2(2, i1, i2) * Btave + &
              Epave * Bpave

      n_aux = real( this%ks * 2._p_double * abs( EdotB ) / sqrt(B2), p_k_part )

    case( p_inj_type_fracEBteo )
      n_aux = real( this%ks * 2._p_double * epar, p_k_part )

    case default
      ! do nothing
      n_aux = 0._p_double

    end select

    ! if ( n_aux .lt. this%electrons%den_min ) then cycle

    q_aux = n_aux * q_norm

    ! create particles
    do i=1, ppcell

      ix(1) = i1
      ix(2) = i2
      x(1) = -0.5
      call random_number(x(2))
      x(2) = x(2) - 0.5

      ! build particle position and momentum
      r = g%r%f1(1, i1)
      t = g%t%f1(2, i2) + g%dt * x(2)

      ! added to test gamma max with no current deposition
      ! should be commented otherwise 
      !if ( t > 0.2 ) cycle                    

      ct = cos( t )
      st = sin( t )
      
      call random_number(ph)
      ph = ph * 2*pi
      
      cph = cos( ph )
      sph = sin( ph )

      x(3) = r * st * cph
      x(4) = r * st * sph
      x(5) = r * ct

      gamma = 1._p_double / sqrt( 1._p_double - g%omega**2 * r**2 * st**2 - this%vs**2 )
      pr  = real( this%vs * gamma * abs( ct / sqrt( ct**2 + 0.25_p_double * st**2 ) ) , p_k_part)
      pt  = real( this%vs * gamma * abs( 0.5_p_double * st / sqrt( ct**2 + 0.25_p_double * st**2 ) ) * ct / abs(ct) , p_k_part)
      pph = real( g%omega * r * st * gamma , p_k_part )

      p(1) = st * cph * pr + ct * cph * pt - sph * pph
      p(2) = st * sph * pr + ct * sph * pt + cph * pph
      p(3) = ct * pr - st * pt

      tmin = max(t - .5 * this%electrons%dx(2), 0.0_p_double)
      tmax = min(t + .5 * this%electrons%dx(2), pi)
      
      q_inj = q_aux * r**3 * abs( cos(tmin) - cos(tmax) )  / ppcell ! * st

      ! create electron
      q = this%electrons%rqm * q_inj
      call this%electrons%create_particle( ix, x, p, q )

      ! create positron/ion
      q = this%positrons%rqm * q_inj
      call this%positrons%create_particle( ix, x, p, q )

    enddo

  enddo

  call freemem( ppos_cell )

end subroutine inject_surface_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine inject_surface_group_gr( this, b, e, dt, no_co )

  use m_vdf_define, only : t_vdf
  use m_node_conf, only : t_node_conf
  implicit none

  ! dummy variables
  class( t_species_group ), intent(inout) :: this
  class( t_vdf ), intent(in) :: b, e
  real( p_double ) :: dt
  class( t_node_conf ), intent(in) :: no_co

  ! local variables
  real(p_k_part) :: r, t, ct, ph !, st

  integer :: i, i1, i2, ppcell

  integer, parameter :: p_max_x_dim = 3
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2
  real(p_k_part), dimension(:,:), pointer :: ppos_cell
  real(p_k_part)  :: x1

  ! attributes for the new particles
  integer          , dimension(p_x_dim  ) :: ix
  real(p_k_part)   , dimension(p_x_dim+1) :: x
  real(p_k_part)   , dimension(p_p_dim  ) :: p
  real(p_k_part) :: q_norm, q_aux, q_inj, n_aux
  real(p_k_part) :: q, gamma, omegagr, bamp
  real(p_k_part) :: braux, btaux, eraux,logaux !, injaux
  real(p_k_part) :: rminus, rplus, aplus,aminus, intr
  real(p_k_part) :: aveeleden, aveposden, epar
  real(p_k_part) :: Erave, Epave, Btave, Bpave
  real(p_k_part) :: nGJpole,nGJlocal, EdotB, B2
  real(p_k_part) :: tmin, tmax
  integer        :: i1max
  logical        :: inject
  real(p_double) :: shift, tc, tc2, tc3, tc4, tc5
  real(p_double) :: tc6, tc7, tc8, tc9, tc10
  real(p_double) :: ctmin, ctmax
  
  class( t_geometry_gr ), pointer :: g

  ! executable statements
  g => this % geometry

  if ( no_co%my_ngp(1) .ne. 1 )  return

  do i = 1, p_x_dim
    ! get half distance between particles
    dxp_2(i) = 0.5_p_k_part/this%electrons%num_par_x(i)
  enddo  

  ! find total number of particles per cell
  ppcell = this%electrons%num_par_x(1)
  do i = 2, p_x_dim
    ppcell = ppcell * this%electrons%num_par_x(i)
  enddo

  ! initialize temp buffers
  call alloc( ppos_cell, (/p_x_dim, ppcell/) )

  ! get particle positions inside the cell depending on the number of particles per cell
  i = 0
  do i1 = 1, this%electrons%num_par_x(1)
    x1 = (2*i1 - 1 - this%electrons%num_par_x(1)) * dxp_2(1)
    do i2 = 1, this%electrons%num_par_x(2)
      i = i + 1
      ppos_cell(1,i) = x1
      ppos_cell(2,i) = (2*i2 - 1 - this%electrons%num_par_x(2)) * dxp_2(2)
    enddo
  enddo

  ! find normalization factor (normalized by the CFL condition)
  q_norm = real( 2*pi * dt / g%dtcfl, p_k_part )

  ! correction of GR to stellar rotation
  omegagr = g%omega - g%beta0

  ! used to evaluate the polynomial approx of trig functions
  shift = pi_2

  ! All quantities are evaluated at r=Rstar
  i1 = 1

  ! find the index interval for local average in density - fracEBteolim injection
  i1max = int(g%nr * log(1.02) / log(g%rmax) + 1)

  ! Each field is composed of an amplitude, an aux function and the theta dependence
  ! Example: Br = bamp * braux * cos(theta)
  ! Auxiliar log function (appears several times)
  logaux = log(1._p_double - g%rs)
  ! B field amplitude
  bamp = 3._p_double * g%bs / (2._p_double * g%rs**3)
  ! Br auxiliar function on the stellar surface (w/o theta dependence and amplitude)
  braux = - ( 2._p_double * logaux + g%rs * (2._p_double + g%rs) )
  ! Bt auxiliar function on the stellar surface (w/o theta dependence and amplitude)
  btaux = ( 2._p_double * (1._p_double - g%rs) * logaux + g%rs * (2._p_double - g%rs) ) / g%alpha%f1(1,i1)
  ! Er auxiliar function on the stellar surface (w/o theta dependence and amplitude)
  eraux = - g%beta0 * (logaux + g%rs) +&
        ((3._p_double - 4._p_double / g%rs) * logaux + g%rs * (1._p_double + g%rs/6._p_double) - 4._p_double) *&
        ( omegagr * (logaux + g%rs) + g%omega * g%rs**2 / 2._p_double ) /&
        (3._p_double * ((1._p_double - 2._p_double / g%rs) * (1._p_double - g%rs) * logaux -&
                         2._p_double * (1._p_double - g%rs) - g%rs**2 / 6._p_double))  
  
  do i2 = 1, b%nx_(2)
    Btave = (b%f2(2, i1, i2) + b%f2(2, i1-1, i2) + b%f2(2, i1, i2+1) +&
             b%f2(2, i1-1, i2+1))/4._p_double
    Bpave = (b%f2(3, i1, i2) + b%f2(3, i1-1, i2))/2._p_double

    ! total magnetic field squared (centered)
    B2 = b%f2(1, i1, i2)**2 + Btave**2 + Bpave**2 + 1._p_double

    ! local GJ density
    select case( g%bprofile )
    case( p_bprofile_monopole )
      nGJlocal = 0._p_double
      epar     = 0._p_double
      nGJpole  = 0._p_double

    case( p_bprofile_dipole )
      !nGJlocal is valid only on the stellar surface
      nGJlocal = 2._p_double * bamp * abs( (1._p_double - 3._p_double * g%ct%f1(2, i2)**2) * eraux )
      !parallel E field is valid only on the stellar surface (exact)
      epar =  bamp * braux *&
                 abs(g%ct%f1(2, i2) * ( (1._p_double - 3._p_double * g%ct%f1(2, i2)**2) * eraux -&
                     omegagr * btaux * g%st%f1(2, i2)**2 / g%alpha%f1(1,i1) )) /&
                 sqrt(braux**2 * g%ct%f1(2, i2)**2 + btaux**2 * g%st%f1(2, i2)**2)
      !nGJpole is valid only on the polar surface (exact)
      nGJpole  = 2._p_double * bamp * abs(- 2._p_double * eraux)

    case default 
      nGJlocal = 0._p_double
      epar     = 0._p_double
      nGJpole  = 0._p_double

    end select

    ! determine if particles should be created according to injection criterion
    select case( this%inj_limit )
    case( p_inj_limit_none )
      inject = .true.

    case( p_inj_limit_fracnGJpole )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden + aveeleden ) < (this%k_limit * nGJpole )

    case( p_inj_limit_fracnGJpolealt )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden - aveeleden ) < (this%k_limit * nGJpole )
    
    case( p_inj_limit_fracnGJlocal )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden + aveeleden ) < (this%k_limit * nGJlocal )

    case( p_inj_limit_fracnGJlocalalt )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden - aveeleden ) < (this%k_limit * nGJlocal )
      
    case( p_inj_limit_fracEB )
      Erave = (e%f2(1, i1, i2) + e%f2(1, i1-1, i2) + e%f2(1, i1, i2+1) + e%f2(1, i1-1, i2+1))/4._p_double
      Epave = (e%f2(3, i1, i2) + e%f2(3, i1, i2+1))/2._p_double
      
      EdotB = Erave * b%f2(1, i1, i2) + &
              e%f2(2, i1, i2) * Btave + &
              Epave * Bpave

      inject = abs( EdotB ) / sqrt(B2) > this%k_limit

    case( p_inj_limit_fracEBteolim )
      aveposden = (sum(this%positrons%charge%f2(1, i1:i1max , i2))+sum(this%positrons%charge%f2(1, i1:i1max  , i2+1))) / i1max / 2._p_double
      aveeleden = (sum(this%electrons%charge%f2(1, i1:i1max , i2))+sum(this%electrons%charge%f2(1, i1:i1max  , i2+1))) / i1max / 2._p_double
      inject = abs( aveposden - aveeleden ) < (abs(this%k_limit - (this%k_limit - 1) * g%st%f1(2, i2)**5)* 2._p_double * epar )

    case( p_inj_limit_fracEBteo )
      aveposden = (sum(this%positrons%charge%f2(1, i1:i1max , i2))+sum(this%positrons%charge%f2(1, i1:i1max  , i2+1))) / i1max / 2._p_double
      aveeleden = (sum(this%electrons%charge%f2(1, i1:i1max , i2))+sum(this%electrons%charge%f2(1, i1:i1max  , i2+1))) / i1max / 2._p_double
      !aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      !aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = abs( aveposden - aveeleden ) < (this%k_limit * 2._p_double * epar )

    case( p_inj_limit_sigma )
      aveposden = 0.5_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1))
      aveeleden = 0.5_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1))
      inject = real( B2 / ( abs( aveposden - aveeleden ) + 1._p_double ) , p_k_part) > this%sig_limit

    case default
      inject = .false.

    end select


    if ( .not. inject ) cycle

    ! determine charge density to inject depending on chosen injection mechanism
    select case( this%inj_type )
    case( p_inj_type_const )
      n_aux = real( this%ks, p_k_part ) 

    case( p_inj_type_fracnGJ )
      n_aux = real( this%ks * nGJlocal, p_k_part ) 

    case( p_inj_type_fracnGJpole )
      n_aux = real( this%ks * nGJpole, p_k_part ) 

    case( p_inj_type_fracEB )
      Erave = (e%f2(1, i1, i2) + e%f2(1, i1-1, i2) + e%f2(1, i1, i2+1) + e%f2(1, i1-1, i2+1))/4._p_double
      Epave = (e%f2(3, i1, i2) + e%f2(3, i1, i2+1))/2._p_double

      EdotB = Erave * b%f2(1, i1, i2) + &
              e%f2(2, i1, i2) * Btave + &
              Epave * Bpave

      n_aux = real( this%ks * 2._p_double * abs( EdotB ) / sqrt(B2), p_k_part ) 

    case( p_inj_type_fracEBteo )
      n_aux = real( this%ks * 2._p_double * epar, p_k_part )

    case default
      ! do nothing
      n_aux = 0._p_double

    end select

    ! if ( n_aux .lt. this%electrons%den_min ) then cycle

    q_aux = n_aux * q_norm 

    ! create particles
    do i=1, ppcell

      ix(1) = i1
      ix(2) = i2
      x(1) = -0.5
      
      call random_number(x(2))
      x(2) = x(2) - 0.5

      ! build particle position and momentum
      r = g%r%f1(1, i1)
      t = g%t%f1(2, i2) + g%dt * x(2) 

      ! added to test gamma max with no current deposition
      ! should be commented otherwise 
      ! if ( t > 0.2 ) cycle

      ! evaluate the sin and cos of theta around pi/2
      tc    = t - shift
      tc2   = tc * tc
      tc3   = tc2 * tc
      tc4   = tc3 * tc
      tc5   = tc4 * tc
      tc6   = tc5 * tc
      tc7   = tc6 * tc
      tc8   = tc7 * tc
      tc9   = tc8 * tc
      tc10  = tc9 * tc

      ! required if injection is sine dependent (commented below)
      !st = 1.0_p_double - 0.500000000000000_p_double * tc2 +&
      !                    0.041666666666667_p_double * tc4 - 0.001388888888889_p_double * tc6  +&
      !                    0.000024801587302_p_double * tc8 - 0.000000275573192_p_double * tc10

      ct = - tc + 0.166666666666667_p_double * tc3 - 0.008333333333333_p_double * tc5 +&
                  0.000198412698413_p_double * tc7 - 0.000002755731922_p_double * tc9
      
      call random_number(ph)
      ph = pi * (2*ph - 1)

      x(3) = ph

      !gamma = 1._p_double / sqrt( 1._p_double - omegagr**2 * r**2 * st**2 / &
      !       (1._p_double - g%rs/r) - this%vs**2 )
      !injaux = sqrt( ct**2 * braux**2 + btaux**2 * st**2 )

      !p(1)  = real( this%vs * gamma * abs(braux * ct / injaux)                 , p_k_part)
      !p(2)  = real( this%vs * gamma * abs(btaux * st / injaux) * ct / abs(ct)  , p_k_part)
      !p(3)  = real( omegagr * r * st * gamma / g%alpha%f1(1,i1) , p_k_part )
      
      gamma = 1._p_double / sqrt( 1._p_double - this%vs**2 ) 
      p(1)  = real( this%vs * gamma * abs(b%f2(1, i1, i2))/sqrt(B2)                , p_k_part)
      p(2)  = real( this%vs * gamma * abs(     Btave     )/sqrt(B2) * ct / abs(ct) , p_k_part)
      !p(3)  = real( this%vs * gamma * (Bpave)/sqrt(B2) , p_k_part)
      !p(3)  = real( this%vs * gamma * abs(Bpave)/sqrt(B2) , p_k_part)
      p(3)  = 0._p_double
      
      ! evaluate cos of tmin
      tmin = max(t - .5_p_double * this%electrons%dx(2), 0.0_p_double)
      tc    = tmin - shift
      tc2   = tc * tc
      tc3   = tc2 * tc
      tc5   = tc2 * tc3
      tc7   = tc2 * tc5
      tc9   = tc2 * tc7

      ctmin = - tc + 0.166666666666667_p_double * tc3 - 0.008333333333333_p_double * tc5 +&
                     0.000198412698413_p_double * tc7 - 0.000002755731922_p_double * tc9

      ! evaluate cos of tmax
      tmax = min(t + .5_p_double * this%electrons%dx(2), pi)
      tc    = tmax - shift
      tc2   = tc * tc
      tc3   = tc2 * tc
      tc5   = tc2 * tc3
      tc7   = tc2 * tc5
      tc9   = tc2 * tc7

      ctmax = - tc + 0.166666666666667_p_double * tc3 - 0.008333333333333_p_double * tc5 +&
                     0.000198412698413_p_double * tc7 - 0.000002755731922_p_double * tc9

      rplus  = 2._p_double * r / (1._p_double + 1._p_double / g%delta) !particle r+
      rminus = 2._p_double * r / (1._p_double + g%delta)               !particle r-
      aplus  = sqrt( 1._p_double - g%rs / rplus )                      !alpha(r+)
      aminus = sqrt( 1._p_double - g%rs / rminus)                      !alpha(r-)
      intr = 0.0416666666666667_p_double *&                            !1/24
          (rplus  * aplus  * (8._p_double * rplus**2  +&
              10._p_double * rplus  * g%rs + 15._p_double * g%rs**2)  -&
            rminus * aminus * (8._p_double * rminus**2 +&
              10._p_double * rminus * g%rs + 15._p_double * g%rs**2)) +&
          0.3125_p_double * g%rs**3 * log(((1._p_double+aplus)*(1._p_double-aminus))/&
                                                  ((1._p_double-aplus)*(1._p_double+aminus)))
      
      q_inj = q_aux * intr * abs( ctmin - ctmax ) / ppcell ! * st

      ! create electron
      q = this%electrons%rqm * q_inj
      call this%electrons%create_particle( ix, x, p, q )

      ! create positron/ion
      q = this%positrons%rqm * q_inj
      call this%positrons%create_particle( ix, x, p, q )

    enddo

  enddo

  call freemem( ppos_cell )

end subroutine inject_surface_group_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine inject_volume_group( this, b, e, dt, no_co )

  use m_vdf_define, only : t_vdf
  use m_node_conf, only : t_node_conf
  implicit none

  ! dummy variables
  class( t_species_group ), intent(inout) :: this
  class( t_vdf ), intent(in) :: b, e
  real( p_double ) :: dt
  class( t_node_conf ), intent(in) :: no_co

  ! local variables
  real(p_k_part) :: r, t, ct, st, ph, cph, sph

  integer :: i, i1, i2, ppcell

  integer, parameter :: p_max_x_dim = 3
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2
  real(p_k_part), dimension(:,:), pointer :: ppos_cell
  real(p_k_part)  :: x1

  ! attributes for the new particles
  integer          , dimension(p_x_dim  ) :: ix
  real(p_k_part)   , dimension(p_x_dim+3) :: x
  real(p_k_part)   , dimension(p_p_dim  ) :: p
  real(p_k_part) :: q_norm, q_aux, q_inj, n_aux
  real(p_k_part) :: q, aveposden, aveeleden, epar
  real(p_k_part) :: nGJpole, nGJlocal, EdotB, B2
  real(p_k_part) :: Erave, Etave, Epave, Brave, Btave, rexp
  real(p_k_part) :: tmin, tmax
  logical        :: inject
  
  class( t_geometry_gr ), pointer :: g

  ! executable statements
  g => this % geometry

  do i = 1, p_x_dim
    ! get half distance between particles
    dxp_2(i) = 0.5_p_k_part/this%electrons%num_par_x(i)
  enddo  

  ! find total number of particles per cell
  ppcell = this%electrons%num_par_x(1)
  do i = 2, p_x_dim
    ppcell = ppcell * this%electrons%num_par_x(i)
  enddo

  ! initialize temp buffers
  call alloc( ppos_cell, (/p_x_dim, ppcell/) )

  ! get particle positions inside the cell depending on the number of particles per cell
  i = 0
  do i1 = 1, this%electrons%num_par_x(1)
    x1 = (2*i1 - 1 - this%electrons%num_par_x(1)) * dxp_2(1)
    do i2 = 1, this%electrons%num_par_x(2)
      i = i + 1
      ppos_cell(1,i) = x1
      ppos_cell(2,i) = (2*i2 - 1 - this%electrons%num_par_x(2)) * dxp_2(2)
    enddo
  enddo

  ! find normalization factor (normalized by the CFL condition)
  q_norm = real( ( (1.+1./g%A)**3 - (1.-1./g%A)**3 ) * 4.188790204786391 * dt / g%dtcfl, p_k_part ) ! 4/3 * pi

  do i1 = 1, b%nx_(1)
    ! r dependent quantities
    rexp = g%r%f1(2, i1)**3

    do i2 = 1, b%nx_(2)
      
      Brave = 0.5_p_double * (b%f2(1, i1, i2) + b%f2(1, i1+1, i2))
      Btave = 0.5_p_double * (b%f2(2, i1, i2) + b%f2(2, i1, i2+1))

      ! total magnetic field squared (centered)
      B2 = Brave**2 + Btave**2 + b%f2(3, i1, i2)**2 + 1.

      ! local GJ density
      select case( g%bprofile )
      case( p_bprofile_monopole )
        nGJlocal = 0._p_double
        epar     = 0._p_double
        nGJpole  = 0._p_double

      case( p_bprofile_dipole )
        ! Quantities would decay with 1/r**5 but here are put to 1/r**3 (efield goes w/ 1/r**4 * div goes w/ 1/r)
        ! local value of the GJ density
        nGJlocal = g%bs * abs( g%omega ) * abs( 1._p_double - 3._p_double * g%ct%f1(2, i2)**2 ) / rexp
        ! local value of the parallel electric field
        epar     = g%bs * abs( g%omega ) * abs( 2._p_double * g%ct%f1(2, i2)**3 ) /&
                           (rexp * sqrt(3._p_double * g%ct%f1(2, i2)**2 + 1._p_double))
        ! value of the GJ density at the polar axis
        nGJpole  = g%bs * abs( g%omega ) * 2._p_double / rexp

      case default 
        nGJlocal = 0._p_double
        epar     = 0._p_double
        nGJpole  = 0._p_double

      end select

      ! determine if particles should be created according to injection criterion
      select case( this%inj_limit )
      case( p_inj_limit_none )
        inject = .true.

      case( p_inj_limit_fracnGJpole )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden + aveeleden ) < (this%k_limit * nGJpole )

      case( p_inj_limit_fracnGJpolealt )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden - aveeleden ) < (this%k_limit * nGJpole )

      case( p_inj_limit_fracnGJlocal )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden + aveeleden ) < (this%k_limit * nGJlocal )

      case( p_inj_limit_fracnGJlocalalt )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden - aveeleden ) < (this%k_limit * nGJlocal )

      case( p_inj_limit_fracEB )
        Erave = 0.5_p_double * (e%f2(1, i1, i2) + e%f2(1, i1, i2+1))
        Etave = 0.5_p_double * (e%f2(2, i1, i2) + e%f2(2, i1+1, i2))
        Epave = 0.25_p_double * (e%f2(3, i1, i2) + e%f2(3, i1+1, i2) + e%f2(3, i1, i2+1) + e%f2(3, i1+1, i2+1))

        EdotB = Erave * Brave + Etave * Btave + Epave * b%f2(3, i1, i2)

        inject = abs( EdotB ) / sqrt(B2) > this%k_limit

      case( p_inj_limit_fracEBteo )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden - aveeleden ) < (this%k_limit * 2._p_double * epar )

      case( p_inj_limit_sigma )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = real( B2 / ( abs( aveposden - aveeleden ) + 1._p_double) , p_k_part) > this%sig_limit

      case default
        inject = .false.

      end select


      if ( .not. inject ) cycle

      ! determine charge density to inject depending on chosen injection mechanism
      select case( this%inj_type )
      case( p_inj_type_const )
        n_aux = real( this%ks, p_k_part )

      case( p_inj_type_fracnGJ )
        n_aux = real( this%ks * nGJlocal, p_k_part )

      case( p_inj_type_fracnGJpole )
        n_aux = real( this%ks * nGJpole, p_k_part )

      case( p_inj_type_fracEB )
        Erave = 0.5_p_double * (e%f2(1, i1, i2) + e%f2(1, i1, i2+1))
        Etave = 0.5_p_double * (e%f2(2, i1, i2) + e%f2(2, i1+1, i2))
        Epave = 0.25_p_double * (e%f2(3, i1, i2) + e%f2(3, i1+1, i2) + e%f2(3, i1, i2+1) + e%f2(3, i1+1, i2+1))

        EdotB = Erave * Brave + Etave * Btave + Epave * b%f2(3, i1, i2)

        n_aux = real( this%ks * 2._p_double * abs( EdotB ) / sqrt(B2), p_k_part )

      case( p_inj_type_fracEBteo )
        n_aux = real( this%ks * 2._p_double * epar, p_k_part )

      case default
        ! do nothing
        n_aux = 0._p_double

      end select

      q_aux = n_aux * q_norm 

      ! create particles
      do i=1, ppcell

        ix(1) = i1
        ix(2) = i2
        call random_number(x(1))
        x(1) = x(1) - 0.5
        call random_number(x(2))
        x(2) = x(2) - 0.5

        ! build particle position and momentum
        r = g%r%f1(2, i1) + g%dr%f1(1, i1) * x(1)
        t = g%t%f1(2, i2) + g%dt * x(2)                    

        ct = cos( t )
        st = sin( t )
        
        call random_number(ph)
        ph = ph * 2*pi
        cph = cos( ph )
        sph = sin( ph )

        x(3) = r * st * cph
        x(4) = r * st * sph
        x(5) = r * ct

        p = 0

        tmin = max(t - .5 * this%electrons%dx(2), 0.0_p_double)
        tmax = min(t + .5 * this%electrons%dx(2), pi)
        
        q_inj = q_aux * r**3 * abs( cos(tmin) - cos(tmax) )  / ppcell ! * st

        ! create electron
        q = this%electrons%rqm * q_inj
        call this%electrons%create_particle( ix, x, p, q )

        ! create positron/ion
        q = this%positrons%rqm * q_inj
        call this%positrons%create_particle( ix, x, p, q )

      enddo

    enddo
  enddo

  call freemem( ppos_cell )

end subroutine inject_volume_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine inject_volume_group_gr( this, b, e, dt, no_co )

  use m_vdf_define, only : t_vdf
  use m_node_conf, only : t_node_conf
  implicit none

  ! dummy variables
  class( t_species_group ), intent(inout) :: this
  class( t_vdf ), intent(in) :: b, e
  real( p_double ) :: dt
  class( t_node_conf ), intent(in) :: no_co

  ! local variables
  real(p_k_part) :: r, t, ct, st, ph

  integer :: i, i1, i2, ppcell

  integer, parameter :: p_max_x_dim = 3
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2
  real(p_k_part), dimension(:,:), pointer :: ppos_cell
  real(p_k_part)  :: x1

  ! attributes for the new particles
  integer          , dimension(p_x_dim  ) :: ix
  real(p_k_part)   , dimension(p_x_dim+1) :: x
  real(p_k_part)   , dimension(p_p_dim  ) :: p
  real(p_k_part) :: q_norm, q_aux, q_inj, n_aux
  real(p_k_part) :: q, omegagr, rexp
  real(p_k_part) :: aveeleden, aveposden, epar
  real(p_k_part) :: bamp, braux, btaux, eraux, logaux
  real(p_k_part) :: rminus, rplus, aplus,aminus, intr
  real(p_k_part) :: nGJpole, nGJlocal, EdotB, B2
  real(p_k_part) :: Erave, Etave, Epave, Brave, Btave
  real(p_k_part) :: tmin, tmax
  logical        :: inject
  
  class( t_geometry_gr ), pointer :: g

  ! executable statements
  g => this % geometry

  do i = 1, p_x_dim
    ! get half distance between particles
    dxp_2(i) = 0.5_p_k_part/this%electrons%num_par_x(i)
  enddo  

  ! find total number of particles per cell
  ppcell = this%electrons%num_par_x(1)
  do i = 2, p_x_dim
    ppcell = ppcell * this%electrons%num_par_x(i)
  enddo

  ! initialize temp buffers
  call alloc( ppos_cell, (/p_x_dim, ppcell/) )

  ! get particle positions inside the cell depending on the number of particles per cell
  i = 0
  do i1 = 1, this%electrons%num_par_x(1)
    x1 = (2*i1 - 1 - this%electrons%num_par_x(1)) * dxp_2(1)
    do i2 = 1, this%electrons%num_par_x(2)
      i = i + 1
      ppos_cell(1,i) = x1
      ppos_cell(2,i) = (2*i2 - 1 - this%electrons%num_par_x(2)) * dxp_2(2)
    enddo
  enddo

  ! find normalization factor (normalized by the CFL condition)
  q_norm = real( 2*pi * dt / g%dtcfl, p_k_part )

  ! correction of GR to stellar rotation
  omegagr = g%omega - g%beta0

  !field amplitude
  bamp = 3._p_double * g%bs / (2._p_double * g%rs**3)
  ! Auxiliar log function (appears several times)
  logaux = log(1._p_double - g%rs)
  ! Br auxiliar function on the stellar surface (w/o theta dependence and amplitude)
  braux = - ( 2._p_double * logaux + g%rs * (2._p_double + g%rs) )
  ! Bt auxiliar function on the stellar surface (w/o theta dependence and amplitude)
  btaux = ( 2._p_double * (1._p_double - g%rs) * logaux +&
           g%rs * (2._p_double - g%rs) ) / g%alpha%f1(1,1)
  ! Er auxiliar function on the stellar surface (w/o theta dependence and amplitude)
  eraux = - g%beta0 * (logaux + g%rs) +&
        ((3._p_double - 4._p_double / g%rs) * logaux + g%rs * (1._p_double + g%rs/6._p_double) - 4._p_double) *&
        ( omegagr * (logaux + g%rs) + g%omega * g%rs**2 / 2._p_double ) /&
        (3._p_double * ((1._p_double - 2._p_double / g%rs) * (1._p_double - g%rs) * logaux -&
                         2._p_double * (1._p_double - g%rs) - g%rs**2 / 6._p_double))  
  

  do i1 = 1, b%nx_(1)
    ! evaluate r dependent quantities
    rexp = g%r%f1(2, i1)**3
    
    do i2 = 1, b%nx_(2)
      
      Brave = 0.5_p_double * (b%f2(1, i1, i2) + b%f2(1, i1+1, i2))
      Btave = 0.5_p_double * (b%f2(2, i1, i2) + b%f2(2, i1, i2+1))

      ! total magnetic field squared (centered)
      B2 = Brave**2 + Btave**2 + b%f2(3, i1, i2)**2 + 1.

      ! local GJ density
      select case( g%bprofile )
      case( p_bprofile_monopole )
        !not implemented
        nGJlocal = 0._p_double
        epar     = 0._p_double
        nGJpole  = 0._p_double

      case( p_bprofile_dipole )
        ! quantities decay w/ 1/r**3 instead of 1/r**5
        ! complex formula for epar is approx by the surface local value times the radial decay (good approx)
        ! approx. nGJ local value
        nGJlocal = 2._p_double * bamp * abs( (1._p_double - 3._p_double * g%ct%f1(2, i2)**2) * eraux ) / rexp
        ! approx. parallel local electric field value 
        epar = bamp * braux *&
                 abs(g%ct%f1(2, i2) * ( (1._p_double - 3._p_double * g%ct%f1(2, i2)**2) * eraux -&
                     omegagr * btaux * g%st%f1(2, i2)**2 / g%alpha%f1(1,1) )) /&
                 sqrt(braux**2 * g%ct%f1(2, i2)**2 + btaux**2 * g%st%f1(2, i2)**2) / rexp

        ! approx. nGJ value along the polar axis
        nGJpole = 2._p_double * bamp * abs(- 2._p_double * eraux) / rexp

      case default 
        nGJlocal = 0._p_double
        epar     = 0._p_double
        nGJpole  = 0._p_double

      end select

      ! determine if particles should be created according to injection criterion
      select case( this%inj_limit )
      case( p_inj_limit_none )
        inject = .true.

      case( p_inj_limit_fracnGJpole )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden + aveeleden ) < (this%k_limit * nGJpole )

      case( p_inj_limit_fracnGJpolealt )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden - aveeleden ) < (this%k_limit * nGJpole )

      case( p_inj_limit_fracnGJlocal )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden + aveeleden ) < (this%k_limit * nGJlocal )

      case( p_inj_limit_fracnGJlocalalt )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden - aveeleden ) < (this%k_limit * nGJlocal )

      case( p_inj_limit_fracEB )
        Erave = 0.5_p_double * (e%f2(1, i1, i2) + e%f2(1, i1, i2+1))
        Etave = 0.5_p_double * (e%f2(2, i1, i2) + e%f2(2, i1+1, i2))
        Epave = 0.25_p_double * (e%f2(3, i1, i2) + e%f2(3, i1+1, i2) + e%f2(3, i1, i2+1) + e%f2(3, i1+1, i2+1))

        EdotB = Erave * Brave + Etave * Btave + Epave * b%f2(3, i1, i2)

        inject = abs( EdotB ) / sqrt(B2) > this%k_limit

      case( p_inj_limit_fracEBteo )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = abs( aveposden - aveeleden ) < (this%k_limit * 2._p_double * epar )

      case( p_inj_limit_sigma )
        aveposden = 0.25_p_double * (this%positrons%charge%f2(1, i1  , i2)+this%positrons%charge%f2(1, i1  , i2+1)+&
                                     this%positrons%charge%f2(1, i1+1, i2)+this%positrons%charge%f2(1, i1+1, i2+1))
        aveeleden = 0.25_p_double * (this%electrons%charge%f2(1, i1  , i2)+this%electrons%charge%f2(1, i1  , i2+1)+&
                                     this%electrons%charge%f2(1, i1+1, i2)+this%electrons%charge%f2(1, i1+1, i2+1))
        inject = real( B2 / ( abs( aveposden - aveeleden ) + 1._p_double) , p_k_part) > this%sig_limit

      case default
        inject = .false.

      end select


      if ( .not. inject ) cycle

      ! determine charge density to inject depending on chosen injection mechanism
      select case( this%inj_type )
      case( p_inj_type_const )
        n_aux = real( this%ks, p_k_part )

      case( p_inj_type_fracnGJ )
        n_aux = real( this%ks * nGJlocal, p_k_part )

      case( p_inj_type_fracnGJpole )
        n_aux = real( this%ks * nGJpole, p_k_part )

      case( p_inj_type_fracEB )
        Erave = 0.5_p_double * (e%f2(1, i1, i2) + e%f2(1, i1, i2+1))
        Etave = 0.5_p_double * (e%f2(2, i1, i2) + e%f2(2, i1+1, i2))
        Epave = 0.25_p_double * (e%f2(3, i1, i2) + e%f2(3, i1+1, i2) + e%f2(3, i1, i2+1) + e%f2(3, i1+1, i2+1))

        EdotB = Erave * Brave + Etave * Btave + Epave * b%f2(3, i1, i2)

        n_aux = real( this%ks * 2._p_double * abs( EdotB ) / sqrt(B2), p_k_part )

      case( p_inj_type_fracEBteo )
        n_aux = real( this%ks * 2._p_double * epar, p_k_part )

      case default
        ! do nothing
        n_aux = 0._p_double

      end select

      q_aux = n_aux * q_norm

      ! create particles
      do i=1, ppcell

        ix(1) = i1
        ix(2) = i2
        call random_number(x(1))
        x(1) = x(1) - 0.5
        call random_number(x(2))
        x(2) = x(2) - 0.5

        ! build particle position and momentum
        r = g%r%f1(2, i1) + g%dr%f1(1, i1) * x(1)
        t = g%t%f1(2, i2) + g%dt * x(2)                    

        ct = cos( t )
        st = sin( t )
        
        call random_number(ph)
        ph = pi * (2*ph - 1)
        
        x(3) = ph

        !particles injected with 0 velocity in volume
        p = 0

        tmin = max(t - .5 * this%electrons%dx(2), 0.0_p_double)
        tmax = min(t + .5 * this%electrons%dx(2), pi)
        
        rplus  = 2._p_double * r / (1._p_double + 1._p_double / g%delta) !particle r+
        rminus = 2._p_double * r / (1._p_double + g%delta)               !particle r-
        aplus  = sqrt( 1._p_double - g%rs / rplus )                      !alpha(r+)
        aminus = sqrt( 1._p_double - g%rs / rminus)                      !alpha(r-)
        intr = 0.0416666666666667_p_double *&                            !1/24
            (rplus  * aplus  * (8._p_double * rplus**2  +&
                10._p_double * rplus  * g%rs + 15._p_double * g%rs**2)  -&
              rminus * aminus * (8._p_double * rminus**2 +&
                10._p_double * rminus * g%rs + 15._p_double * g%rs**2)) +&
            0.3125_p_double * g%rs**3 * log(((1._p_double+aplus)*(1._p_double-aminus))/&
                                                    ((1._p_double-aplus)*(1._p_double+aminus)))
      
        q_inj = q_aux * intr * abs( cos(tmin) - cos(tmax) )  / ppcell ! * st

        ! create electron
        q = this%electrons%rqm * q_inj
        call this%electrons%create_particle( ix, x, p, q )

        ! create positron/ion
        q = this%positrons%rqm * q_inj
        call this%positrons%create_particle( ix, x, p, q )

      enddo

    enddo
  enddo

  call freemem( ppos_cell )

end subroutine inject_volume_group_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine report_group(this, g_space, grid, no_co, tstep, t)
  
  use m_space, only : t_space
  use m_grid_define, only : t_grid
  use m_node_conf, only : t_node_conf
  use m_time_step, only : t_time_step, test_if_report
  use m_species_diagnostics_gr, only : report_nemit_gr

  implicit none

  class( t_species_group ), intent(inout) :: this
  type( t_space ),      intent(in) :: g_space
  class( t_grid ),      intent(in) :: grid
  class( t_node_conf ),  intent(in) :: no_co
  type( t_time_step ),  intent(in) :: tstep
  real( p_double ) :: t

  if (test_if_report( tstep, this%diag%ndump_fac_nemit )) then
    call report_nemit_gr( this%electrons, g_space, grid, no_co, tstep, t )
    call report_nemit_gr( this%positrons, g_space, grid, no_co, tstep, t )

    if (this%cooling_type == p_cooling_qed) then
      call report_nemit_gr( this%photons, g_space, grid, no_co, tstep, t )
    endif
  endif

end subroutine report_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine cleanup_group( this )


  implicit none

  class( t_species_group ), intent(inout) :: this
  
  if (this % diag % ndump_fac_nemit > 0) then
    call this % electrons % nemit % cleanup()
    call this % positrons % nemit % cleanup()

    if (this%cooling_type == p_cooling_qed) then
      call this % photons % nemit % cleanup()
    endif
  endif

  deallocate( this % diag )

  if ( this % inj_limit == p_inj_limit_fracnGJpole .or. &
       this % inj_limit == p_inj_limit_fracnGJpolealt .or. &
       this % inj_limit == p_inj_limit_fracnGJlocal .or. &
       this % inj_limit == p_inj_limit_fracnGJlocalalt .or. &
       this % inj_limit == p_inj_limit_fracEBteo .or. &
       this % inj_limit == p_inj_limit_fracEBteolim .or. &
       this % inj_limit == p_inj_limit_sigma ) then
    call this % electrons % charge % cleanup()
    call this % positrons % charge % cleanup()
  endif

end subroutine cleanup_group
!-------------------------------------------------------------------------------

end module m_species_group_gr