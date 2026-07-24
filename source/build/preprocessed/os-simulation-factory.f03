# 1 "os-simulation-factory.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "os-simulation-factory.f03" 2
!
! Creating a new simluation type involves:
!
! 1. a 'use' statement to include the new code
! 2. a new entry in the 'case' statement within the 'create_core_simulation' function.
!
module m_simulation_factory

  use m_system
  use m_simulation


    use m_simulation_pgc



    use m_simulation_qed



    use m_simulation_cyl_modes



    use m_simulation_shear



    use m_simulation_rad



    use m_simulation_group



    use m_simulation_overdense



    use m_simulation_overdense_cyl



    use m_simulation_neutral_spin



    use m_simulation_xxfel



    use m_simulation_qedcyl



    use m_simulation_rad_cyl



    use m_simulation_gr



    use m_simulation_pkt






  implicit none

  private

  interface create
    module procedure create_simulation
  end interface

  public :: create, get_algorithm_from_string

  contains

!---------------------------------------------------------------------------------------------------
! Externally called routine that ultimatly passes back a t_simulation pointer.
!---------------------------------------------------------------------------------------------------
  function create_simulation( opts )
!---------------------------------------------------------------------------------------------------
    use m_input_file
    use m_parameters
    implicit none

    type( t_options ), intent(in) :: opts
    class( t_simulation ), pointer :: create_simulation

    class( t_simulation ), pointer :: sim

    sim => create_core_simulation( opts )

    create_simulation => sim
  end function
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Utility function (called from os-main.f03) that maps simulation type strings
! to the internally used integer parameters
!---------------------------------------------------------------------------------------------------
function get_algorithm_from_string(simulation_label)
  use m_parameters
  implicit none
  character(len=20) :: simulation_label
  integer :: get_algorithm_from_string

  select case( trim( simulation_label ) )
    case( "standard" )
      get_algorithm_from_string = p_standard
    case( "pgc" )
      get_algorithm_from_string = p_sim_pgc
    case( "qed" )
      get_algorithm_from_string = p_sim_qed
    case( "shear" )
      get_algorithm_from_string = p_sim_shear
    case( "quasi-3D" )
      get_algorithm_from_string = p_sim_cyl_modes
    case( "radiative", "radiat" )
      get_algorithm_from_string = p_sim_rad
    case( "tiles" )
      get_algorithm_from_string = p_sim_tiles
    case( "split-damp", "overdense" )
      get_algorithm_from_string = p_sim_overdense
    case( "split-damp-q3d", "overdense-cyl" )
      get_algorithm_from_string = p_sim_overdense_cyl
    case( "neutral-spin" )
      get_algorithm_from_string = p_sim_neut_spin
    case( "exact-ext-emf", "xxfel" )
      get_algorithm_from_string = p_sim_xxfel
    case( "qed-q3d", "qed-cyl" )
      get_algorithm_from_string = p_sim_qed_cyl
    case( "radiative-q3d", "radiat-cyl" )
      get_algorithm_from_string = p_sim_rad_cyl
    case( "gr" )
      get_algorithm_from_string = p_sim_gr
    case( "pkt" )
      get_algorithm_from_string = p_sim_pkt
    case( "cuda" )
      get_algorithm_from_string = p_sim_cuda
    case default
      call error_and_die( "Error reading global simulation parameters"//NEW_LINE('A')//&
                         &"algorithm: \'"//trim( simulation_label )//"\' is invalid or&
                         & the code is not compiled to support it. Aborting.")
  end select
end function

!---------------------------------------------------------------------------------------------------
! Core that creates the basic simulation types.
! This is the code that used to be in 'os-main%initialize'
!---------------------------------------------------------------------------------------------------
function create_core_simulation( opts )
!---------------------------------------------------------------------------------------------------
  use m_input_file
  use m_parameters
  implicit none

  type( t_options ), intent(in) :: opts
  class( t_simulation ), pointer :: create_core_simulation

  class( t_simulation ), pointer :: sim

    ! Allocate simulation object based on algorithm
  select case ( opts%algorithm )
    case( p_standard )
      ! Standard EM-PIC simulation
      allocate( t_simulation :: sim )
    case( p_sim_pgc )
      ! Ponderomotive guiding center simulation

        allocate( t_simulation_pgc :: sim )



    case( p_sim_qed )
      ! QED simulation

        allocate( t_simulation_qed :: sim )



    case( p_sim_shear )
      ! shear framework

        allocate( t_simulation_shear :: sim )



    case( p_sim_cyl_modes )
      ! quasi-3D (cylindrical modes) simulation

        allocate( t_simulation_cyl_modes :: sim )



    case( p_sim_rad )

        allocate( t_simulation_rad:: sim )




    case( p_sim_tiles )

        allocate( t_simulation_group :: sim )




    case( p_sim_overdense )

        allocate( t_simulation_overdense :: sim )




    case( p_sim_overdense_cyl )

        allocate( t_simulation_overdense_cyl :: sim )




    case( p_sim_neut_spin )

        allocate( t_simulation_neutral_spin:: sim )




    case( p_sim_xxfel )

        allocate( t_simulation_xxfel :: sim )




    case( p_sim_qed_cyl )

        allocate( t_simulation_qedcyl :: sim )




    case( p_sim_rad_cyl )

        allocate( t_simulation_rad_cyl:: sim )




    case( p_sim_gr )

        allocate( t_simulation_gr :: sim )




    case( p_sim_pkt )

        allocate( t_simulation_pkt :: sim )




    case( p_sim_cuda )



        call error_and_die('(*error*) You requested a simulation using the cuda &
                            & algorithm but the code is not compiled with cuda support')

    case default
      call error_and_die( '(*error*) Invalid simulation mode, only standard, pgc, qed, shear, quasi-3D,&
                          & radiative, tiles, overdense, overdense-q3d, neutral-spin, xxfel, qed-q3d, gr,&
                          & radiative-q3d and cuda are currently implemented.')
  end select
  create_core_simulation => sim
end function
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Utility function to print out fatal error conditions.
!---------------------------------------------------------------------------------------------------
subroutine error_and_die( message )
!---------------------------------------------------------------------------------------------------
    use m_parameters
    implicit none
    character(len=*) :: message

    if ( mpi_node() == 0 ) then
      write(0,*) message
      call abort_program( p_err_invalid )
    endif

  end subroutine
!---------------------------------------------------------------------------------------------------
end module


!---------------------------------------------------------------------------------------------------
! Version of the simulation factory that rests outside a module
! to avoid circular depencies in cases where a simulation wants to,
! inside it's own code, use the factory to create other simulation objects
! ( an example of this is the 'tandem' simulation class )
!---------------------------------------------------------------------------------------------------
! function global_simulation_factory( opts )

! use m_system, only: t_options
! use m_simulation, only: t_simulation
! use m_simulation_factory

! implicit none

! !integer, intent(in) :: algorithm
! type( t_options ), intent(inout) :: opts

! class( t_simulation ), pointer :: global_simulation_factory


! global_simulation_factory => create( opts )

! end function
