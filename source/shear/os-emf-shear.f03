!  update_shear_part - updates grad |a|^2 (n), grad |a|^2 (n+1/2), |a|^2 (n-1/2)
!  Must be called after advance a, before calling this routine

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_shear

#include "memory/memory.h"

use m_system
use m_parameters

use m_emf_define, only : t_emf
use m_vdf_define, only : t_vdf
use m_vdf_memory, only : alloc, freemem
use m_tridiag,    only : t_tridiag
use m_vdf_smooth, only : t_smooth

use m_logprof

implicit none

! IBM XL compilers
! This has to be explicitly made public so that you can access the superclass of t_emf_shear
public :: t_emf

private

! string to id restart data
character(len=*), parameter :: p_emf_shear_rst_id = "emf_shear rst data - 0x0001"

! parameter for timing event

!-------------------------------------------------------------------------------
! t_emf_shear class definition
!-------------------------------------------------------------------------------

type, extends( t_emf ) :: t_emf_shear

! include new type structure for SHEAR algorithm
   real( p_double ) :: alpha

contains

    procedure :: read_input => read_input_shear
    procedure :: init       => init_shear
    procedure :: advance    => advance_shear

    procedure :: allocate_objs => allocate_objs_shear

end type t_emf_shear

! exported symbols
public :: t_emf_shear

!-------------------------------------------------------------------------------
contains

!-----------------------------------------------------------------------------------------
subroutine allocate_objs_shear( this )
!-----------------------------------------------------------------------------------------
!       Allocate any objects contained within the EMF (t_emf) object
!-----------------------------------------------------------------------------------------

    implicit none

    class( t_emf_shear ), intent(inout) :: this

    ! Allocate superclass objects
    call this % t_emf % allocate_objs()

end subroutine

!-------------------------------------------------------------------------------
! Read input file data for SHEAR algorithm
!   - this routine also calls the superclass (t_emf) read_nml routine to read
!     that section first
!-------------------------------------------------------------------------------
subroutine read_input_shear( this, input_file, periodic, if_move, grid, dx, dt, gamma )

  use m_input_file
  use m_grid_define
  use m_vdf_smooth, only : read_nml

  implicit none

  class( t_emf_shear ), intent( inout )  :: this
  class( t_input_file ), intent(inout)  :: input_file
  logical, dimension(:), intent(in)    :: periodic, if_move
  class( t_grid ), intent(in)           :: grid
  real(p_double), dimension(:), intent(in) :: dx
  real(p_double), intent(in) :: dt
  real(p_double), intent(in) :: gamma 

  ! namelist parameters
  real( p_double )             :: alpha

  integer :: ierr

  namelist /nl_shear/ alpha

  ! -------------- Allocate any sub-objects needed --------------------
  call this%allocate_objs()

  ! first read the superclass (emf) section
  call this%t_emf%read_input(input_file, periodic, if_move, grid, dx, dt, gamma)

  ! Now read this class input section

  !> default values <
  ! general shear parameter
  alpha          = 0.0

  ! read values from file
  call get_namelist( input_file, "nl_shear", ierr )

  if (ierr == 0) then
    read (input_file%nml_text, nml = nl_shear, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading shear parameters"
      print *, "aborting..."
      stop
    endif
  else
    SCR_ROOT(" - no shear parameters specified")
    stop
  endif

  ! store namelist parameters
  this%alpha       = alpha

end subroutine read_input_shear
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Initialize emf_shear object
!   - this routine also calls the superclass (t_emf) setup routine first
!-------------------------------------------------------------------------------
subroutine init_shear( this, part_grid_center, part_interpolation, &
                     g_space, grid, gc_min, dx, tstep, tmin, tmax, no_co, &
                     send_msg, recv_msg, &
                     restart,  restart_handle, sim_options )

  use m_space, only : t_space, xmin
  use m_node_conf, only : t_node_conf
  use m_restart, only : t_restart_handle
  use m_grid_define, only : t_grid
  use m_vdf_smooth, only : if_smooth, setup
  use m_vdf_comm, only : t_vdf_msg
  use m_time_step, only : t_time_step

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  class( t_emf_shear ), intent( inout ), target  ::  this

  logical, intent(in) :: part_grid_center
  integer, intent(in) :: part_interpolation
  type( t_space ),     intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  integer, dimension(:,:), intent(in) :: gc_min
  real( p_double ), dimension(:), intent(in) :: dx
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax
  class( t_node_conf ), intent(in), target :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type( t_options ), intent(in) :: sim_options

  real(p_k_fld), dimension(p_max_dim) ::  g_xmin

  ! setup superclass data
  ! (this will also setup the diagnostics)
  call this % t_emf % init( part_grid_center, part_interpolation, &
                            g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                            no_co, send_msg, recv_msg, &
                            &restart, restart_handle, sim_options )

  !executable statements
  g_xmin(1:p_x_dim) = xmin(g_space)

end subroutine init_shear
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Advance EMF fields for SHEAR algorithm
!-------------------------------------------------------------------------------
subroutine advance_shear( this, jay, dt, send_msg, recv_msg )

  use m_node_conf, only : t_node_conf
  use m_vdf_comm, only : update_boundary, t_vdf_msg
  use m_current_define, only : t_current

  use m_vdf_smooth, only : smooth

  implicit none

  class( t_emf_shear ), intent( inout )  ::  this
  class( t_current ),   intent(inout) :: jay
  real(p_double),  intent(in) :: dt
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  !local variables

  real(p_double) :: dt_b, dt_e, dt_e2

  type( t_vdf ) :: etmp, btmp

  dt_b = dt / 2.0_p_double
  dt_e = dt

  dt_e2 = dt / 2.0_p_double

  select case( p_x_dim )
  case( 1 )

    ! Advance B half time step
    call dbdt_RK2_1d( this%b, this%e, dt_b, this%alpha )

  case( 2 )

    ! Advance B half time step
    call dbdt_RK2_2d( this%b, this%e, dt_b, this%alpha )

  case( 3 )

    ! Advance B half time step
    call dbdt_shear_3d( this%b, this%e, dt_b )

  case default

    print *, "error: shear B1 solver not loaded!"
    print *, "aborting ..."
    stop

  end select

  ! Save a temporary E and B field
  call etmp %new( this%e, copy=.true., f_dim = 3)
  call btmp %new( this%b, copy=.true., f_dim = 3)

  select case( p_x_dim )
  case( 1 )

    ! Advance E half time step for the correction term
    call dedt_RK2_1d( btmp, etmp, jay % pf(1), dt_e2, this%alpha )

    ! Advance E one full time step
    call dedt_shear_1d( this%e, etmp, this%b, jay % pf(1), dt, this%alpha )

  case( 2 )

    ! Advance E half time step for the correction term
    call dedt_RK2_2d( btmp, etmp, jay % pf(1), dt_e2, this%alpha )

    ! Advance E one full time step
    call dedt_shear_2d( this%e, etmp, this%b, jay % pf(1), dt, this%alpha )

  case( 3 )

    ! Advance E one full time step
    call dedt_shear_3d( this%e, this%b, jay % pf(1), dt )

  case default

    print *, "error: shear E solver not loaded!"
    print *, "aborting ..."
    stop

  end select

  call btmp%cleanup()
  call etmp%cleanup()

  select case( p_x_dim )
  case( 1 )

    ! Advance B another half time step
    call dbdt_RK2_1d( this%b, this%e, dt_b, this%alpha )

  case( 2 )

    ! Advance B another half time step
    call dbdt_RK2_2d( this%b, this%e, dt_b, this%alpha )

  case( 3 )

    ! Advance B another half time step
    call dbdt_shear_3d( this%b, this%e, dt_b )

  case default

    print *, "error: shear B2 solver not loaded!"
    print *, "aborting ..."
    stop

  end select
  ! TODO: remove: this%n_current_emf = this%n_current_emf + 1
  
end subroutine advance_shear
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dbdt_RK2_1d( b, e, dt, alpha )
!-------------------------------------------------------------------------------
!    advance time step magnetic field using RK method at the second order
!-------------------------------------------------------------------------------

implicit none
  integer, parameter :: rank = 1

!       dummy variables
  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt, alpha

!       local variables
  real(p_k_fld) :: dtdx1, dx1
  integer :: i1
  type(t_vdf) :: k1b, k2b

!       executable statements
  call k1b%new( b, copy = .false., f_dim = 3)
  call k2b%new( b, copy = .false., f_dim = 3)

  dtdx1 = real(dt/b%dx_(1), p_k_fld) !*c*c !c=1
  dx1 = real(1/b%dx_(1), p_k_fld)

  do i1 = 0, b%nx_(1)+1
    k1b%f1(1, i1) = 0.0
    k1b%f1(2, i1) = 0.0
    k1b%f1(3, i1) = - 1.5_p_double * alpha * b%f1( 2, i1 )

  enddo

  do i1 = 0, b%nx_(1)+1
    k2b%f1(1, i1) = 0.0
    k2b%f1(2, i1) = 0.0
    k2b%f1(3, i1) =  - 1.5_p_double * alpha * b%f1( 2, i1 ) &
                  - 1.5_p_double * alpha * (2.0_p_double/3.0_p_double) &
                  * dt * k1b%f1(3, i1)

  enddo

  do i1 = 0, b%nx_(1)+1

    b%f1(2, i1) = b%f1(2, i1) + dtdx1 * ( e%f1( 3, i1+1 ) - e%f1( 3, i1 ))
    b%f1(3, i1) = b%f1(3, i1) - dtdx1 * ( e%f1( 2, i1+1 ) - e%f1( 2, i1 )) &

             + 0.5d0*(k1b%f1(3, i1) - 1.5_p_double * alpha * b%f1( 2, i1 )) * dt
  enddo

  call k2b%cleanup()
  call k1b%cleanup()

end subroutine dbdt_RK2_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dbdt_RK2_2d( b, e, dt, alpha )
!-------------------------------------------------------------------------------
!    advance time step magnetic field using RK method
!    at the second order
!    the framework is rotated in order to have the normal poloidal plane
!-------------------------------------------------------------------------------

implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt, alpha

!       local variables
  real(p_k_fld) :: dtdx1, dtdx2, dx1, dx2, s1, s2, dt2
  integer :: i1, i2
  type(t_vdf) :: k1b, k2b

!       executable statements
  call k1b % new( b, copy = .false., f_dim = 3)
  call k2b % new( b, copy = .false., f_dim = 3)

  dtdx1 = real(dt/b%dx_(1), p_k_fld) !*c*c !c=1
  dtdx2 = real(dt/b%dx_(2), p_k_fld) !*c*c !c=1
  dx1 = real(1/b%dx_(1), p_k_fld)

  dx2 = real(1/b%dx_(2), p_k_fld)
  s1 = real(1.5_p_double * alpha * 0.5d0, p_k_fld)
  s2 = real(1.5_p_double * alpha * (2.0_p_double/3.0_p_double) * dt, p_k_fld)
  dt2 = real(0.5d0 * dt, p_k_fld)

  do i2 = 0, b%nx_(2)+1

    do i1 = 0, b%nx_(1)+1

    k1b%f2(1, i1, i2) = 0.0

    k1b%f2(2, i1, i2) = 0.0

    k1b%f2(3, i1, i2) = - s1 * (b%f2( 1, i1, i2 ) + b%f2( 1, i1, i2+1 ))

    enddo
  enddo

  do i2 = 0, b%nx_(2)+1

    do i1 = 0, b%nx_(1)+1

    k2b%f2(1, i1, i2) = 0.0

    k2b%f2(2, i1, i2) = 0.0

    k2b%f2(3, i1, i2) = - s1 * (b%f2( 1, i1, i2 ) + b%f2( 1, i1, i2+1 )) &
                   - s2 * k1b%f2(3, i1, i2)

    enddo
  enddo

  ! write( *, * ) ">> k2b2 = ", k2b%f2(2, 0, 0)

  do i2 = 0, b%nx_(2)+1

    do i1 = 0, b%nx_(1)+1

    b%f2(1, i1, i2) = b%f2(1, i1, i2) - dtdx2 * ( e%f2( 3, i1, i2+1 ) - e%f2( 3, i1, i2 ))
    b%f2(2, i1, i2) = b%f2(2, i1, i2) + dtdx1 * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1, i2 ))
    b%f2(3, i1, i2) = b%f2(3, i1, i2) - dtdx1 * ( e%f2( 2, i1+1, i2 ) - e%f2( 2, i1, i2 )) &
                + dtdx2 * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 )) &
                + dt2 * (k1b%f2( 3, i1, i2 ) &
                - s1 * (b%f2( 1, i1, i2+1 ) + b%f2( 1, i1, i2 )))

    enddo
  enddo

    call k2b%cleanup()
    call k1b%cleanup()

end subroutine dbdt_RK2_2d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dedt_RK2_1d( b, e, jay, dt, alpha )
!-------------------------------------------------------------------------------
!    advance time step electric field using RK method at the second order
!-------------------------------------------------------------------------------

implicit none
  integer, parameter :: rank = 1

!       dummy variables
  type( t_vdf ), intent( in ) :: b, jay
  type( t_vdf ),    intent(inout) :: e
  real(p_double),               intent(in) :: dt, alpha

!       local variables
  real(p_k_fld) :: dtdx1, dx1, dtif
  integer :: i1
  type(t_vdf) :: k1e, k2e

!       executable statements

  dtdx1 = real(dt/e%dx_(1), p_k_fld) !*c*c !c=1
  dx1 = real(1/e%dx_(1), p_k_fld)

  dtif = real( dt, p_k_fld )

  call k1e % new( b, copy = .false., f_dim = 3)
  call k2e % new( b, copy = .false., f_dim = 3)

  do i1 = 0, e%nx_(1)+1
    k1e%f1(1, i1) = 0.0
    k1e%f1(2, i1) = 0.0
    k1e%f1(3, i1) = - 1.5_p_double * alpha * e%f1( 2, i1 )
  enddo

  do i1 = 0, e%nx_(1)+1
    k2e%f1(1, i1) = 0.0
    k2e%f1(2, i1) = 0.0
    k2e%f1(3, i1) = - 1.5_p_double * alpha * e%f1( 2, i1 ) &
                 - 1.5_p_double * alpha * ((2.0_p_double/3.0_p_double) &
                 * dt * k1e%f1(3, i1))
  enddo

  do i1 = 1, e%nx_(1)+1

    e%f1(1, i1) = e%f1(1, i1) - dtif * jay%f1(1, i1)
    e%f1(2, i1) = e%f1(2, i1) - dtdx1 * ( b%f1(3, i1) - b%f1(3, i1-1) ) &

              - dtif * jay%f1(2, i1)
    e%f1(3, i1) = e%f1(3, i1) + dtdx1 * ( b%f1(2, i1) - b%f1(2, i1-1) ) &

              + 0.5d0*(k1e%f1(3, i1)- 1.5_p_double * alpha * e%f1( 2, i1 )) * dt &
              - dtif * jay%f1(3, i1)
  enddo

  call k2e%cleanup()
  call k1e%cleanup()

end subroutine dedt_RK2_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dedt_RK2_2d( b, e, jay, dt, alpha )
!-------------------------------------------------------------------------------
!    advance time step electric field using RK method
!    at the second order
!-------------------------------------------------------------------------------

implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent( in ) :: b, jay
  type( t_vdf ),    intent(inout) :: e
  real(p_double),               intent(in) :: dt, alpha

!       local variables
  real(p_k_fld) :: dtdx1, dtdx2, dx1, dx2, dtif, s1, s2, dt2
  integer :: i1, i2
  type(t_vdf) :: k1e, k2e

!       executable statements

  dtdx1 = real(dt/e%dx_(1), p_k_fld) !*c*c !c=1
  dtdx2 = real(dt/e%dx_(2), p_k_fld) !*c*c !c=1
  dx1 = real(1/e%dx_(1), p_k_fld)

  dx2 = real(1/e%dx_(2), p_k_fld)
  dtif = real( dt, p_k_fld )
  s1 = real(1.5_p_double * alpha * 0.5d0, p_k_fld)
  s2 = real(1.5_p_double * alpha * (2.0_p_double/3.0_p_double) * dt, p_k_fld)
  dt2 = real(0.5d0 * dt, p_k_fld)

  call k1e % new( b, copy = .false., f_dim = 3)
  call k2e % new( b, copy = .false., f_dim = 3)

  do i2 = 1, e%nx_(2)+1

    do i1 = 1, e%nx_(1)+1

    k1e%f2(1, i1, i2) = 0.0

    k1e%f2(2, i1, i2) = 0.0

    k1e%f2(3, i1, i2) = - s1 * (e%f2( 1, i1, i2 ) + e%f2( 1, i1-1, i2))

    enddo
  enddo

  ! write( *, * ) ">> k1e3 = ", k1e%f2(3, 0, 0)

  do i2 = 1, e%nx_(2)+1

    do i1 = 1, e%nx_(1)+1

    k2e%f2(1, i1, i2) = 0.0

    k2e%f2(2, i1, i2) = 0.0

    k2e%f2(3, i1, i2) = - s1 * (e%f2( 1, i1, i2 ) + e%f2( 1, i1-1, i2)) &
                   - s2 * k1e%f2(3, i1, i2)

    enddo
  enddo

  do i2 = 1, e%nx_(2)+1

    do i1 = 1, e%nx_(1)+1

    e%f2(1, i1, i2) = e%f2(1, i1, i2) + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) ) &
              - dtif * jay%f2(1, i1, i2)
    e%f2(2, i1, i2) = e%f2(2, i1, i2) - dtdx1 * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) ) &
                - dtif * jay%f2(2, i1, i2)
    e%f2(3, i1, i2) = e%f2(3, i1, i2) + dtdx1 * ( b%f2(2, i1, i2) - b%f2(2, i1-1, i2) ) &
              - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1) ) &
              + dt2 * ( k1e%f2(3, i1, i2) - s1 &
              * (e%f2( 1, i1, i2 ) + e%f2( 1, i1, i2-1 ))) &
              - dtif * jay%f2(3, i1, i2)

    enddo
  enddo

    call k2e%cleanup()
    call k1e%cleanup()

end subroutine dedt_RK2_2d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dedt_shear_1d( e, etmp, b, jay, dt, alpha )
!-------------------------------------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!-------------------------------------------------------------------------------

implicit none
  integer, parameter :: rank = 1
!       dummy variables
  type( t_vdf ),    intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay, etmp
  real(p_double),   intent(in) :: dt, alpha

!       local variables
  real(p_k_fld) :: dtdx1, dtif
  integer :: i1
  real(p_k_fld) :: dump_coef

!       executable statements
  dtdx1 = real( dt/e%dx_(1), p_k_fld )
  dtif = real( dt, p_k_fld )

  dump_coef = real(1.5_p_double * alpha * dt )

  !$omp parallel do
  do i1 = 1, e%nx_(1)+1
    ! E1
    e%f1(1, i1) = e%f1(1, i1) &
            - dtif * jay%f1(1, i1)

    ! E2
    e%f1(2, i1) = e%f1(2, i1) &
                      - dtif * jay%f1(2, i1) &
                      - dtdx1 * ( b%f1(3, i1) - b%f1(3, i1-1) )

    ! E3
    e%f1(3, i1) = e%f1(3, i1) &
                      - dtif * jay%f1(3, i1) &
                      + dtdx1 * ( b%f1(2, i1) - b%f1(2, i1-1)) &
                      - dump_coef * etmp%f1( 2, i1 )
  enddo
  !$omp end parallel do

end subroutine dedt_shear_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dedt_shear_2d( e, etemp, b, jay, dt, alpha )
!-------------------------------------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!       This routine will calculate
!       E1( 1:nx(1)  , 1:nx(2)+1 )
!       E2( 1:nx(1)+1, 1:nx(2)   )
!       E3( 1:nx(1)+1, 1:nx(2)+1 )
!
!   MODIFIED FOR SHEAR TERMS
!-------------------------------------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay, etemp

  real(p_double),               intent(in) :: dt, alpha

!       local variables
!  integer :: nx1m, nx2m, nx1b, nx2b
  integer :: i1, i2
  real(p_k_fld) :: dtdx1, dtdx2, dtif
  real(p_k_fld) :: dump_coef

!       executable statements
!  nx1m = nx(e,1)
!  nx2m = nx(e,2)
!  nx1b = nx1m+1
!  nx2b = nx2m+1

  dtdx1 = real( dt/e%dx_(1), p_k_fld )
  dtdx2 = real( dt/e%dx_(2), p_k_fld )
  dtif = real( dt, p_k_fld )

  dump_coef = real(1.5_p_double * alpha * dt / 2_p_double)

  ! version 1, advance 1 cell at a time

  !$omp parallel do private(i1)
  do i2 = 1, e%nx_(2)+1
    do i1 = 1, e%nx_(1)+1
    ! E1
    e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                      - dtif * jay%f2(1, i1, i2) &
                      + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )

    ! E2
    e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                      - dtif * jay%f2(2, i1, i2) &
                      - dtdx1 * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) )

    ! E3
    e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                      - dtif * jay%f2(3, i1, i2) &
                      + dtdx1 * ( b%f2(2, i1, i2) - b%f2(2, i1-1, i2)) &
                      - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1)) &
                      - dump_coef * (etemp%f2( 1, i1, i2 ) + etemp%f2( 1, i1-1, i2))
  enddo
  enddo
  !$omp end parallel do

end subroutine dedt_shear_2d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dbdt_shear_3d( b, e, dt )
!-------------------------------------------------------------------------------
!       advances the magnetic field
!-------------------------------------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

  ! dummy variables
  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt

  ! local variables
  real(p_k_fld) ::  dtdx1,  dtdx2,  dtdx3
  integer :: i1, i2, i3

  ! executable statements

  dtdx1 = real( dt/b%dx_(1), p_k_fld )

  dtdx2 = real( dt/b%dx_(2), p_k_fld )
  dtdx3 = real( dt/b%dx_(3), p_k_fld )

  ! version 1, advance 1 cell at a time
  ! loop through i3 first

  !$omp parallel do private(i2,i1)
  do i3 = 0, b%nx_(3)+1

    do i2 = 0, b%nx_(2)+1

      do i1 = 0, b%nx_(1)+1

        !B1
        b%f3( 1, i1, i2, i3 ) = b%f3( 1, i1, i2, i3 ) &
                            - dtdx2 * ( e%f3( 3, i1, i2+1, i3 ) - e%f3( 3, i1, i2, i3 )) &
                            + dtdx3 * ( e%f3( 2, i1, i2, i3+1 ) - e%f3( 2, i1, i2, i3 ))

        !B2
        b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, i2, i3 ) &
                            + dtdx1 * ( e%f3( 3, i1+1, i2, i3 ) - e%f3( 3, i1, i2, i3 )) &
                            - dtdx3 * ( e%f3( 1, i1, i2, i3+1 ) - e%f3( 1, i1, i2, i3 ))
        !B3
        b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, i3 ) &
                            - dtdx1 * ( e%f3( 2, i1+1, i2, i3 ) - e%f3( 2, i1, i2, i3 )) &
                            + dtdx2 * ( e%f3( 1, i1, i2+1, i3 ) - e%f3( 1, i1, i2, i3 ))

      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_shear_3d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dedt_shear_3d( e, b, jay, dt )
!-------------------------------------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!-------------------------------------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

!       dummy variables

  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

  ! local variables
  real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
  integer :: i1, i2, i3

  ! executable statements

  dtdx1 = real( dt/e%dx_(1), p_k_fld )
  dtdx2 = real( dt/e%dx_(2), p_k_fld )
  dtdx3 = real( dt/e%dx_(3), p_k_fld )
  dtif = real( dt, p_k_fld )

  !$omp parallel do private(i2,i1)
  do i3 = 1, e%nx_(3)+1
    do i2 = 1, e%nx_(2)+1
      do i1 = 1, e%nx_(1)+1

      e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
        - dtif  * jay%f3( 1, i1, i2, i3 )   &
        + dtdx2 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1, i2-1, i3 ) ) &
        - dtdx3 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1, i2, i3-1 ) )

      e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
        - dtif  * jay%f3( 2, i1, i2, i3 ) &
        + dtdx3 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2, i3-1 ) ) &
        - dtdx1 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1-1, i2, i3 ) )

      e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
        - dtif  * jay%f3( 3, i1, i2, i3 )   &
        + dtdx1 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1-1, i2, i3 ) ) &
        - dtdx2 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2-1, i3 ) )

      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine dedt_shear_3d
!-------------------------------------------------------------------------------

end module m_emf_shear
