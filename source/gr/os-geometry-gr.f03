#include "os-config.h"
#include "os-preprocess.fpp"

module m_geometry_gr

#include "memory/memory.h"

use m_system
use m_parameters

use m_input_file,  only : t_input_file

use m_vdf_define, only : t_vdf
use m_math, only : pi, pi_2

implicit none

private

integer, parameter :: p_rc_dim = 1, &
                      p_tc_dim = 2, &
                      p_pc_dim = 3

! type of metric
integer, parameter :: p_geometry_minkowski     = 0, &
                      p_geometry_schwarzschild = 1, &
                      p_geometry_kerr_slow     = 2, &
                      p_geometry_kerr          = 3

! type of b profile
integer, parameter :: p_bprofile_none = 0, &
                      p_bprofile_monopole = 1, &
                      p_bprofile_dipole   = 2

! type of spline rule
integer, parameter :: p_spline_rectangular = 0, &
                      p_spline_integral    = 1

! string to id restart data
character(len=*), parameter, public :: p_geometry_gr_rst_id = "geometry_gr rst data - 0x0006"                 
                      
!-------------------------------------------------------------------------------
! t_geometry_gr class definition
!-------------------------------------------------------------------------------

type :: t_geometry_gr

  ! type of metric to use
  integer :: metric

  ! grid positions
  type( t_vdf ) :: r
  type( t_vdf ) :: dr
  type( t_vdf ) :: r3

  type( t_vdf ) :: t
  type( t_vdf ) :: st
  type( t_vdf ) :: ct

  real( p_double ) :: dt

  ! CFL condition dt
  real( p_double ) :: dtcfl 

  ! line and area elements
  type( t_vdf ) :: lr
  type( t_vdf ) :: lt
  type( t_vdf ) :: lp
  type( t_vdf ) :: ar
  type( t_vdf ) :: at
  type( t_vdf ) :: ap

  ! volume elements
  type( t_vdf ) :: vr
  type( t_vdf ) :: vt

  ! grid/normalization constants
  real(p_k_part) :: delta
  real(p_k_part) :: A
  real(p_k_part) :: nr
  real(p_k_part) :: rmax

  ! metric dependent quantities
  real(p_double) :: rs
  real(p_double) :: beta0

  type( t_vdf ) :: alpha
  type( t_vdf ) :: beta

  ! grid polar axes flags
  logical :: north_axis
  logical :: south_axis
  
  ! physical quantities
  integer :: bprofile
  real(p_double) :: bs
  real(p_double) :: omega ! time evolving amplitude
  real(p_double) :: omega0 ! final amplitude
  real(p_double) :: trise

  ! classical rad. reaction 
  logical :: classical_rr
  real(p_double) :: b_schwinger, rr_norm

  ! type of spline rule to use
  integer :: spline_rule

  ! gravitational force switch
  real(p_double) :: grav_switch

contains

  procedure :: read_input    => read_input_gr
  procedure :: init          => init_gr
  procedure :: cleanup       => cleanup_gr
  procedure :: advance       => advance_gr

  procedure :: spline_r_rectangular    => spline_r_rect_gr
  procedure :: spline_r_integral       => spline_r_int_gr
  procedure :: spline_r_rectangular_gr => spline_r_rect_curved_gr
  procedure :: spline_r_integral_gr    => spline_r_int_curved_gr
  procedure :: spline_t                => spline_t_gr
  procedure :: spline_rh               => spline_rh_gr
  procedure :: spline_th               => spline_th_gr
  
  procedure :: write_checkpoint => write_checkpoint_gr
  procedure :: restart_read     => restart_read_gr

end type t_geometry_gr

! exported symbols
public :: t_geometry_gr
public :: p_rc_dim, p_tc_dim, p_pc_dim
public :: p_geometry_minkowski, p_geometry_schwarzschild
public :: p_geometry_kerr_slow, p_geometry_kerr
public :: p_bprofile_none, p_bprofile_monopole, p_bprofile_dipole
public :: p_spline_rectangular, p_spline_integral

contains

!-------------------------------------------------------------------------------
! Read input file data for gr geometry object
!-------------------------------------------------------------------------------
subroutine read_input_gr( this, input_file )

  use m_input_file

  implicit none

  class( t_geometry_gr ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  ! local variables
  real( p_double ) :: rs, omega, omega0, bs, trise
  character(20) :: metric, bprofile, spline_rule
  logical :: classical_rr, gravity
  real( p_double ) :: dtcfl, b_schwinger

  namelist /nl_geometry/ metric, rs, omega, trise, bprofile, bs, &
    classical_rr, b_schwinger, spline_rule, gravity

  integer :: ierr

  ! executable statements
  metric = "minkowski"
  rs = 0
  omega0 = 0.2
  bprofile = "none"
  bs = 100
  trise = 10

  ! CFL condition
  dtcfl = 1.0

  ! classical Rad. Reaction
  classical_rr = .false.
  b_schwinger = 1.0d6

  ! gravity switch
  gravity = .true.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_geometry", ierr )

  if ( ierr /= 0 ) then
    if (ierr < 0) then
      print *, "Error reading geometry parameters"
    else
      print *, "Error: geometry parameters missing"
    endif
    print *, "aborting..."
    stop
  endif


  read (input_file%nml_text, nml = nl_geometry, iostat = ierr)
  if (ierr /= 0) then
    print *, "Error reading geometry parameters"
    print *, "aborting..."
    stop
  endif

  select case ( trim( metric ) )
  case ( "minkowski" )
    this%metric = p_geometry_minkowski
  case ( "schwarzschild" )
    this%metric = p_geometry_schwarzschild
  case ( "kerr_slow" )
    this%metric = p_geometry_kerr_slow
  case ( "kerr" )
    this%metric = p_geometry_kerr
  case default
    this%metric = p_geometry_minkowski
  end select

  this%rs = rs
  this%omega0 = omega
  this%omega = 0
  this%beta0 = 0
  this%trise = trise

  if (this%metric /= p_geometry_minkowski .and. this%rs < 0.0001 ) then
    print *, "Error reading geometry parameters"
    print *, "please select Minkowski metric"
    print *, "if rs < 0.0001"
    print *, "aborting..."
    stop
  endif

  if (this%metric == p_geometry_minkowski .and. this%rs /= 0. ) then
    print *, "Error reading geometry parameters"
    print *, "please select non-Minkowski metric"
    print *, "if rs is not zero!"
    print *, "aborting..."
    stop
  endif

  select case ( trim( bprofile ) )
  case ( "monopole" )
    this%bprofile = p_bprofile_monopole
  case ( "dipole" )
    this%bprofile = p_bprofile_dipole
  case default
    this%bprofile = p_bprofile_none
  end select

  this%bs = bs

  ! Reads if classical radiation reaction is used
  this%classical_rr = classical_rr
  this%b_schwinger = b_schwinger
  ! evaluates the norm factor for radiation reaction (used for GR-OSIRIS module)
  ! 2./3./alpha , with alpha = 1/137.
  this%rr_norm = 0.004866180049 / this%b_schwinger

  select case ( trim( spline_rule ) )
  case ( "rectangular" )
    this%spline_rule = p_spline_rectangular
  case ( "integral" )
    this%spline_rule = p_spline_integral
  case default
    this%spline_rule = p_spline_integral
  end select

  ! gravitational force switch
  select case ( gravity )
  case ( .true. )
    this%grav_switch = 1.0
  case ( .false. )
    this%grav_switch = 0.0
  end select

end subroutine

!-------------------------------------------------------------------------------
! Initialize gr geometry object
!-------------------------------------------------------------------------------
subroutine init_gr( this, g_space, grid )

  use m_space, only : t_space
  use m_grid_define, only : t_grid

  implicit none

  class( t_geometry_gr ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  type( t_grid ), intent(in) :: grid

  ! executable statements
  this % north_axis = grid % my_nx( p_lower, p_tc_dim ) - 2 < 0
  this % south_axis = grid % g_nx( p_tc_dim ) - grid % my_nx( p_upper, p_tc_dim ) - 1 < 0

  select case (this % metric)
    case ( p_geometry_minkowski )
      call init_geometry_kerrslow( this, g_space, grid )

    case ( p_geometry_schwarzschild )
      call init_geometry_kerrslow( this, g_space, grid )

    case ( p_geometry_kerr_slow )
      call init_geometry_kerrslow( this, g_space, grid )

    case ( p_geometry_kerr )
      ! call init_metric_kerr( this, g_space, grid )
      ! not yet implemented
  end select
  
end subroutine

!-------------------------------------------------------------------------------
! Initializes slow rotation limit of Kerr metric
!-------------------------------------------------------------------------------
subroutine init_geometry_kerrslow( this, g_space, grid )

  use m_space, only : t_space, xmin, xmax
  use m_grid_define, only : t_grid

  implicit none

  ! dummy variables
  class( t_geometry_gr ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  type( t_grid ), intent(in) :: grid

  ! local variables
  real( p_double ), dimension(p_max_dim) ::  g_xmin, g_xmax
  real( p_double ) :: r_min, r_max, exp, delta, &
                     t_min, t_max, dt, t_eval, &
                     dp

  integer :: i1, i2, i1_tot, i2_tot
  real( p_double ), dimension(2) :: dx_vec

  integer, dimension(2, 1) :: gc_num_f1r, gc_num_f1t
  integer, dimension(2, 2) :: gc_num_f2

  real( p_double ) :: rs, beta0
  real( p_k_fld ), pointer :: r(:,:), dr(:,:), r3(:,:), t(:,:), st(:,:), ct(:,:), &
                              lr(:,:), lt(:,:), ap(:,:), vr(:,:), vt(:,:), alpha(:,:)
  real( p_k_fld ), pointer :: lp(:,:,:), ar(:,:,:), at(:,:,:), beta(:,:,:)

  ! executable statements

  ! get grid boundaries
  g_xmin(1:p_x_dim) = xmin(g_space)
  g_xmax(1:p_x_dim) = xmax(g_space)

  ! get grid properties
  r_min = g_xmin(p_rc_dim)
  r_max = g_xmax(p_rc_dim)
  
  this % nr = grid % g_nx( p_rc_dim )
  this % rmax = r_max
  exp = real(1._p_double / grid % g_nx( p_rc_dim ),p_double)
  delta = (r_max/r_min)**exp

  this % delta = delta
  this % A = real( ( 1._p_double + 1._p_double/delta ) / ( 1._p_double - 1._p_double/delta ), p_double )
  
  t_min = g_xmin(p_tc_dim)
  t_max = g_xmax(p_tc_dim)
  dt  = (t_max-t_min) / grid % g_nx( p_tc_dim )
  this % dt = dt 

  ! change this in 3D
  dp = 1._p_double

  ! dummy array to create vdfs
  dx_vec = (/0,0/)

  ! guard cells for radial f1s
  gc_num_f1r(p_lower, 1) = 2
  gc_num_f1r(p_upper, 1) = 2

  ! guard cells for polar f1s
  gc_num_f1t(p_lower, 1) = 1
  gc_num_f1t(p_upper, 1) = 2

  ! guard cells for f2s
  gc_num_f2(p_lower, 1) = 3 !2
  gc_num_f2(p_upper, 1) = 2
  gc_num_f2(p_lower, 2) = 2 !1
  gc_num_f2(p_upper, 2) = 2

  ! radial position (1: i , 2: i+1/2)
  call this%r%new( 1, 2, (/grid % my_nx(3, p_rc_dim )/), &
            gc_num_f1r + 1, (/dx_vec(1)/), .true. )

  ! radial increment (1: i)
  call this%dr%new( 1, 1, (/grid % my_nx(3, p_rc_dim )/), &
            gc_num_f1r + 1, (/dx_vec(1)/), .true. )

  ! cubed radial position (1: i , 2: i+1/2)
  call this%r3%new( 1, 2, (/grid % my_nx(3, p_rc_dim )/), &
            gc_num_f1r + 1, (/dx_vec(1)/), .true. )

  ! gr alpha lapse function (1: i, 2: i+1/2)
  call this%alpha%new( 1, 2, (/grid % my_nx(3, p_rc_dim )/), &
            gc_num_f1r + 1, (/dx_vec(1)/), .true. )

  ! gr beta shift function (1: (i+1/2,j), 2: (i,j+1/2))
  call this%beta%new( 2, 2, (/grid % my_nx(3, p_rc_dim ), &
            grid % my_nx(3, p_tc_dim )/), gc_num_f2, dx_vec, .true. )

  ! polar position (1: j, 2: j+1/2)
  call this%t%new( 1, 2, (/grid % my_nx(3, p_tc_dim )/), &
            gc_num_f1t + 1, (/dx_vec(2)/), .true. )

  ! sine at the polar position (1: j, 2: j+1/2)
  call this%st%new( 1, 2, (/grid % my_nx(3, p_tc_dim )/), &
            gc_num_f1t + 1, (/dx_vec(2)/), .true. )

  ! cosine at the polar position (1: j, 2: j+1/2)
  call this%ct%new( 1, 2, (/grid % my_nx(3, p_tc_dim )/), &
            gc_num_f1t + 1, (/dx_vec(2)/), .true. )


  ! r component of line coeficient (1: i, 2: i+1/2)
  call this%lr%new( 1, 2, (/grid % my_nx(3, p_rc_dim )/), &
            gc_num_f1r, (/dx_vec(1)/), .true. )

  ! theta component of line coeficient (1: i, 2: i+1/2)
  call this%lt%new( 1, 2, (/grid % my_nx(3, p_rc_dim )/), &
            gc_num_f1r, (/dx_vec(1)/), .true. )

  ! phi component of line coeficient (1: (i,j), 2: (i+1/2,j+1/2))
  call this%lp%new( 2, 2, (/grid % my_nx(3, p_rc_dim ), &
            grid % my_nx(3, p_tc_dim )/), gc_num_f2, dx_vec, .true. )

  ! r component of areal coeficient (1: (i+1/2,j), 2: (i,j+1/2))
  call this%ar%new( 2, 2, (/grid % my_nx(3, p_rc_dim ), &
            grid % my_nx(3, p_tc_dim )/), gc_num_f2, dx_vec, .true. )

  ! theta component of areal coeficient (1: (i+1/2,j), 2: (i,j+1/2))
  call this%at%new( 2, 2, (/grid % my_nx(3, p_rc_dim ), &
            grid % my_nx(3, p_tc_dim )/), gc_num_f2, dx_vec, .true. )

  ! phi component of areal coeficient (1: i, 2: i+1/2)
  call this%ap%new( 1, 2, (/grid % my_nx(3, p_rc_dim )/), &
            gc_num_f1r, (/dx_vec(1)/), .true. )

  ! radial component of volume element (1: i)
  call this%vr%new( 1, 2, (/grid % my_nx(3, p_rc_dim )/), &
            gc_num_f1r + 1, (/dx_vec(1)/), .true. )

  ! polar component of volume element (1: i)
  call this%vt%new( 1, 2, (/grid % my_nx(3, p_tc_dim )/), &
            gc_num_f1t + 1, (/dx_vec(2)/), .true. )

  ! setup pointers to improve readibility
  r     => this%r%f1
  dr    => this%dr%f1
  r3    => this%r3%f1

  t     => this%t%f1
  st    => this%st%f1
  ct    => this%ct%f1

  lr    => this%lr%f1
  lt    => this%lt%f1
  lp    => this%lp%f2

  ar    => this%ar%f2
  at    => this%at%f2
  ap    => this%ap%f1

  vr    => this%vr%f1
  vt    => this%vt%f1

  alpha => this%alpha%f1
  beta  => this%beta%f2

  rs = this%rs
  beta0 = this%beta0

  do i1 = - gc_num_f1r(p_lower, 1), grid % my_nx(3, p_rc_dim) + gc_num_f1r(p_upper, 1) + 1
    ! get total grid index
    i1_tot = grid % my_nx( p_lower, p_rc_dim ) - 1 + i1

    ! compute radial positions
    r(1, i1)  = r_min * delta**(i1_tot-1)
    dr(1, i1) = r_min * delta**i1_tot - r_min * delta**(i1_tot-1)
    r(2, i1)  = r(1, i1) + 0.5_p_double * dr(1, i1)

    r3(1, i1) = r(1, i1)**3
    r3(2, i1) = r(2, i1)**3

    ! compute alpha
    alpha(1, i1) = ( 1 - rs / r(1,i1) )**0.5
    alpha(2, i1) = ( 1 - rs / r(2,i1) )**0.5

  enddo

  ! CFL condition
  this%dtcfl = 1.0_p_double / sqrt( (1._p_double/dr(1,1))**2 + (1._p_double/dt)**2 )

  do i2 = - gc_num_f1t(p_lower, 1), grid % my_nx(3, p_tc_dim ) + gc_num_f1t(p_upper, 1) + 1
    ! get total grid index
    i2_tot = grid % my_nx( p_lower, p_tc_dim ) - 1 + i2 - 1

    ! compute polar positions
    t_eval = t_min + i2_tot * dt
    t(1, i2) = t_eval
    t(2, i2) = t_eval + 0.5_p_double * dt

    ! compute complex functions of polar positions
    st(1, i2) = sin( t(1, i2) )
    st(2, i2) = sin( t(2, i2) )
    ct(1, i2) = cos( t(1, i2) )
    ct(2, i2) = cos( t(2, i2) )

  enddo

  ! correct sines and cosines at polar boundaries
  if ( this % north_axis ) then
    st(1, 1) = st(2, 1)
    st(2, 0) = 0.
    ct(2, 0) = 1.
  endif

  if ( this % south_axis ) then
    st(1, grid % my_nx(3, p_tc_dim )+1) = st(2, grid % my_nx(3, p_tc_dim ))
    st(2, grid % my_nx(3, p_tc_dim )+1) = 0.
    ct(2, grid % my_nx(3, p_tc_dim )+1) = -1.
  endif

  do i1 = 1 - gc_num_f1r(p_lower, 1), grid % my_nx(3, p_rc_dim) + gc_num_f1r(p_upper, 1)
    ! compute line and area elements
    lr(1, i1)  = r(2, i1) * alpha(2, i1) - r(2, i1-1) * alpha(2, i1-1) + &
                 0.5_p_double * rs * log( (2._p_double * r(2, i1  ) * (alpha(2, i1  ) + 1._p_double) - rs)/&
                 (2._p_double * r(2, i1-1) * (alpha(2, i1-1) + 1._p_double) - rs) )

    lr(2, i1)  = r(1, i1+1) * alpha(1, i1+1) - r(1, i1 ) * alpha(1, i1) + &
                 0.5_p_double * rs * log( (2._p_double * r(1, i1+1) * (alpha(1, i1+1) + 1._p_double) - rs)/&
                 (2._p_double * r(1, i1  ) * (alpha(1, i1  ) + 1._p_double) - rs) )

    lt(1, i1) = r(1, i1) * dt
    lt(2, i1) = r(2, i1) * dt

    ap(1, i1) = dt * (0.375_p_double * rs**2 * log( (2._p_double * r(2, i1  ) * (alpha(2, i1  ) + 1._p_double) - rs) /&
                                                    (2._p_double * r(2, i1-1) * (alpha(2, i1-1) + 1._p_double) - rs) ) + &
                      0.25_p_double * ( r(2, i1  ) * alpha(2, i1  ) * (2._p_double * r(2, i1  ) + 3._p_double * rs) -&
                                        r(2, i1-1) * alpha(2, i1-1) * (2._p_double * r(2, i1-1) + 3._p_double * rs) ) )

    ap(2, i1) = dt * (0.375_p_double * rs**2 * log( (2._p_double * r(1, i1+1) * (alpha(1, i1+1) + 1._p_double) - rs) /&
                                                    (2._p_double * r(1, i1  ) * (alpha(1, i1  ) + 1._p_double) - rs) ) + &
                      0.25_p_double * ( r(1, i1+1) * alpha(1, i1+1) * (2._p_double * r(1, i1+1) + 3._p_double * rs) -&
                                        r(1, i1  ) * alpha(1, i1  ) * (2._p_double * r(1, i1  ) + 3._p_double * rs) ) )
  enddo

  select case ( this%metric )
  case ( p_geometry_minkowski )
    do i1 = 1 - gc_num_f1r(p_lower, 1), grid % my_nx(3, p_rc_dim) + gc_num_f1r(p_upper, 1)
    ! compute radial component of volume element
    vr(1, i1) = r3(2, i1  ) - r3(2, i1-1)
    vr(2, i1) = r3(1, i1+1) - r3(1, i1  )
    enddo
  case default
    do i1 = 1 - gc_num_f1r(p_lower, 1), grid % my_nx(3, p_rc_dim) + gc_num_f1r(p_upper, 1)

      vr(1, i1) = 0.3125_p_double * rs**3 * log(((1._p_double + alpha(2, i1))*(1._p_double - alpha(2, i1-1))) /&
                                                ((1._p_double - alpha(2, i1))*(1._p_double + alpha(2, i1-1)))) +&
                  r(2, i1  ) * alpha(2, i1  ) * (8._p_double*r(2, i1  )**2 + 10._p_double*r(2, i1  )*rs + 15._p_double*rs**2) / 24._p_double -&
                  r(2, i1-1) * alpha(2, i1-1) * (8._p_double*r(2, i1-1)**2 + 10._p_double*r(2, i1-1)*rs + 15._p_double*rs**2) / 24._p_double

      vr(2, i1) = 0.3125_p_double * rs**3 * log(((1._p_double + alpha(1, i1+1))*(1._p_double - alpha(1, i1))) /&
                                                ((1._p_double - alpha(1, i1+1))*(1._p_double + alpha(1, i1)))) +&
                  r(1, i1+1) * alpha(1, i1+1) * (8._p_double*r(1, i1+1)**2 + 10._p_double*r(1, i1+1)*rs + 15._p_double*rs**2) / 24._p_double -&
                  r(1, i1  ) * alpha(1, i1  ) * (8._p_double*r(1, i1  )**2 + 10._p_double*r(1, i1  )*rs + 15._p_double*rs**2) / 24._p_double
    enddo
  end select

  do i2 = 1 - gc_num_f1t(p_lower, 1), grid % my_nx(3, p_tc_dim ) + gc_num_f1t(p_upper, 1)
    ! compute polar component of volume element
    vt(1, i2) = ct(2, i2-1) - ct(2, i2  )
    vt(2, i2) = ct(1, i2  ) - ct(1, i2+1)

  enddo

  do i1 = 1 - gc_num_f1r(p_lower, 1), grid % my_nx(3, p_rc_dim ) + gc_num_f1r(p_upper, 1)
    do i2 = 1 - gc_num_f1t(p_lower, 1), grid % my_nx(3, p_tc_dim ) + gc_num_f1t(p_upper, 1)

      ! compute beta
      beta(1, i1, i2) = - st(1, i2) / r(2, i1)**2
      beta(2, i1, i2) = - st(2, i2) / r(1, i1)**2

      ! compute line and area elements
      lp(1, i1, i2) = r(1, i1) * st(1, i2) * dp
      lp(2, i1, i2) = r(2, i1) * st(2, i2) * dp

      ar(1, i1, i2) = r(2, i1)**2 * (ct(2, i2-1) - ct(2, i2  )) * dp
      ar(2, i1, i2) = r(1, i1)**2 * (ct(1, i2  ) - ct(1, i2+1)) * dp

      at(1, i1, i2) = st(1, i2) * ap(2, i1) / dt * dp
      at(2, i1, i2) = st(2, i2) * ap(1, i1) / dt * dp
    enddo
  enddo

  ! correct line and area elements at polar boundaries
  if ( this % south_axis ) then
    at(2, :, grid % my_nx(3, p_tc_dim )+1) = at(2, :, grid % my_nx(3, p_tc_dim ))
  endif

  ! clear pointers
  nullify(r)
  nullify(dr)
  nullify(r3)

  nullify(t)
  nullify(st)
  nullify(ct)

  nullify(lr)
  nullify(lt)
  nullify(lp)

  nullify(ar)
  nullify(at)
  nullify(ap)

  nullify(vr)
  nullify(vt)

  nullify(alpha)
  nullify(beta)

end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Cleanup geometry object
!-------------------------------------------------------------------------------
subroutine cleanup_gr( this )

  implicit none

  class(t_geometry_gr) , intent(inout) :: this

  ! cleanup coordinate data structures
  call this%r%cleanup()
  call this%dr%cleanup()
  call this%r3%cleanup()

  call this%t%cleanup()
  call this%st%cleanup()
  call this%ct%cleanup()

  ! cleanup integration coeficients
  call this%lr%cleanup()
  call this%lt%cleanup()
  call this%lp%cleanup()

  call this%ar%cleanup()
  call this%at%cleanup()
  call this%ap%cleanup()

  call this%vr%cleanup()
  call this%vt%cleanup()

  call this%alpha%cleanup()
  call this%beta%cleanup()

end subroutine cleanup_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Advance geometry object
!-------------------------------------------------------------------------------
subroutine advance_gr( this, t )

  implicit none

  class(t_geometry_gr) , intent(inout) :: this
  real( p_double ), intent( in ) :: t

  if (t .le. this%trise) then
    ! update star angular frequency
    this%omega = this%omega0 * sin(pi_2 * t / this%trise)**2
  else
    this%omega = this%omega0
  endif

  ! update magnitude of beta vector
  select case ( this%metric )
  case ( p_geometry_kerr_slow )
    this%beta0 = 0.21_p_double * this%omega * this%rs / (1._p_double - this%rs)
  end select
  

end subroutine advance_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  First-order radial spline function - rectangular rule - Minkowski
!-------------------------------------------------------------------------------
subroutine spline_r_rect_gr( this, ix, x, s )

  implicit none

  class( t_geometry_gr ), intent(in) :: this
  integer, intent(in) :: ix
  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:1), intent(out) :: s

  ! local variables
  real(p_k_fld) :: r, rmin, rmax, r0min, r0max
  logical :: cross

  ! get particle position
  r = real( this%r%f1(2, ix)+x*this%dr%f1(1,ix), p_k_fld )

  ! get particle boundaries
  rmin = real( r*(1.-1./this%A), p_k_fld )
  rmax = real( r*(1.+1./this%A), p_k_fld )

  ! get integration limits in cell ix
  r0min = max(rmin, this%r%f1(2, ix-1));
  r0max = min(rmax, this%r%f1(2, ix  ));

  ! compute spline in cell ix
  s(0) = (r0max**3 - r0min**3) / r**3 / this%vr%f1(1, ix)

  ! get crossing criterion
  cross = x .ge. real( -0.25_p_double*(1.+1./this%delta), p_k_fld )

  if (cross) then
      ! particle deposits in cell ix+1
      s(-1) = 0
      s( 1) = (rmax**3 - this%r3%f1(2, ix)) / r**3 / this%vr%f1(1, ix+1)
  else
      ! particle deposits in cell ix-1
      s(-1) = (this%r3%f1(2, ix-1) - rmin**3) / r**3 / this%vr%f1(1, ix-1)
      s( 1) = 0
  endif

end subroutine spline_r_rect_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  First-order radial spline function - integral rule - Minkowski
!-------------------------------------------------------------------------------
subroutine spline_r_int_gr( this, ix, x, s )

  implicit none

  class( t_geometry_gr ), intent(in) :: this
  integer, intent(in) :: ix
  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:1), intent(out) :: s

  ! local variables
  real(p_k_fld) :: r, rmin, rmax, r0min, r0max
  logical :: cross

  ! get particle position
  r = real( this%r%f1(2, ix)+x*this%dr%f1(1,ix), p_k_fld )

  ! get particle boundaries
  rmin = real( r*(1.-1./this%A), p_k_fld )
  rmax = real( r*(1.+1./this%A), p_k_fld )

  ! get integration limits in cell ix
  r0min = max(rmin, this%r%f1(2, ix-1));
  r0max = min(rmax, this%r%f1(2, ix  ));

  ! compute spline in cell ix
  s(0) = (r0max**3 - r0min**3) / r**3 / this%vr%f1(1, ix)

  ! get crossing criterion
  cross = rmin .ge. this%r%f1(2, ix-1)

  if (cross) then
      ! particle deposits in cell ix+1
      s(-1) = 0
      s( 1) = (rmax**3 - this%r3%f1(2, ix)) / r**3 / this%vr%f1(1, ix+1)
  else
      ! particle deposits in cell ix-1
      s(-1) = (this%r3%f1(2, ix-1) - rmin**3) / r**3 / this%vr%f1(1, ix-1)
      s( 1) = 0
  endif

end subroutine spline_r_int_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  First-order radial spline function - rectangular rule - Schwarzschild
!-------------------------------------------------------------------------------
subroutine spline_r_rect_curved_gr( this, ix, x, s )

  implicit none

  class( t_geometry_gr ), intent(in) :: this
  integer, intent(in) :: ix
  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:1), intent(out) :: s

  ! local variables
  real(p_k_fld) :: r, rmin, rmax, r0min, r0max
  logical :: cross

  ! get particle position
  r = real( this%r%f1(2, ix)+x*this%dr%f1(1,ix), p_k_fld )

  ! get particle boundaries
  rmin = real( r*(1.-1./this%A), p_k_fld )
  rmax = real( r*(1.+1./this%A), p_k_fld )

  ! get integration limits in cell ix
  r0min = max(rmin, this%r%f1(2, ix-1));
  r0max = min(rmax, this%r%f1(2, ix  ));

  ! compute spline in cell ix
  s(0) =  ( 0.3125_p_double * this%rs**3 * &
              log(((1._p_double + sqrt(1._p_double-this%rs/r0max )) * (1._p_double - sqrt(1._p_double-this%rs/r0min ))) / &
                  ((1._p_double - sqrt(1._p_double-this%rs/r0max )) * (1._p_double + sqrt(1._p_double-this%rs/r0min )))) +&
              r0max * sqrt(1._p_double-this%rs/r0max ) * (8._p_double * r0max**2 + 10._p_double * r0max * this%rs +&
                                                          15._p_double * this%rs**2 ) / 24._p_double -&
              r0min * sqrt(1._p_double-this%rs/r0min ) * (8._p_double * r0min**2 + 10._p_double * r0min * this%rs +&
                                                          15._p_double * this%rs**2 ) / 24._p_double  &
          ) * sqrt(1._p_double-this%rs/r) / r**3 / this%vr%f1(1, ix)
              
  ! get crossing criterion
  cross = rmin .ge. this%r%f1(2, ix-1)

  if (cross) then
      ! particle deposits in cell ix+1
      s( 1) = ( 0.3125_p_double * this%rs**3 * &
                log(((1._p_double + sqrt(1._p_double - this%rs/rmax )) * (1._p_double - this%alpha%f1(2, ix))) /&
                    ((1._p_double - sqrt(1._p_double - this%rs/rmax )) * (1._p_double + this%alpha%f1(2, ix)))) +&
                rmax * sqrt(1._p_double-this%rs/rmax) * (8._p_double * rmax**2 + 10._p_double * rmax * this%rs +&
                                                         15._p_double* this%rs**2 ) / 24._p_double -&
                this%r%f1(2, ix) * this%alpha%f1(2, ix) * (8._p_double * this%r%f1(2, ix)**2 + 10._p_double * this%r%f1(2, ix) * this%rs +&
                                                           15._p_double * this%rs**2) / 24._p_double &
              ) * sqrt(1._p_double-this%rs/r) / r**3 / this%vr%f1(1, ix+1)

      s(-1) = 0

  else
      ! particle deposits in cell ix-1
      s(-1) = ( 0.3125_p_double * this%rs**3 *&
                log(((1._p_double + this%alpha%f1(2, ix-1)) * (1._p_double - sqrt(1._p_double-this%rs/rmin))) /&
                    ((1._p_double - this%alpha%f1(2, ix-1)) * (1._p_double + sqrt(1._p_double-this%rs/rmin)))) +&
                this%r%f1(2, ix-1) * this%alpha%f1(2, ix-1) * (8._p_double * this%r%f1(2, ix-1)**2 + 10._p_double * this%r%f1(2, ix-1) * this%rs +&
                                                               15._p_double * this%rs**2) / 24._p_double -&
                rmin * sqrt(1._p_double-this%rs/rmin) * (8._p_double * rmin**2 + 10._p_double * rmin * this%rs +&
                                                               15._p_double * this%rs**2) / 24._p_double &
              ) * sqrt(1._p_double-this%rs/r) / r**3 / this%vr%f1(1, ix-1)
      
      s( 1) = 0
  endif

end subroutine spline_r_rect_curved_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  First-order radial spline function - integral rule - Schwarzschild
!-------------------------------------------------------------------------------
subroutine spline_r_int_curved_gr( this, ix, x, s )

  implicit none

  class( t_geometry_gr ), intent(in) :: this
  integer, intent(in) :: ix
  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:1), intent(out) :: s

  ! local variables
  real(p_k_fld) :: r, rmin, rmax, r0min, r0max
  real(p_k_fld) :: rplus, rminus, aplus, aminus, intr
  logical :: cross

  ! get particle position
  r = real( this%r%f1(2, ix)+x*this%dr%f1(1,ix), p_k_fld )

  ! the splines are modified in using the rectangular or the integral rule
  ! integral rule quantities
  rplus  = 2._p_double * r / (1._p_double + 1._p_double / this%delta) !ri+1
  rminus = 2._p_double * r / (1._p_double + this%delta)               !ri
  aplus  = sqrt( 1._p_double - this%rs / rplus)                       !alpha(ri+1)
  aminus = sqrt( 1._p_double - this%rs / rminus)                      !alpha(ri)
  intr = 0.0416666666666667_p_double *&
          (rplus  * aplus  * (8._p_double * rplus**2  + 10._p_double * rplus  * this%rs + 15._p_double * this%rs**2)  -&
           rminus * aminus * (8._p_double * rminus**2 + 10._p_double * rminus * this%rs + 15._p_double * this%rs**2)) +&
         0.3125_p_double * this%rs**3 * log(((1._p_double+aplus)*(1._p_double-aminus))/&
                                            ((1._p_double-aplus)*(1._p_double+aminus)))

  ! get particle boundaries
  rmin = real( r*(1.-1./this%A), p_k_fld )
  rmax = real( r*(1.+1./this%A), p_k_fld )

  ! get integration limits in cell ix
  r0min = max(rmin, this%r%f1(2, ix-1));
  r0max = min(rmax, this%r%f1(2, ix  ));

  ! compute spline in cell ix
  s(0) =  ( 0.3125_p_double * this%rs**3 * &
              log(((1._p_double + sqrt(1._p_double-this%rs/r0max )) * (1._p_double - sqrt(1._p_double-this%rs/r0min ))) / &
                  ((1._p_double - sqrt(1._p_double-this%rs/r0max )) * (1._p_double + sqrt(1._p_double-this%rs/r0min )))) +&
              r0max * sqrt(1._p_double-this%rs/r0max ) * (8._p_double * r0max**2 + 10._p_double * r0max * this%rs +&
                                                          15._p_double * this%rs**2 ) / 24._p_double -&
              r0min * sqrt(1._p_double-this%rs/r0min ) * (8._p_double * r0min**2 + 10._p_double * r0min * this%rs +&
                                                          15._p_double * this%rs**2 ) / 24._p_double  &
          ) / intr / this%vr%f1(1, ix)
              
  ! get crossing criterion
  cross = rmin .ge. this%r%f1(2, ix-1)

  if (cross) then
      ! particle deposits in cell ix+1
      s( 1) = ( 0.3125_p_double * this%rs**3 * &
                log(((1._p_double + sqrt(1._p_double - this%rs/rmax )) * (1._p_double - this%alpha%f1(2, ix))) /&
                    ((1._p_double - sqrt(1._p_double - this%rs/rmax )) * (1._p_double + this%alpha%f1(2, ix)))) +&
                rmax * sqrt(1._p_double-this%rs/rmax) * (8._p_double * rmax**2 + 10._p_double * rmax * this%rs +&
                                                          15._p_double* this%rs**2 ) / 24._p_double -&
                this%r%f1(2, ix) * this%alpha%f1(2, ix) * (8._p_double * this%r%f1(2, ix)**2 + 10._p_double * this%r%f1(2, ix) * this%rs +&
                                                            15._p_double * this%rs**2) / 24._p_double &
              ) / intr / this%vr%f1(1, ix+1)

      s(-1) = 0
  else
      ! particle deposits in cell ix-1
      s(-1) = ( 0.3125_p_double * this%rs**3 *&
                log(((1._p_double + this%alpha%f1(2, ix-1)) * (1._p_double - sqrt(1._p_double-this%rs/rmin))) /&
                    ((1._p_double - this%alpha%f1(2, ix-1)) * (1._p_double + sqrt(1._p_double-this%rs/rmin)))) +&
                this%r%f1(2, ix-1) * this%alpha%f1(2, ix-1) * (8._p_double * this%r%f1(2, ix-1)**2 + 10._p_double * this%r%f1(2, ix-1) * this%rs +&
                                                                15._p_double * this%rs**2) / 24._p_double -&
                rmin * sqrt(1._p_double-this%rs/rmin) * (8._p_double * rmin**2 + 10._p_double * rmin * this%rs +&
                                                                15._p_double * this%rs**2) / 24._p_double &
              ) / intr / this%vr%f1(1, ix-1)
      
      s( 1) = 0
  endif

end subroutine spline_r_int_curved_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  First-order polar spline function
!-------------------------------------------------------------------------------
subroutine spline_t_gr( this, jx, y, s )

  implicit none

  class( t_geometry_gr ), intent(in) :: this
  integer, intent(in) :: jx
  real( p_k_fld ), intent(in) :: y
  real( p_k_fld ), dimension(0:1), intent(out) :: s

  ! local variables
  ! real(p_k_fld) :: t, tmin, tmax, ctmin, ctmax
  real(p_k_fld) :: dt, aux_0, aux_1
  real(p_k_fld) :: ys
  
  dt = this%dt

  !! exact expression
  ! get particle position
  ! t = real( this%t%f1(2, jx)+y*dt, p_k_fld )

  ! dt = 0.5_p_double * dt

  ! get particle boundaries
  ! tmin = max(t - dt, 0.0_p_double)
  ! tmax = min(t + dt, pi)

  ! get complex functions
  ! ctmin = cos( tmin )
  ! ctmax = cos( tmax )
  ! ctmin = min(this%ct%f1(2, jx) + this%st%f1(2, jx) * tmin, 1.)
  ! ctmax = max(this%ct%f1(2, jx) + this%st%f1(2, jx) * tmax, -1.)

  ! get weight functions
  ! aux_0 = (ctmin - this%ct%f1(2, jx)) / (ctmin - ctmax)
  ! aux_1 = (this%ct%f1(2, jx) - ctmax) / (ctmin - ctmax)
  
  !! second order taylor expansion
  ! get auxiliary variable
  ys = y - 0.5_p_double

  ! get weight functions
  aux_0 = (- ys * this%st%f1(2, jx) - 0.5_p_double * dt * ys * ys * this%ct%f1(2, jx)) / \
          (this%st%f1(2, jx) + dt * y * this%ct%f1(2, jx))
  aux_1 = 1.0_p_double - aux_0

  ! compute spline
  s(0) = aux_0 / this%vt%f1(1, jx)
  s(1) = aux_1 / this%vt%f1(1, jx+1)

end subroutine spline_t_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  First-order radial spline function at half grid points
!-------------------------------------------------------------------------------
subroutine spline_rh_gr( this, ix, x, s )

  implicit none

  class( t_geometry_gr ), intent(in) :: this
  integer, intent(in) :: ix
  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:1), intent(out) :: s

  ! local variables
  real(p_k_fld) :: r

  ! get particle position
  r = real( this%r%f1(2, ix)+x*this%dr%f1(1,ix), p_k_fld )

  ! compute spline
  s = 0

  if (x < 0) then
    s(-1) = -x
  else
    s( 1) = +x
  endif

  s( 0) = 1._p_double - s(-1) - s( 1)

  s = s * sqrt(1._p_double - this%rs / r) / r**3

end subroutine spline_rh_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  First-order polar spline function at half grid points
!-------------------------------------------------------------------------------
subroutine spline_th_gr( this, jx, y, s )

  implicit none

  class( t_geometry_gr ), intent(in) :: this
  integer, intent(in) :: jx
  real( p_k_fld ), intent(in) :: y
  real( p_k_fld ), dimension(-1:1), intent(out) :: s

  ! local variables
  real(p_k_fld) :: t, tmin, tmax, ctmin, ctmax, dt
  
  dt = this%dt

  ! get particle position
  t = real( this%t%f1(2, jx)+y*dt, p_k_fld )

  dt = 0.5 * dt

  ! get particle boundaries
  tmin = max(t - dt, 0.0_p_double)
  tmax = min(t + dt, pi)

  ! get complex functions
  ctmin = cos( tmin )
  ctmax = cos( tmax )

  ! compute spline
  s = 0

  if (y < 0) then
    s(-1) = -y
  else
    s( 1) = +y
  endif

  s( 0) = 1._p_double - s(-1) - s( 1)

  s = s / (ctmin - ctmax)

end subroutine spline_th_gr
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_gr( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------
  use m_restart
  use m_parameters
  
  implicit none

  class( t_geometry_gr ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for geometry_gr object.'
  integer :: ierr

  restart_io_wr( p_geometry_gr_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  !do nothing for now
  

end subroutine write_checkpoint_gr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_gr( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file
!-----------------------------------------------------------------------------------------
  use m_restart
  use m_parameters
  
  implicit none

  class( t_geometry_gr ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for geometry_gr object.'
  character(len=len(p_geometry_gr_rst_id)) :: rst_id
  integer :: ierr

  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  !do nothing for now

end subroutine restart_read_gr
!-----------------------------------------------------------------------------------------


end module m_geometry_gr