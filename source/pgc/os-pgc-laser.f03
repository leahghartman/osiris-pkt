#include "os-config.h"
#include "os-preprocess.fpp"

!-------------------------------------------------------------------------------
! Module with different PGC laser types
!-------------------------------------------------------------------------------
module m_pgc_laser

use m_system
use m_parameters

use m_emf_pgc_define

implicit none

private

interface calculate_per_plane
  module procedure calculate_per_plane_pgc
end interface

interface calculate_per_gaussian
  module procedure calculate_per_gaussian_pgc
end interface

interface calculate_per_gaussian_asym
  module procedure calculate_per_gaussian_asym_pgc
end interface

interface calculate_lon_polynomial
  module procedure calculate_lon_polynomial_pgc
end interface

interface calculate_lon_sin2
  module procedure calculate_lon_sin2_pgc
end interface

interface calculate_lon_gaussian
  module procedure calculate_lon_gaussian_pgc
end interface calculate_lon_gaussian

interface launch_laser
  module procedure launch_laser_pgc
end interface launch_laser

public :: launch_laser

contains

!-------------------------------------------------------------------------------
! Interface to launch the pgc laser
!-------------------------------------------------------------------------------
subroutine launch_laser_pgc(this , g_xmin , n_x_pmin , dx , dt)

  use m_system
  use m_parameters

  use m_emf_pgc_define

  implicit none

  class( t_emf_pgc ), intent(inout)  ::  this
  real(p_k_fld), dimension(:), intent(in) ::  g_xmin
  real(p_double) , dimension(:), intent(in) :: dx
  integer, dimension(:), intent(in) :: n_x_pmin
  real( p_double ),   intent(in) :: dt

  ! initialize max of envelope
  select case(p_x_dim)
  case(2)
    this%env%z2(:, 1:this%env%nx_(1), 1:this%env%nx_(2)) = this%a0
  case(3)
    this%env%z3(:, 1:this%env%nx_(1), 1:this%env%nx_(2), 1:this%env%nx_(3)) = this%a0
  end select

  ! calculate the perpendicular profile
  select case ( this%per_type )
  case ( p_pgc_per_plane )
    call calculate_per_plane( this )

  case ( p_pgc_per_gaussian )
    call calculate_per_gaussian( this , g_xmin , n_x_pmin , dx , dt )

  case ( p_pgc_per_gaussian_asym )
    call calculate_per_gaussian_asym( this , g_xmin , n_x_pmin , dx , dt )

  case default
    print * , "[pgc :: error] perpendicular laser pulse profile."
    print * , "aborting ..."
    stop

  end select

  ! calculate the longitudinal profile
  select case ( this%lon_type )
  case ( p_pgc_lon_gaussian )
    call calculate_lon_gaussian( this , g_xmin , n_x_pmin , dx , dt )

  case ( p_pgc_lon_polynomial )
    call calculate_lon_polynomial( this , g_xmin , n_x_pmin , dx , dt )

  case ( p_pgc_lon_sin2 )
    call calculate_lon_sin2( this , g_xmin , n_x_pmin , dx , dt )

  case default
    print * , "[pgc :: error] longitudinal laser pulse profile."
    print * , "aborting ..."
    stop

  end select

end subroutine launch_laser_pgc
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine calculate_per_plane_pgc( this )

  implicit none

  ! input/output variables
  class( t_emf_pgc ), intent(inout), target  ::  this

  ! aux variables definition
  integer :: i , j , k
  integer :: nx_max , ny_max , nz_max , nx_min , ny_min , nz_min

  ! only cartesian
  if ( this%coordinates /= p_cartesian ) then
    print * , "[pgc :: error] plane wave only for cartesian coordinates"
    print * , "aborting ..."
    stop
  endif

  select case (p_x_dim)
  case (2)
    nx_min = 1
    ny_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)

    do j = ny_min , ny_max
      do i = nx_min, nx_max
        this%env%z2(p_env_nm,i,j) = &
          this%env%z2(p_env_nm,i,j) * per_env_plane_pgc()
        this%env%z2(p_env_n,i,j) = &
          this%env%z2(p_env_n,i,j) * per_env_plane_pgc()
        this%env%z2(p_env_np,i,j) = &
          this%env%z2(p_env_np,i,j) * per_env_plane_pgc()
      enddo
    enddo

  case (3)
    nx_min = 1
    ny_min = 1
    nz_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)
    nz_max = this%env%nx_(3)

    do k = nz_min , nz_max
      do j = ny_min , ny_max
        do i = nx_min, nx_max
          this%env%z3(p_env_nm,i,j,k) = &
            this%env%z3(p_env_nm,i,j,k) * per_env_plane_pgc()
          this%env%z3(p_env_n,i,j,k) = &
            this%env%z3(p_env_n,i,j,k) * per_env_plane_pgc()
          this%env%z3(p_env_np,i,j,k) = &
            this%env%z3(p_env_np,i,j,k) * per_env_plane_pgc()
        enddo
      enddo
    enddo
  end select

end subroutine calculate_per_plane_pgc
!-------------------------------------------------------------------------------
subroutine calculate_per_gaussian_pgc( this , g_xmin , n_x_pmin , dx , dt )

  implicit none

  ! input/output variables
  class( t_emf_pgc ), intent(inout), target  ::  this
  real( p_k_fld) , dimension(:), intent(in) ::  g_xmin
  integer, dimension(:), intent(in) :: n_x_pmin
  real( p_double ) , dimension(:), intent(in) :: dx
  real( p_double ),   intent(in) :: dt

  ! aux variables definition
  real(p_k_fld) :: ximin , yimin , zimin
  real(p_k_fld) :: xi , yi , zi
  real(p_k_fld) :: dxi , dyi , dzi
  integer :: i , j , k , imin , jmin , kmin
  integer :: nx_max , ny_max , nz_max , nx_min , ny_min , nz_min
  complex(p_k_fld) :: per_nm, per_n, per_np
  real(p_k_fld) :: y0, z0

  select case (p_x_dim)
  case (2)
    ximin = g_xmin(1)
    yimin = g_xmin(2)

    dxi = real(dx(1), p_k_fld)
    dyi = real(dx(2), p_k_fld)

    imin = n_x_pmin(1)
    jmin = n_x_pmin(2)

    nx_min = 1
    ny_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)

    ! extract the shift of the optical axis
    y0 = this%per_center(1)
    z0 = 0.0_p_k_fld

    ! no dependence on z-coordinate
    zi = 0.0_p_k_fld

    do j = ny_min , ny_max
      yi = ( jmin + j - 2 ) * dyi + yimin
      do i = nx_min, nx_max
        xi = ( imin + i - 2 ) * dxi + ximin

        per_nm = per_env_gaussian_pgc( &
          xi-0.5_p_k_fld*real(dt, p_k_fld), yi, zi, this%w0, this%omega, &
          this%per_focus, y0, z0 )
        per_n  = per_env_gaussian_pgc( &
          xi+0.5_p_k_fld*real(dt, p_k_fld), yi, zi, this%w0, this%omega, this%per_focus, y0, z0 )
        per_np = per_env_gaussian_pgc( &
          xi+1.5_p_k_fld*real(dt, p_k_fld), yi, zi, this%w0, this%omega, &
          this%per_focus, y0,  z0 )

        this%env%z2(p_env_nm,i,j) = this%env%z2(p_env_nm,i,j) * per_nm
        this%env%z2(p_env_n,i,j) = this%env%z2(p_env_n,i,j) * per_n
        this%env%z2(p_env_np,i,j) = this%env%z2(p_env_np,i,j) * per_np
      enddo
    enddo

  case (3)
    ximin = g_xmin(1)
    yimin = g_xmin(2)
    zimin = g_xmin(3)

    dxi = real(dx(1), p_k_fld)
    dyi = real(dx(2), p_k_fld)
    dzi = real(dx(3), p_k_fld)

    imin = n_x_pmin(1)
    jmin = n_x_pmin(2)
    kmin = n_x_pmin(3)

    nx_min = 1
    ny_min = 1
    nz_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)
    nz_max = this%env%nx_(3)

    ! extract the shift of the optical axis
    y0 = this%per_center(1)
    z0 = this%per_center(2)

    do k = nz_min , nz_max
      zi = ( kmin + k - 2 ) * dzi + zimin

      do j = ny_min , ny_max
        yi = ( jmin + j - 2 ) * dyi + yimin

        do i = nx_min, nx_max
          xi = ( imin + i - 2 ) * dxi + ximin

          per_nm = per_env_gaussian_pgc( &
            xi-0.5_p_k_fld*real(dt, p_k_fld), yi, zi, this%w0, this%omega, &
            this%per_focus, y0, z0 )
          per_n = per_env_gaussian_pgc( &
            xi, yi, zi, this%w0, this%omega, this%per_focus, y0,  z0 )
          per_np = per_env_gaussian_pgc( &
            xi+0.5_p_k_fld*real(dt, p_k_fld), yi, zi, this%w0, this%omega, &
            this%per_focus, y0, z0 )

          this%env%z3(p_env_nm,i,j,k) = this%env%z3(p_env_nm,i,j,k) * per_nm
          this%env%z3(p_env_n,i,j,k) = this%env%z3(p_env_n,i,j,k) * per_n
          this%env%z3(p_env_np,i,j,k) = this%env%z3(p_env_np,i,j,k) * per_np
        enddo
      enddo
    enddo

  end select

end subroutine calculate_per_gaussian_pgc
!-------------------------------------------------------------------------------
subroutine calculate_per_gaussian_asym_pgc( this , g_xmin , n_x_pmin , dx , dt )

  implicit none

  ! input/output variables
  class( t_emf_pgc ), intent(inout), target  ::  this
  real( p_k_fld) , dimension(:), intent(in) ::  g_xmin
  integer, dimension(:), intent(in) :: n_x_pmin
  real( p_double ) , dimension(:), intent(in) :: dx
  real( p_double ),   intent(in) :: dt

  ! aux variables definition
  real(p_k_fld) :: ximin , yimin , zimin
  real(p_k_fld) :: xi , yi , zi
  real(p_k_fld) :: dxi , dyi , dzi
  integer :: i , j , k , imin , jmin , kmin
  integer :: nx_max , ny_max , nz_max , nx_min , ny_min , nz_min
  complex(p_k_fld) :: per_nm, per_n, per_np
  real(p_k_fld) :: y0 , z0 , w0_y , w0_z , focus_y , focus_z , omega

  select case (p_x_dim)
  case (2)
    print * , "[pgc :: error] asymmetric gaussian profile not available as"
    print * , "[pgc :: error] perpendicular type for 2d simulations"
    print * , "aborting ..."
    stop

  case (3)
    ximin = g_xmin(1)
    yimin = g_xmin(2)
    zimin = g_xmin(3)

    dxi = real(dx(1), p_k_fld)
    dyi = real(dx(2), p_k_fld)
    dzi = real(dx(3), p_k_fld)

    imin = n_x_pmin(1)
    jmin = n_x_pmin(2)
    kmin = n_x_pmin(3)

    nx_min = 1
    ny_min = 1
    nz_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)
    nz_max = this%env%nx_(3)

    ! extract shift of the optical axis
    y0 = this%per_center(1)
    z0 = this%per_center(2)

    ! extract asymmetric beam waist
    w0_y = this%w0_asym(1)
    w0_z = this%w0_asym(2)

    ! extract position of focal plane
    focus_y = this%per_focus_asym(1)
    focus_z = this%per_focus_asym(2)

    ! extract omega
    omega = this%omega

    do k = nz_min , nz_max
      zi = ( kmin + k - 2 ) * dzi + zimin

      do j = ny_min , ny_max
        yi = ( jmin + j - 2 ) * dyi + yimin

        do i = nx_min, nx_max
          xi = ( imin + i - 2 ) * dxi + ximin

          per_nm = per_env_gaussian_asymm_pgc( &
            xi - 0.5_p_k_fld * real(dt, p_k_fld), yi, zi, w0_y, w0_z, omega, &
            focus_y, focus_z, y0,  z0 )
          per_n = per_env_gaussian_asymm_pgc( &
            xi, yi, zi, w0_y, w0_z, omega, focus_y, focus_z, y0,  z0 )
          per_np = per_env_gaussian_asymm_pgc( &
            xi + 0.5_p_k_fld * real(dt, p_k_fld), yi, zi, w0_y, w0_z, omega, &
            focus_y, focus_z, y0,  z0 )

          this%env%z3(p_env_nm,i,j,k) = this%env%z3(p_env_nm,i,j,k) * per_nm
          this%env%z3(p_env_n,i,j,k) = this%env%z3(p_env_n,i,j,k) * per_n
          this%env%z3(p_env_np,i,j,k) = this%env%z3(p_env_np,i,j,k) * per_np
        enddo
      enddo
    enddo

  end select

end subroutine calculate_per_gaussian_asym_pgc
!-------------------------------------------------------------------------------
subroutine calculate_lon_gaussian_pgc( this , g_xmin , n_x_pmin , dx , dt )

  implicit none

  ! input/output variables
  class( t_emf_pgc ), intent(inout), target  ::  this
  real( p_k_fld) , dimension(:), intent(in) ::  g_xmin
  integer, dimension(:), intent(in) :: n_x_pmin
  real( p_double ) , dimension(:), intent(in) :: dx
  real( p_double ),   intent(in) :: dt

  ! aux variables
  real(p_k_fld) :: lon_nm, lon_n, lon_np
  integer :: i , imin , j , k
  integer :: nx_max, ny_max, nz_max, ny_min, nx_min, nz_min   ! grid numbers
  real(p_k_fld) :: ximin, xi, dxi
  real(p_k_fld) :: tau , lon_center , lon_duration
  real( p_k_fld ) , parameter :: t = 0.0_p_k_fld

  select case (p_x_dim)
  case (2)
    ximin = g_xmin(1)

    imin = n_x_pmin(1)

    dxi = real(dx(1), p_k_fld)

    nx_min = 1
    ny_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)

    ! extract envelope parameters
    lon_duration = this%lon_range
    lon_center   = this%lon_center
    tau          = this%lon_duration

    do i = nx_min, nx_max
      xi = ( imin + i - 2 ) * dxi + ximin

      lon_nm = lon_env_gauss_pgc( &
        xi, t - 0.5_p_k_fld*real(dt, p_k_fld), tau, lon_center, lon_duration )
      lon_n = lon_env_gauss_pgc( &
        xi, t + 0.5_p_k_fld*real(dt, p_k_fld), tau, lon_center, lon_duration )
      lon_np = lon_env_gauss_pgc( &
        xi, t + 1.5_p_k_fld*real(dt, p_k_fld), tau, lon_center, lon_duration )

      do j = ny_min , ny_max
        this%env%z2(p_env_nm,i,j) = this%env%z2(p_env_nm,i,j) * lon_nm
        this%env%z2(p_env_n,i,j) = this%env%z2(p_env_n,i,j) * lon_n
        this%env%z2(p_env_np,i,j) = this%env%z2(p_env_np,i,j) * lon_np
      enddo

    enddo

  case (3)
    ximin = g_xmin(1)

    imin = n_x_pmin(1)

    dxi = real(dx(1), p_k_fld)

    nx_min = 1
    ny_min = 1
    nz_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)
    nz_max = this%env%nx_(3)

    ! extract envelope parameters
    lon_duration = this%lon_range
    lon_center   = this%lon_center
    tau          = this%lon_duration

    do i = nx_min, nx_max
      xi = ( imin + i - 2 ) * dxi + ximin

      lon_nm = lon_env_gauss_pgc( &
        xi, t - 0.5_p_k_fld*real(dt, p_k_fld), tau, lon_center, lon_duration )
      lon_n  = lon_env_gauss_pgc( xi, t, tau, lon_center, lon_duration )
      lon_np = lon_env_gauss_pgc( &
        xi, t + 0.5_p_k_fld*real(dt, p_k_fld), tau, lon_center, lon_duration )

      do k = nz_min , nz_max
        do j = ny_min , ny_max
          this%env%z3(p_env_nm,i,j,k) = this%env%z3(p_env_nm,i,j,k) * lon_nm
          this%env%z3(p_env_n,i,j,k) = this%env%z3(p_env_n,i,j,k) * lon_n
          this%env%z3(p_env_np,i,j,k) = this%env%z3(p_env_np,i,j,k) * lon_np
        enddo
      enddo
    enddo

  end select

end subroutine calculate_lon_gaussian_pgc
!-------------------------------------------------------------------------------
subroutine calculate_lon_polynomial_pgc( this , g_xmin , n_x_pmin , dx , dt )

  implicit none

  ! input/output variables
  class( t_emf_pgc ), intent(inout), target  ::  this
  real( p_k_fld) , dimension(:), intent(in) ::  g_xmin
  integer, dimension(:), intent(in) :: n_x_pmin
  real( p_double ) , dimension(:), intent(in) :: dx
  real( p_double ),   intent(in) :: dt

  ! aux variables
  real(p_k_fld) :: lon_nm, lon_n, lon_np
  integer :: i , imin , j , k
  integer :: nx_max, ny_max, nz_max, ny_min, nx_min, nz_min   ! grid numbers
  real(p_k_fld) :: ximin, xi, dxi
  real(p_k_fld) :: lon_start , lon_rise , lon_fall , lon_flat
  real( p_k_fld ) , parameter :: t = 0.0_p_k_fld

  select case (p_x_dim)
  case (2)
    ximin = g_xmin(1)

    imin = n_x_pmin(1)

    dxi = real(dx(1), p_k_fld)

    nx_min = 1
    ny_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)

    ! extract envelope parameters
    lon_start = this%lon_start
    lon_rise = this%lon_rise
    lon_fall = this%lon_fall
    lon_flat = this%lon_flat

    do i = nx_min, nx_max
      xi = ( imin + i - 2 ) * dxi + ximin

      lon_nm = lon_env_poly_pgc( &
        xi, t - 0.5_p_k_fld*real(dt, p_k_fld), &
        lon_start, lon_rise, lon_fall, lon_flat )
      lon_n = lon_env_poly_pgc( xi, t + 0.5_p_k_fld*real(dt, p_k_fld), &
        lon_start, lon_rise, lon_fall, lon_flat )
      lon_np = lon_env_poly_pgc( xi, t + 1.5_p_k_fld*real(dt, p_k_fld), &
        lon_start, lon_rise, lon_fall, lon_flat )

      do j = ny_min , ny_max
        this%env%z2(p_env_nm,i,j) = this%env%z2(p_env_nm,i,j) * lon_nm
        this%env%z2(p_env_n,i,j) = this%env%z2(p_env_n,i,j) * lon_n
        this%env%z2(p_env_np,i,j) = this%env%z2(p_env_np,i,j) * lon_np
      enddo
    enddo

  case (3)
    ximin = g_xmin(1)

    imin = n_x_pmin(1)

    dxi = real(dx(1), p_k_fld)

    nx_min = 1
    ny_min = 1
    nz_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)
    nz_max = this%env%nx_(3)

    ! extract envelope parameters
    lon_start = this%lon_start
    lon_rise = this%lon_rise
    lon_fall = this%lon_fall
    lon_flat = this%lon_flat

    do i = nx_min, nx_max
      xi = ( imin + i - 2 ) * dxi + ximin

      lon_nm = lon_env_poly_pgc( &
        xi, t-0.5_p_k_fld*real(dt, p_k_fld), &
        lon_start, lon_rise, lon_fall, lon_flat )
      lon_n = lon_env_poly_pgc( xi, t, lon_start, lon_rise, lon_fall, lon_flat )
      lon_np = lon_env_poly_pgc( &
        xi, t+0.5_p_k_fld*real(dt, p_k_fld), &
        lon_start, lon_rise, lon_fall, lon_flat )

      do k = nz_min , nz_max
        do j = ny_min , ny_max
          this%env%z3(p_env_nm,i,j,k) = this%env%z3(p_env_nm,i,j,k) * lon_nm
          this%env%z3(p_env_n,i,j,k) = this%env%z3(p_env_n,i,j,k) * lon_n
          this%env%z3(p_env_np,i,j,k) = this%env%z3(p_env_np,i,j,k) * lon_np
        enddo
      enddo
    enddo
  end select

end subroutine calculate_lon_polynomial_pgc
!-------------------------------------------------------------------------------
subroutine calculate_lon_sin2_pgc( this , g_xmin , n_x_pmin , dx , dt )

  implicit none

  ! input/output variables
  class( t_emf_pgc ), intent(inout), target  ::  this
  real( p_k_fld) , dimension(:), intent(in) ::  g_xmin
  integer, dimension(:), intent(in) :: n_x_pmin
  real( p_double ) , dimension(:), intent(in) :: dx
  real( p_double ),   intent(in) :: dt

  ! aux variables
  real(p_k_fld) :: lon_nm, lon_n, lon_np
  integer :: i , imin , j , k
  integer :: nx_max, ny_max, nz_max, ny_min, nx_min, nz_min   ! grid numbers
  real(p_k_fld) :: ximin, xi, dxi
  real(p_k_fld) :: lon_center , lon_rise , lon_fall , lon_flat
  real( p_k_fld ) , parameter :: t = 0.0_p_k_fld


  ximin = g_xmin(1)
  dxi = real(dx(1), p_k_fld)
  imin = n_x_pmin(1)

  select case (p_x_dim)
  case (2)
    nx_min = 1
    ny_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)

    lon_center = this%lon_center
    lon_rise = this%lon_rise
    lon_fall = this%lon_fall
    lon_flat = this%lon_flat

    do i = nx_min, nx_max

      xi = ( imin + i - 2 ) * dxi + ximin

      lon_nm = lon_env_sin2_pgc( xi, t-0.5_p_k_fld*real(dt, p_k_fld), &
        lon_center, lon_rise, lon_fall, lon_flat )
      lon_n = lon_env_sin2_pgc(xi, t+0.5_p_k_fld*real(dt, p_k_fld), &
        lon_center, lon_rise, lon_fall, lon_flat)
      lon_np = lon_env_sin2_pgc( xi, t+1.5_p_k_fld*real(dt, p_k_fld), &
        lon_center, lon_rise, lon_fall, lon_flat )

      do j = ny_min , ny_max
        this%env%z2(p_env_nm,i,j) = this%env%z2(p_env_nm,i,j) * lon_nm
        this%env%z2(p_env_n,i,j) = this%env%z2(p_env_n,i,j) * lon_n
        this%env%z2(p_env_np,i,j) = this%env%z2(p_env_np,i,j) * lon_np
      enddo
    enddo

  case (3)
    nx_min = 1
    ny_min = 1
    nz_min = 1

    nx_max = this%env%nx_(1)
    ny_max = this%env%nx_(2)
    nz_max = this%env%nx_(3)

    lon_center = this%lon_center
    lon_rise = this%lon_rise
    lon_fall = this%lon_fall
    lon_flat = this%lon_flat

    do i = nx_min, nx_max

      xi = ( imin + i - 2 ) * dxi + ximin

      lon_nm = lon_env_sin2_pgc( &
        xi, t-0.5_p_k_fld*real(dt, p_k_fld), &
        lon_center, lon_rise, lon_fall, lon_flat )
      lon_n = lon_env_sin2_pgc(xi, t, lon_center, lon_rise, lon_fall, lon_flat)
      lon_np = lon_env_sin2_pgc( &
        xi, t+0.5_p_k_fld*real(dt, p_k_fld), &
        lon_center, lon_rise, lon_fall, lon_flat )

      do k = nz_min , nz_max
        do j = ny_min , ny_max
          this%env%z3(p_env_nm,i,j,k) = this%env%z3(p_env_nm,i,j,k) * lon_nm
          this%env%z3(p_env_n,i,j,k) = this%env%z3(p_env_n,i,j,k) * lon_n
          this%env%z3(p_env_np,i,j,k) = this%env%z3(p_env_np,i,j,k) * lon_np
        enddo
      enddo
    enddo

  end select

end subroutine calculate_lon_sin2_pgc
!-------------------------------------------------------------------------------
function per_env_plane_pgc(  )

  implicit none

  ! function output
  complex(p_k_fld) :: per_env_plane_pgc

  per_env_plane_pgc = cmplx( 1.0_p_k_fld , 0.0_p_k_fld, kind=p_k_fld )

end function per_env_plane_pgc
!-------------------------------------------------------------------------------
function per_env_gaussian_pgc( xi, yi , zi , w0 , omega , per_focus , y0 , z0 )

  implicit none

  ! input parameters
  real(p_k_fld), intent(in) :: xi , yi , zi , w0 , omega , per_focus , y0 , z0

  ! auxiliary calculations parameter
  real(p_k_fld) :: xr , rWl2 , curv , gouy_shift , x , rho , Amp

  ! function output
  complex(p_k_fld) :: per_env_gaussian_pgc

  ! Rayleigh length
  xr = 0.5_p_k_fld * omega * (w0**2.0_p_k_fld)

  ! radial coordinate for gaussian beam
  rho = sqrt( (yi - y0)**2.0_p_k_fld + (zi - z0)**2.0_p_k_fld )

  ! Default curvature and Gouy shift
  rWl2 = 1.0_p_k_fld
  curv = 0.0_p_k_fld
  gouy_shift = 0.0_p_k_fld

  ! longitudinal position
  x = xi - per_focus

  if ( x /= 0.0_p_k_fld ) then
     rWl2       = xr**2.0_p_k_fld / ( xr**2.0_p_k_fld + x**2.0_p_k_fld )
     curv       = &
      0.5_p_k_fld * rho**2.0_p_k_fld * x / ( x**2.0_p_k_fld + xr**2.0_p_k_fld )
     gouy_shift = atan2( x , xr )
  endif

  select case (p_x_dim)
  case (2)
    Amp = rWl2**0.25_p_k_fld
  case (3)
    Amp = rWl2**0.5_p_k_fld
  end select

  per_env_gaussian_pgc = &
    Amp * exp( -rho**2.0_p_k_fld * rWl2 / (w0**2.0_p_k_fld) ) * &
                cmplx( cos( -omega * curv + gouy_shift ) , &
                       sin( -omega * curv + gouy_shift ) , kind=p_k_fld)
end function per_env_gaussian_pgc
!-------------------------------------------------------------------------------
function per_env_gaussian_asymm_pgc(xi, yi, zi, w0y, w0z, omega, fy, fz, y0, z0)

  implicit none

  ! input parameters
  real(p_k_fld), intent(in) :: xi, yi, zi, w0y, w0z, omega, fy, fz, y0, z0

  ! auxiliary calculations parameter
  real(p_k_fld) :: x, y, z, gouy_shift
  real(p_k_fld) :: xr_y, rWl2_y, curv_y, gouy_shift_y
  real(p_k_fld) :: xr_z, rWl2_z, curv_z, gouy_shift_z

  ! function output
  complex(p_k_fld) :: per_env_gaussian_asymm_pgc

  ! Rayleigh length
  xr_y = 0.5_p_k_fld * omega * (w0y**2.0_p_k_fld)
  xr_z = 0.5_p_k_fld * omega * (w0z**2.0_p_k_fld)

  ! shift of the optical axis
  y = yi - y0
  z = zi - z0

  ! Default curvature and Gouy shift
  rWl2_y = 1.0_p_k_fld
  curv_y = 0.0_p_k_fld
  gouy_shift_y = 0.0_p_k_fld

  rWl2_z = 1.0_p_k_fld
  curv_z = 0.0_p_k_fld
  gouy_shift_z = 0.0_p_k_fld

  ! beam parameters for y-dependence
  x = xi - fy
  if ( x /= 0.0_p_k_fld ) then
     rWl2_y       = xr_y**2.0_p_k_fld / ( xr_y**2.0_p_k_fld + x**2.0_p_k_fld )
     curv_y       = &
      0.5_p_k_fld * y**2.0_p_k_fld * x / ( x**2.0_p_k_fld + xr_y**2.0_p_k_fld )
     gouy_shift_y = atan2( x , xr_y )
  endif

  ! beam parameters for z-dependence
  x = xi - fz
  if ( x /= 0.0_p_k_fld ) then
     rWl2_z       = xr_z**2.0_p_k_fld / ( xr_z**2.0_p_k_fld + x**2.0_p_k_fld )
     curv_z       = &
      0.5_p_k_fld * z**2.0_p_k_fld * x / ( x**2.0_p_k_fld + xr_z**2.0_p_k_fld )
     gouy_shift_z = atan2( x , xr_z )
  endif

  gouy_shift = 0.5_p_k_fld * ( gouy_shift_y + gouy_shift_z )

  per_env_gaussian_asymm_pgc = &
    sqrt( sqrt(rWl2_y) * sqrt(rWl2_z) ) * &
    exp( - y**2.0_p_k_fld * rWl2_y / (w0y**2.0_p_k_fld) ) * &
    exp( - z**2.0_p_k_fld * rWl2_z / (w0z**2.0_p_k_fld) ) * &
    cmplx( cos( -omega*(curv_y+curv_z) + gouy_shift ) , &
           sin( -omega*(curv_y+curv_z) + gouy_shift ) , kind=p_k_fld)
end function per_env_gaussian_asymm_pgc
!-------------------------------------------------------------------------------
function lon_env_gauss_pgc( xi , t , tau , lon_center , lon_duration )

  implicit none

  real(p_k_fld), intent(in) :: xi , t , tau , lon_center , lon_duration
  real(p_k_fld) :: lon_env_gauss_pgc , center

  center = lon_center - t

  if( abs(xi - center) < lon_duration/2.0_p_k_fld ) then
    lon_env_gauss_pgc = exp( - 2.0_p_k_fld*((xi-center)/tau)**2.0_p_k_fld )
  else
    lon_env_gauss_pgc = 0.0_p_k_fld
  endif

end function lon_env_gauss_pgc
!-------------------------------------------------------------------------------
function lon_env_poly_pgc( xi, t, lon_start, lon_rise, lon_fall, lon_flat )

  implicit none

  real(p_k_fld), intent(in) :: xi , t
  real(p_k_fld), intent(in) :: lon_start , lon_rise , lon_fall , lon_flat
  real(p_k_fld) :: lon_env_poly_pgc, center

  center = (lon_start - lon_rise - lon_flat/2.0_p_k_fld) - t

  if ( xi < (center - lon_flat/2.0_p_k_fld - lon_fall)) then
    lon_env_poly_pgc = 0.0_p_k_fld
  elseif ( xi < (center - lon_flat/2.0_p_k_fld) ) then
    lon_env_poly_pgc = &
      fenv(( xi - (center - lon_flat/2.0_p_k_fld - lon_fall))/lon_fall )
  elseif ( xi < (center + lon_flat/2.0_p_k_fld) ) then
    lon_env_poly_pgc = 1.0_p_k_fld
  elseif ( xi < (center + lon_flat/2.0_p_k_fld + lon_rise) ) then
    lon_env_poly_pgc = &
      fenv(-( xi - (center + lon_flat/2.0_p_k_fld + lon_rise))/lon_rise )
  else
    lon_env_poly_pgc = 0.0_p_k_fld
  endif

end function lon_env_poly_pgc

! function declaration of polynomial profile
function fenv(tt)

  implicit none

  real(p_k_fld) :: fenv
  real(p_k_fld), intent(in) :: tt

  fenv = &
    10.0_p_k_fld * tt**3.0_p_k_fld &
    -15.0_p_k_fld * tt**4.0_p_k_fld &
    +6.0_p_k_fld * tt**5.0_p_k_fld

end function fenv
!-------------------------------------------------------------------------------
function lon_env_sin2_pgc( xi , t , lon_center, lon_rise, lon_fall, lon_flat )

  implicit none

  real(p_k_fld), intent(in) :: xi , t
  real(p_k_fld), intent(in) :: lon_center , lon_rise , lon_fall , lon_flat
  real(p_k_fld) :: lon_env_sin2_pgc , center
  real(p_k_fld) , parameter :: pi = 3.141598

   center = lon_center - t

   if ( xi < (center - lon_flat/2.0_p_k_fld - lon_fall)) then
     lon_env_sin2_pgc = 0.0_p_k_fld
   elseif ( xi < (center - lon_flat/2.0_p_k_fld) ) then
     lon_env_sin2_pgc = &
      sin(pi * (xi-center)/(2.0_p_k_fld*lon_fall) - pi*0.5_p_k_fld)**2.0_p_k_fld
   elseif ( xi < (center + lon_flat/2.0_p_k_fld) ) then
     lon_env_sin2_pgc = 1.0_p_k_fld
   elseif ( xi < (center + lon_flat/2.0_p_k_fld + lon_rise) ) then
     lon_env_sin2_pgc = &
      sin(pi * (xi-center)/(2.0_p_k_fld*lon_rise) - pi*0.5_p_k_fld)**2.0_p_k_fld
   else
     lon_env_sin2_pgc = 0.0_p_k_fld
   endif

end function lon_env_sin2_pgc
!-------------------------------------------------------------------------------
end module m_pgc_laser
