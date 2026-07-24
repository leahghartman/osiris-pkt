#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_phasespace_cyl_modes

#include "memory/memory.h"

use m_parameters
use m_species_define, only : t_species, t_phasespace_list, t_phasespace_diagnostics
use m_species_define, only : t_phasespace_params, t_phasespace, p_max_phasespace_dims
use m_species_phasespace
use m_space, only : t_space, xmin_initial, xmax_initial, if_move, total_xmoved
use m_diagnostic_utilities,  only : p_diag_prec

private

! phasespace parameters
integer, parameter, public :: p_phase_xc = 7

type, extends( t_phasespace_list ) :: t_phasespace_list_cyl_modes
  contains
  procedure :: add_to_list => add_phasespace_to_list_cyl_modes
end type t_phasespace_list_cyl_modes

type, extends( t_phasespace_diagnostics ) :: t_phasespace_diag_cyl_modes

  ! physical range for phasespace data dumps
  real(p_diag_prec), dimension(p_max_dim) :: xcmin, xcmax

  ! resolutions for phasespace data dumps
  integer, dimension(p_max_dim) :: nxc, nxc_3D

  contains

  procedure :: allocate_objs => allocate_objs_pspace_diag_cyl_modes
  procedure :: init_list => init_pspace_list_cyl_modes
  procedure :: get_params => get_pspace_params_cyl_modes
  procedure :: size => size_phasespace_cyl_modes

end type t_phasespace_diag_cyl_modes

public :: t_phasespace_diag_cyl_modes

contains

!-----------------------------------------------------------------------------------------
subroutine parse_phasespace_cyl_modes( pha_name, ndims, o_x_or_p, o_xp_dim, ps_type, ierr )
!-----------------------------------------------------------------------------------------
! parses a phasespace name
!-----------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: pha_name
  integer, intent(out) :: ndims, ierr
  integer, dimension(p_max_phasespace_dims), intent(out) :: o_x_or_p, o_xp_dim
  integer, intent(out) :: ps_type

  integer:: i, pos
  integer, dimension(p_max_phasespace_dims) :: x_or_p, xp_dim
  logical:: read_dim = .false.

  ierr = 0

  ! parse the string

  ps_type = p_pstype_normal
  ndims = 0
  pos = 1
  do

    if (pos > len(pha_name)) exit

    ! parse phasespace type

    if (pha_name(pos:pos) == '_') then
      if ((ndims == 0) .or. (pos >=len(pha_name))) then
        ierr = pos
        return
      endif

      select case (trim(pha_name(pos+1:)))

      case( p_psext_normal )
        ps_type = p_pstype_normal

      case( p_psext_mass )
        ps_type = p_pstype_mass

      case( p_psext_ene )
        ps_type = p_pstype_ene

      case( p_psext_q1 )
        ps_type = p_pstype_q1

      case( p_psext_q2 )
        ps_type = p_pstype_q2

      case( p_psext_q3 )
        ps_type = p_pstype_q3

      case( p_psext_abs )
        ps_type = p_pstype_abs

      case( p_psext_j1 )
        ps_type = p_pstype_j1

      case( p_psext_j2 )
        ps_type = p_pstype_j2

      case( p_psext_j3 )
        ps_type = p_pstype_j3
      case default
        ierr = -2
        return

      end select

      exit
    endif

    ! check dimensions

    if (ndims == 3) then
      ierr = pos
      return
    endif

    ! parse phasespace axis

    select case (pha_name(pos:pos))
    case ("x")
      x_or_p(ndims+1) = p_phase_x
      ! only addition here, allow for x-cartesian
      if (pos < len(pha_name)) then
        if (pha_name(pos+1:pos+1) == "c") then
          x_or_p(ndims+1) = p_phase_xc
          pos = pos+1
        endif
      endif
      read_dim = .true.
    case ("p")
      x_or_p(ndims+1) = p_phase_p
      read_dim = .true.
    case ("g")
      x_or_p(ndims+1) = p_phase_g
      ! Check to see if phasespace is log gamma, "gl"
      if (pos < len(pha_name)) then
        if (pha_name(pos+1:pos+1) == "l") then
          ! Check to see if that "l" actually corresponds to angular momentum
          if (pos+1 < len(pha_name)) then
            if ( pha_name(pos+2:pos+2) == "1" .or. &
                 pha_name(pos+2:pos+2) == "2" .or. &
                 pha_name(pos+2:pos+2) == "3" ) then
              ! "l" corresponds to angular momentum, move on
            else
              x_or_p(ndims+1) = p_phase_gl
              pos = pos+1
            endif
          else
            x_or_p(ndims+1) = p_phase_gl
            pos = pos+1
          endif
        endif
      endif
      xp_dim(ndims+1) = 0
      read_dim = .false.
    case ("l")
      x_or_p(ndims+1) = p_phase_l
      read_dim = .true.
    case ("k")
      ! Next two characters must be "el", corresponding to phasespace "kel",
      ! or log kinetic energy
      if (pos+2 <= len(pha_name)) then
        if (pha_name(pos+1:pos+2) /= "el") then
          ierr = pos
        endif
      else
        ierr = pos
      endif
      x_or_p(ndims+1) = p_phase_kel
      pos = pos+2
      xp_dim(ndims+1) = 0
      read_dim = .false.
    case default
      ierr = pos
      return
    end select
    pos = pos+1

    ! read dimension for phasespace axis if necessary

    if (read_dim) then
      if (pos >len(pha_name)) then
        ierr = pos
        return
      endif

      ! parse the dimension
      select case (pha_name(pos:pos))
      case ("1")
        xp_dim(ndims+1) = 1
      case ("2")
        xp_dim(ndims+1) = 2
      case ("3")
        xp_dim(ndims+1) = 3
      case default
        ierr = pos
        return
      end select
      pos = pos+1

    endif

    ndims = ndims + 1

  enddo

  ! check if phasespace is valid and save results
  do pos =1, ndims

    ! check for x_dim > p_x_dim
    if ((x_or_p(pos) == p_phase_x) .and. (xp_dim(pos) > p_x_dim)) then
      ierr = -1
      return
    endif

    ! check for x_dim > p_max_dim for x-cartesian
    if ((x_or_p(pos) == p_phase_xc) .and. (xp_dim(pos) > p_max_dim)) then
      ierr = -1
      return
    endif

    ! no x-cartesian for 1st dimension, default to regular
    if ((x_or_p(pos) == p_phase_xc) .and. (xp_dim(pos) == 1)) then
      x_or_p(pos) = p_phase_x
      SCR_ROOT("   - NOTE: xc phasespace only available in dims 2 and 3,")
      SCR_ROOT("     xc1 defaulting to x1.")
    endif

    ! check for repeated axis
    do i = 1, pos-1

      if ((x_or_p(pos) == x_or_p(i)) .and. (xp_dim(pos)==xp_dim(i))) then
        ierr = -3
        return
      endif
    enddo

    ! for historical reasons these are in reverse order
    o_x_or_p(ndims - pos + 1) = x_or_p(pos)
    o_xp_dim(ndims - pos + 1) = xp_dim(pos)

  enddo

end subroutine parse_phasespace_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine allocate_objs_pspace_diag_cyl_modes( this )

  implicit none

  class( t_phasespace_diag_cyl_modes ), intent(inout) :: this

  if ( .not. associated( this%phasespace_list ) ) then

    allocate( t_phasespace_list_cyl_modes :: this % phasespace_list )
    allocate( t_phasespace_list_cyl_modes :: this % pha_ene_bin_list )
    allocate( t_phasespace_list_cyl_modes :: this % pha_cell_avg_list )
    allocate( t_phasespace_list_cyl_modes :: this % pha_time_avg_list )

  endif

  ! If needed in the future
  ! call this % phasespace_list % allocate_objs()
  ! call this % pha_ene_bin_list % allocate_objs()
  ! call this % pha_cell_avg_list % allocate_objs()
  ! call this % pha_time_avg_list % allocate_objs()

end subroutine allocate_objs_pspace_diag_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initializes a phasepace list from the supplied string array
!-----------------------------------------------------------------------------------------
subroutine init_pspace_list_cyl_modes( this, list, phasespaces, msg, time_average )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_phasespace_diag_cyl_modes ), intent(in) :: this
  class( t_phasespace_list ), intent(inout) :: list
  character(len=*), dimension(:), intent(in) :: phasespaces
  character(len=*), intent(in) :: msg
  logical, intent(in), optional :: time_average

  ! local variables
  integer :: i, j, npha, ierr
  integer :: ndims
  integer, dimension(p_max_phasespace_dims) :: x_or_p, xp_dim
  integer :: ps_type

  npha = size(phasespaces)
  i = 1
  do

    if (i > npha) exit
    if (trim(phasespaces(i)) == "-") exit

    call parse_phasespace_cyl_modes( trim(phasespaces(i)), ndims, x_or_p, xp_dim, &
                            ps_type, ierr)

    if (ierr/=0) then
      print *, "(*error*) Invalid ",trim(msg)," requested ", trim(phasespaces(i))
      select case (ierr)
      case(-1)

        print *, "(*error*) Invalid x dimension"
      case(-2)
        print *, "(*error*) Invalid deposit parameter"

      case(-3)
        print *, "(*error*) Repeated axis"

      case default
      end select
      stop
    endif

    ! If doing time averages special attention is required
    if ( present(time_average) ) then
      ! check if doing time average of an autoranged axis
      do j = 1, ndims

        if ((x_or_p(j) == p_phase_p) .and. (this%if_p_auto(xp_dim(j)))) then
          print *, "(*error*) time average with an autorange axis"
          print *, "(*error*) Cannot do a time average phasespace with a p", xp_dim(j), " axis"
          print *, "(*error*) with if_ps_p_auto(",xp_dim(j)," = .true. (",&
          trim(phasespaces(i)),")"

          stop
        endif

        if (((x_or_p(j) == p_phase_g) .or. (x_or_p(j) == p_phase_gl)) .and. &
            (this%if_gamma_auto)) then
          print *, "(*error*) time average with an autorange axis"
          print *, "(*error*) Cannot do a time average phasespace with a gamma or log10(gamma) axis"
          print *, "(*error*) with if_ps_gamma_auto = .true. (",trim(phasespaces(i)),")"

          stop
        endif

        if ((x_or_p(j) == p_phase_l) .and. (this%if_l_auto(xp_dim(j)))) then
          print *, "(*error*) time average with an autorange axis"
          print *, "(*error*) Cannot do a time average phasespace with a l", xp_dim(j), " axis"
          print *, "(*error*) with if_l_auto(",xp_dim(j)," = .true. (",&
          trim(phasespaces(i)),")"

          stop
        endif

        if ((x_or_p(j) == p_phase_kel) .and. (this%if_ke_auto)) then
          print *, "(*error*) time average with an autorange axis"
          print *, "(*error*) Cannot do a time average phasespace with a log10(ke) axis"
          print *, "(*error*) with if_ps_ke_auto = .true. (",trim(phasespaces(i)),")"

          stop
        endif

      enddo
    endif

    call list % add_to_list( ndims, x_or_p, xp_dim, ps_type )
    i = i + 1

  enddo

end subroutine init_pspace_list_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_pspace_params_cyl_modes( this, phasespace, g_space, params )
!-----------------------------------------------------------------------------------------
!   determines the parameters for the phasespace requested
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_phasespace_diag_cyl_modes ), intent(in) :: this
  type( t_phasespace ), intent(in) :: phasespace
  type( t_space ), intent(in) :: g_space
  type (t_phasespace_params), intent(out) :: params

  real(p_k_part), dimension(2,p_x_dim) :: ltotal_xmoved
  real(p_k_part), dimension(p_x_dim) :: lxmin, lxmax
  integer :: i
  integer, parameter :: izero = ichar('0')
  integer :: nxax, npax, ngax, nlax, nxcax, nkeax

  ! get phasespace name and boundaries

  params%ndims = phasespace%ndims
  params%name = phasespace%name

  params%long_name = ''
  params%offset_x = (/ 0.5_p_double, 0.5_p_double, 0.5_p_double /)

  nxax = 0
  npax = 0
  ngax = 0
  nlax = 0
  nxcax = 0
  nkeax = 0

  do i = 1, phasespace%ndims

    select case (phasespace%x_or_p(i))
    case (p_phase_x) ! x

      ! check if supplied limits are ok
      if (this%xmin(phasespace%xp_dim(i)) < this%xmax(phasespace%xp_dim(i))) then
        params%min(i) = this%xmin(phasespace%xp_dim(i))
        params%max(i) = this%xmax(phasespace%xp_dim(i))
      else
        ! use full simulation box size
        lxmin = real( xmin_initial(g_space), p_k_part )
        lxmax = real( xmax_initial(g_space), p_k_part )
        params%min(i) = lxmin(phasespace%xp_dim(i))
        params%max(i) = lxmax(phasespace%xp_dim(i))
      endif

      ! correct for moving window
      if (if_move(g_space,phasespace%xp_dim(i))) then
        ltotal_xmoved = real( total_xmoved(g_space), p_k_part )
        params%min(i) = params%min(i) + ltotal_xmoved(1,phasespace%xp_dim(i))
        params%max(i) = params%max(i) + ltotal_xmoved(2,phasespace%xp_dim(i))
        params%move_c(i) = .true.
      else
        params%move_c(i) = .false.
      endif

      params%n(i) = this%nx(phasespace%xp_dim(i))
      params%long_name = 'x_'//char(izero+phasespace%xp_dim(i))//trim(params%long_name)
      params%xname(i)  = 'x'//char(izero + phasespace%xp_dim(i))
      params%xlabel(i) = 'x_'//char(izero + phasespace%xp_dim(i))
      params%xunits(i) = 'c / \omega_p'
      params%offset_t_ax(i) = 0.0_p_double

      nxax = nxax + 1

    case (p_phase_p) ! p
      if ( this%pmin(phasespace%xp_dim(i)) < this%pmax(phasespace%xp_dim(i)) ) then
        params%min(i) = this%pmin(phasespace%xp_dim(i))
        params%max(i) = this%pmax(phasespace%xp_dim(i))
      else
        ! set to some reasonable (but arbitrary ) value

        ! if not set in the input file
        params%min(i) = -30
        params%max(i) =  30
      endif
      params%move_c(i) = .false.

      params%n(i) = this%np(phasespace%xp_dim(i))
      params%long_name = 'p_'//char(izero+phasespace%xp_dim(i))//trim(params%long_name)
      params%xname(i)  = 'p'//char(izero + phasespace%xp_dim(i))
      params%xlabel(i) = 'p_'//char(izero + phasespace%xp_dim(i))
      params%xunits(i) = 'm_e c'
      params%offset_t_ax(i) = -0.5_p_double

      npax = npax + 1

    case (p_phase_g) ! gamma
      params%long_name = '\gamma '//trim(params%long_name)
      params%move_c(i) = .false.
      params%n(i) = this%ngamma
      params%min(i) = this%gammamin
      params%max(i) = this%gammamax
      params%xname(i)  = 'gamma'
      params%xlabel(i) = '\gamma'
      params%xunits(i) = ''
      params%offset_t_ax(i) = -0.5_p_double

      ngax = ngax + 1

    case (p_phase_gl) ! log(gamma)
      params%long_name = 'log_{10}(\gamma)'//trim(params%long_name)
      params%move_c(i) = .false.
      params%n(i) = this%ngamma
      params%min(i) = log10(this%gammamin)
      params%max(i) = log10(this%gammamax)
      params%xname(i)  = 'gamma'
      params%xlabel(i) = 'log_{10}(\gamma)'
      params%xunits(i) = ''
      params%offset_t_ax(i) = -0.5_p_double

      ngax = ngax + 1

    case (p_phase_l) ! l (angular momentum)
      if ( this%lmin(phasespace%xp_dim(i)) < this%lmax(phasespace%xp_dim(i)) ) then
        params%min(i) = this%lmin(phasespace%xp_dim(i))
        params%max(i) = this%lmax(phasespace%xp_dim(i))
      else
        ! set to some reasonable (but arbitrary ) value

        ! if not set in the input file
        params%min(i) = -30
        params%max(i) =  30
      endif
      params%move_c(i) = .false.

      params%n(i) = this%nl(phasespace%xp_dim(i))
      params%long_name = 'l'//char(izero+phasespace%xp_dim(i))//trim(params%long_name)
      params%xname(i)  = 'l'//char(izero + phasespace%xp_dim(i))
      params%xlabel(i) = 'l_'//char(izero + phasespace%xp_dim(i))
      params%xunits(i) = 'm_e c^2 \omega_p^{-1}'
      ! This is not well-defined, since it uses both x and p
      params%offset_t_ax(i) = -0.25_p_double

      nlax = nlax + 1

    case (p_phase_kel) ! log(ke)
      params%long_name = 'log_{10}(KE)'//trim(params%long_name)
      params%move_c(i) = .false.
      params%n(i) = this%nke
      params%min(i) = log10(this%kemin)
      params%max(i) = log10(this%kemax)
      params%xname(i)  = 'KE'
      params%xlabel(i) = 'log_{10}(KE)'
      params%xunits(i) = ''
      params%offset_t_ax(i) = -0.5_p_double

      nkeax = nkeax + 1

    case (p_phase_xc) ! x-cartesian

      ! check if supplied limits are ok
      if (this%xcmin(phasespace%xp_dim(i)) < this%xcmax(phasespace%xp_dim(i))) then
        params%min(i) = this%xcmin(phasespace%xp_dim(i))
        params%max(i) = this%xcmax(phasespace%xp_dim(i))
      else
        ! use full simulation box size
        lxmin = real( xmin_initial(g_space), p_k_part )
        lxmax = real( xmax_initial(g_space), p_k_part )
        ! use "reflected" lower boundary since this is cylindrical symmetry
        ! (only dims 2 and 3 allowed for x-cartesian)
        params%min(i) = - lxmax(2)
        ! use standard upper limit
        params%max(i) = lxmax(2)
      endif

      ! no moving window since we only allow for dims 2 and 3
      params%move_c(i) = .false.

      params%n(i) = this%nxc(phasespace%xp_dim(i))
      params%long_name = 'xcart_'//char(izero+phasespace%xp_dim(i))//trim(params%long_name)
      params%xname(i)  = 'xcart'//char(izero + phasespace%xp_dim(i))
      params%xlabel(i) = 'xcart_'//char(izero + phasespace%xp_dim(i))
      params%xunits(i) = 'c / \omega_p'
      params%offset_t_ax(i) = 0.0_p_double

      nxcax = nxcax + 1

    end select

  enddo

  select case( phasespace%ps_type )

  case(p_pstype_normal)
    params%offset_t = 0.0_p_double
    continue

  case(p_pstype_mass)
    params%long_name = trim(params%long_name)//" ("//p_psext_mass//")"
    params%offset_t = 0.0_p_double

  case(p_pstype_abs)
    params%long_name = trim(params%long_name)//" ("//p_psext_abs//")"
    params%offset_t = 0.0_p_double

  case(p_pstype_ene)
    params%long_name = trim(params%long_name)//" ("//p_psext_ene//")"
    params%offset_t = -0.5_p_double

  case(p_pstype_q1)
    params%long_name = trim(params%long_name)//" ("//p_psext_q1//")"
    params%offset_t = -0.5_p_double

  case(p_pstype_q2)
    params%long_name = trim(params%long_name)//" ("//p_psext_q2//")"
    params%offset_t = -0.5_p_double

  case(p_pstype_q3)
    params%long_name = trim(params%long_name)//" ("//p_psext_q3//")"
    params%offset_t = -0.5_p_double

  case(p_pstype_j1)
    params%long_name = trim(params%long_name)//" ("//p_psext_j1//")"
    params%offset_t = -0.5_p_double

  case(p_pstype_j2)
    params%long_name = trim(params%long_name)//" ("//p_psext_j2//")"
    params%offset_t = -0.5_p_double

  case(p_pstype_j3)
    params%long_name = trim(params%long_name)//" ("//p_psext_j3//")"
    params%offset_t = -0.5_p_double

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "(*error*) Unsupported phasespace type ", __FILE__, ":", __LINE__
      write(0,*) "(*error*) aborting..."
    endif
    call abort_program(p_err_invalid)

  end select

  ! temporary until writing a routine that calculates
  ! proper units

  params%unit = "a.u."

end subroutine get_pspace_params_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Determine phasespace size
!-----------------------------------------------------------------------------------------
function size_phasespace_cyl_modes( this, phasespace )
!-----------------------------------------------------------------------------------------
  implicit none

  class( t_phasespace_diag_cyl_modes ), intent(in) :: this
  type( t_phasespace ), intent(in) :: phasespace
  integer :: size_phasespace_cyl_modes

  integer :: i

  size_phasespace_cyl_modes = this % t_phasespace_diagnostics % size( phasespace )

  do i = 1, phasespace%ndims
    select case (phasespace%x_or_p(i))
    case (p_phase_xc) !x-cartesian

      size_phasespace_cyl_modes = size_phasespace_cyl_modes * this%nxc(phasespace%xp_dim(i))
    end select
  enddo

end function size_phasespace_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine add_phasespace_to_list_cyl_modes( list, ndims, x_or_p, xp_dim, ps_type )
!-----------------------------------------------------------------------------------------
! adds a phasespace to the phasespace list
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_phasespace_list_cyl_modes ), intent( inout )  ::  list
  integer, intent(in) :: ndims
  integer, dimension(:), intent(in) :: x_or_p, xp_dim
  integer, intent(in) :: ps_type

  character( len=32 ) :: name = ""
  integer, parameter  :: izero =  ichar('0')
  integer :: i

  ! call superclass method, which creates the phasespace

  call list % t_phasespace_list % add_to_list( ndims, x_or_p, xp_dim, ps_type )

  ! modify the name, allowing for the additional case of x-cartesian

  ASSERT(ndims>0)

  name = ""
  do i=1, ndims
    if (xp_dim(i) /= 0) then
      name = char(izero + xp_dim(i)) // trim(name)

    endif
    select case (x_or_p(i))
    case (p_phase_x)

      name =  "x"//trim(name)

    case (p_phase_p)
      name =  "p"//trim(name)
    case (p_phase_g)
      name =  "gamma"//trim(name)
    case (p_phase_gl)
      name =  "log_gamma"//trim(name)
    case (p_phase_l)
      name =  "l"//trim(name)
    case (p_phase_kel)
      name =  "log_ke"//trim(name)
    case (p_phase_xc)
      name =  "xc"//trim(name)
    end select
  enddo

  select case( ps_type )

  case(p_pstype_normal)
    continue

  case(p_pstype_mass)
    name = trim(name)//"_"//p_psext_mass

  case(p_pstype_abs)
    name = trim(name)//"_"//p_psext_abs

  case(p_pstype_ene)
    name = trim(name)//"_"//p_psext_ene

  case(p_pstype_q1)
    name = trim(name)//"_"//p_psext_q1

  case(p_pstype_q2)
    name = trim(name)//"_"//p_psext_q2

  case(p_pstype_q3)
    name = trim(name)//"_"//p_psext_q3

  case(p_pstype_j1)
    name = trim(name)//"_"//p_psext_j1

  case(p_pstype_j2)
    name = trim(name)//"_"//p_psext_j2

  case(p_pstype_j3)
    name = trim(name)//"_"//p_psext_j3

  case default
    ERROR("Invalid phasespace type, or ")
    ERROR("phasespace type not implemented")
    call abort_program(p_err_invalid)

  end select

  list%tail%name  = name

end subroutine add_phasespace_to_list_cyl_modes
!-----------------------------------------------------------------------------------------

end module m_species_phasespace_cyl_modes

!-----------------------------------------------------------------------------------------
subroutine get_phasespace_axis_cyl_modes( spec, xp, l, lp, x_or_p, xp_dim )
!-----------------------------------------------------------------------------------------

  use m_parameters
  use m_diagnostic_utilities, only : p_diag_prec
  use m_species_cyl_modes_define, only : t_species_cyl_modes
  use m_species_phasespace_cyl_modes, only : p_phase_xc

  implicit none

  class(t_species_cyl_modes), intent(in) :: spec
  real(p_diag_prec), dimension(:), intent(out) :: xp
  integer, intent(in) :: l, lp
  integer, intent(in) :: x_or_p
  integer, intent(in) :: xp_dim

  select case(x_or_p)
  case (p_phase_xc) ! x-cartesian

    ! Get correct coordinate for 2nd/3rd cartesian dimensions
    select case(xp_dim)
    case(2)
      call spec % get_position( 3, l, lp, xp )
    case(3)
      call spec % get_position( 4, l, lp, xp )
    case default
      ERROR("Should not get to this point, 1 not allowed")
      call abort_program(p_err_invalid)
    end select

  case default
    call spec % t_species % get_phasespace_axis( xp, l, lp, x_or_p, xp_dim )
  end select

end subroutine get_phasespace_axis_cyl_modes
!-----------------------------------------------------------------------------------------
