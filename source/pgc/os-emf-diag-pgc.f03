#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_pgc_diag

use m_emf_diag

private

! PGC specific diagnostics
character(len=10), dimension(5), parameter :: p_report_quants_pgc = &
   (/  'a_mod     ', 'fp1       ', 'fp2       ', 'fp3       ', 'chi       '/)

integer, parameter :: p_a_mod = 1, p_fp1 = 2, p_fp2 = 3, p_fp3 = 4, p_chi = 5

!-------------------------------------------------------------------------------
! parameters for timing events
!-------------------------------------------------------------------------------
! longitudinal profiles
integer :: pgc_advance_env  = 0      ! event for advancing envelope
integer :: pgc_advance_part          ! event for whole advancing

type, extends( t_diag_emf ) :: t_diag_emf_pgc

contains

  procedure :: avail_report_quants => avail_report_quants_emf_pgc
  procedure :: init_report_quants  => init_report_quants_emf_pgc
  procedure :: init                => init_diag_emf_pgc
  procedure :: quant_offset        => quant_offset_emf_pgc

end type t_diag_emf_pgc

public :: t_diag_emf_pgc
public :: p_a_mod, p_fp1, p_fp2, p_fp3, p_chi
public :: pgc_advance_env, pgc_advance_part

contains

! Returns the integer value that must be subtracted from the diagnostic
! quantity to match the local quantity definitions

function quant_offset_emf_pgc( this )

  class( t_diag_emf_pgc ), intent(in) :: this
  integer :: quant_offset_emf_pgc

  ! Offset by the number of quantities defined by the superclass
  quant_offset_emf_pgc = this % t_diag_emf % avail_report_quants( )

end function quant_offset_emf_pgc


function avail_report_quants_emf_pgc( this )

  class( t_diag_emf_pgc ), intent(in) :: this
  integer :: avail_report_quants_emf_pgc

  integer :: n

  ! get the reports available on the superclass
  n = this % t_diag_emf % avail_report_quants( )

  ! add the additional local reports
  n = n + size( p_report_quants_pgc )

  avail_report_quants_emf_pgc = n

end function avail_report_quants_emf_pgc


subroutine init_report_quants_emf_pgc( this )

  class( t_diag_emf_pgc ), intent(inout) :: this
  integer :: n

  ! print *, '[diag-pgc] In init_report_quants_emf_pgc'

  if ( .not. associated( this % report_quants ) ) then
    allocate( this % report_quants( this % avail_report_quants( ) ) )
  endif

  n = this % t_diag_emf % avail_report_quants(  )

  ! Initialize superclass report quantities
  call this % t_diag_emf % init_report_quants()

  ! add the additional reports
  this % report_quants( n+1 : n+size(p_report_quants_pgc) ) = p_report_quants_pgc

end subroutine init_report_quants_emf_pgc

!-------------------------------------------------------------------------------
subroutine init_diag_emf_pgc( this, no_ext_fld, part_fld_alloc, interpolation )

  use m_logprof

  use m_vdf_define, only : t_vdf_report

  implicit none

  class( t_diag_emf_pgc ), intent(inout) :: this
  logical, intent(in) :: no_ext_fld
  logical, intent(in) :: part_fld_alloc
  integer, intent(in) :: interpolation

  type(t_vdf_report), pointer :: report

  integer :: quant

  ! Initialize superclass reports
  call this % t_diag_emf % init( no_ext_fld, part_fld_alloc, interpolation )

  ! Normal reports
  report => this%reports
  do
    if ( .not. associated( report ) ) exit

    quant = report%quant - this % quant_offset( )

  select case ( quant )

    case ( p_a_mod )
    report%label    = '|a_0|'
    report%units    = 'm_e c e^{-1}'

    case ( p_fp1 )
    report%label    = '\nabla_{x_1} a^2'
    report%units    = 'm_e^2 c^3 e^{-2} / \omega_p '

    case ( p_fp2 )
    report%label    = '\nabla_{x_2} a^2'
    report%units    = 'm_e^2 c^3 e^{-2} / \omega_p '

    case ( p_fp3 )
    report%label    = '\nabla_{x_3} a^2'
    report%units    = 'm_e^2 c^3 e^{-2} / \omega_p '

    case ( p_chi )
    report%label    = '\chi'
    report%units    = 'adimensional'

    case default
    ! must be a superclass diagnostic, ignore
    ! print *, 'In setup_diag_pgc, unknown report = ', report % name
    continue

  end select

    ! process next report
    report => report % next
  enddo

  if (pgc_advance_env==0) then
    pgc_advance_env  = create_event('PGC advance envelope ')
    pgc_advance_part = create_event('PGC advance particles')
  endif

end subroutine init_diag_emf_pgc
!-------------------------------------------------------------------------------

end module m_emf_pgc_diag
