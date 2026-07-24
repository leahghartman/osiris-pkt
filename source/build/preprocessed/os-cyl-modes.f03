# 1 "cyl_modes/os-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-cyl-modes.f03" 2
# 1 "./os-config.h" 1
! Configuration file for osiris

! ----------------------------------------------------------------------------------------
! Algorithm options
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! System options
! ----------------------------------------------------------------------------------------

! MPI supports MPI_IN_PLACE in global operations
!#define 1

! Use OpenMP
!#define 1

! Use SIMD optimized code
!#define SIMD
!#define SIMD_SSE
!#define SIMD_AVX
!#define SIMD_BGQ
!#define SIMD_MIC

! Use SION for checkpointing (not available on all systems)
!#define __RST_IO__ = __RST_SION__

! Use log files
!#define __USE_LOG__

! Use MPE for logging and profiling
!#define __USE_MPE__

! Compiler does not fully support the sizeof intrinsic which is a Fortran 2003 feature
!#define __NO_SIZEOF__

! Use the PAPI library for profiling
!#define __USE_PAPI__


! ----------------------------------------------------------------------------------------
! Distribution Options
! ----------------------------------------------------------------------------------------
!
! Optional modules to be removed for distribution
! These are all turned on by default unless the __DISTRO__ preprocessor macro is defined
!
! ----------------------------------------------------------------------------------------



! Include Ionization module


! Include binary collisions module


! Include particle tracking module


! Include perfectly matched layers boundary conditions for EMF


! Include spin advance module which is disabled by default
! #define __HAS_SPIN__

! Include a debug flag for cylindrical modes simulations
! #define __CYL_MODES_DEBUG__

! Compile functions for reporting B fields in the RaDiO module
! #define __HAS_RAD_BFLD__
# 2 "cyl_modes/os-cyl-modes.f03" 2
# 1 "./os-preprocess.fpp" 1
!
! File: os-preprocess.fpp
!
! A set of preprocessing macros for osiris, and the fpp/cpp preprocessor
!
!



! Macros for the IBM CPP/GNU preprocessor




! Assertion macro



! Debug functions
# 33 "./os-preprocess.fpp"
! SCR hostname functions







! LOG functions




! ERROR functions






! WARNING functions




! functions for restart io
# 3 "cyl_modes/os-cyl-modes.f03" 2

module m_cyl_modes
# 1 "./memory/memory.h" 1
! Include file for the memory module
! Files must include this file using #include "memory/memory.h" rather than just using the module
! the #include must be placed where the module use statement would usually be i.e.
!
! module module2
!
! use module1
! #include "memory/memory.h"
!
# 19 "./memory/memory.h"
use memory
# 6 "cyl_modes/os-cyl-modes.f03" 2

use m_parameters
use m_vdf_define
use m_vdf
use m_vdf_comm
use m_vdf_report
use m_space
use m_node_conf
use m_vdf_memory
use m_grid_define
use m_time_step
use m_vdf_average
use m_vdf_math, only : reduce, add, mult, pow, add_pow
use m_vdf_smooth, only: smooth, t_smooth
use stringutil
use m_restart

implicit none

private

! constants defined here
integer, public, parameter:: p_max_cyl_modes = 15
integer, public, parameter:: p_real = 0
integer, public, parameter:: p_imag = 1


type :: t_cyl_modes

  type( t_vdf ), dimension(:), pointer :: pf_re => null()
  type( t_vdf ), dimension(:), pointer :: pf_im => null()

end type t_cyl_modes

interface setup
  module procedure setup_cyl_modes
end interface

interface new_copy
  module procedure new_copy_cyl_modes
end interface

interface copy_component
  module procedure copy_component_cyl_modes
end interface

interface cleanup
  module procedure cleanup_cyl_modes
end interface

interface restart_write
  module procedure restart_write_cyl_modes
end interface

interface restart_read
  module procedure restart_read_cyl_modes
end interface

interface move_window
  module procedure move_window_cyl_modes
end interface

interface update_boundary
  module procedure update_boundary_cyl_modes
end interface

interface report_cyl_modes
  module procedure report_cyl_modes_
end interface

interface reshape_nocopy
  module procedure reshape_cyl_modes_nocopy
  module procedure reshape_all_cyl_modes_nocopy
end interface

interface reshape_copy
  module procedure reshape_cyl_modes_copy
  module procedure reshape_all_cyl_modes_copy
end interface

interface reduce
  module procedure reduce_cyl_modes
end interface

interface check_nan
  module procedure check_nan_cyl_modes
end interface

interface total
  module procedure total_cyl_modes
end interface

interface add
  module procedure add_cyl_modes
end interface

interface mult
  module procedure mult_cyl_modes
end interface

interface pow
  module procedure pow_cyl_modes
end interface

interface add_pow
  module procedure add_pow_cyl_modes
end interface

interface div
  module procedure div_cyl_modes
end interface

interface smooth
  module procedure smooth_cyl_modes
end interface

public :: t_cyl_modes, setup, new_copy, cleanup, restart_write, restart_read, move_window
public :: update_boundary, report_cyl_modes, reshape_nocopy, reshape_copy
public :: reduce, check_nan, total, add, mult, pow, add_pow, div, smooth, copy_component

contains

!-----------------------------------------------------------------------------------------
! Link pf_re(0) to source, then allocate remainder of arrays
!-----------------------------------------------------------------------------------------
subroutine setup_cyl_modes( this, n_modes, source )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  integer, intent(in) :: n_modes
  type( t_vdf ), intent(in) :: source

  integer :: i

  call alloc(this%pf_re, (/ 0 /), (/ n_modes /),"cyl_modes/os-cyl-modes.f03",141)
  if (n_modes>0) call alloc(this%pf_im, (/ 1 /), (/ n_modes /),"cyl_modes/os-cyl-modes.f03",142)

  ! mode 0 has special treatment because the real part is
  ! just a link to the main source vdf
  call link_vdf( this%pf_re(0), source )
  ! No imaginary mode 0

  ! allocate remaining modes
  do i = 1, n_modes
    call this % pf_re(i) % new( source, zero = .true. )
    call this % pf_im(i) % new( source, zero = .true. )
  enddo

end subroutine setup_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Add capability for initializing all arrays (including mode 0) and copying data
!-----------------------------------------------------------------------------------------
subroutine new_copy_cyl_modes( this, source, f_dim, zero, copy, fc )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(in) :: source
  integer, intent(in), optional :: f_dim
  logical, intent(in), optional :: zero
  logical, intent(in), optional :: copy
  integer, intent(in), optional :: fc

  integer :: f_dim_, n_modes, i
  logical :: zero_, copy_

  if (present(f_dim)) then
    if (f_dim < 1) then
      write(err_buf__,*) "Field components (f_dim) must be >= 1";call err__("cyl_modes/os-cyl-modes.f03",177)
      call abort_program(p_err_invalid)
    endif
    f_dim_ = f_dim
  else
    f_dim_ = source%pf_re(0)%f_dim()
  endif

  if ( present(zero) ) then
    zero_ = zero
  else
    zero_ = .false.
  endif

  if ( present( copy ) ) then
    copy_ = copy
    if ( copy_ .and. (.not. present( fc )) ) then
      if ( source%pf_re(0)%f_dim() /= f_dim_) then
        write(0,*) "(*error*) new_copy_cyl_modes was called with copy = .true. but target and source"
        write(0,*) "(*error*) field dimensions do not match (may need to set the fc parameter)."
        call abort_program(p_err_invalid)
      endif
    endif
  else
    copy_ = .false.
  endif

  n_modes = ubound( source%pf_re, 1 )

  call alloc(this%pf_re, (/ 0 /), (/ n_modes /),"cyl_modes/os-cyl-modes.f03",206)
  if (n_modes>0) call alloc(this%pf_im, (/ 1 /), (/ n_modes /),"cyl_modes/os-cyl-modes.f03",207)

  if ( present(fc) ) then

    call this % pf_re(0) % new( source%pf_re(0), f_dim=f_dim_, zero=zero_, copy=copy_, fc=fc )
    do i = 1, n_modes
      call this % pf_re(i) % new( source%pf_re(i), f_dim=f_dim_, zero=zero_, copy=copy_, fc=fc )
      call this % pf_im(i) % new( source%pf_im(i), f_dim=f_dim_, zero=zero_, copy=copy_, fc=fc )
    enddo

  else

    call this % pf_re(0) % new( source%pf_re(0), f_dim=f_dim_, zero=zero_, copy=copy_ )
    do i = 1, n_modes
      call this % pf_re(i) % new( source%pf_re(i), f_dim=f_dim_, zero=zero_, copy=copy_ )
      call this % pf_im(i) % new( source%pf_im(i), f_dim=f_dim_, zero=zero_, copy=copy_ )
    enddo

  endif

end subroutine new_copy_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! copy one component of all the vdfs in source into this.
! assume that the arrays all have the same dimensions.
!-----------------------------------------------------------------------------------------
subroutine copy_component_cyl_modes( this, source, source_component, dest_component )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(in) :: source
  integer, intent(in) :: source_component, dest_component

  integer :: n_modes, i

  n_modes = ubound( source%pf_re, 1 )

  call this % pf_re(0) % copy( source % pf_re(0), source_component, dest_component )

  do i = 1, n_modes

    call this % pf_re(i) % copy( source % pf_re(i), source_component, dest_component )

    call this % pf_im(i) % copy( source % pf_re(i), source_component, dest_component )

  enddo

end subroutine copy_component_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_cyl_modes( this, cleanup_0 )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  logical, intent(in), optional :: cleanup_0

  integer :: n_modes, i

  if (.not. associated( this%pf_re ) ) then
    return
  endif

  n_modes = ubound( this%pf_re, 1 )

  ! Mode 0 has special treatment because the real part is
  ! usually just a link to the main source vdf.
  ! call this % pf_re(0) % cleanup() ! usually not needed
  ! No imaginary mode 0

  ! But cleanup mode 0 if asked for
  if ( present(cleanup_0) ) then
    if (cleanup_0) then
      call this % pf_re(0) % cleanup()
    endif
  endif

  ! cleanup remaining modes
  do i = 1, n_modes
    call this % pf_re(i) % cleanup()
    call this % pf_im(i) % cleanup()
  enddo

  ! free local vdf array
  call freemem(this%pf_re,"cyl_modes/os-cyl-modes.f03",294)
  call freemem(this%pf_im,"cyl_modes/os-cyl-modes.f03",295)

end subroutine cleanup_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_write_cyl_modes( this, restart_handle )
!-----------------------------------------------------------------------------------------
! write object information into a restart file
!-----------------------------------------------------------------------------------------
  implicit none

  type( t_cyl_modes ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  integer :: mode, n_modes
  integer :: ierr
  character(len=*), parameter :: err_msg = 'error writing restart data for cylindrical mode object.'

  n_modes = ubound(this%pf_re,1)

  call restart_io_write("n_modes", n_modes, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-cyl-modes.f03",317)

  do mode = 1, n_modes
    call this % pf_re(mode) % write_checkpoint( restart_handle )
    call this % pf_im(mode) % write_checkpoint( restart_handle )
  enddo

end subroutine restart_write_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_cyl_modes( this, source, restart_handle )
!-----------------------------------------------------------------------------------------
! write object information into a restart file
!-----------------------------------------------------------------------------------------
  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: source
  type( t_restart_handle ), intent(in) :: restart_handle

  integer :: mode, n_modes
  integer :: ierr
  character(len=*), parameter :: err_msg = 'error reading restart data for cylindrical mode object.'

  call restart_io_read(n_modes, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-cyl-modes.f03",343)

  call alloc(this%pf_re, (/ 0 /), (/ n_modes /),"cyl_modes/os-cyl-modes.f03",345)
  if (n_modes>0) call alloc(this%pf_im, (/ 1 /), (/ n_modes /),"cyl_modes/os-cyl-modes.f03",346)

  call link_vdf( this%pf_re(0), source )

  do mode=1, n_modes
    call this % pf_re(mode) % read_checkpoint( restart_handle )
    call this % pf_im(mode) % read_checkpoint( restart_handle )
  enddo

end subroutine restart_read_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine move_window_cyl_modes( this, g_space )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space

  integer :: i, n_modes

  n_modes = ubound( this%pf_re, 1 )

  ! This is taken care of by normal routines
  ! call move_window( this%pf_re(0), g_space )

  do i = 1, n_modes
    call move_window( this%pf_re(i), g_space )
    call move_window( this%pf_im(i), g_space )
  enddo

end subroutine move_window_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_boundary_cyl_modes( this, update_type, no_co, send_msg, recv_msg )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  integer, intent(in) :: update_type
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, n_modes

  n_modes = ubound( this%pf_re, 1 )

  ! This is taken care of by normal routines
  ! call update_boundary( this%pf_re(0), update_type, no_co, send_msg, recv_msg )





  do i = 1, n_modes
    call update_boundary( this%pf_re(i), update_type, no_co, send_msg, recv_msg )
    call update_boundary( this%pf_im(i), update_type, no_co, send_msg, recv_msg )
  enddo





end subroutine update_boundary_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Process all reports for the quantity specified in report
!-----------------------------------------------------------------------------------------
subroutine report_cyl_modes_( report, cyl_modes_data, fc, g_space, grid, no_co, tstep, t )
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_vdf_report), pointer :: report
  type(t_cyl_modes), intent(in), target :: cyl_modes_data
  integer, intent(in) :: fc
  type(t_space), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: t

  integer :: n_modes, mode

  character(len=256) :: orig_filename, orig_path, prefix

  orig_filename = ''; orig_path = ''; prefix = ''

  n_modes = ubound(cyl_modes_data%pf_re,1)

  ! accumulate data for time averaged reports
  !if ( if_add_tavg_data( report, n(tstep) ) ) then
  ! call add_tavg_data( report, cyl_modes_data%pf_re(mode), fc, n(tstep) )
  !endif

  orig_filename = trim(report%fileLabel)
  if ( trim(orig_filename) /= '' ) prefix = '-'
  orig_path = trim(report%basePath)

  do mode = 0, n_modes

    write( report%fileLabel, '(A,I0,A)' ) trim(orig_filename)//trim(prefix), mode, "-re"
    write( report%basePath, '(A,I0,A)' ) trim(orig_path)//"MODE-", mode, "-RE"
    call report_vdf( report, cyl_modes_data%pf_re(mode), fc, g_space, grid, no_co, tstep, t )

    if (mode > 0) then
      write( report%fileLabel, '(A,I0,A)' ) trim(orig_filename)//trim(prefix), mode, "-im"
      write( report%basePath, '(A,I0,A)' ) trim(orig_path)//"MODE-", mode, "-IM"
      call report_vdf( report, cyl_modes_data%pf_im(mode), fc, g_space, grid, no_co, tstep, t )
    endif

  enddo ! mode

  report%fileLabel = trim(orig_filename)
  report%basePath = trim(orig_path)

end subroutine report_cyl_modes_
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Reshape cyl_modes object when node grids change, but do not copy data.
! This routine assumes that mode 0 has already been reshaped.
!-----------------------------------------------------------------------------------------
subroutine reshape_cyl_modes_nocopy( this, new_grid, source, dbg )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  class( t_grid ), intent(in) :: new_grid
  type( t_vdf ), intent(in) :: source

  ! Initialize the new array with NaN values
  logical, optional, intent(in) :: dbg

  integer :: i, n_modes

  n_modes = ubound( this%pf_re, 1 )

  ! Link mode 0 to the source vdf
  call link_vdf( this%pf_re(0), source )

  if ( present(dbg) ) then

    do i = 1, n_modes
      call reshape_nocopy( this%pf_re(i), new_grid, dbg )
      call reshape_nocopy( this%pf_im(i), new_grid, dbg )
    enddo

  else

    do i = 1, n_modes
      call reshape_nocopy( this%pf_re(i), new_grid )
      call reshape_nocopy( this%pf_im(i), new_grid )
    enddo

  endif

  if (n_modes > 0) then
    if ( this%pf_re(0)%nx_(1) /= this%pf_re(1)%nx_(1) .or. &
         this%pf_re(0)%nx_(2) /= this%pf_re(1)%nx_(2) ) then
      write(err_buf__,*) 'In reshape_cyl_modes_nocopy, the mode 0 vdf was not reshaped beforehand.';call err__("cyl_modes/os-cyl-modes.f03",511)
      call abort_program(p_err_invalid)
    endif
  endif

end subroutine reshape_cyl_modes_nocopy
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Reshape cyl_modes object when node grids change and redistribute the data through all
! nodes. Waits for all messages to arive and then unpacks them.
! This routine assumes that mode 0 has already been reshaped.
!-----------------------------------------------------------------------------------------
subroutine reshape_cyl_modes_copy( this, old_lb, new_grid, source, no_co, send_msg, &
                                   recv_msg )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  class( t_grid ), intent(in) :: old_lb, new_grid
  type( t_vdf ), intent(in) :: source
  class( t_node_conf ), intent(in) :: no_co
  type( t_vdf_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, n_modes

  n_modes = ubound( this%pf_re, 1 )

  ! Link mode 0 to the source vdf
  call link_vdf( this%pf_re(0), source )

  do i = 1, n_modes
    call reshape_copy( this%pf_re(i), old_lb, new_grid, no_co, send_msg, recv_msg )
    call reshape_copy( this%pf_im(i), old_lb, new_grid, no_co, send_msg, recv_msg )
  enddo

  if (n_modes > 0) then
    if ( this%pf_re(0)%nx_(1) /= this%pf_re(1)%nx_(1) .or. &
         this%pf_re(0)%nx_(2) /= this%pf_re(1)%nx_(2) ) then
      write(err_buf__,*) 'In reshape_cyl_modes_copy, the mode 0 vdf was not reshaped beforehand.';call err__("cyl_modes/os-cyl-modes.f03",550)
      call abort_program(p_err_invalid)
    endif
  endif

end subroutine reshape_cyl_modes_copy
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Reshape cyl_modes object when node grids change, but do not copy data.
!-----------------------------------------------------------------------------------------
subroutine reshape_all_cyl_modes_nocopy( this, new_grid, dbg )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  class( t_grid ), intent(in) :: new_grid

  ! Initialize the new array with NaN values
  logical, optional, intent(in) :: dbg

  integer :: i, n_modes

  n_modes = ubound( this%pf_re, 1 )

  if ( present(dbg) ) then

    call reshape_nocopy( this%pf_re(0), new_grid, dbg )

    do i = 1, n_modes
      call reshape_nocopy( this%pf_re(i), new_grid, dbg )
      call reshape_nocopy( this%pf_im(i), new_grid, dbg )
    enddo

  else

    call reshape_nocopy( this%pf_re(0), new_grid )

    do i = 1, n_modes
      call reshape_nocopy( this%pf_re(i), new_grid )
      call reshape_nocopy( this%pf_im(i), new_grid )
    enddo

  endif

end subroutine reshape_all_cyl_modes_nocopy
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Reshape cyl_modes object when node grids change and redistribute the data through all
! nodes. Waits for all messages to arive and then unpacks them.
!-----------------------------------------------------------------------------------------
subroutine reshape_all_cyl_modes_copy( this, old_lb, new_grid, no_co, send_msg, recv_msg )

  implicit none

  type( t_cyl_modes ), intent(inout) :: this
  class( t_grid ), intent(in) :: old_lb, new_grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_vdf_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, n_modes

  n_modes = ubound( this%pf_re, 1 )

  call reshape_copy( this%pf_re(0), old_lb, new_grid, no_co, send_msg, recv_msg )

  do i = 1, n_modes
    call reshape_copy( this%pf_re(i), old_lb, new_grid, no_co, send_msg, recv_msg )
    call reshape_copy( this%pf_im(i), old_lb, new_grid, no_co, send_msg, recv_msg )
  enddo

end subroutine reshape_all_cyl_modes_copy
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Reduce all pf arrays for an array of t_cyl_modes
!-----------------------------------------------------------------------------------------
subroutine reduce_cyl_modes( cyl_modes )
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_cyl_modes), dimension(:), intent(inout) :: cyl_modes

  type(t_vdf_arr), dimension(1:size(cyl_modes)) :: pf
  integer :: nt, n_modes, i, j

  nt = size(cyl_modes)
  n_modes = ubound( cyl_modes(1)%pf_re, 1 )

  do i = 0, n_modes

    do j = 1, nt
      pf(j)%v => cyl_modes(j)%pf_re(i)
      ! Point buffer directly to avoid errors with some Intel compilers
      pf(j)%v%buffer => cyl_modes(j)%pf_re(i)%buffer
    enddo
    call reduce(pf)

    if (i>0) then
      do j = 1, nt
        pf(j)%v => cyl_modes(j)%pf_im(i)
        ! Point buffer directly to avoid errors with some Intel compilers
        pf(j)%v%buffer => cyl_modes(j)%pf_im(i)%buffer
      enddo
      call reduce(pf)
    endif

  enddo

end subroutine reduce_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Check for nans in t_cyl_modes
!-----------------------------------------------------------------------------------------
subroutine check_nan_cyl_modes( cyl_modes, msg_in )
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_cyl_modes), intent(in) :: cyl_modes
  character( len = * ), intent(in), optional :: msg_in

  character(len=4) :: mode_string
  integer :: n_modes, mode

  n_modes = ubound( cyl_modes%pf_re, 1 )

  mode = 0
  write(mode_string,'(I4)') mode
  call cyl_modes % pf_re(mode) % check_nan(msg_in//' real mode = '//mode_string)

  do mode = 1, n_modes

    write(mode_string,'(I4)') mode
    call cyl_modes % pf_re(mode) % check_nan(msg_in//' real mode = '//mode_string)
    call cyl_modes % pf_im(mode) % check_nan(msg_in//' imag mode = '//mode_string)

  enddo

end subroutine check_nan_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! returns sum over all elements of r^r_pow * vdf_a()^pow
! multiplies by 1/2 for modes with m > 0 due to sin/cos factors in integral
! note, we do not multiply by 2*pi here
! this calculation is always performed in double precision
!-----------------------------------------------------------------------------------------
subroutine total_cyl_modes( cyl_a, total, gir_pos, pow, r_pow, include_guard_cells, off )

  implicit none

  type(t_cyl_modes), intent(in) :: cyl_a
  real(p_double), dimension(:), intent(out) :: total
  integer, intent(in) :: gir_pos
  integer, intent(in), optional :: pow, r_pow
  logical, intent(in), optional :: include_guard_cells
  real(p_double), dimension(:), intent(in), optional :: off

  integer, dimension(p_max_dim) :: vdf_lbound, vdf_ubound
  logical :: inc_guard_cells
  real(p_double), dimension(cyl_a%pf_re(0)%f_dim_) :: offs, r
  integer :: x_dim, f_dim, i, lpow, rpow, i1, i2, n_modes, mode, shift_i2, re_im, lb
  type(t_vdf), pointer :: vdf_a

  vdf_a => cyl_a%pf_re(0)
  x_dim = vdf_a%x_dim_
  f_dim = vdf_a%f_dim_
  n_modes = ubound( cyl_a%pf_re, 1 )

  ! process optional parameters
  ! field power
  if (present(pow)) then
    lpow = pow
  else
    lpow = 1
  endif

  ! radius power (can be 0 if r not desired)
  if (present(r_pow)) then
    rpow = r_pow
  else
    rpow = 1
  endif

  if (present(include_guard_cells)) then
    inc_guard_cells = include_guard_cells
  else
    inc_guard_cells = .false.
  endif

  ! offset of each field component from the lower-left cell corner
  ! note, cell (1,1) has lower-left corner with radial position -dr/2
  ! e field should be (0,0.5,0) and b field should be (0.5,0,0.5)
  if (present(off)) then
    offs = off(1:f_dim)
  else
    offs = 0.0_p_double
  endif

  ! get array size and guard cells
  if (inc_guard_cells) then
    vdf_lbound(1:x_dim) = 1 - vdf_a%gc_num_(1,1:x_dim)
    vdf_ubound(1:x_dim) = vdf_a%nx_(1:x_dim) + vdf_a%gc_num_(2,1:x_dim)
  else
    vdf_lbound(1:x_dim) = 1
    vdf_ubound(1:x_dim) = vdf_a%nx_(1:x_dim)
  endif

  ! don't count field values with r <= 0
  shift_i2 = gir_pos - 2
  if ( shift_i2 < 0 ) then
    vdf_lbound(2) = 2
  endif

  do i = 1, f_dim*(2*n_modes+1)
    total(i) = 0.0_p_double
  enddo

  do mode = 0, n_modes

    do re_im = p_real, p_imag

      if ( re_im == p_real ) then
        vdf_a => cyl_a%pf_re(mode)
      elseif ( re_im == p_imag .and. mode > 0 ) then
        vdf_a => cyl_a%pf_im(mode)
      else
        vdf_a => null()
        cycle
      endif

      ! one less than lower bound of total array to be used in the loop
      if ( mode == 0 ) then
        lb = 0
      else
        lb = f_dim * ( 2*mode - 1 + re_im )
      endif

      ! special case of 3 field components, unrolls inner loop
      ! also multiplies by radius
      if ( f_dim == 3 ) then

        do i2 = vdf_lbound(2), vdf_ubound(2)
          r = ( ( i2 + shift_i2 + offs - 0.5_p_double ) * vdf_a%dx_(2) ) ** rpow
          do i1 = vdf_lbound(1), vdf_ubound(1)
            total(lb+1) = total(lb+1) + r(1) * vdf_a%f2(1,i1,i2)**lpow
            total(lb+2) = total(lb+2) + r(2) * vdf_a%f2(2,i1,i2)**lpow
            total(lb+3) = total(lb+3) + r(3) * vdf_a%f2(3,i1,i2)**lpow
          enddo
        enddo

      else

        do i2 = vdf_lbound(2), vdf_ubound(2)
          r = ( ( i2 + shift_i2 + offs - 0.5_p_double ) * vdf_a%dx_(2) ) ** rpow
          do i1 = vdf_lbound(1), vdf_ubound(1)
            do i = 1, f_dim
              total(lb+i) = total(lb+i) + r(i) * vdf_a%f2(i,i1,i2)**lpow
            enddo
          enddo
        enddo

      endif

    enddo

  enddo

  ! divide by 2 for sin(theta) and cos(theta) in integrals
  do i = f_dim+1, f_dim*(2*n_modes+1)
    total(i) = total(i) * 0.5_p_double
  enddo

end subroutine total_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Add cyl_modes vdf data and stores the values in cyl_modes_a
!-----------------------------------------------------------------------------------------
subroutine add_cyl_modes( cyl_modes_a, cyl_modes_b )
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_cyl_modes), intent(inout) :: cyl_modes_a
  type(t_cyl_modes), intent(in) :: cyl_modes_b

  integer :: n_modes, mode

  n_modes = ubound( cyl_modes_a%pf_re, 1 )

  call add( cyl_modes_a%pf_re(0), cyl_modes_b%pf_re(0) )

  do mode = 1, n_modes

    call add( cyl_modes_a%pf_re(mode), cyl_modes_b%pf_re(mode) )
    call add( cyl_modes_a%pf_im(mode), cyl_modes_b%pf_im(mode) )

  enddo

end subroutine add_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! a = a * b
!-----------------------------------------------------------------------------------------
subroutine mult_cyl_modes( cyl_modes_a, cyl_modes_b )
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_cyl_modes), intent(inout) :: cyl_modes_a
  type(t_cyl_modes), intent(in) :: cyl_modes_b

  integer :: n_modes, mode

  n_modes = ubound( cyl_modes_a%pf_re, 1 )

  call mult( cyl_modes_a%pf_re(0), cyl_modes_b%pf_re(0) )

  do mode = 1, n_modes

    call mult( cyl_modes_a%pf_re(mode), cyl_modes_b%pf_re(mode) )
    call mult( cyl_modes_a%pf_im(mode), cyl_modes_b%pf_im(mode) )

  enddo

end subroutine mult_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! a = a ** int
!-----------------------------------------------------------------------------------------
subroutine pow_cyl_modes( cyl_modes_a, int_exp )
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_cyl_modes), intent(inout) :: cyl_modes_a
  integer, intent(in) :: int_exp

  integer :: n_modes, mode

  n_modes = ubound( cyl_modes_a%pf_re, 1 )

  call pow( cyl_modes_a%pf_re(0), int_exp )

  do mode = 1, n_modes

    call pow( cyl_modes_a%pf_re(mode), int_exp )
    call pow( cyl_modes_a%pf_im(mode), int_exp )

  enddo

end subroutine pow_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! cyl_modes_a(1) = cyl_modes_a(1) + cyl_modes_b(fc)^int_exp
!-----------------------------------------------------------------------------------------
subroutine add_pow_cyl_modes( cyl_modes_a, cyl_modes_b, int_exp, fc )
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_cyl_modes), intent(inout) :: cyl_modes_a
  type(t_cyl_modes), intent(in) :: cyl_modes_b
  integer, intent(in) :: int_exp, fc

  integer :: n_modes, mode

  n_modes = ubound( cyl_modes_a%pf_re, 1 )

  call add_pow( cyl_modes_a%pf_re(0), cyl_modes_b%pf_re(0), int_exp, fc )

  do mode = 1, n_modes

    call add_pow( cyl_modes_a%pf_re(mode), cyl_modes_b%pf_re(mode), int_exp, fc )
    call add_pow( cyl_modes_a%pf_im(mode), cyl_modes_b%pf_im(mode), int_exp, fc )

  enddo

end subroutine add_pow_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Calculate the divergence of a vector field. Note that this routine assumes
! that the quantities are scattered in a Yee like mesh (which is true for all
! vector fields in Osiris) meaning that all spatial derivatives will fall in the
! same point.
!
! Here we also assume cylindrical coordinates and calculate the divergence accordingly.
!-----------------------------------------------------------------------------------------
subroutine div_cyl_modes( field, div, gir_pos, phi_at_axis )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2
  integer, parameter :: f_dim = 3

  type(t_cyl_modes), intent(in) :: field
  type(t_cyl_modes), intent(inout) :: div
  integer, intent(in) :: gir_pos
  logical, intent(in) :: phi_at_axis

  integer :: oz, or, i1, i2, n_modes, mode, re_im, shift_i2
  integer, dimension(rank) :: vdf_lbound, vdf_ubound
  real(p_k_fld), dimension(rank) :: rdx
  real(p_k_fld) :: rr_2
  real(p_double) :: dr, off
  real(p_k_fld), dimension(:,:,:), pointer :: fld, fld_opp, dv

  n_modes = ubound( field%pf_re, 1 )

  rdx(1:rank) = real( 1.0_p_double/field%pf_re(0)%dx_(1:rank), p_k_fld )
  dr = field%pf_re(0)%dx_(p_r_dim)

  if ( field%pf_re(0)%f_dim_ /= f_dim ) then
    write(err_buf__,*) 'In div_cyl_modes, only written for f_dim = 3';call err__("cyl_modes/os-cyl-modes.f03",973)
    call abort_program(p_err_notimplemented)
  endif

  ! Phi component is located at the axis for B field, off-axis for E field.
  ! The off parameter is r offset of the z/phi field components of the first cell compared
  ! to r = 0. For E it is -0.5 and for B it is 0.0.
  ! We also define integer offsets to make divergence valid
  if ( phi_at_axis ) then
    ! Like B field
    off = 0.0_p_double
    oz = 1
    or = 1
  else
    ! Like E field
    off = -0.5_p_double
    oz = 0
    or = 0
  endif

  ! get array size and guard cells
  vdf_lbound(1:rank) = 1
  vdf_ubound(1:rank) = field%pf_re(0)%nx_(1:rank)

  ! don't calculate divergence for cell at r <= 0
  ! Maybe employ the divergence theorem to calculate div B at r = 0?
  shift_i2 = gir_pos - 2
  if ( shift_i2 < 0 ) then
    vdf_lbound(2) = 2
  endif

  ! Divergence is calculated as
  ! div F = Fr/r + dEr/dr + dFz/dz + i*m*Fphi/r

  ! mode 0
  fld => field%pf_re(0)%f2
  dv => div%pf_re(0)%f2
  do i2 = vdf_lbound(2), vdf_ubound(2)
    ! 1/(2*r)
    rr_2 = real( 1 / (( i2 + shift_i2 + off ) * dr * 2.0_p_double), p_k_fld )
    do i1 = vdf_lbound(1), vdf_ubound(1)
      dv(1,i1,i2) = (fld(2,i1,i2+or) + fld(2,i1,i2+or-1)) * rr_2 + &
                    (fld(2,i1,i2+or) - fld(2,i1,i2+or-1)) * rdx(2) + &
                    (fld(1,i1+oz,i2) - fld(1,i1+oz-1,i2)) * rdx(1)
    enddo
  enddo

  ! do other modes
  do mode = 1, n_modes
    do re_im = p_real, p_imag
      if ( re_im == p_real ) then
        fld => field%pf_re(mode)%f2
        dv => div%pf_re(mode)%f2
        fld_opp => field%pf_im(mode)%f2
      else
        fld => field%pf_im(mode)%f2
        dv => div%pf_im(mode)%f2
        fld_opp => field%pf_re(mode)%f2
      endif

      do i2 = vdf_lbound(2), vdf_ubound(2)
        ! 1/(2*r)
        rr_2 = real( 1 / (( i2 + shift_i2 + off ) * dr * 2.0_p_double), p_k_fld )
        do i1 = vdf_lbound(1), vdf_ubound(1)
          dv(1,i1,i2) = (fld(2,i1,i2+or) + fld(2,i1,i2+or-1)) * rr_2 + &
                        (fld(2,i1,i2+or) - fld(2,i1,i2+or-1)) * rdx(2) + &
                        (fld(1,i1+oz,i2) - fld(1,i1+oz-1,i2)) * rdx(1) + &
                        (-1) ** (re_im) * mode * 2.0_p_k_fld * rr_2 * fld_opp(3,i1,i2)
        enddo
      enddo

    enddo
  enddo

end subroutine div_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Smooth cyl_modes vdf data
!-----------------------------------------------------------------------------------------
subroutine smooth_cyl_modes( cyl_modes, smoother )
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_cyl_modes), intent(inout) :: cyl_modes
  type( t_smooth ), intent(in) :: smoother

  integer :: n_modes, mode

  n_modes = ubound( cyl_modes%pf_re, 1 )

  call smooth( cyl_modes%pf_re(0), smoother )

  do mode = 1, n_modes

    call smooth( cyl_modes%pf_re(mode), smoother )
    call smooth( cyl_modes%pf_im(mode), smoother )

  enddo

end subroutine smooth_cyl_modes
!-----------------------------------------------------------------------------------------

end module m_cyl_modes
