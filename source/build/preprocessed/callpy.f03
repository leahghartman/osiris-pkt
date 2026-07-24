# 1 "py-util/call_py_fort/callpy.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "py-util/call_py_fort/callpy.f03" 2
!-----------------------------------------------------------------------------------------
! This file is a modified version of the src/callpy_mod.f90 file from the call_py_fort
! repo by nbren12. That software is distributed under the Apache License 2.0, which can
! be found within OSIRIS in the codegen/py-util folder.
!-----------------------------------------------------------------------------------------

! module for calling python from C
module m_callpy
  use, intrinsic :: iso_c_binding
  use m_system
  implicit none

  private
# 34 "py-util/call_py_fort/callpy.f03"
  interface set_state
    module procedure set_state_double_3d
    module procedure set_state_double_2d
    module procedure set_state_double_1d
    module procedure set_state_float_3d
    module procedure set_state_float_2d
    module procedure set_state_float_1d
    module procedure set_state_integer_1d
    module procedure set_state_char
  end interface

  interface get_state
    module procedure get_state_float_3d
    module procedure get_state_float_2d
    module procedure get_state_float_1d
    module procedure get_state_double_3d
    module procedure get_state_double_2d
    module procedure get_state_double_1d
    module procedure get_state_integer_1d
    module procedure get_state_char
  end interface

  public :: if_py_util, get_state, set_state, call_function, set_state_scalar


contains
# 393 "py-util/call_py_fort/callpy.f03"
  ! The below subroutines were added to support compilation without the Python interface

  pure function if_py_util()

    implicit none

    logical :: if_py_util

    if_py_util = .false.

  end function if_py_util

  subroutine print_error( fname )

    implicit none

    character(len=*), intent(in) :: fname

    if ( mpi_node() == 0 ) then
      write(0,*) "   Error: The function ", trim(fname), " was called, but code was"
      write(0,*) "   not compiled with python support. Please see the"
      write(0,*) "   PY_FCOMPILEFLAGS and PY_FLINKFLAGS flags in your config file."
      write(0,*) "   aborting..."
    endif
    stop

  end subroutine print_error

  subroutine call_function(module_name, function_name)
    character(len=*) :: module_name, function_name
    call print_error( "call_function" )
  end subroutine call_function

  subroutine set_state_integer_1d(tag, t)
    character(len=*) :: tag
    integer :: t(:)
    call print_error( "set_state_integer_1d" )
  end subroutine set_state_integer_1d

  subroutine set_state_double_1d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:)
    call print_error( "set_state_double_1d" )
  end subroutine set_state_double_1d

  subroutine set_state_float_1d(tag, t)
    character(len=*) :: tag
    real :: t(:)
    call print_error( "set_state_float_1d" )
  end subroutine set_state_float_1d

  subroutine set_state_double_2d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:,:)
    call print_error( "set_state_double_2d" )
  end subroutine set_state_double_2d

  subroutine set_state_float_2d(tag, t)
    character(len=*) :: tag
    real :: t(:,:)
    call print_error( "set_state_float_2d" )
  end subroutine set_state_float_2d

  subroutine set_state_double_3d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:,:,:)
    call print_error( "set_state_double_3d" )
  end subroutine set_state_double_3d

  subroutine set_state_float_3d(tag, t)
    character(len=*) :: tag
    real :: t(:,:,:)
    call print_error( "set_state_float_3d" )
  end subroutine set_state_float_3d

  subroutine set_state_scalar(tag, t)
    character(len=*) :: tag
    real :: t
    call print_error( "set_state_scalar" )
  end subroutine set_state_scalar

  subroutine get_state_double_3d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:, :, :)
    call print_error( "get_state_double_3d" )
  end subroutine get_state_double_3d

  subroutine get_state_double_2d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:, :)
    call print_error( "get_state_double_2d" )
  end subroutine get_state_double_2d

  subroutine get_state_double_1d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:)
    call print_error( "get_state_double_1d" )
  end subroutine get_state_double_1d

  subroutine get_state_float_3d(tag, t)
    character(len=*) :: tag
    real :: t(:, :, :)
    call print_error( "get_state_float_3d" )
  end subroutine get_state_float_3d

  subroutine get_state_float_2d(tag, t)
    character(len=*) :: tag
    real :: t(:, :)
    call print_error( "get_state_float_2d" )
  end subroutine get_state_float_2d

  subroutine get_state_float_1d(tag, t)
    character(len=*) :: tag
    real :: t(:)
    call print_error( "get_state_float_1d" )
  end subroutine get_state_float_1d

  subroutine get_state_integer_1d(tag, t)
    character(len=*) :: tag
    integer :: t(:)
    call print_error( "get_state_integer_1d" )
  end subroutine get_state_integer_1d

  subroutine set_state_char(tag, chr)
    character(len=*) :: tag, chr
    call print_error( "set_state_char" )
  end subroutine set_state_char

  subroutine get_state_char(tag, chr)
    character(len=*) :: tag, chr
    call print_error( "get_state_char" )
  end subroutine get_state_char



end module m_callpy
