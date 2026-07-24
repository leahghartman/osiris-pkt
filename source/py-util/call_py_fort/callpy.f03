!-----------------------------------------------------------------------------------------
! This file is a modified version of the src/callpy_mod.f90 file from the call_py_fort
! repo by nbren12.  That software is distributed under the Apache License 2.0, which can
! be found within OSIRIS in the codegen/py-util folder.
!-----------------------------------------------------------------------------------------

! module for calling python from C
module m_callpy
  use, intrinsic :: iso_c_binding
  use m_system
  implicit none

  private

#ifdef ENABLE_PY_UTIL
  interface
    function set_state_py(tag, t, nx, ny, nz) result(y) bind(c)
      use iso_c_binding
      character(c_char) :: tag
      integer(c_int) :: nx, ny, nz
      real(c_double) t(nx, ny, nz)
      integer(c_int) :: y
    end function set_state_py

    function get_state_py(tag, t, n) result(y) bind(c)
      use iso_c_binding
      character(c_char) :: tag
      integer(c_int) :: n, y
      real(c_double) t(n)
    end function get_state_py
  end interface
#endif

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

#ifdef ENABLE_PY_UTIL

  ! This function was added to the original file
  pure function if_py_util()

    implicit none

    logical :: if_py_util

    if_py_util = .true.

  end function if_py_util

  subroutine call_function(module_name, function_name)
    interface
      function call_function_py(mod_name_c, fun_name_c) &
        result(y) bind(c, name='call_function')
        use iso_c_binding
        character(kind=c_char) mod_name_c, fun_name_c
        integer(c_int) :: y
      end function call_function_py
    end interface

    character(len=*) :: module_name, function_name
    character(kind=c_char, len=256) :: mod_name_c, fun_name_c

    mod_name_c = trim(module_name)//char(0)
    fun_name_c = trim(function_name)//char(0)

    call check(call_function_py(mod_name_c, fun_name_c))

  end subroutine call_function

  subroutine set_state_integer_1d(tag, t)
    character(len=*) :: tag
    integer :: t(:)
    ! work arrays
    real(c_double) :: tmp(size(t, 1))
    integer(c_int) :: nx, ny, nz
    character(len=256) :: tag_c


    tag_c = trim(tag)//char(0)
    tmp = dble(t)
    nx = size(tmp, 1)
    ny = -1
    nz = -1

    call check(set_state_py(tag_c, tmp, nx, ny, nz))
  end subroutine set_state_integer_1d

  subroutine set_state_double_1d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:)
    ! work arrays
    real(c_double) :: tmp(size(t, 1))
    integer(c_int) :: nx, ny, nz
    character(len=256) :: tag_c


    tag_c = trim(tag)//char(0)

    tmp = t

    nx = size(tmp, 1)
    ny = -1
    nz = -1

    call check(set_state_py(tag_c, tmp, nx, ny, nz))

  end subroutine set_state_double_1d

  subroutine set_state_float_1d(tag, t)
    character(len=*) :: tag
    real :: t(:)
    ! work arrays
    real(c_double) :: tmp(size(t, 1))
    integer(c_int) :: nx, ny, nz
    character(len=256) :: tag_c


    tag_c = trim(tag)//char(0)

    tmp = t

    nx = size(tmp, 1)
    ny = -1
    nz = -1

    call check(set_state_py(tag_c, tmp, nx, ny, nz))

  end subroutine set_state_float_1d

  subroutine set_state_double_2d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:,:)
    ! work arrays
    real(c_double) :: tmp(size(t, 1), size(t, 2))
    integer(c_int) :: nx, ny, nz
    character(len=256) :: tag_c


    tag_c = trim(tag)//char(0)

    tmp = t

    nx = size(tmp, 1)
    ny = size(tmp, 2)
    nz = -1

    call check(set_state_py(tag_c, tmp, nx, ny, nz))

  end subroutine set_state_double_2d

  subroutine set_state_float_2d(tag, t)
    character(len=*) :: tag
    real :: t(:,:)
    ! work arrays
    real(c_double) :: tmp(size(t, 1), size(t, 2))
    integer(c_int) :: nx, ny, nz
    character(len=256) :: tag_c


    tag_c = trim(tag)//char(0)

    tmp = t

    nx = size(tmp, 1)
    ny = size(tmp, 2)
    nz = -1

    call check(set_state_py(tag_c, tmp, nx, ny, nz))

  end subroutine set_state_float_2d

  subroutine set_state_double_3d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:,:,:)
    ! work arrays
    real(c_double) :: tmp(size(t, 1), size(t, 2), size(t, 3))
    integer(c_int) :: nx, ny, nz
    character(len=256) :: tag_c


    tag_c = trim(tag)//char(0)

    tmp = t

    nx = size(tmp, 1)
    ny = size(tmp, 2)
    nz = size(tmp, 3)

    call check(set_state_py(tag_c, tmp, nx, ny, nz))

  end subroutine set_state_double_3d

  subroutine set_state_float_3d(tag, t)
    character(len=*) :: tag
    real :: t(:,:,:)
    ! work arrays
    real(c_double) :: tmp(size(t, 1), size(t, 2), size(t, 3))
    integer(c_int) :: nx, ny, nz
    character(len=256) :: tag_c


    tag_c = trim(tag)//char(0)

    tmp = t

    nx = size(tmp, 1)
    ny = size(tmp, 2)
    nz = size(tmp, 3)

    call check(set_state_py(tag_c, tmp, nx, ny, nz))

  end subroutine set_state_float_3d

  subroutine set_state_scalar(tag, t)
    character(len=*) :: tag
    real :: t
    real(c_double) :: t_
    character(len=256) :: tag_c
    interface
       function set_state_scalar_py(tag, t) result(y)&
            bind(c, name='set_state_scalar')
         use iso_c_binding
         character(c_char) :: tag
         real(c_double) t
         integer(c_int) :: y
       end function set_state_scalar_py
    end interface

    t_ = t
    tag_c = trim(tag)//char(0)
    call check(set_state_scalar_py(tag_c, t_))
  end subroutine set_state_scalar


  subroutine get_state_double_3d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:, :, :)
    real(c_double) :: t_(size(t, 1), size(t, 2), size(t, 3))
    character(len=256) :: tag_c

    integer(c_int) :: n
    n  = size(t)
    tag_c = trim(tag)//char(0)
    call check(get_state_py(tag_c, t_, n))
    t = dble(t_)
  end subroutine get_state_double_3d

  subroutine get_state_double_2d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:, :)
    real(c_double) :: t_(size(t, 1), size(t, 2))
    character(len=256) :: tag_c

    integer(c_int) :: n
    n  = size(t)
    tag_c = trim(tag)//char(0)
    call check(get_state_py(tag_c, t_, n))
    t = dble(t_)
  end subroutine get_state_double_2d

  subroutine get_state_double_1d(tag, t)
    character(len=*) :: tag
    real(8) :: t(:)
    real(c_double) :: t_(size(t, 1))
    character(len=256) :: tag_c

    integer(c_int) :: n
    n  = size(t)
    tag_c = trim(tag)//char(0)
    call check(get_state_py(tag_c, t_, n))
    t = dble(t_)
  end subroutine get_state_double_1d

  subroutine get_state_float_3d(tag, t)
    character(len=*) :: tag
    real :: t(:, :, :)
    real(c_double) :: t_(size(t, 1), size(t, 2), size(t, 3))
    character(len=256) :: tag_c

    integer(c_int) :: n
    n  = size(t)
    tag_c = trim(tag)//char(0)
    call check(get_state_py(tag_c, t_, n))
    t = real(t_)
  end subroutine get_state_float_3d

  subroutine get_state_float_2d(tag, t)
    character(len=*) :: tag
    real :: t(:, :)
    real(c_double) :: t_(size(t, 1), size(t, 2))
    character(len=256) :: tag_c

    integer(c_int) :: n
    n  = size(t)
    tag_c = trim(tag)//char(0)
    call check(get_state_py(tag_c, t_, n))
    t = real(t_)
  end subroutine get_state_float_2d

  subroutine get_state_float_1d(tag, t)
    character(len=*) :: tag
    real :: t(:)
    real(c_double) :: t_(size(t, 1))
    character(len=256) :: tag_c

    integer(c_int) :: n
    n  = size(t)
    tag_c = trim(tag)//char(0)
    call check(get_state_py(tag_c, t_, n))
    t = real(t_)
  end subroutine get_state_float_1d

  subroutine get_state_integer_1d(tag, t)
    character(len=*) :: tag
    integer :: t(:)
    real(c_double) :: t_(size(t, 1))
    character(len=256) :: tag_c

    integer(c_int) :: n
    n  = size(t)
    tag_c = trim(tag)//char(0)
    call check(get_state_py(tag_c, t_, n))
    t = nint(t_)
  end subroutine get_state_integer_1d


  subroutine set_state_char(tag, chr)
    interface
      function set_state_char_py(tag, chr) result(y)&
        bind(c, name='set_state_char')
        use iso_c_binding
        implicit none
        character(c_char) :: tag
        character(c_char) :: chr
        integer(c_int) :: y
      end function set_state_char_py
    end interface
    character(len=*) :: tag, chr
    character(len=256) :: tag_, chr_

    tag_ = trim(tag)//char(0)
    chr_ = trim(chr)//char(0)
    call check(set_state_char_py(tag_, chr_))
  end subroutine set_state_char

  subroutine get_state_char(tag, chr)
    interface
      function get_state_char_py(tag, chr, n) result(y)&
        bind(c, name='get_state_char')
        use iso_c_binding
        implicit none
        character(c_char) :: tag
        character(c_char) :: chr
        integer(c_int) :: y, n
      end function
    end interface
    character(len=*) :: tag, chr
    chr = ""
    call check(get_state_char_py(trim(tag) // char(0), chr, len(chr)))
  end subroutine get_state_char

  subroutine check(ret)
    integer :: ret
    if (ret /= 0) stop -1
  end subroutine check

#else

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

#endif

end module m_callpy
