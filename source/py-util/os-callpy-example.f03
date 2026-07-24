!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     module to add interface to python functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_py_util

use m_system
use m_parameters
use m_callpy

implicit none

private

interface call_py
  module procedure call_py_add_one
end interface

public :: call_py

contains

!---------------------------------------------------------------------------------------------------
! Call your custom python file
!---------------------------------------------------------------------------------------------------
subroutine call_py_add_one( x )

  implicit none

  real(p_k_fld), dimension(:), intent(inout) :: x

  ! Set the STATE dictionary with key "VAR_NAME" and value x (makes a copy of the array)
  call set_state( "VAR_NAME", x )

  ! Call your python function (add_one) defined in a python file named "cool_file"
  call call_function( "cool_file", "add_one" )

  ! You can also call built-in functions of STATE
  call call_function( "builtins", "print" )

  ! Get the value with key "VAR_NAME" from the STATE dictionary and copy it into x
  call get_state( "VAR_NAME", x )

end subroutine call_py_add_one

end module m_py_util
