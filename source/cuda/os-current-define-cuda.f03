#include "os-config.h"
#include "os-preprocess.fpp"

module m_current_cuda

use m_current_define
use m_node_conf, only : t_node_conf

implicit none

private

type, extends( t_current ) :: t_current_cuda

contains

  procedure :: num_threads => num_threads_current_cuda

end type t_current_cuda

public :: t_current_cuda

contains

!---------------------------------------------------------------------------------------------------
! Get the number of threads to initialize current arrays (just 1 when running on GPU)
!---------------------------------------------------------------------------------------------------
function num_threads_current_cuda( this, no_co )

  implicit none

  class( t_current_cuda ), intent(in) :: this
  class( t_node_conf ), intent(in) :: no_co

  integer :: num_threads_current_cuda

  num_threads_current_cuda = 1

end function num_threads_current_cuda

end module m_current_cuda
