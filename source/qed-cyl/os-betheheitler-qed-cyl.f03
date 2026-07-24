#include "os-preprocess.fpp"
#include "os-config.h"

module m_betheheitler_qedcyl

#include "memory/memory.h"
  !use m_random_class

  use m_parameters
  use m_species_define_qed
  use m_system
  use m_node_conf
  use m_random
  use m_vdf_define
  use m_vdf
  use m_species_define
  use m_betheheitler
  use m_cyl_modes
  use m_species_cyl_modes_define


implicit none

private

! -----
! CLASS DEFINITIONS
! -----

type, extends( t_betheheitler_info ) :: t_betheheitler_info_qedcyl

  type( t_cyl_modes ) :: rho_ion_cyl_m
  class( t_species_cyl_modes ), pointer :: ion_cyl_m

  contains

end type

public :: t_betheheitler_info_qedcyl

contains

end module m_betheheitler_qedcyl