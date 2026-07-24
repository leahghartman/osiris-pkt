#include "os-config.h"
#include "os-preprocess.fpp"

module m_photons_define_gr

#include "memory/memory.h"

use m_system
use m_parameters

use stringutil

use m_species_define_gr

implicit none

private

type, extends( t_species_gr ) :: t_photons_gr

contains

end type t_photons_gr

public t_photons_gr

end module m_photons_define_gr