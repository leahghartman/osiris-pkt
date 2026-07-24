#include "os-preprocess.fpp"
#include "os-config.h"

module m_bremsstrahlung

#include "memory/memory.h"
  !use m_random_class

  use m_parameters
  use m_system
  use m_node_conf
  use m_random
  use m_vdf_define
  use m_vdf
  use m_species_define


implicit none

private

! -----
! CLASS DEFINITIONS
! -----

type :: t_bremsstrahlung_info

  logical :: if_bremsstrahlung
  integer :: Z_ion
  type( t_vdf ) :: rho_ion

  ! multiplicative factor applied on the probability
  real( p_double ) :: proba_mult

  ! Type of cross section
  integer :: cs_model
  logical :: energy_damp
  real( p_double ) :: p_emit_cutoff

  ! pointer to the ion species for bremsstrahlung
  class( t_species ), pointer :: ion

  integer, dimension(:), allocatable :: i_ion

  contains

  procedure :: read_input => read_input_bremsstrahlung

end type

public :: t_bremsstrahlung_info

contains

  !-----------------------------------------------------------------------------------------
  ! Reads info for bremsstrahlung
  !-----------------------------------------------------------------------------------------

  subroutine read_input_bremsstrahlung( this, input_file, num_qed)

  use m_input_file, only : t_input_file, get_namelist

  implicit none

  class( t_bremsstrahlung_info ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: num_qed
  
  logical :: if_bremsstrahlung, energy_damp
  integer :: ierr, Z_ion
  integer :: cs_model
  real(p_double) :: proba_mult
  integer, dimension(:), allocatable :: i_ion

  namelist /nl_qed_bremsstrahlung/ if_bremsstrahlung, Z_ion, proba_mult, cs_model, energy_damp, i_ion

  if_bremsstrahlung = .false.
  Z_ion = -1
  proba_mult = 1.0
  cs_model = -1
  energy_damp = .true.
  allocate(i_ion(num_qed))
  i_ion = 0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_qed_bremsstrahlung", ierr )
  
  if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_qed_bremsstrahlung, iostat = ierr)
      if (ierr /= 0) then
          print *, ""
          print *, "   Error reading qed_bremsstrahlung parameters "
          print *, "   aborting..."
          stop
      endif
  else
      SCR_ROOT("   - no particle bremsstrahlung specified")
  endif

  this%if_bremsstrahlung = if_bremsstrahlung
  this%Z_ion = Z_ion
  this%proba_mult = proba_mult
  this%cs_model = cs_model
  this%energy_damp = energy_damp
  this%i_ion = i_ion

  end subroutine read_input_bremsstrahlung

end module m_bremsstrahlung