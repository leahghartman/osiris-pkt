# 1 "gr/os-psource-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-psource-gr.f03" 2
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
# 2 "gr/os-psource-gr.f03" 2
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
# 3 "gr/os-psource-gr.f03" 2

module m_psource_gr

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
# 7 "gr/os-psource-gr.f03" 2

use m_species_define, only :t_psource
use m_psource_std
use m_math, only : pi

implicit none

private

interface new_source
    module procedure new_source
end interface new_source

public :: new_source

type, extends( t_psource_std ) :: t_psource_std_gr

contains

procedure :: inject => inject_std_gr

end type t_psource_std_gr

public t_psource_std_gr

contains

!-------------------------------------------------------------------------------
function new_source( type, input_file, coordinates )

  use m_input_file

  implicit none

  class( t_psource ), pointer :: new_source
  character(len=*), intent(in) :: type
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates
  ! create t_profile object of the selected kind
  select case ( trim(type) )
  case('standard','profile')
    allocate( t_psource_std_gr :: new_source )
  case default
    new_source => null()
  end select

  ! If successfull read input parameters and finish initalization
  if ( associated(new_source) ) then
    call new_source % read_input( input_file, coordinates )
  endif

end function new_source
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function inject_std_gr( this, species, ig_xbnd_inj, jay, no_co, &
  bnd_cross, node_cross, send_msg, recv_msg ) result(num_inj)

  use m_system
  use m_parameters
  use m_species_define, only : t_species, t_part_idx, p_spec_buf_block
  use m_current_define, only : t_current
  use m_node_conf, only : t_node_conf
  use m_species_comm, only : t_spec_msg
  use m_species_udist

  use m_geometry_gr!, only : t_geometry_gr
  use m_species_define_gr, only : t_species_gr

  implicit none

  class(t_psource_std_gr), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  integer, parameter :: p_max_x_dim = 3

  ! cell size (p_x_dim)
  real(p_double), dimension(p_max_x_dim) :: dx
  real(p_double), dimension(p_max_x_dim) :: g_xmin

  ! half distance between particles
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2

  integer :: i, i1, i2, ipart

  real(p_k_part) :: x1

  ! volume that each particle occupies
  real(p_k_part) :: pvol

  ! number of particles per cell
  integer :: ppcell

  ! particle positions (global / inside cell )
  real(p_k_part), dimension(:,:), pointer :: ppos, ppos_cell

  ! particle charges
  real(p_k_part), dimension(:), pointer :: pcharge

  integer :: my_i1
  class( t_geometry_gr ), pointer :: geometry
  real(p_k_part) :: ph

  ! radial integral quantities
  real(p_k_part) :: rplus,rminus,aplus,aminus,intr

  select type( s => species )
  class is ( t_species_gr )
    geometry => s%geometry
  class default
    write(err_buf__,*) 'pointer to geometry must be set with a t_species_geometry object';call err__("gr/os-psource-gr.f03",125)
    call abort_program( p_err_invalid )
  end select

  ! Set total number of injected particles to 0
  num_inj = 0

  ! check if num_par_x is > 0 for all directions
  ! if not return silently
  do i=1, p_x_dim
    if (species%num_par_x(i) <= 0) return
  enddo

  do i = 1, p_x_dim
    ! get cell size
    dx(i) = species%dx(i)

    ! get global mininum (this is shifted by +0.5 cells from global simulation values)
    g_xmin( i ) = species%g_box( p_lower , i )

    ! get half distance between particles
    dxp_2(i) = 0.5_p_k_part/species%num_par_x(i)
  enddo

  ! find total number particles per cell
  ppcell = species%num_par_x(1)
  do i = 2, p_x_dim
    ppcell = ppcell * species%num_par_x(i)
  enddo

  ! initialize temp buffers
  call alloc(ppos, (/p_x_dim, ppcell/),"gr/os-psource-gr.f03",156)
  call alloc(ppos_cell, (/p_x_dim, ppcell/),"gr/os-psource-gr.f03",157)
  call alloc(pcharge, (/ppcell/),"gr/os-psource-gr.f03",158)

  ! find normalization factor
  pvol = sign( 1.0_p_k_part/ppcell, species%rqm )
  select case ( geometry%metric )
    case ( p_geometry_minkowski )
      pvol = real(pvol * ( (1.+1./geometry%A)**3 - (1.-1./geometry%A)**3 ) * 2.0943951023931955_p_double, p_k_part)!2pi/3
    case default !curved case
      pvol = real(pvol * 2*pi, p_k_part)
  end select

  ! Get position of particles inside the cell
  ! The position will always be in the range [-0.5, +0.5 [ regardless of interpolation type

  i = 0
  select case ( p_x_dim )
  case (2)
    do i1 = 1, species%num_par_x(1)
      x1 = (2*i1 - 1 - species%num_par_x(1)) * dxp_2(1)
      do i2 = 1, species%num_par_x(2)
        i = i + 1
        ppos_cell(1,i) = x1
        ppos_cell(2,i) = (2*i2 - 1 - species%num_par_x(2)) * dxp_2(2)
      enddo
    enddo
  end select

  ! Inject particles
  ipart = species%num_par + 1

  select case ( p_x_dim )
  case (2)
  ! loop through all cells and count number of particles to inject
  do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)

    my_i1 = i1 - species%my_nx_p( p_lower, 1 ) + 1
    ppos(1,:) = real( geometry%r%f1(2, my_i1) + ppos_cell(1,:) * geometry%dr%f1(1, my_i1), p_k_part )


    do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
      ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

      call this%get_den_value( ppos, ppcell, pcharge )

      do i = 1, ppcell
        if (pcharge(i) > this % den_min ) then
          num_inj = num_inj + 1
        endif
      enddo

    enddo
  enddo

  if ( num_inj > 0 ) then

    ! check if the buffer size is sufficient and grow it if necessary
    if ( num_inj > species%num_par_max - species%num_par ) then
      call species%grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
    endif


    ! loop through all the injection cells and
    ! inject particles, normalizing charge
    do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)

      my_i1 = i1 - species%my_nx_p( p_lower, 1 ) + 1
      ppos(1,:) = real( geometry%r%f1(2, my_i1) + ppos_cell(1,:) * geometry%dr%f1(1, my_i1), p_k_part )

      do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
        ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

        call this%get_den_value( ppos, ppcell, pcharge )

        do i=1, ppcell
          if (pcharge(i) > this % den_min) then
            ! add particle
            species%x(1, ipart) = ppos_cell(1,i)
            species%x(2, ipart) = ppos_cell(2,i)
            species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
            species%ix(2, ipart) = i2 - species%my_nx_p( p_lower, 2 ) + 1

            ! add cartesian data
            call random_number(ph)
            ph = ph * 2*pi
            ! ph = 0

            ! the radial contribution is incorporated in the next select
            species%q(ipart) = pcharge(i) * pvol * &
              ( cos(ppos(2,i) - 0.5_p_k_part*dx(2)) -&
                cos(ppos(2,i) + 0.5_p_k_part*dx(2)) )

            select case (geometry%metric)
            case (p_geometry_minkowski)
              species%q(ipart) = species%q(ipart) * ppos(1,i)**3
              species%x(3, ipart) = ppos(1,i) * sin(ppos(2,i)) * cos(ph)
              species%x(4, ipart) = ppos(1,i) * sin(ppos(2,i)) * sin(ph)
              species%x(5, ipart) = ppos(1,i) * cos(ppos(2,i))
            case default !curved spacetime
              rplus = 2._p_double * ppos(1,i) / (1._p_double + 1._p_double / geometry%delta) !particle r+
              rminus = 2._p_double * ppos(1,i) / (1._p_double + geometry%delta) !particle r-
              aplus = sqrt( 1._p_double - geometry%rs / rplus ) !alpha(r+)
              aminus = sqrt( 1._p_double - geometry%rs / rminus) !alpha(r-)
              intr = 0.0416666666666667_p_double *& !1/24
                  (rplus * aplus * (8._p_double * rplus**2 +&
                     10._p_double * rplus * geometry%rs + 15._p_double * geometry%rs**2) -&
                   rminus * aminus * (8._p_double * rminus**2 +&
                     10._p_double * rminus * geometry%rs + 15._p_double * geometry%rs**2)) +&
                  0.3125_p_double * geometry%rs**3 * log(((1._p_double+aplus)*(1._p_double-aminus))/&
                                                         ((1._p_double-aplus)*(1._p_double+aminus)))

              species%q(ipart) = species%q(ipart) * intr
              species%x(3, ipart) = ph
            end select
            ipart = ipart + 1
          endif
        enddo
      enddo
    enddo
  endif
  end select

  ! free temporary memory
  call freemem(ppos,"gr/os-psource-gr.f03",280)
  call freemem(ppos_cell,"gr/os-psource-gr.f03",281)
  call freemem(pcharge,"gr/os-psource-gr.f03",282)

  ! Set momentum of injected particles
  if ( num_inj > 0 ) then
    ! Initialize particle momentum
    call set_momentum( species, species%num_par+1, species%num_par + num_inj )
  endif

end function inject_std_gr
!-------------------------------------------------------------------------------

end module m_psource_gr
