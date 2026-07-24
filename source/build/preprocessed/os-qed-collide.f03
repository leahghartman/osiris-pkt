# 1 "qed/os-qed-collide.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "qed/os-qed-collide.f03" 2
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
# 2 "qed/os-qed-collide.f03" 2
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
# 3 "qed/os-qed-collide.f03" 2

module m_qed_collide

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
# 7 "qed/os-qed-collide.f03" 2
 !use m_random_class

  use m_parameters
  use m_species_define_qed
  use m_photons_define
  use m_system
  use m_node_conf
  use m_random
  use m_vdf_define
  use m_vdf



implicit none

private

! -----
! CLASS DEFINITIONS
! -----

! collision type
integer, parameter :: eph = 1 ! Compton scattering electron-photon
integer, parameter :: pph = 2 ! Compton scattering positron-photon

type :: t_qed_collide_info

  logical :: if_qed_collide
  integer, dimension(p_x_dim) :: npc_per_ccell ! number of pic cells per collision cells
  logical :: if_recoil
  logical :: if_split
  real(p_double) :: sigma_coeff ! 1.0 Thomson/Compton cross section, this is only a multiplicativce coefficient

  contains

  procedure :: read_input => read_input_qed_collide_info

end type

type :: t_qed_clist

! index list of space-binned particles
integer, dimension(:), pointer :: lst_e => null() ! electrons
integer, dimension(:), pointer :: lst_p => null() ! positrons
integer, dimension(:), pointer :: lst_ph => null() ! photons

integer :: ic_mates, d_buff = p_vecwidth, Ncoll_cells ! number of macro-collisions
integer, dimension(p_x_dim) :: num_coll_cells, npc_per_ccell ! number of collision cells, number of pic cells per collision cell

! indexes in collision cell
integer, dimension(:), pointer :: index_e => null()
integer, dimension(:), pointer :: index_p => null()
integer, dimension(:), pointer :: index_ph => null()

! indexes for shuffle
integer, dimension(:), pointer :: Shuf_ind_ep => null()
integer, dimension(:), pointer :: Shuf_ind_ph => null()

! Mask
integer, dimension(:), pointer :: Msk_ep => null()

integer, dimension(:,:), pointer :: idx_c_mates => null() ! index list of colliding couples by index (N,3) where the third column indicates the collision type

contains
  procedure :: allocate_objs => allocate_obj_qed_clist
  procedure :: init => init_qed_clist
  procedure :: cleanup => cleanup_qed_clist
  procedure :: grow_buffer => grow_buffer_qed_clist

end type t_qed_clist

public :: t_qed_collide_info, t_qed_clist


interface get_coll_listNTC
  module procedure get_coll_listNTC
end interface

interface Compton_scatter
  module procedure Compton_scatter
end interface

interface angle_Thomson_gen
  module procedure angle_Thomson_gen
end interface

!interface cdf_Compton
! module procedure cdf_Compton
!end interface

public :: eph, pph, Compton_scatter, angle_Compton_gen, cdf_Compton, angle_Thomson_gen, get_coll_listNTC

contains

  subroutine allocate_obj_qed_clist( this, if_sort, if_list )
    implicit none

    class( t_qed_clist ), intent(inout) :: this
    logical, intent(in) :: if_sort, if_list

    integer, dimension(p_x_dim) :: num_coll_cells
    integer :: Ncoll_cells, i

    Ncoll_cells = this%Ncoll_cells
    if (if_sort) then
      call alloc(this%lst_e, (/ Ncoll_cells /),"qed/os-qed-collide.f03",112)
      call alloc(this%lst_p, (/ Ncoll_cells /),"qed/os-qed-collide.f03",113)
      call alloc(this%lst_ph, (/ Ncoll_cells /),"qed/os-qed-collide.f03",114)
    endif
    if (if_list) then
      call alloc(this%Shuf_ind_ep, (/ this%d_buff /),"qed/os-qed-collide.f03",117)
      call alloc(this%Shuf_ind_ph, (/ this%d_buff /),"qed/os-qed-collide.f03",118)
      call alloc(this%idx_c_mates, (/Ncoll_cells*this%d_buff, 3 /),"qed/os-qed-collide.f03",119)
      call alloc(this%index_e, (/ this%d_buff /),"qed/os-qed-collide.f03",120)
      call alloc(this%index_p, (/ this%d_buff /),"qed/os-qed-collide.f03",121)
      call alloc(this%index_ph, (/ this%d_buff /),"qed/os-qed-collide.f03",122)
      call alloc(this%Msk_ep, (/ this%d_buff /),"qed/os-qed-collide.f03",123)
    endif

    this%ic_mates = 0

  end subroutine allocate_obj_qed_clist

  subroutine init_qed_clist( this, if_sort, if_list, if_mates )
    implicit none

    class( t_qed_clist ), intent(inout) :: this
    logical, intent(in) :: if_sort, if_list, if_mates

    integer :: i
    if (if_sort) then
      do i = 1, this%Ncoll_cells
        this%lst_e(i) = 0
        this%lst_p(i) = 0
        this%lst_ph(i) = 0
      enddo
    endif
    if (if_list) then
      do i = 1, this%d_buff
        this%idx_c_mates(i,:) = (/ 0, 0, 0 /)
        this%index_e(i) = 0
        this%index_p(i) = 0
        this%index_ph(i) = 0
        this%Shuf_ind_ep = i
        this%Shuf_ind_ph = i
        this%Msk_ep(i) = 0
      enddo
    endif

    if (if_mates) this%ic_mates = 0

  end subroutine init_qed_clist

  subroutine cleanup_qed_clist( this, if_sort, if_list, if_mates )
    implicit none

    class( t_qed_clist ), intent(inout) :: this
    logical, intent(in) :: if_sort, if_list, if_mates

    if (if_sort) then
      call freemem(this%lst_e,"qed/os-qed-collide.f03",167)
      call freemem(this%lst_p,"qed/os-qed-collide.f03",168)
      call freemem(this%lst_ph,"qed/os-qed-collide.f03",169)
    endif
    if (if_list) then
      call freemem(this%index_e,"qed/os-qed-collide.f03",172)
      call freemem(this%index_p,"qed/os-qed-collide.f03",173)
      call freemem(this%index_ph,"qed/os-qed-collide.f03",174)
      call freemem(this%Shuf_ind_ep,"qed/os-qed-collide.f03",175)
      call freemem(this%Shuf_ind_ph,"qed/os-qed-collide.f03",176)
      call freemem(this%Msk_ep,"qed/os-qed-collide.f03",177)
    endif

    if (if_mates) call freemem(this%idx_c_mates,"qed/os-qed-collide.f03",180)

  end subroutine cleanup_qed_clist

  subroutine grow_buffer_qed_clist( this, num_par_req )
    implicit none

    class( t_qed_clist ), intent(inout) :: this
    integer, intent(in) :: num_par_req

    integer :: num_par_new

    ! The buffer size must always be a multiple of vector width in size because of SIMD code
    num_par_new = (( num_par_req + p_vecwidth - 1 ) / p_vecwidth) * p_vecwidth
    print*,'resizing qed colliison buffer',this%d_buff,' --> ',num_par_new

    ! cleanup list
    call this%cleanup( .false. , .true. , .true. ) ! dont clean sorting

    this%d_buff = num_par_new

    ! allocate list
    call this%allocate_objs( .false., .true.) ! dont allocate sorting

    ! init list
    call this%init( .false., .true., .true. ) ! dont init sorting

  end subroutine grow_buffer_qed_clist

  !-----------------------------------------------------------------------------------------
  ! Reads info for qed collisions
  !-----------------------------------------------------------------------------------------

  subroutine read_input_qed_collide_info( this, input_file, num_qed )

  use m_input_file, only : t_input_file, get_namelist

  implicit none

  class( t_qed_collide_info ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: num_qed

  logical :: if_qed_collide, if_recoil, if_split
  integer, dimension(p_x_dim) :: npc_per_ccell
  real(p_double) :: sigma_coeff

  integer :: i, ierr

  namelist /nl_qed_collide_info/ if_qed_collide, npc_per_ccell, if_recoil, if_split, sigma_coeff

  if_qed_collide = .false.
  if_recoil = .false.
  if_split = .true.
  sigma_coeff = 1.0 ! default to Thomson cross section
  do i=1,p_x_dim
      npc_per_ccell(i) = 1
  enddo

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_qed_collide_info", ierr )
  if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_qed_collide_info, iostat = ierr)
      if (ierr /= 0) then
          print *, ""
          print *, "   Error reading qed_collision_info parameters "
          print *, "   aborting..."
          stop
      endif
  else
      if (mpi_node()==0) print *,"   - no particle qed_collide specified"
  endif

  this%if_qed_collide = if_qed_collide
  this%if_recoil = if_recoil
  this%if_split = if_split
  this%sigma_coeff = sigma_coeff
  do i=1,p_x_dim
      this%npc_per_ccell(i) = npc_per_ccell(i)
  enddo

  if ( (if_qed_collide) .and. (num_qed /= 1) ) then
    if (mpi_node()==0) print *,"(*error*) num_qed \= 1 not implemented for qed collisions"
    if (mpi_node()==0) print *,"(*error*) aborting..."
    stop
  endif

  end subroutine read_input_qed_collide_info




  ! --------------------------------------------------------------------------------------------------------------------
  ! ------------------------------------------------- COLLISION LIST ------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------



subroutine Knuth_shuffle( VV, il, iu )

  implicit none
  integer, intent(inout) :: VV(:)
  integer, intent(in) :: il, iu ! lower and upper indexex, range of array to shuffle
  integer :: N, Idice, i, hat

  N = size(VV)
  ! sanity checke
  if ( (il>=1) .and. (iu<=N) .and. (iu>=il) ) then ! il, iu are ordered indexes of the array
    do i = iu, il, -1
      Idice = ceiling( rng % genrand_real3() * i ) ! generates a rnd from 1:i
      hat = VV( Idice )
      VV( Idice ) = VV( i )
      VV( i ) = hat
    enddo
  else
    print*,'error on knuth shuffle'
  endif

end subroutine Knuth_shuffle


function Pnu_C( w1, gama, w0, Coeff, q )

  implicit none
  real(p_double) :: Pnu_C ! collision probability
  real(p_double) :: w1, gama, w0, q, Coeff, w1fw1

  if (w1 > 3E-4) then ! Klein Nishina
    w1fw1 = 3.0_p_double/8.0_p_double * ( (1.0_p_double-2.0_p_double/w1-2.0_p_double/w1**2)*log( 1.0_p_double+2.0_p_double*w1 ) + &
                              0.5_p_double * ( 1.0_p_double-1.0_p_double/(1.0_p_double+2.0_p_double*w1)**2 ) + 4.0_p_double/w1 )
  else ! Thomson
    w1fw1 = w1
  endif

  ! compute the probability of scattering in the lab frame
    Pnu_C = w1fw1 /gama/w0 * q * Coeff

end function Pnu_C

subroutine get_coll_listNTC( this,electrons,positrons,photons,dt,sigma_coeff ) ! based on the No-Time-Counter (NTC) collision model

  implicit none

  class( t_qed_clist ), intent(inout) :: this
  class(t_species_qed), intent(in) :: electrons, positrons
  class(t_photons), intent(in) :: photons
  real(p_double), intent(in) :: dt, sigma_coeff

  ! dummy variables
  real(p_double) :: w0, w1, ti_list, te_list, P_max, Coeff, gama, q, weight_max, PP, weight_e,weight_p,weight_ph,dumNe,dumNp
  real(p_double), dimension(3) :: p,k

  integer :: i, i_cell, NC_max, Ne, Np, Nph, num_c, ind, ind_ep, ind_ph, Msk, i_buff
  integer :: ie_l, ie_u, ip_l, ip_u, iph_l, iph_u ! lower and upper indexes in the collision cell

  real, parameter :: r_e = 2.8179403227E-13, s_c = 2.99792458E10! classical electron radius and speed of light in cgs

  !call cpu_time(ti_list)
  num_c = this%Ncoll_cells
  Coeff = 2.0_p_double/3.0_p_double * r_e * electrons%omega_p0/s_c * dt *sigma_coeff

  ! --- initialize the list size for the first collision cell

  Ne = this%lst_e(1)
  Np = this%lst_p(1)
  Nph = this%lst_ph(1)

  dumNp = 0.0_p_double
  dumNe = 0.0_p_double

  do i_cell = 1, num_c

    ! initialise the maximum probability
    weight_e=0.0_p_double
    weight_p=0.0_p_double
    weight_ph=0.0_p_double
    P_max = 0.0_p_double
    if ( ( (Ne > 0) .or. (Np > 0) ) .and. (Nph > 0) ) then ! check if there are scattering macro-particles

      ! get macro-particle indexes
      if (Ne>0) then
        do i=1, Ne
          this%index_e(i) = this%lst_e(i_cell) - Ne+i
        enddo
        weight_e = maxval(abs( electrons%q( this%index_e(1):this%index_e(Ne) ) ))
      endif
      if (Np>0) then
        do i=1, Np
          this%index_p(i) = this%lst_p(i_cell) - Np+i
        enddo
        weight_p = maxval( positrons%q( this%index_p(1):this%index_p(Np) ) )
      endif
      !print*, this%index_p

      do i=1, Nph
        this%index_ph(i) = this%lst_ph(i_cell) - Nph+i
      enddo
      weight_ph = maxval( photons%q( this%index_ph(1):this%index_ph(Nph) ) )


      ! get maximum weight
      weight_max = maxval( (/weight_e,weight_p,weight_ph /) )

      ! compute the maximum probability of interaction with Thomson cross section
      ! the factor 2 acounts for relativistic dynamics in the collision :
      ! sigma_L = sigma_0 w_0/gamma_L/w_L where max(w_0) = 2*gamma_L*w_L
      P_max = 2*Pnu_C( 1D-6, 1.0_p_double, 1D-6, Coeff, weight_max )

      ! Maximum number of scattering macro-particles
      NC_max = floor(P_max* (Ne+Np) * Nph)
      if (rng % genrand_real1()<= (P_max* (Ne+Np) * Nph-NC_max)) NC_max=NC_max+1 ! random rounding preserving the probability
      if ( NC_max>minval( (/Ne+Np,Nph /) ) ) then
        print*,'   ***   - warning -   ***   reduce timestep: collision frequency too high'
        NC_max = minval( (/Ne+Np,Nph /) )
      endif

        !write (48,*) NC_max


      ! Generate the collision list
      if (NC_max>0) then

        ! Initialise mask before shuffle Msk(electrons, positrons)
        do i=1, (Ne+Np)
          if (i<=Ne) then
            this%Msk_ep(i) = eph
          else
            this%Msk_ep(i) = pph
          endif
        enddo

        ! --- Reset shuffle Indexes
        do i_buff=1,this%d_buff
          this%Shuf_ind_ep(i_buff) = i_buff
          this%Shuf_ind_ph(i_buff) = i_buff
        enddo

        ! Shuffle the indexes for electrons and positrons
        call Knuth_shuffle(this%Shuf_ind_ep,1,Ne+Np)
        call Knuth_shuffle(this%Shuf_ind_ph,1,Nph)
      endif

      ! Generate collision list
      do i = 1, NC_max
        ind_ph = this%index_ph( this%Shuf_ind_ph(i) )
        Msk = this%Msk_ep( this%Shuf_ind_ep(i) )

        weight_max = 0.0_p_double
        k = photons%p(:, ind_ph )
        weight_max = photons%q(ind_ph)
        select case (Msk)
        case (eph)
            ind_ep = this%index_e( this%Shuf_ind_ep(i) )
            p = electrons%p(:, ind_ep )
            q = abs( electrons%q(ind_ep) )
            weight_max = maxval( (/ q, weight_max /) )
          case (pph)
            if((Ne+Np>Ne) .or. (Ne+Np>Np)) this%Shuf_ind_ep(i) = this%Shuf_ind_ep(i)-Ne
            ind_ep = this%index_p( this%Shuf_ind_ep(i) )
            p = positrons%p(:, ind_ep )
            q = abs( positrons%q(ind_ep) )
            weight_max = maxval( (/ q, weight_max /) )
          case (0)
            print*,'wrong mask'
        end select

        ! compute probability
        gama = sqrt( p(1)**2 + p(2)**2 + p(3)**2 + 1.0_p_double )
        w0 = sqrt( k(1)**2 + k(2)**2 + k(3)**2 )
        w1 = gama*w0 - p(1)*k(1) - p(2)*k(2) - p(3)*k(3)

        PP = Pnu_C( w1, gama, w0, Coeff, weight_max )
        !write(50,*) PP/P_max, gama
        if (rng % genrand_real1()<= PP/P_max ) then ! add couple to the list
          ! increase the number of scattering couples by 1
          this%ic_mates = this%ic_mates + 1
          ! update the list of scattering couples
          select case (Msk)
          case (eph)
              this%idx_c_mates(this%ic_mates,:) = (/ ind_ep, ind_ph, eph/)
              dumNe = dumNe + 1.0_p_double
            case (pph)
              this%idx_c_mates(this%ic_mates,:) = (/ ind_ep, ind_ph, pph/)
              dumNp = dumNp + 1.0_p_double

          end select
        endif
      enddo
      !write(49,*) this%ic_mates
    endif

    ! --- Update range of macroparticle indexes for next collision cell
    if(i_cell<num_c) then
      Ne = this%lst_e(i_cell+1)-this%lst_e(i_cell)
      Np = this%lst_p(i_cell+1)-this%lst_p(i_cell)
      Nph = this%lst_ph(i_cell+1)-this%lst_ph(i_cell)

    endif
  enddo

  !call cpu_time(te_list)
  !write(48,*) te_list-ti_list

end subroutine get_coll_listNTC



! --------------------------------------------------------------------------------------------------------------------
! ------------------------------------------------- COLLIDE PARTICLES ------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------

subroutine Ltz_boost(Ltz,p4)

  implicit none
  real(p_k_part), dimension(4), intent(in) :: p4
  real(p_k_part), dimension(4,4), intent(out) :: Ltz

  ! executive statements
  Ltz(1,1:4) = (/ p4(1) , - p4(2) , - p4(3) , - p4(4) /)
  Ltz(2,1:4) = (/ -p4(2) , 1.0_p_k_part + p4(2)**2/(1.0_p_k_part+p4(1)) , p4(2)*p4(3)/(1.0_p_k_part+p4(1)) , p4(2)*p4(4)/(1.0_p_k_part+p4(1)) /)
  Ltz(3,1:4) = (/ -p4(3) , p4(3)*p4(2)/(1.0_p_k_part+p4(1)) , 1.0_p_k_part + p4(3)**2/(1.0_p_k_part+p4(1)) , p4(3)*p4(4)/(1.0_p_k_part+p4(1)) /)
  Ltz(4,1:4) = (/ -p4(4) , p4(4)*p4(2)/(1.0_p_k_part+p4(1)) , p4(4)*p4(3)/(1.0_p_k_part+p4(1)) , 1.0_p_k_part + p4(4)**2/(1.0_p_k_part+p4(1)) /)

end subroutine Ltz_boost

function cdf_Compton( mu, rw1 )

  implicit none
  real(p_double) :: cdf_Compton ! cumulative probability between 0 and 1
  real(p_double) :: mu
  real(p_k_part) :: rw1


  cdf_Compton =( (1.0_p_double-2.0_p_double/rw1-2.0_p_double/rw1**2) * log( (1.0_p_double+2.0_p_double*rw1) / (1.0_p_double+rw1*(1.0_p_double-mu)) ) + &
          0.5_p_double * ( 1.0_p_double/(1.0_p_double+rw1*(1.0_p_double-mu))**2 - 1.0_p_double/(1.0_p_double+2.0_p_double*rw1)**2 ) + &
          1.0_p_double/rw1**2 * (2.0_p_double+rw1*(1.0_p_double-mu)) * ( (1.0_p_double+2.0_p_double*rw1) / (1.0_p_double+rw1*(1.0_p_double-mu)) - 1.0_p_double ) ) / &
          ( (1.0_p_double-2.0_p_double/rw1-2.0_p_double/rw1**2)*log( (1.0_p_double+2.0_p_double*rw1) )+0.5_p_double * ( 1.0_p_double-1.0_p_double/(1.0_p_double+2.0_p_double*rw1)**2 ) + 4.0_p_double/rw1 )

end function cdf_Compton

function pdf_Compton( mu, rw1 )

  implicit none
  real(p_double) :: pdf_Compton ! probability density function
  real(p_double) :: mu
  real(p_k_part) :: rw1

  pdf_Compton = rw1/(1.0_p_double+rw1*(1.0_p_double-mu))**2 * (1.0_p_double/(1.0_p_double+rw1*(1.0_p_double-mu))+rw1*(1.0_p_double-mu) + mu**2 ) / &
          ( (1.0_p_double-2.0_p_double/rw1-2.0_p_double/rw1**2)*log( 1.0_p_double+2.0_p_double*rw1 )+0.5_p_double * ( 1.0_p_double-1.0_p_double/(1.0_p_double+2.0_p_double*rw1)**2 ) + 4.0_p_double/rw1 )

end function pdf_Compton

subroutine angle_Compton_gen( fi, mu, w1)

  implicit none
  real(p_double), intent(out) :: fi, mu
  real(p_k_part), intent(in) :: w1

  ! dummy variables
  real(p_double) :: mu_old, mu_new, err, toll, bi, bs, dice
  integer :: n, n_max
  real(p_double), parameter :: pi = 3.14159265358979323846264_p_double ! pi


  ! executive statements
  fi = 2.0_p_double*pi*rng % genrand_real1()
  dice = rng % genrand_real1()


  ! Get the root of the nonlinear equation cdf_Compton(mu)-dice = 0

  err = 1.0_p_double

  n_max = 1000
  n = 0

  toll = 1E-12 !1.0_p_double/10.0_p_double**4
  bi = -1.0_p_double
  bs = 1.0_p_double

  do while( (err>toll) .and. (n<=n_max) ) ! Bisection method
    n = n + 1
    mu_new = (bi + bs ) / 2.0_p_double

    if ( (cdf_Compton(mu_new+toll/2, w1)-dice)*(cdf_Compton(mu_new-toll/2, w1)-dice)<0 ) then

      err = toll/2

    elseif ( (cdf_Compton(mu_new, w1)-dice)*(cdf_Compton(bi, w1)-dice)<0 ) then
      bs = mu_new
      err = abs( (bs-bi)/2.0_p_double )
    else ! if ( (cdf_Compton(mu_new, w1)-dice)*(cdf_Compton(bs, w1)-dice)<0 )
      bi = mu_new
      err = abs( (bs-bi)/2.0_p_double )
    endif

  enddo
  if (n==n_max) print*,' ** warning **   reached max iter for convergence on mu (compton), err: ', err

  mu = mu_new

  if (isnan(mu)) print*,"warning mu is nan"
  if (abs(mu) .gt. 1 ) print*,"warning mu is out of bounds"

end subroutine angle_Compton_gen

subroutine angle_Thomson_gen(fi,mu)

  implicit none

  real(p_double), intent(out) :: mu
  real(p_double), intent(out) :: fi

  real(p_double), parameter :: pi = 3.14159265358979323846264_p_double
  real(p_double) :: c

  fi = 2.0_p_double*pi*rng % genrand_real1()
  c = rng % genrand_real1()

  mu = (4.0_p_double*c-2.0_p_double+sqrt(5.0_p_double-16.0_p_double*c+16.0_p_double*c**2))**(1.0_p_double/3.0_p_double)- &
        1.0_p_double/(4.0_p_double*c-2.0_p_double+sqrt(5.0_p_double-16.0_p_double*c+16.0_p_double*c**2))**(1.0_p_double/3.0_p_double)

end subroutine angle_Thomson_gen

subroutine Compton_scatter(p,k,if_recoil)

  implicit none

  real(p_k_part), dimension(3), intent(inout) :: p
  real(p_k_part), dimension(3), intent(inout) :: k
  logical, intent(in) :: if_recoil

  ! dummy variables
  real(p_k_part), dimension(4) :: p4, k4, k4_erf, k4_new, p4_new! four vectors
  real(p_k_part), dimension(4,4) :: Lf, Lb ! forward and backward boosting matricies
  real(p_k_part), dimension(3) :: ek1, ek2, ek3, ek_perp ! versors of the basis in the electron rest frame (ERF)
  real(p_k_part), dimension(3) :: k_par, k_perp
  real(p_k_part) :: ekm, w2, err
  real(p_double) :: fi, mu
  integer :: i, j, im

  ! executive statements

  ! ----- preliminary quantities -----
  p4(1) = sqrt( p(1)**2 + p(2)**2 + p(3)**2 + 1.0_p_double ) ! particle gamma/energy
  p4(2:4) = p
  k4(1) = sqrt( k(1)**2 + k(2)**2 + k(3)**2 ) ! photon energy in the simulation frame of reference
  k4(2:4) = k

  ! forward boost matrix
  call Ltz_boost(Lf,p4)

  ! backward boost matrix (p -> -p)
  call Ltz_boost(Lb,(/p4(1), -p4(2), -p4(3) , -p4(4)/) )

  ! ----- STEPS -----

  ! 1 ) - boost forward the photon ( to the electron rest frame )
  do i=1,4
   k4_erf(i) = Lf(i,1)*k4(1) + Lf(i,2)*k4(2) + Lf(i,3)*k4(3) + Lf(i,4)*k4(4)
  enddo

  ! 2 ) - generate a base in the ERF [ek1, ek2, ek3]
  ek1 = k4_erf(2:4)/sqrt(k4_erf(2)**2+k4_erf(3)**2+k4_erf(4)**2)
  ekm = maxval(abs(ek1)) ! get the largest component of ek1
  im = 0

  do i=1,3
   if (abs(ek1(i)) == ekm) im = i
  enddo

  select case(im) ! generate a vector perp tp ek
   case (1)
    ek2 = (/ -(ek1(2)+ek1(3))/ek1(1), 1.0_p_k_part , 1.0_p_k_part /)
   case (2)
    ek2 = (/ 1.0_p_k_part, -(ek1(1)+ek1(3))/ek1(2), 1.0_p_k_part /)
   case (3)
    ek2 = (/ 1.0_p_k_part, 1.0_p_k_part, -(ek1(1)+ek1(2))/ek1(3) /)
  end select
  ek2 = ek2 / sqrt(ek2(1)**2 + ek2(2)**2 + ek2(3)**2 ) ! normalization of the versor
  ek3(1) = ek1(2)*ek2(3) - ek1(3)*ek2(2)
  ek3(2) = ek1(3)*ek2(1) - ek1(1)*ek2(3)
  ek3(3) = ek1(1)*ek2(2) - ek1(2)*ek2(1)
  ek3 = ek3 / sqrt(ek3(1)**2 + ek3(2)**2 + ek3(3)**2 )

  ! 3 ) - generate the scattering angles
  if (k4_erf(1)>3.0E-4) then
    call angle_Compton_gen( fi, mu, k4_erf(1)) ! fi[0;2pi], mu[-1;1] = cos(theta) theta[0;pi]
  else ! Thomson regime
    call angle_Thomson_gen( fi, mu )
  endif

  ek_perp = cos(fi)*ek2 + sin(fi)*ek3
  ek_perp = ek_perp / sqrt( ek_perp(1)**2 + ek_perp(2)**2 + ek_perp(3)**2 )


  ! 4 ) - four vector of outcoming photon in the ERF
  w2 = k4_erf(1)/( 1.0_p_k_part + k4_erf(1) * (1.0_p_double-mu) )

  k_par = w2*mu*ek1 ! component of the new k vector along ek1 (initial propagation direction of the incoming photon)
  k_perp = w2*sqrt(1.0_p_double-mu**2) * ek_perp

  k4_new(2:4) = (/ k_par(1)+k_perp(1) , k_par(2)+k_perp(2) , k_par(3)+k_perp(3) /)
  k4_new(1) = sqrt( k4_new(2)**2 + k4_new(3)**2 + k4_new(4)**2 )

  ! 5 ) - opt - recoil on the charged particle
  if (if_recoil) then
   p4_new = k4_erf-k4_new+(/ 1.0_p_k_part , 0.0_p_k_part , 0.0_p_k_part , 0.0_p_k_part /) ! new electron 4-momentum

   ! 6 ) - boost back to the simulation frame
   do i=1,4
     k4(i) = Lb(i,1)*k4_new(1) + Lb(i,2)*k4_new(2) + Lb(i,3)*k4_new(3) + Lb(i,4)*k4_new(4)
     p4(i) = Lb(i,1)*p4_new(1) + Lb(i,2)*p4_new(2) + Lb(i,3)*p4_new(3) + Lb(i,4)*p4_new(4)
   enddo

   k = k4(2:4)/sqrt( k4(2)**2 + k4(3)**2 + k4(4)**2 ) * k4(1)

   p = p4(2:4) /sqrt( p4(2)**2 + p4(3)**2 + p4(4)**2 ) * sqrt(p4(1)**2-1)

  else
   ! 6 ) - boost back to the simulation frame
   do i=1,4
     k4(i) = Lb(i,1)*k4_new(1) + Lb(i,2)*k4_new(2) + Lb(i,3)*k4_new(3) + Lb(i,4)*k4_new(4)
   enddo
    k = k4(2:4)

  endif

end subroutine Compton_scatter

subroutine Compton_scatter_new(p,k,if_recoil)

  use m_utilities, only : lorentz_transform

  implicit none

  real(p_k_part), dimension(3), intent(inout) :: p
  real(p_k_part), dimension(3), intent(inout) :: k
  logical, intent(in) :: if_recoil

  ! dummy variables
  real(p_k_part), dimension(3) :: u, p1, kERF, kERF_1, k1
  real(p_k_part) :: g, w, wERF, wERF_1, Diff
  real(p_k_part), dimension(3) :: ek1, ek2, ek3, ek_perp ! versors of the basis in the electron rest frame (ERF)
  real(p_k_part), dimension(3) :: k_par, k_perp
  real(p_k_part) :: ekm
  real(p_double) :: fi, mu
  integer :: i, j, im

  ! executive statements


  ! 1) ---- Boost to ERF
  g = sqrt( p(1)**2 + p(2)**2 + p(3)**2 + 1.0_p_double ) ! particle gamma/energy
  u = p/g ! particle velocity v/c
  w = sqrt( k(1)**2 + k(2)**2 + k(3)**2 ) ! photon energy in the simulation frame of reference
  call lorentz_transform( u, g, k, w, kERF )
  wERF = sqrt( kERF(1)**2 + kERF(2)**2 + kERF(3)**2 ) ! photon energy in the electron frame of reference

  ! 2 ) - generate a base in the ERF [ek1, ek2, ek3]
  ek1 = kERF/wERF
  ekm = maxval(abs(ek1)) ! get the largest component of ek1
  im = 0

  do i=1,3
   if (abs(ek1(i)) == ekm) im = i
  enddo

  select case(im) ! generate a vector perp tp ek
   case (1)
    ek2 = (/ -(ek1(2)+ek1(3))/ek1(1), 1.0_p_k_part , 1.0_p_k_part /)
   case (2)
    ek2 = (/ 1.0_p_k_part, -(ek1(1)+ek1(3))/ek1(2), 1.0_p_k_part /)
   case (3)
    ek2 = (/ 1.0_p_k_part, 1.0_p_k_part, -(ek1(1)+ek1(2))/ek1(3) /)
  end select
  ek2 = ek2 / sqrt(ek2(1)**2 + ek2(2)**2 + ek2(3)**2 ) ! normalization of the versor
  ek3(1) = ek1(2)*ek2(3) - ek1(3)*ek2(2)
  ek3(2) = ek1(3)*ek2(1) - ek1(1)*ek2(3)
  ek3(3) = ek1(1)*ek2(2) - ek1(2)*ek2(1)
  ek3 = ek3 / sqrt(ek3(1)**2 + ek3(2)**2 + ek3(3)**2 )

  ! 3 ) - generate the scattering angles
  if (wERF>3.0E-4) then
    call angle_Compton_gen( fi, mu, wERF) ! fi[0;2pi], mu[-1;1] = cos(theta) theta[0;pi]
  else ! Thomson regime
    call angle_Thomson_gen( fi, mu )
  endif

  ek_perp = cos(fi)*ek2 + sin(fi)*ek3
  ek_perp = ek_perp / sqrt( ek_perp(1)**2 + ek_perp(2)**2 + ek_perp(3)**2 )


  ! 4 ) - four vector of scattered photon in the ERF
  wERF_1 = wERF/( 1.0_p_k_part + wERF * (1.0_p_double-mu) )

  k_par = wERF_1*mu*ek1 ! component of the new k vector along ek1 (initial propagation direction of the incoming photon)
  k_perp = wERF_1*sqrt(1.0_p_double-mu**2) * ek_perp

  kERF_1 = (/ k_par(1)+k_perp(1) , k_par(2)+k_perp(2) , k_par(3)+k_perp(3) /)
  ! check if energy is correct
  Diff = sqrt(abs(wERF_1**2-kERF_1(1)**2 - kERF_1(2)**2 - kERF_1(3)**2))/wERF_1 ! relative error
  if ( Diff > 1.0E-6) then
    print*,'wrong energy', Diff
  endif

  ! 5 ) - boost back to the simulation frame
  call lorentz_transform( -u, g, kERF_1, wERF_1, k1 )

  ! 6 ) - opt - recoil on the charged particle
  if (if_recoil) then
   p1 = p + k - k1 ! new electron momentum
  endif

  k = k1
  p = p1

end subroutine Compton_scatter_new


end module m_qed_collide
