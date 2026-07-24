#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_radiat_cyl_define

#include "memory/memory.h"
!-----------------------------------------------------------------------------------------

use m_system
use m_parameters
use m_space, only : t_space
use m_grid_define, only : t_grid
use m_node_conf, only : t_node_conf, n_threads
use m_restart, only : t_restart_handle
use m_time_step, only : t_time_step, dt, n
use m_input_file, only : t_input_file, get_namelist
use m_current_define, only : t_current
use m_emf_define, only : t_emf
use m_vdf_define, only : t_vdf
use m_vdf_comm, only : t_vdf_msg
use detector_cls
use m_detlist
use m_diagnostic_utilities, only : p_diag_prec
use m_fparser, only : t_fparser, setup, eval, p_k_fparse, p_max_expr_len
use m_species_define, only : t_species, t_part_idx, t_spec_msg
use m_species_udist, only : use_accelerate, use_q_incr
use m_random, only : t_random, rng, new_random
use m_species_cyl_modes_define,  only: t_species_cyl_modes
use m_species_cyl_modes
use m_emf_cyl_modes ,only: t_emf_cyl_modes
use m_cyl_modes,only: t_cyl_modes
implicit none

private

!-----------------------------------------------------------------------------------------
! t_species_radiat_cyl
!  - Charged particle species with radiat pusher and merging
!-----------------------------------------------------------------------------------------
type, extends( t_species_cyl_modes ) :: t_species_radiat_cyl
  class( t_detector_list ), pointer ::  detector_list => null()
  real(p_k_part) :: cvol
  integer :: rad_select
  real(p_diag_prec) :: rad_gamma_limit
  real(p_diag_prec) :: rad_fraction
  character(len = p_max_expr_len) :: rad_math_expr = ""
  type(t_fparser) :: rad_func
  logical :: rad_tags

  ! particle data from previous time step
  real(p_k_part),dimension(:,:), pointer :: x_prev => null()
  real(p_k_part), dimension(:,:), pointer :: p_prev => null()
  integer, dimension(:,:), pointer :: ix_prev => null()

contains

  procedure :: create_particle_single_cell_p  => create_particle_single_cell_p_radiat_cyl
  procedure :: create_particle_single_cell    => create_particle_single_cell_radiat_cyl
  procedure :: read_input  => read_input_species_radiat_cyl
  procedure :: cleanup     => cleanup_species_radiat_cyl
  procedure :: init_buffer => init_buffer_spec_radiat_cyl
  procedure :: grow_buffer => grow_buffer_spec_radiat_cyl
  procedure :: copy_alld   => copy_alld_spec_radiat_cyl
  procedure :: rad_calc    => rad_calc_wrap_spec_radiat_cyl
  procedure :: rad_calc_inner => rad_calc_inner_spec_radiat_cyl
  procedure :: init        => init_species_radiat_cyl
  procedure :: push_cyl        => push_species_radiat_cyl
  procedure :: report      => report_species_radiat_cyl
  procedure :: restart_write => restart_write_species_radiat_cyl
  procedure :: restart_read => restart_read_species_radiat_cyl
  !!positions in 3D Space
  procedure :: position_3d_range => position_3d_cyl_range
  procedure :: position_3d => position_3d_cyl
  procedure :: prev_position_3d_range => prev_position_3d_cyl_range
  procedure :: prev_position_3d => prev_position_3d_cyl
  generic   :: get_prev_position_3d => prev_position_3d_range, prev_position_3d
  generic   :: get_position_3d=> position_3d_range, position_3d

end type t_species_radiat_cyl
!-----------------------------------------------------------------------------------------

public :: t_species_radiat_cyl

! Add interface to superclass init function to avoid type casting
interface
  subroutine create_particle_single_cell( this, ix, x, q )

    import t_species, p_k_part

    class(t_species), intent(inout) :: this
    integer, dimension(:), intent(in) :: ix
    real(p_k_part), dimension(:), intent(in) :: x
    real(p_k_part),               intent(in) :: q

  end subroutine
end interface

interface
  subroutine create_particle_single_cell_p( this, ix, x, p, q )

    import t_species, p_k_part

    class(t_species), intent(inout) :: this
    integer, dimension(:), intent(in) :: ix
    real(p_k_part), dimension(:), intent(in) :: x
    real(p_k_part), dimension(:), intent(in) :: p
    real(p_k_part),               intent(in) :: q

  end subroutine
end interface

interface
  subroutine create_particle_single_cell_cyl( this, ix, x, q )

  import t_species_cyl_modes , p_k_part

  class(t_species_cyl_modes), intent(inout) :: this
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x
  real(p_k_part),               intent(in) :: q

  end subroutine
end interface

interface
  subroutine create_particle_single_cell_p_cyl( this, ix, x, p, q )

  import t_species_cyl_modes, p_k_part

  class(t_species_cyl_modes), intent(inout) :: this
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x
  real(p_k_part), dimension(:), intent(in) :: p
  real(p_k_part),               intent(in) :: q

  end subroutine
end interface

contains

!-----------------------------------------------------------------------------------------
subroutine read_input_species_radiat_cyl( this, input_file, def_name, periodic, if_move, &
                                      grid, dt, read_prof, sim_options )

  implicit none

  class( t_species_radiat_cyl ),  intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  character(len = *), intent(in) :: def_name
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  real(p_double), intent(in) :: dt
  logical, intent(in) :: read_prof
  type(t_options), intent(in) :: sim_options

  integer :: ierr

  ! parameters for rad particle selection
  real(p_single)  :: rad_gamma_limit
  real(p_single)  :: rad_fraction
  character(len = p_max_expr_len) :: rad_math_expr
  logical :: rad_tags
  character(len = p_max_filename_len) :: rad_file_tags

  namelist /nl_selection/ rad_gamma_limit, rad_fraction, rad_math_expr, rad_tags, &
                          rad_file_tags

  ! call superclass read_input
  call this % t_species_cyl_modes % read_input( input_file, def_name, periodic, if_move, grid, &
                                      dt, read_prof, sim_options )

  !radiative particle selection section
  rad_gamma_limit = 0.0_p_k_part
  rad_fraction = 1.0_p_k_part

  ! default values for math function
  rad_math_expr = ''
  rad_file_tags = ''
  rad_tags = .False.

  call get_namelist( input_file, "nl_selection", ierr )
  if (ierr == 0) then
    read (input_file%nml_text, nml = nl_selection, iostat = ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "   Error reading selection parameters"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading selection parameters"
        write(0,*) "aborting..."
        stop
      endif
    endif
  endif

  ! Rad diagnostics data
  this%rad_gamma_limit = rad_gamma_limit
  this%rad_fraction = rad_fraction

  ! set math func variables
  this%rad_math_expr = trim(rad_math_expr)

  if (trim(rad_math_expr) /= '') then
    this%rad_select = 2
  elseif (rad_gamma_limit /= 0.0_p_k_part .or. rad_fraction < 1.0_p_k_part) then
    this%rad_select = 1
  else
    this%rad_select = 0
  endif

  this%rad_tags=rad_tags
  if(this%rad_tags) then
    if(this%diag%ndump_fac_tracks<1) then
      !to avoid dumping unwanted track files
      this%diag%ndump_fac_tracks=huge(this%diag%ndump_fac_tracks)
      this%diag%tracks%niter=huge(this%diag%ndump_fac_tracks)
    endif
    if(rad_file_tags/='') this%diag%tracks%file_tags = rad_file_tags
    if ( .not. file_exists(this%diag%tracks%file_tags)) then
      print *, "(*error*) ndump_fac_tracks is set but tags file '" // &
               trim(rad_file_tags) // "' cannot be found."
      stop
    endif
  endif

  allocate( this%detector_list )

  call this%detector_list%read_input(input_file,dt)

end subroutine read_input_species_radiat_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_species_radiat_cyl( this, sp_id, interpolation, grid_center, grid, g_space, &
                                emf, jay, no_co, send_vdf, recv_vdf, bnd_cross, &
                                node_cross, send_spec, recv_spec, ndump_fac, restart, &
                                restart_handle, sim_options, tstep, tmin, tmax )

  implicit none

  class( t_species_radiat_cyl ), intent(inout) :: this

  integer, intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in) :: grid
  type( t_space ), intent(in) :: g_space
  class( t_emf ), intent(inout) :: emf
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_vdf, recv_vdf
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_spec, recv_spec
  integer, intent(in) :: ndump_fac
  class( t_current ), intent(inout) :: jay
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type(t_options), intent(in) :: sim_options
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax

  integer :: n_x_dim, ierr

  call this % t_species_cyl_modes % init( sp_id, interpolation, grid_center, grid, g_space, emf, &
                                jay, no_co, send_vdf, recv_vdf, bnd_cross, node_cross, &
                                send_spec, recv_spec, ndump_fac, restart, restart_handle,&
                                sim_options, tstep, tmin, tmax )


  n_x_dim = this%get_n_x_dims()
  call alloc(this%x_prev,  (/n_x_dim,max(this%num_par_max,1) /) )
  call alloc(this%ix_prev, (/p_x_dim,max(this%num_par_max,1) /) )
  call alloc(this%p_prev,  (/p_p_dim,max(this%num_par_max,1) /) )

  ! set variable list for rad_math_expr
  if (this%rad_math_expr /= '') then
    select case (p_x_dim)
    case (1)
      call setup(this%rad_func, trim(this%rad_math_expr), &
                 (/'x1', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)

    case (2)
      call setup(this%rad_func, trim(this%rad_math_expr), &
                 (/'x1', 'x2', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)

    case (3)
      call setup(this%rad_func, trim(this%rad_math_expr), &
                 (/'x1', 'x2', 'x3', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)
    end select

    ! check if function compiled ok
    if (ierr /= 0) then
      ERROR("Error compiling supplied function :")
      ERROR(trim(this%rad_math_expr))
      call abort_program(-1)
    endif
  endif

  ! cell volume
  this%cvol = product( this%dx(1:p_x_dim) )

  call this%detector_list%init( this%name, no_co, restart, restart_handle )

end subroutine
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! This routine dispatches the code to the appropriate push.
! Currently the following of push are available:
!  i)   Standard PIC push (Fortran or Hardware Optimized)
!  ii)  Accelerate pusher (used for initializning beams )
!  iii)  Radiation cooling pusher
!-----------------------------------------------------------------------------------------
subroutine push_species_radiat_cyl( this, emf, jay_cyl_m, t, tstep, tid, n_threads, options )

  implicit none

  class( t_species_radiat_cyl ), intent(inout) :: this
  class( t_emf_cyl_modes ), intent( inout )  ::  emf
  type( t_cyl_modes ), intent(inout) :: jay_cyl_m

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  type( t_options ), intent(in) :: options

  integer, intent(in) :: tid        ! local thread id
  integer, intent(in) :: n_threads  ! total number of threads

  ! if before push start time return silently
  if ( t >= this%push_start_time ) then

    call this%copy_alld(tstep,tid,n_threads)


    call this%t_species_cyl_modes%push_cyl( emf, jay_cyl_m, t, tstep, tid, n_threads, options )

  endif

end subroutine push_species_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! cleanup all dynamic memory allocated by this object
!-----------------------------------------------------------------------------------------
subroutine cleanup_species_radiat_cyl( this )

  implicit none

  class( t_species_radiat_cyl ), intent(inout) :: this

  ! cleanup memory used by the particle buffers
  call freemem(this%x_prev)
  call freemem(this%ix_prev)
  call freemem(this%p_prev)
  call this%detector_list%cleanup()

  call this % t_species % cleanup()

end subroutine cleanup_species_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine rad_calc_wrap_spec_radiat_cyl(this,t,tstep,no_co)

  implicit none

  class( t_species_radiat_cyl ), intent(inout) :: this
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  class( t_node_conf ), intent(in) :: no_co

  integer :: nt, tid

  ! if before push start time return silently
  if ( t < this%push_start_time ) return

  nt = n_threads(no_co)
  if ( nt > 1 ) then

    !$omp parallel do
    do tid = 1, nt
      call this % rad_calc_inner( t, tstep, tid-1, nt )
    enddo
    !$omp end parallel do

  else

    call this % rad_calc_inner( t, tstep, 0, 1 )

  endif

end subroutine rad_calc_wrap_spec_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine rad_calc_inner_spec_radiat_cyl(this,t,tstep,tid,nt)

  implicit none

  class( t_species_radiat_cyl ), intent(inout) :: this
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  integer, intent(in) :: tid, nt

  real(p_k_part) :: dtcycle, q(p_cache_size)
  real(p_k_part), dimension(3,p_cache_size) :: x3d, x_prev3d
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: beta, betaprev
  integer :: pp, i, idx, n_x_dim, ptrcur, np
  real(p_k_fparse) :: raw_eval
  real(p_k_fparse), dimension(p_p_dim + this%get_n_x_dims() + 2) :: rad_var
  real(p_k_part), dimension(4) :: pos ! was 3
  real(p_k_part):: gamma, rnd, gamma_par_sqr, gamma_lim_sqr
  class(t_random), pointer :: rng_loc

  dtcycle = real( dt(tstep), p_k_part )

  ! In this construction, we do not support n_x_dim > 3
  n_x_dim = this%get_n_x_dims()

  ! In case running with 1 or 2 dimensions
  x3d = 0.0_p_k_part
  x_prev3d = 0.0_p_k_part

  select case ( this%coordinates )
  !case ( p_cylindrical_b )
    ! Currently not implemented
  !  continue
  case default

    ! If tags are used, other selection methods are ignored
    if (this%rad_tags) then

      if ( this%diag%tracks%npresent > 0 ) then

        do ptrcur = 1, this%diag%tracks%npresent, p_cache_size
          ! check if last copy of table and set np
          if( ptrcur + p_cache_size > this%diag%tracks%npresent ) then
            np = this%diag%tracks%npresent - ptrcur + 1
          else
            np = p_cache_size
          endif

          do i = 1, np

            idx=this%diag%tracks%present(ptrcur+i-1)
            pp=this%diag%tracks%tracks(idx)%part_idx

            call this%get_position_3d(pp,x3d(:,i))
            call this%get_prev_position_3d(pp,x_prev3d(:,i))

            beta(1:p_p_dim,i) = this%p(:,pp)/sqrt(this%p(1,pp)**2&
                               +this%p(2,pp)**2+this%p(3,pp)**2+1.0_p_k_part)

            betaprev(1:p_p_dim,i) = this%p_prev(:,pp)/sqrt(&
                                    this%p_prev(1,pp)**2&
                                   +this%p_prev(2,pp)**2&
                                   +this%p_prev(3,pp)**2+1.0_p_k_part)

            q(i) = this%q(pp) * this%cvol

          enddo

          call this%detector_list%fill_emf(x3d,beta,t,x_prev3d,betaprev,q,np,tid)

        enddo

      endif

    else

      select case( this%rad_select )
      case(0)

        ! Select all particles
        do ptrcur = 1, this%num_par, p_cache_size
          pp = ptrcur
          ! check if last copy of table and set np
          if( ptrcur + p_cache_size > this%num_par ) then
            np = this%num_par - ptrcur + 1
          else
            np = p_cache_size
          endif

          call this%get_position_3d(pp,pp+np-1,x3d(:,1:np))
          call this%get_prev_position_3d(pp,pp+np-1,x_prev3d(:,1:np))

          do i = 1, np
            betaprev(1:p_p_dim,i)=this%p_prev(:,pp)/sqrt(&
                                  this%p_prev(1,pp)**2&
                                 +this%p_prev(2,pp)**2&
                                 +this%p_prev(3,pp)**2+1.0_p_k_part)

            beta(1:p_p_dim,i)=this%p(:,pp)/sqrt(this%p(1,pp)**2&
                             +this%p(2,pp)**2+this%p(3,pp)**2+1.0_p_k_part)

            q(i) =  this%cvol*this%q(pp)
            pp = pp + 1
          enddo

          call this%detector_list%fill_emf(x3d,beta,t,x_prev3d,betaprev,q,np,tid)

        enddo

      case(1)

        if ( nt > 1 ) then
          rng_loc => new_random()
          ! Provide a different seed each time step, but give same seed for all threads
          call rng_loc % init_genrand_scalar( 1 + n(tstep) )
        else
          rng_loc => rng
        endif

        ! Select particles based on gamma limit
        gamma_lim_sqr = this%rad_gamma_limit**2
        pp = 1
        np = 0 ! number of particles we have found
        do
          if ( pp > this%num_par ) then
            if ( np > 0 ) then
              call this%detector_list%fill_emf(x3d,beta,t,x_prev3d,betaprev,q,np,tid)
            endif
            exit
          endif

          if ( np == p_cache_size ) then
            call this%detector_list%fill_emf(x3d,beta,t,x_prev3d,betaprev,q,np,tid)
            np = 0
          endif

          gamma_par_sqr = this%p(1,pp)**2 + this%p(2,pp)**2 &
                        + this%p(3,pp)**2 + 1.0_p_k_part

          if ( gamma_par_sqr >= gamma_lim_sqr ) then

            call rng_loc % harvest_real2( rnd, this%ix(:,pp) )
            if ( rnd <= this%rad_fraction ) then

              np = np + 1

              call this%get_position_3d(pp,x3d(:,np))
              call this%get_prev_position_3d(pp,x_prev3d(:,np))

              betaprev(1:p_p_dim,np)=this%p_prev(:,pp)/sqrt(&
                                      this%p_prev(1,pp)**2&
                                     +this%p_prev(2,pp)**2&
                                     +this%p_prev(3,pp)**2+1.0_p_k_part)

              beta(1:p_p_dim,np)=this%p(:,pp)/sqrt(gamma_par_sqr)

              q(np) = this%q(pp) * this%cvol

            endif

          endif

          pp = pp + 1
        enddo

        if ( nt > 1 ) then
          deallocate(rng_loc)
        endif

      case(2)

        if ( nt > 1 ) then
          rng_loc => new_random()
          ! Provide a different seed each time step, but give same seed for all threads
          call rng_loc % init_genrand_scalar( 1 + n(tstep) )
        else
          rng_loc => rng
        endif

        ! Select particles based on math expression
        pp = 1
        np = 0 ! number of particles we have found
        do
          if ( pp > this%num_par ) then
            if ( np > 0 ) then
              call this%detector_list%fill_emf(x3d,beta,t,x_prev3d,betaprev,q,np,tid)
            endif
            exit
          endif

          if ( np == p_cache_size ) then
            call this%detector_list%fill_emf(x3d,beta,t,x_prev3d,betaprev,q,np,tid)
            np = 0
          endif

          call this%get_position(pp,pos)
          rad_var(1:n_x_dim) = real( pos(1:n_x_dim), p_k_fparse )
          rad_var(n_x_dim+1) = real( this%p(1,pp), p_k_fparse )
          rad_var(n_x_dim+2) = real( this%p(2,pp), p_k_fparse )
          rad_var(n_x_dim+3) = real( this%p(3,pp), p_k_fparse )
          gamma=sqrt(this%p(1,pp)**2+&
                     this%p(2,pp)**2+&
                     this%p(3,pp)**2+1.0_p_k_part)
          rad_var(n_x_dim + p_p_dim + 1) = real( gamma, p_k_fparse )
          rad_var(n_x_dim + p_p_dim + 2) = real( t, p_k_fparse )

          ! evaluate
          raw_eval = eval( this%rad_func, rad_var )

          if (raw_eval > 0) then

            call rng % harvest_real2( rnd, this%ix(:,pp) )
            if ( rnd <= this%rad_fraction ) then

              np = np + 1
              call this%get_position_3d(pp,x3d(:,np))
              call this%get_prev_position_3d(pp,x_prev3d(:,np))
              betaprev(1:p_p_dim,np)=this%p_prev(:,pp)/sqrt(&
                                     this%p_prev(1,pp)**2&
                                    +this%p_prev(2,pp)**2&
                                    +this%p_prev(3,pp)**2+1.0_p_k_part)

              beta(1:p_p_dim,np)=this%p(:,pp)/gamma

              q(np) = this%q(pp) * this%cvol

            endif

          endif

          pp = pp + 1
        enddo

        if ( nt > 1 ) then
          deallocate(rng_loc)
        endif

      end select

    endif

  end select


end subroutine rad_calc_inner_spec_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine create_particle_single_cell_radiat_cyl( this, ix, x, q )

  implicit none

  class(t_species_radiat_cyl), intent(inout) :: this
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x
  real(p_k_part),               intent(in) :: q

  integer :: n_x_dim

  call create_particle_single_cell_cyl( this, ix, x, q )

  n_x_dim = this%get_n_x_dims()
  this%x_prev(1:n_x_dim,this%num_par)=this%x(1:n_x_dim,this%num_par)
  this%ix_prev(1:p_x_dim,this%num_par)=this%ix(1:p_x_dim,this%num_par)

end subroutine create_particle_single_cell_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine create_particle_single_cell_p_radiat_cyl( this, ix, x, p, q )

  implicit none

  class(t_species_radiat_cyl), intent(inout) :: this
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x
  real(p_k_part), dimension(:), intent(in) :: p
  real(p_k_part),               intent(in) :: q

  integer :: n_x_dim

  call create_particle_single_cell_p_cyl( this, ix, x, p, q )

  n_x_dim = this%get_n_x_dims()
  this%x_prev(1:n_x_dim,this%num_par)  = x(1:n_x_dim)
  this%ix_prev(1:p_x_dim,this%num_par) = ix(1:p_x_dim)
  this%p_prev(1:p_p_dim,this%num_par)  = p(1:p_p_dim)

end subroutine create_particle_single_cell_p_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine init_buffer_spec_radiat_cyl( this, num_par_req )

  implicit none

  class(t_species_radiat_cyl), intent(inout) :: this
  integer, intent(in) :: num_par_req

  integer :: n_x_dim

  call this%t_species_cyl_modes%init_buffer(num_par_req)

  ! setup position buffer
  n_x_dim = this%get_n_x_dims()
  call freemem( this%x_prev )
  call alloc(  this%x_prev, (/ n_x_dim, this%num_par_max /))

  ! initialize particle cell information
  call freemem(this%ix_prev )
  call alloc( this%ix_prev, (/ p_x_dim, this%num_par_max /) )

  ! setup momenta buffer
  call freemem( this%p_prev )
  call alloc( this%p_prev, (/ p_p_dim, this%num_par_max /) )

end subroutine init_buffer_spec_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Grow the particle buffers
!-----------------------------------------------------------------------------------------
subroutine grow_buffer_spec_radiat_cyl( this, num_par_req )

  implicit none

  class(t_species_radiat_cyl), intent(inout) :: this
  integer, intent(in) :: num_par_req

  integer :: num_par_old, num_par_new, n_x_dim

  !call this%t_species_cyl_modes%grow_buffer(num_par_req)
  call this%t_species_cyl_modes%grow_buffer(num_par_req)
  ! The buffer size must always be a multiple of vector width in size because of SIMD code
  num_par_new = (( num_par_req + p_vecwidth - 1 ) / p_vecwidth) * p_vecwidth

  if ( this%num_par > 0 ) then

    num_par_old = this%num_par

    write(0,'(A,I0,A,A)') '[', mpi_node(), &
                          '] (* warning *) resizing particle buffers for this ', &
                          trim(this%name)
    write(0,'(A,I0,A,I0,A,I0)') '[', mpi_node(), '] (* warning *) Buffer size: ', &
                                this%num_par_max, ' -> ', num_par_new
    write(0,'(A,I0,A,I0)') '[', mpi_node(), &
                           '] (* warning *) Number of particles currently in buffer: ', &
                           this%num_par

    if ( num_par_new <= num_par_old ) then
      ERROR('Invalid size for new buffer')
      call abort_program( p_err_invalid )
    endif

    ! particle positions (may not be p_x_dim)
    n_x_dim = this%get_n_x_dims()
    ! setup position buffer
    call freemem( this%x_prev )
    call alloc(  this%x_prev, (/ n_x_dim, num_par_new /))

    ! initialize particle cell information
    call freemem(this%ix_prev )
    call alloc( this%ix_prev, (/ p_x_dim, num_par_new /) )

    ! setup momenta buffer
    call freemem( this%p_prev )
    call alloc( this%p_prev, (/ p_p_dim, num_par_new /) )

  else

    ! no particles in buffer, simply reallocate the buffers
    call this % init_buffer( num_par_new )

  endif

end subroutine grow_buffer_spec_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine copy_alld_spec_radiat_cyl( this, tstep, tid, n_threads )

  implicit none

  class( t_species_radiat_cyl ), intent(inout) :: this
  type( t_time_step ), intent(in) :: tstep
  integer, intent(in) :: tid        ! local thread id
  integer, intent(in) :: n_threads  ! total number of threads

  integer :: pp, np, ptrcur, chunk, i0, i1, push_type

  push_type = this%push_type
  if ( use_accelerate( this%udist, n(tstep) ) ) push_type = p_beam_accel
  if ( use_q_incr( this%udist, n(tstep) ) ) push_type = p_beam_accel

  ! Chunk size may be different for vectorized pusher
  if ( push_type == p_simd ) then

    ! range of particles for each thread
    chunk = ( this%num_par + n_threads - 1 ) / n_threads
    ! round up to the nearest multiple of vector width
    chunk = ( ( chunk + p_vecwidth - 1 ) / p_vecwidth ) * p_vecwidth
    i0    = tid * chunk + 1
    i1    = min( (tid+1) * chunk, this%num_par )

  else

    ! range of particles for each thread
    chunk = ( this%num_par + n_threads - 1 ) / n_threads
    i0    = tid * chunk + 1
    i1    = min( (tid+1) * chunk, this%num_par )

  endif

  ! loop through all particles
  do ptrcur = i0, i1, p_cache_size
    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif
    pp = ptrcur
    call memcpy(this%x_prev(:,pp:pp+np-1),this%x(:,pp:pp+np-1),np*this%get_n_x_dims())
    call memcpy(this%ix_prev(:,pp:pp+np-1),this%ix(:,pp:pp+np-1),np*p_x_dim)
    call memcpy(this%p_prev(:,pp:pp+np-1),this%p(:,pp:pp+np-1),np*p_p_dim)
  enddo

end subroutine copy_alld_spec_radiat_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine report_species_radiat_cyl( this, emf, g_space, grid, no_co, tstep, t, tmin, &
                                  send_msg, recv_msg )

  implicit none

  class( t_species_radiat_cyl ), intent(inout) :: this
  class( t_emf ),      intent(inout) :: emf
  type( t_space ),       intent(in) :: g_space
  class( t_grid ),        intent(in) :: grid
  class( t_node_conf ),   intent(in) :: no_co
  type( t_time_step ),   intent(in) :: tstep
  real(p_double),        intent(in) :: t, tmin
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  call this%t_species_cyl_modes%report(emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg)

  call this%detector_list%report(t,tstep,no_co)

end subroutine report_species_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine restart_write_species_radiat_cyl( this, restart_handle )

  implicit none

  class( t_species_radiat_cyl ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  call this%t_species_cyl_modes%restart_write( restart_handle )

  ! Extra species_radiat quantities
  call this%detector_list%restart_write( restart_handle )

end subroutine restart_write_species_radiat_cyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine restart_read_species_radiat_cyl( this, restart_handle )

  implicit none

  class( t_species_radiat_cyl ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  call this%t_species_cyl_modes%restart_read( restart_handle )

  ! Extra species_radiat quantities
  call this%detector_list%restart_read( restart_handle )

end subroutine restart_read_species_radiat_cyl
!-----------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------
subroutine position_3d_cyl( this, idx, pos )


  implicit none

  class( t_species_radiat_cyl ), intent(in) :: this
  integer, intent(in) :: idx
  real( p_k_part ), dimension(:), intent(out) :: pos

  integer :: j,n_x_dim

  n_x_dim = this%get_n_x_dims()

  pos(1) = real( ( (this%ix( 1, idx ) + this%my_nx_p(p_lower, 1) - 2) + &
            this%x( 1, idx ) ) * this%dx( 1 ) + &
            this%g_box( p_lower, 1 ), p_k_part )
  !pos is only 3d
  do j = 3, n_x_dim
      pos(j-1) = real(this%x(j,idx),p_k_part)
  enddo

end subroutine position_3d_cyl
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine prev_position_3d_cyl( this, idx, pos )


  implicit none

  class( t_species_radiat_cyl ), intent(in) :: this
  integer, intent(in) :: idx
  real( p_k_part ), dimension(:), intent(out) :: pos

  integer :: j,n_x_dim

  n_x_dim = this%get_n_x_dims()

  pos(1) = real( ( (this%ix_prev( 1, idx ) + this%my_nx_p(p_lower, 1) - 2) + &
            this%x_prev( 1, idx ) ) * this%dx( 1 ) + &
            this%g_box( p_lower, 1 ), p_k_part )

  do j = 3, n_x_dim
      pos(j-1) = real(this%x_prev(j,idx),p_k_part)
  enddo

end subroutine prev_position_3d_cyl
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine position_3d_cyl_range( this, idx0, idx1, pos )

  implicit none

  class( t_species_radiat_cyl ), intent(in) :: this
  integer, intent(in) :: idx0, idx1
  real( p_k_part ), dimension(:,:), intent(out) :: pos

  integer :: j,k
  integer ixmin(p_x_dim)

  real( p_k_part ) :: xmin( p_x_dim ), dx( p_x_dim )
  integer :: n_x_dim

  n_x_dim = this%get_n_x_dims()

  xmin(1:p_x_dim)  = real( this%g_box( p_lower, 1:p_x_dim ), p_k_part )
  ixmin(1:p_x_dim) = this%my_nx_p( 1, 1:p_x_dim ) - 2
  dx(1:p_x_dim) = real( this%dx( 1:p_x_dim ), p_k_part )


  do j= idx0, idx1
    pos(1,j-idx0+1) = xmin(1) + dx(1) * (real(( this%ix( 1, j ) + ixmin(1) ) + &
          this%x( 1, j ), p_k_part ) )
    do k = 3, n_x_dim
        pos(k-1,j-idx0+1) = real(this%x(k,j),p_k_part)
    enddo
  enddo

end subroutine position_3d_cyl_range
!---------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------
subroutine prev_position_3d_cyl_range( this, idx0, idx1, pos )

  implicit none

  class( t_species_radiat_cyl ), intent(in) :: this
  integer, intent(in) :: idx0, idx1
  real( p_k_part ), dimension(:,:), intent(out) :: pos

  integer :: j,k
  integer ixmin(p_x_dim)

  real( p_k_part ) :: xmin( p_x_dim ), dx( p_x_dim )
  integer :: n_x_dim

  n_x_dim = this%get_n_x_dims()


  xmin(1:p_x_dim)  = real( this%g_box( p_lower, 1:p_x_dim ), p_k_part )
  ixmin(1:p_x_dim) = this%my_nx_p( 1, 1:p_x_dim ) - 2
  dx(1:p_x_dim) = real( this%dx( 1:p_x_dim ), p_k_part )


  do j= idx0, idx1
    pos(1,j-idx0+1) = xmin(1) + dx(1) * (real(( this%ix_prev( 1, j ) + ixmin(1) ) + &
          this%x_prev( 1, j ), p_k_part ) )
    do k = 3, n_x_dim
        pos(k-1,j-idx0+1) = real(this%x_prev(k,j),p_k_part)
    enddo
  enddo

end subroutine prev_position_3d_cyl_range
!---------------------------------------------------------------------------------------




end module m_species_radiat_cyl_define
