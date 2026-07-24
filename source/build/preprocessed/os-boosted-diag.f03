# 1 "boost/os-boosted-diag.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "boost/os-boosted-diag.f03" 2
! this module provides a data structure which abstracts
! boosted grid diagnostics. specifically it handles the conversion
! of lab and boosted coordinates


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
# 7 "boost/os-boosted-diag.f03" 2
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
# 8 "boost/os-boosted-diag.f03" 2

module m_boosted_diag

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
# 12 "boost/os-boosted-diag.f03" 2

use m_system
use m_parameters

use stringutil
use m_vdf_define
use m_vdf_memory
use m_diagnostic_utilities
use m_grid_define
use m_grid
use m_space
use m_time_step
use m_utilities
use m_boosted_grid_mgr
use m_input_file
use m_parameters
use m_restart
use m_diagfile, only : p_time_length
use m_fparser

implicit none


character(len=*), parameter :: p_boosted_diag_rst_id = "boosted diag rst data - 0x0000"

! configurations for boosted diagnostics to work with
! these change how the data is stored, updated, and lorentz-transformed
integer, parameter, public :: p_boosted_fields = 0
integer, parameter, public :: p_boosted_4vector = 1
integer, parameter :: p_slice_max = 20 ! maximum number of the slices we can obtain

! Mirrors the p_max_spname_len parameter in os-spec-define.f03
integer, parameter, private :: p_max_spname_len_ = 64

! Minimal block size to allocate and grow the raw particle buffers
integer, parameter :: p_raw_buf_block = 131072 ! 128k

! Chunk size for raw particle diagnostic
integer, parameter :: p_raw_chunk_size = 1000


! note: in all the code below, the 1 and 2 indexing makes life a lot easier.
! just remember that the simulation frame is frame 1, and the lab frame
! which we're transforming to is frame 2.


!-------------------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------------------
type :: t_boosted_diag

    ! shared object to keep track of coordinates for each calculation
    type(t_boosted_grid_mgr), pointer :: grid_mgr => null()

    ! spatial dimensionality of the data
    integer :: x_dim_

    ! current options: p_boosted_4vector, p_boosted_fields
    integer :: diag_type
    character(len=p_max_spname_len_) :: spec_name

    integer :: ndump_fac ! ndump rate in simulation time t1
    integer :: ndump_fac_raw ! ndump rate in simulation time t1 for raw data

    ! previous and current data
    ! the implementation could be generalized by instead keeping a slot
    ! for generic data here (e.g. type(t_vdf) :: cur_data, prev_data),
    ! but the current implementation makes life easier since there aren't any other boosted
    ! quantities we would want to look at.

    ! type(t_vdf), pointer :: cur_e, cur_b
    type(t_vdf), pointer :: prev_e => null(), prev_b => null()
    type(t_vdf), pointer :: prev_4vector => null()

    ! buffers storing partial data of a time slice in frame 2
    ! for each timeslice
    type(t_vdf), pointer :: buffered_e => null(), buffered_b => null()
    type(t_vdf), pointer :: buffered_4vector => null()

    ! buffers storing the raw particle data of a time slice in frame 2 for each timeslice
    real(p_k_part), dimension(:,:), pointer :: buffered_part => null()
    integer, dimension(:,:), pointer :: buffered_tags => null()
    integer, dimension(:), pointer :: idx => null(), buffered_diag_idx => null(), &
                                      num_par_diag => null()
    integer :: spec_dim, num_par, num_par_max
    logical :: add_tag
    real(p_k_part) :: raw_gamma_limit, raw_fraction
    character(len = p_max_expr_len) :: raw_math_expr = " "
    type(t_fparser) :: raw_func
    ! the below option currently does nothing, but could be implemented in the future
    logical, dimension(p_x_dim) :: raw_if_pos_ref_box

    ! current and previous leftmost global cell index (needed for interp)
    real(p_double) :: xmin1_prev

    ! size of data buffer, calculated to fit max amount of data given
    integer :: data_buf_size

    ! number of data currently stored in buffer
    integer :: num_data_in_buf

    ! index of grid % boosted_cell_coords that buffer index 0 corresponds to
    integer :: buf_start_idx

    ! directory where this data will be stored
    character(len=1024) :: basepath

    ! number of components to write to disk
    integer :: num_diag_components
    ! number of slice components to write to disk
    integer :: num_diag_components_slice

    ! indices of the specific components to write to disk (size = ndiags)
    ! integer, dimension(6) :: diag_components

    integer, dimension(6) :: diag_components
    integer, dimension(p_slice_max) :: diag_components_slice ! e1_boost_slice, etc.
    integer, dimension(p_slice_max) :: direction
    integer, dimension(p_slice_max) :: gipos

    integer :: num_diag_components_to_buffer
    integer, dimension(6) :: diag_components_to_buffer

    ! names of all components
    character(len=128), dimension(6) :: component_names
    character(len=128), dimension(p_slice_max) :: component_names_slice ! slice name

    contains

    procedure :: restart_write => restart_write_boosted_diag
    procedure :: restart_read => restart_read_boosted_diag

    procedure :: init => init_boosted_diag
    procedure :: compute_data_buf_size => compute_data_buf_size_boosted_diag
    procedure :: init_grid_mgr => init_grid_mgr_boosted_diag
    procedure :: cleanup => cleanup_boosted_diag
    procedure :: read_input => read_input_boosted_diag

    procedure :: parse_reports_string => parse_reports_string_boosted_diag
    procedure :: init_components_to_buffer => init_components_to_buffer_boosted_diag
    procedure :: get_dir_fname_label => get_dir_fname_label_boosted_diag
    procedure :: get_dir_fname_label_slice => get_dir_fname_label_boosted_diag_slice

    procedure :: create_files => create_files_boosted_diag
    procedure :: create_files_component => create_files_component_boosted_diag
    procedure :: create_files_component_slice => create_files_component_boosted_diag_slice
    procedure :: create_files_raw => create_files_raw_boosted_diag
    procedure :: populate_raw_diagFile => populate_raw_diagFile_boosted_diag
    procedure :: init_raw_func => init_raw_func_boosted_diag

    procedure :: report => report_boosted_diag
    procedure :: report_component => report_component_boosted_diag
    procedure :: report_component_slice => report_component_boosted_diag_slice

    procedure :: extract => extract_boosted_diag

    procedure :: grow_buffer => grow_buffer_boosted_diag
    procedure :: push_back => push_back_boosted_diag

    procedure :: allocate_objs => allocate_objs_boosted_diag

end type
!-------------------------------------------------------------------------------------------------

interface generate_sort_idx
  module procedure generate_sort_idx_
end interface

public :: t_boosted_diag, generate_sort_idx


contains



!-------------------------------------------------------------------------------------------------
subroutine restart_write_boosted_diag( this, restart_handle )

    use m_restart

    implicit none

    class( t_boosted_diag ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for boosted_diag object.'
    integer :: ierr

    call restart_io_write("p_boosted_diag_rst_id", p_boosted_diag_rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",201)

    call restart_io_write("this%xmin1_prev", this%xmin1_prev, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",204)

    call restart_io_write("this%buf_start_idx", this%buf_start_idx, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",207)

    call restart_io_write("this%num_data_in_buf", this%num_data_in_buf, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",210)

    select case( this%diag_type )

    case( p_boosted_4vector )

        if ( this%ndump_fac > 0 ) then

            call this % buffered_4vector % write_checkpoint( restart_handle )
            call this % prev_4vector % write_checkpoint( restart_handle )

        endif

        call restart_io_write("this%num_par", this%num_par, restart_handle, ierr)
        call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",224)

        call restart_io_write("this%num_par_max", this%num_par_max, restart_handle, ierr)
        call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",227)

        if ( this%num_par > 0 ) then
            call restart_io_write("this%buffered_part(:,1:this%num_par)", this%buffered_part(:,1:this%num_par), restart_handle, ierr)
            call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",231)

            call restart_io_write("this%buffered_diag_idx(1:this%num_par)", this%buffered_diag_idx(1:this%num_par), restart_handle, ierr)
            call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",234)

            call restart_io_write("this%num_par_diag", this%num_par_diag, restart_handle, ierr)
            call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",237)

            if ( this%add_tag ) then
                call restart_io_write("this%buffered_tags(:,1:this%num_par)", this%buffered_tags(:,1:this%num_par), restart_handle, ierr)
                call check_error(ierr,err_msg,p_err_rstwrt,"boost/os-boosted-diag.f03",241)
            endif
        endif

    case( p_boosted_fields )

        call this % buffered_e % write_checkpoint( restart_handle )
        call this % buffered_b % write_checkpoint( restart_handle )

        call this % prev_e % write_checkpoint( restart_handle )
        call this % prev_b % write_checkpoint( restart_handle )

    end select

end subroutine restart_write_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine restart_read_boosted_diag( this, restart_handle )

    use m_restart

    implicit none

    class( t_boosted_diag ), intent(inout) :: this
    type( t_restart_handle ), intent(in) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for boosted_diag object.'
    integer :: ierr
    character(len=len(p_boosted_diag_rst_id)) :: rst_id

    integer, dimension(3) :: data_buf_nx


    call restart_io_read(rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",278)

    ! check if restart file is compatible
    if ( rst_id /= p_boosted_diag_rst_id ) then
        write(err_buf__,*) 'Corrupted restart file, or restart file from incompatible binary (emf)';call err__("boost/os-boosted-diag.f03",282)
        call abort_program(p_err_rstrd)
    endif

    call restart_io_read(this%xmin1_prev, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",287)

    call restart_io_read(this%buf_start_idx, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",290)

    call restart_io_read(this%num_data_in_buf, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",293)

    select case( this%diag_type )

    case( p_boosted_4vector )

        if ( this%ndump_fac > 0 ) then

            call this % buffered_4vector % read_checkpoint( restart_handle )

            data_buf_nx = this%buffered_4vector%nx_
            data_buf_nx(1) = this%data_buf_size

            call reshape_array_vdf( this % buffered_4vector, data_buf_nx )

            call this % prev_4vector % read_checkpoint( restart_handle )

        endif

        call restart_io_read(this%num_par, restart_handle, ierr)
        call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",313)

        call restart_io_read(this%num_par_max, restart_handle, ierr)
        call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",316)

        if ( this%ndump_fac_raw > 0 ) then
            call alloc(this%buffered_part, (/ this%spec_dim + p_p_dim + 2, this%num_par_max /),"boost/os-boosted-diag.f03",319)
            call alloc(this%buffered_diag_idx, (/ this%num_par_max /),"boost/os-boosted-diag.f03",320)
            call alloc(this%num_par_diag, (/ this%grid_mgr%nt2 /),"boost/os-boosted-diag.f03",321)
            if ( this%add_tag ) then
                call alloc(this%buffered_tags, (/ 2, this%num_par_max /),"boost/os-boosted-diag.f03",323)
            endif
        endif

        if ( this%num_par > 0 ) then
            call restart_io_read(this%buffered_part(:,1:this%num_par), restart_handle, ierr)
            call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",329)

            call restart_io_read(this%buffered_diag_idx(1:this%num_par), restart_handle, ierr)
            call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",332)

            call restart_io_read(this%num_par_diag, restart_handle, ierr)
            call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",335)

            if ( this%add_tag ) then
                call restart_io_read(this%buffered_tags(:,1:this%num_par), restart_handle, ierr)
                call check_error(ierr,err_msg,p_err_rstrd,"boost/os-boosted-diag.f03",339)
            endif
        endif

    case( p_boosted_fields )

        call this % buffered_e % read_checkpoint( restart_handle )
        call this % buffered_b % read_checkpoint( restart_handle )

        data_buf_nx = this%buffered_e%nx_
        data_buf_nx(1) = this%data_buf_size

        call reshape_array_vdf( this % buffered_e, data_buf_nx )
        call reshape_array_vdf( this % buffered_b, data_buf_nx )

        call this % prev_e % read_checkpoint( restart_handle )
        call this % prev_b % read_checkpoint( restart_handle )

    end select

end subroutine restart_read_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine init_grid_mgr_boosted_diag( this, grid, g_space, dt1, tmin1, tmax1 )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    class( t_grid ), intent(in) :: grid
    type( t_space ), intent(in) :: g_space
    real(p_double), intent(in) :: dt1, tmin1, tmax1

    call this % grid_mgr % set_space1_params( grid, g_space, dt1, tmin1, tmax1 )
    call this % grid_mgr % init( this%ndump_fac, this%ndump_fac_raw )

end subroutine init_grid_mgr_boosted_diag
!-------------------------------------------------------------------------------------------------


! after reading the input file, allocate buffered_e/b. For slice reports, the data buffer
! will have great space left.
!-------------------------------------------------------------------------------------------------
subroutine init_boosted_diag( this, &
                    grid, g_space, dt1, tmin1, tmax1, &
                    ndump_fac_global, basepath, ref_vdf, &
                    gc_num, restart, restart_handle, init_grid_mgr_, &
                    spec_name_, spec_dim_, add_tag_ )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    class( t_grid ), intent(in) :: grid
    type( t_space ), intent(in) :: g_space
    real(p_double), intent(in) :: dt1, tmin1, tmax1
    integer, intent(in) :: ndump_fac_global
    character(len=*), intent(in) :: basepath
    class(t_vdf), intent(in) :: ref_vdf
    integer, dimension (2, p_x_dim), intent(in) :: gc_num
    logical, intent(in) :: restart
    type( t_restart_handle ), intent(in) :: restart_handle
    logical, intent(in), optional :: init_grid_mgr_
    character(len=*), intent(in), optional :: spec_name_
    integer, intent(in), optional :: spec_dim_
    logical, intent(in), optional :: add_tag_

    integer :: ndump_total, spec_dim
    integer, dimension(3) :: data_buf_nx
    logical :: init_grid_mgr, add_tag
    character(len=p_max_spname_len_) :: spec_name

    if( present( init_grid_mgr_ )) then
        init_grid_mgr = init_grid_mgr_
    else
        init_grid_mgr = .true.
    endif

    if( present( spec_name_ )) then
        spec_name = spec_name_
    else
        spec_name = ""
    endif

    if( present( spec_dim_ )) then
        spec_dim = spec_dim_
    else
        spec_dim = 0
    endif

    if( present( add_tag_ )) then
        add_tag = add_tag_
    else
        add_tag = .false.
    endif

    ndump_total = this%ndump_fac * ndump_fac_global

    ! this is an input because in quasi-3D we only call it for one of the 2m+1 diags
    if( init_grid_mgr ) then
        call this % init_grid_mgr( grid, g_space, dt1, tmin1, tmax1 )

        ! calculate the size to be used for the data buffer
        ! this is the most data that a local node would ever need to store
        call this % compute_data_buf_size( ndump_total )
    endif

    ! allocate the appropriate vdfs -- needs to happen before either restart
    ! or initialization later in this function
    call this % allocate_objs()

    ! non-restart initialization
    this%basepath = basepath
    this%spec_name = spec_name
    this%spec_dim = spec_dim
    this%add_tag = add_tag
    this%x_dim_ = ref_vdf%x_dim_
    call this % init_components_to_buffer()

    ! the rest of the initialization doesn't apply if we're restarting
    if( restart ) then

        call this % restart_read( restart_handle )

    else

        this%xmin1_prev = this%grid_mgr%xmin1

        ! init vars not passed in
        this % num_data_in_buf = 0
        this % buf_start_idx = 1

        data_buf_nx(1) = this%data_buf_size
        data_buf_nx(2:3) = ref_vdf%nx_(2:3)

        select case( this%diag_type )

        case( p_boosted_4vector )

            if ( this%ndump_fac > 0 ) then
                call this%prev_4vector %new( ref_vdf%x_dim_, 4, ref_vdf%nx_, gc_num, ref_vdf%dx_, .true. )
                call this%buffered_4vector%new( ref_vdf%x_dim_, 4, data_buf_nx, gc_num, ref_vdf%dx_, .true. )
            endif

            this%num_par = 0
            this%num_par_max = p_raw_buf_block

            if ( this%ndump_fac_raw > 0 ) then
                ! Include, positions, momenta, charge and energy
                call alloc(this%buffered_part, (/ spec_dim + p_p_dim + 2, p_raw_buf_block /),"boost/os-boosted-diag.f03",491)
                call alloc(this%buffered_diag_idx, (/ p_raw_buf_block /),"boost/os-boosted-diag.f03",492)
                call alloc(this%num_par_diag, (/ this%grid_mgr%nt2 /),"boost/os-boosted-diag.f03",493)
                this%num_par_diag = 0
                if ( this%add_tag ) then
                    call alloc(this%buffered_tags, (/ 2, p_raw_buf_block /),"boost/os-boosted-diag.f03",496)
                endif

                call this % init_raw_func()
            endif

        case( p_boosted_fields )

            call this%prev_e%new( ref_vdf%x_dim_, 3, ref_vdf%nx_, gc_num, ref_vdf%dx_, .true. )
            call this%prev_b%new( ref_vdf%x_dim_, 3, ref_vdf%nx_, gc_num, ref_vdf%dx_, .true. )

            call this%buffered_e%new( ref_vdf%x_dim_, 3, data_buf_nx, gc_num, ref_vdf%dx_, .true. )
            call this%buffered_b%new( ref_vdf%x_dim_, 3, data_buf_nx, gc_num, ref_vdf%dx_, .true. )

        case default
            if (mpi_node()==0) print *,"ERROR in os-boosted-diag: unrecognized boosted diagnostic type request.d"
            stop

        end select

    endif

end subroutine init_boosted_diag
!-------------------------------------------------------------------------------------------------


!max t2 data num between t1 and t1+ndump_total*dt1
!-------------------------------------------------------------------------------------------------
subroutine compute_data_buf_size_boosted_diag( this, ndump_total )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    integer, intent(in) :: ndump_total

    integer :: i1, tmp_buf_size, start_idx, final_idx

    this%data_buf_size = 1

    do i1 = 1, ceiling( this%grid_mgr%nt1 * 1.0 / ndump_total )

        start_idx = ( i1 - 1 ) * ndump_total + 1
        final_idx = min( i1 * ndump_total, this%grid_mgr%nt1 )

        tmp_buf_size = sum( this%grid_mgr%num_t2_points_at_t1(start_idx : final_idx))
        this%data_buf_size = max( this%data_buf_size, tmp_buf_size )
    enddo

end subroutine
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! allocate all objects except our grid manager
! the grid manager is allocated separately because in quasi-3D we only
! want to allocate grid manager which is shared by the different diags.
!-------------------------------------------------------------------------------------------------
subroutine allocate_objs_boosted_diag( this )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    select case( this%diag_type )

    case( p_boosted_4vector )

        call alloc(this%prev_4vector,"boost/os-boosted-diag.f03",564)
        call alloc(this%buffered_4vector,"boost/os-boosted-diag.f03",565)

    case( p_boosted_fields )

        call alloc(this%prev_e,"boost/os-boosted-diag.f03",569)
        call alloc(this%prev_b,"boost/os-boosted-diag.f03",570)
        call alloc(this%buffered_e,"boost/os-boosted-diag.f03",571)
        call alloc(this%buffered_b,"boost/os-boosted-diag.f03",572)

    end select

end subroutine allocate_objs_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine cleanup_boosted_diag( this )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    if( associated(this%grid_mgr) ) then
        call this%grid_mgr%cleanup()
        deallocate( this%grid_mgr )
    endif

    select case( this%diag_type )
    case( p_boosted_4vector )

        if ( associated(this%prev_4vector) ) then
            call this%prev_4vector%cleanup()
            call this%buffered_4vector%cleanup()
        endif

        call freemem(this%prev_4vector,"boost/os-boosted-diag.f03",601)
        call freemem(this%buffered_4vector,"boost/os-boosted-diag.f03",602)

        call freemem(this%buffered_part,"boost/os-boosted-diag.f03",604)
        call freemem(this%buffered_diag_idx,"boost/os-boosted-diag.f03",605)
        call freemem(this%num_par_diag,"boost/os-boosted-diag.f03",606)
        call freemem(this%buffered_tags,"boost/os-boosted-diag.f03",607)
        call freemem(this%idx,"boost/os-boosted-diag.f03",608)

    case( p_boosted_fields )

        if ( associated(this%prev_e) ) then
            call this%prev_e%cleanup()
            call this%prev_b%cleanup()
            call this%buffered_e%cleanup()
            call this%buffered_b%cleanup()
        endif

        call freemem(this%prev_e,"boost/os-boosted-diag.f03",619)
        call freemem(this%prev_b,"boost/os-boosted-diag.f03",620)
        call freemem(this%buffered_e,"boost/os-boosted-diag.f03",621)
        call freemem(this%buffered_b,"boost/os-boosted-diag.f03",622)

    end select

end subroutine cleanup_boosted_diag
!-------------------------------------------------------------------------------------------------


! read parameters, parse reports
!-------------------------------------------------------------------------------------------------
subroutine read_input_boosted_diag( this, input_file, quant_list, sim_gamma, diag_type, &
                                    if_use_boosted_diag, if_use_boosted_raw, x_dim )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    character(len=*), dimension(:), intent(in) :: quant_list
    real(p_double), intent(in) :: sim_gamma
    integer, intent(in) :: diag_type
    integer, intent(in) :: x_dim !
    logical, intent(out) :: if_use_boosted_diag, if_use_boosted_raw

    integer :: ndump_fac, ndump_fac_raw
    real(p_double) :: boost_signed_gamma
    real(p_double) :: boost_tmin, boost_tmax, boost_dt
    real(p_double) :: boost_xmin, boost_xmax
    integer :: boost_nx
    logical :: boost_if_move
    real(p_double) :: boost_move_vel
    real(p_k_part) :: raw_gamma_limit, raw_fraction
    character(len = p_max_expr_len) :: raw_math_expr
    logical :: raw_if_pos_ref_box

    integer :: ierr

    character( len = p_max_reports_len ), dimension( 12 ) :: boosted_reports

    namelist /nl_boosted_diag/ ndump_fac, ndump_fac_raw, boost_signed_gamma, &
                           boost_tmin, boost_tmax, boost_dt, &
                           boost_xmin, boost_xmax, boost_nx, &
                           boost_if_move, boost_move_vel, &
                           boosted_reports, raw_gamma_limit, raw_fraction, &
                           raw_math_expr, raw_if_pos_ref_box

    ndump_fac = 0
    ndump_fac_raw = 0
    boost_signed_gamma = -sim_gamma
    boost_tmin = 0.
    boost_tmax = 0.
    boost_dt = 0.
    boost_xmin = 0.
    boost_xmax = 0.
    boost_nx = 0
    boost_if_move = .true.
    boost_move_vel = 1.
    boosted_reports = "-"
    raw_gamma_limit = 0
    raw_fraction = 1.0_p_k_part
    raw_math_expr = ''
    ! the below option currently does nothing
    raw_if_pos_ref_box = .false.

    this%diag_type = diag_type

    this%x_dim_ = x_dim
    ! disable boosted diags by default
    this%ndump_fac = ndump_fac
    this%ndump_fac_raw = ndump_fac_raw

    ! return and say we didn't find a boosted diag if there's no boosted_diag namelist
    if_use_boosted_diag = .false.
    if_use_boosted_raw = .false.
    if (trim(input_file % get_section_name()) /= 'nl_boosted_diag') return

    ! Get namelist text from input file
    call get_namelist( input_file, "nl_boosted_diag", ierr )

    if (ierr == 0) then
        read (input_file%nml_text, nml = nl_boosted_diag, iostat = ierr)
        if (ierr /= 0) then
            if (mpi_node()==0) print *,"Error reading boosted_diag parameters"
            if (mpi_node()==0) print *,"aborting..."
            stop
        endif
    else
        if (ierr < 0) then
            if (mpi_node()==0) print *,"Error reading boosted_diag parameters"
            if (mpi_node()==0) print *,"aborting..."
            stop
        else
            if (disp_out(input_file)) then
                if (mpi_node()==0) print *," - no diagnostics specified"
            endif
        endif
    endif

    this%ndump_fac = ndump_fac
    if_use_boosted_diag = ( ndump_fac > 0 )

    if ( diag_type == p_boosted_4vector ) then
        this%ndump_fac_raw = ndump_fac_raw
        if_use_boosted_raw = ( ndump_fac_raw > 0 )
    else
        if_use_boosted_raw = .false.
    endif

    this%raw_gamma_limit = raw_gamma_limit
    this%raw_fraction = raw_fraction
    this%raw_math_expr = trim(raw_math_expr)
    this%raw_if_pos_ref_box = raw_if_pos_ref_box

    if (.not. (if_use_boosted_diag .or. if_use_boosted_raw)) return

    if( .not. associated( this % grid_mgr ) ) then
        allocate( this%grid_mgr )
    else
        if (mpi_node()==0) print *,'ERROR in os-boosted-diag: grid manager was already allocated in read_input.'
        if (mpi_node()==0) print *,'This should never happen.'
        stop
    endif


    call this % parse_reports_string( boosted_reports, quant_list )

    call this % grid_mgr % set_space2_params( boost_signed_gamma, &
                                         boost_tmin, boost_tmax, boost_dt, &
                                         boost_xmin, boost_xmax, boost_nx, &
                                         boost_if_move, boost_move_vel )

end subroutine read_input_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! boosted_reports_str: list of reports to track. only these quantities will be
! buffered. e.g. ( "e1", "b3" )
! quant_list: list of ALL options for components for this diagnostic. example:
! ( "e1", "e2", "e3", "b1", "b2", "b3" )
!-------------------------------------------------------------------------------------------------
subroutine parse_reports_string_boosted_diag( this, boosted_reports_requested, quant_list )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: boosted_reports_requested
    character(len=*), dimension(:), intent(in) :: quant_list

    integer :: component !index of the quant list
    character(len=64) :: requested_report

    integer :: i, j

    character( len = 64) :: substr !report_items: boost_e1, boost_e2,...
    integer, dimension( 16 ) :: splitIdx
    integer :: nsplit, direction
    character(len=64) :: name
    integer :: gipos !slice parameter
    integer :: ierr

    this%num_diag_components = 0
    this%num_diag_components_slice = 0
    ! loop thru all requested quants
    do i = 1, size( boosted_reports_requested )

        requested_report = boosted_reports_requested(i)

        if ( trim(requested_report ) == '-' ) exit

        !scanall return splitIdx: the index of each comma, nsplit: how many comma
        call scanall( requested_report, ',', splitIdx, nsplit )
        if ( nsplit == 0 ) splitIdx(1) = len(requested_report)+1
        if ( splitIdx(1) < 2 ) then
            ! invalid quantity
            ierr = -1
            return
        endif
        substr = adjustl(requested_report(1 : splitIdx(1) - 1))
        component = -1

        ! loop through all available components
        ! make sure we match one, and if so store the component
        do j = 1, size(quant_list)

            if ( trim(substr) == trim(quant_list(j)) ) then
                component = j
                exit
            endif

        enddo

        if ( component == -1 ) then
            if (mpi_node()==0) print *,"ERROR in os-boosted-diag: invalid boosted quantity requested: ", substr
            stop
        endif

        name = quant_list(component)
        ! pfull
        if (nsplit == 0) then
            this%num_diag_components = this%num_diag_components + 1
            this%diag_components(this%num_diag_components) = component
            this%component_names(this%num_diag_components) = name
        ! pslice
        else
            if ( nsplit == 1 ) splitIdx(2) = len(requested_report)+1
            substr = adjustl(requested_report(splitIdx(1)+1:splitIdx(2)-1))
            ! name = trim(name) // '-' // trim(substr)

            select case (trim(substr))
            case( 'savg' )
                ! invalid type
                print *, 'savg has not been written"', substr, '"'
                ierr = -2
                stop
            case( 'senv' )
                ! invalid type
                print *, 'senv has not been written"', substr, '"'
                ierr = -2
                stop
            case( 'line' )
                ! invalid type
                print *, 'line has not been written"', substr, '"'
                ierr = -2
                stop
            case( 'slice' )
                this%num_diag_components_slice = this%num_diag_components_slice + 1
                if ( this%num_diag_components_slice > p_slice_max ) then
                    ! we only have less than p_slice_max slices
                    print *, 'the number of slices should be less than "', p_slice_max, '"'
                    ierr = -3
                    stop
                endif
                this%diag_components_slice(this%num_diag_components_slice) = component
                if ( this%x_dim_ /= 3 ) then
                    ! slices are only available in 3D
                    ierr = -3
                    stop
                endif
                if ( nsplit /= 3 ) then
                    ! wrong parameter count
                    ierr = -2
                    stop
                endif
                ! default values
                gipos = -1
                direction = -1
                ! direction
                substr = adjustl(requested_report(splitIdx(2)+1:splitIdx(3) - 1))
                select case (trim(substr))
                    case ('x1')
                    direction = 1
                    print *, 'x1 direction does not need boost diag '
                    ierr = -5
                    stop
                    case ('x2')
                    direction = 2
                    case ('x3')
                    direction = 3
                    case default
                    ! invalid direction
                    ierr = -5
                    stop
                end select
                if ( direction > this%x_dim_ ) then
                    ! invalid direction
                    ierr = -5
                    stop
                endif
                this%direction(this%num_diag_components_slice) = direction

                    !strtoint: read the slice position
                if ( nsplit == 3 ) then
                    gipos = strtoint( requested_report( splitIdx(3)+1: ), ierr )
                    if ( ierr /= 0 ) then
                    ! invalid position
                    ierr = -6
                    stop
                    endif
                else
                    print *, 'Invalid slice input, we only have one parameter for slice now ', trim(substr)
                    ierr = -6
                    stop
                    ! gipos(1) = strtoint( input( splitIdx(3)+1:splitIdx(4)-1 ), ierr )
                    ! if ( ierr /= 0 ) then
                    ! ! invalid position
                    ! ierr = -6
                    ! return
                    ! endif
                    ! gipos(2) = strtoint( input( splitIdx(4)+1: ), ierr )
                    ! if ( ierr /= 0 ) then
                    ! ! invalid position
                    ! ierr = -6
                    ! return
                    ! endif
                endif
                this%gipos(this%num_diag_components_slice) = gipos
                this%component_names_slice(this%num_diag_components_slice) = name
            case default
                ! invalid report type
                print *, 'Invalid report type, ', trim(substr)
                ierr = -4
                stop
            end select
        endif
    enddo

    !this%component_names(1:size(quant_list)) = quant_list

    ! if (mpi_node()==0) print *,'requested: ', boosted_reports_requested
    ! if (mpi_node()==0) print *,'quant_list: ', quant_list

end subroutine parse_reports_string_boosted_diag
!-------------------------------------------------------------------------------------------------


! used in os-spec-diagnostics.f03
!-------------------------------------------------------------------------------------------------
subroutine init_components_to_buffer_boosted_diag( this )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    integer :: i

    this%num_diag_components_to_buffer = 0

    ! compute necessary diag components here
    select case( this%diag_type )

    case( p_boosted_4vector)

        i = 1

        if( any( this%diag_components(1:this%num_diag_components) .eq. 1 ) &
            .or. any( this%diag_components(1:this%num_diag_components) .eq. 2 ) .or. &
            any(this%diag_components_slice(1:this%num_diag_components_slice) .eq. 1) &
            .or. any( this%diag_components_slice(1:this%num_diag_components_slice) .eq. 2) ) then

            this%diag_components_to_buffer(i) = 1
            i = i + 1

            this%diag_components_to_buffer(i) = 2
            i = i + 1

            this%num_diag_components_to_buffer = this%num_diag_components_to_buffer + 2

        endif

        if( any( this%diag_components(1:this%num_diag_components) .eq. 3 ) &
            .or. any(this%diag_components_slice(1:this%num_diag_components_slice) .eq. 3)) then
            this%diag_components_to_buffer(i) = 3
            this%num_diag_components_to_buffer = this%num_diag_components_to_buffer + 1
            i = i + 1
        endif

        if( any( this%diag_components(1:this%num_diag_components) .eq. 4 ) &
            .or. any(this%diag_components_slice(1:this%num_diag_components_slice) .eq. 4)) then
            this%diag_components_to_buffer(i) = 4
            this%num_diag_components_to_buffer = this%num_diag_components_to_buffer + 1
            i = i + 1
        endif

    end select

end subroutine init_components_to_buffer_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine get_dir_fname_label_boosted_diag( this, component, t2_idx, dir, fname, label )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    integer, intent(in) :: t2_idx
    integer, intent(in) :: component
    character(len=1024 ), intent(out) :: dir
    character(len=128 ), intent(out) :: fname
    character(len=128 ), intent(out) :: label

    character(len=128) :: name
    character(len=p_max_spname_len_+1) :: spec_name

    if ( this%diag_type == p_boosted_4vector ) then
        spec_name = replace_blanks(trim(this%spec_name)) // "-"
    else
        spec_name = ""
    endif

    name = this%component_names( component )

    dir = trim(this%basepath) // trim(name) // p_dir_sep

    label = trim(name)

    fname = trim( name ) // "-" // trim(spec_name) // idx_string( t2_idx, p_time_length )

end subroutine get_dir_fname_label_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine get_dir_fname_label_boosted_diag_slice( this, component, t2_idx, dir, fname, label )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    integer, intent(in) :: t2_idx
    integer, intent(in) :: component
    character(len=1024 ), intent(out) :: dir
    character(len=128 ), intent(out) :: fname
    character(len=128 ), intent(out) :: label

    character(len=128) :: name
    character(len=p_max_spname_len_+1) :: spec_name

    if ( this%diag_type == p_boosted_4vector ) then
        spec_name = trim(this%spec_name) // "-"
    else
        spec_name = ""
    endif

    name = this%component_names_slice( component )

    dir = trim(this%basepath) // trim(name) // "-slice" // p_dir_sep

    select case(this%direction(component))
    case(2)
    label = trim(name) // "_x2_slice"

    fname = trim( name ) // "-slice-" // trim(spec_name) // "x2-" // &
            trim(tostring_int(this%gipos(component))) // "-" // idx_string( t2_idx, p_time_length )
    case(3)
    label = trim(name) // "_x3_slice"

    fname = trim( name ) // "-slice-" // trim(spec_name) // "x3-" // &
            trim(tostring_int(this%gipos(component))) // "-" // idx_string( t2_idx, p_time_length )
    end select
end subroutine get_dir_fname_label_boosted_diag_slice
!-------------------------------------------------------------------------------------------------




!-------------------------------------------------------------------------------------------------
subroutine create_files_boosted_diag( this, grid )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    class( t_grid ), intent(in) :: grid

    integer :: i

    if ( this%ndump_fac > 0 ) then

        do i = 1, this%num_diag_components
            call this % create_files_component( this%diag_components(i), grid, i )
        enddo

        do i = 1, this%num_diag_components_slice
            call this % create_files_component_slice( this%diag_components_slice(i), grid, i )
        enddo

    endif

    ! The creation of raw files is handled in the write_boosted_raw_ function, which calls
    ! the create_files_raw function.

end subroutine create_files_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine create_files_component_boosted_diag( this, component, grid, component_idx )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    integer, intent(in) :: component
    integer, intent(in) :: component_idx !index of the component
    class( t_grid ), intent(in) :: grid

    character(len=1024) :: dir
    character(len=128 ) :: filename
    character(len=128 ) :: label

    integer :: i2, j
    real(p_double) :: t2, xmin_current2, xmax_current2
    integer, dimension(3) :: dataset_dims

    class( t_diag_file ), allocatable :: diagfile
    type( t_diag_dataset ) :: dset
    character(len=12), dimension(3) :: xnames, xlabels, xunits
    character(len=8) :: label_
    integer, parameter :: izero = ichar('0') !zero in ascii code

    xnames = (/'x1', 'x2', 'x3'/)
    xlabels = (/'x_1', 'x_2', 'x_3'/)
    xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    !number of diagfile: nt2
    do i2 = 0, this%grid_mgr%nt2 - 1
        ! get the time and spatial bounds for this slice
        t2 = this%grid_mgr%tmin2 + i2 * this%grid_mgr%dt2

        if( this%grid_mgr%if_move2 ) then
            xmin_current2 = this%grid_mgr%xmin2 + t2
        else
            xmin_current2 = this%grid_mgr%xmin2
        endif

        xmax_current2 = xmin_current2 + ( this%grid_mgr%xmax2 - this%grid_mgr%xmin2 )

        ! this allocates the diagfile -- must be deallocated later
        call create_diag_file( diagfile )

        call this % get_dir_fname_label( component_idx, i2, dir, filename, label )

        diagfile % ftype = p_diag_grid

        diagfile % filepath = dir
        diagfile % filename = filename
        diagfile % name = label

        ! Iteration info
        diagfile%iter%n = i2
        diagfile%iter%t = t2

        ! Grid info
        diagfile%grid%ndims = this % x_dim_


        select case( this%diag_type )

        case( p_boosted_4vector )

            diagfile%grid%units = 'n_0'

            select case( component )

            case(1)
                label_ = '\rho'

                if ( p_x_dim == 1 ) then
                    diagfile%grid%units = 'e \omega_p / c'
                else
                    diagfile%grid%units = 'e \omega_p^'//char(izero+p_x_dim)// &
                                           '/ c^'//char(izero+p_x_dim)
                endif

            case(2)
                label_ = 'j_1'
            case(3)
                label_ = 'j_2'
            case(4)
                label_ = 'j_3'
            end select

            ! diagfile%units = 'm_e c \omega_p e^{-1}'


        case( p_boosted_fields )

            diagfile%grid%units = 'm_e c \omega_p e^{-1}'

            select case( component )
            case(1)
                label_ = 'E_1'
            case(2)
                label_ = 'E_2'
            case(3)
                label_ = 'E_3'
            case(4)
                label_ = 'B_1'
            case(5)
                label_ = 'B_2'
            case(6)
                label_ = 'B_3'
            end select

        end select

        diagFile%grid%label = trim( label_ )
        diagfile%grid%name = diagfile%name

        ! handle i=1 separately
        diagfile%grid%count(1) = this%grid_mgr%nx2
        diagfile%grid%axis(1)%type = diag_axis_linear
        diagfile%grid%axis(1)%min = xmin_current2
        diagfile%grid%axis(1)%max = xmax_current2

        diagfile%grid%axis(1)%name = xnames(1)
        diagfile%grid%axis(1)%label = xlabels(1)
        diagfile%grid%axis(1)%units = xunits(1)

        ! handle transverse directions
        do j = 2, this % x_dim_
            diagfile%grid%count(j) = this%grid_mgr%gnxperp(j-1)
            diagfile%grid%axis(j)%type = diag_axis_linear
            diagfile%grid%axis(j)%min = this%grid_mgr%gxperp_min(j-1)
            diagfile%grid%axis(j)%max = this%grid_mgr%gxperp_max(j-1)

            diagfile%grid%axis(j)%name = xnames(j)
            diagfile%grid%axis(j)%label = xlabels(j)
            diagfile%grid%axis(j)%units = xunits(j)
        enddo

        call diagFile % open( p_diag_create, grid%io%comm )

        dataset_dims(1) = this%grid_mgr%nx2
        dataset_dims(2:3) = this%grid_mgr%gnxperp(1:2)

        ! don't chunk right now, use contiguous dataset
        call diagFile % start_cdset( label, this%x_dim_, dataset_dims, &
                                     diag_type_real( p_single ), dset )

        call diagFile % close()

        deallocate( diagFile )
    enddo

end subroutine create_files_component_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine create_files_component_boosted_diag_slice( this, component, grid, component_idx )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    integer, intent(in) :: component
    integer, intent(in) :: component_idx !index of the component
    class( t_grid ), intent(in) :: grid

    character(len=1024) :: dir
    character(len=128 ) :: filename
    character(len=128 ) :: label

    integer :: i2
    real(p_double) :: t2, xmin_current2, xmax_current2
    integer, dimension(2) :: dataset_dims_slice

    class( t_diag_file ), allocatable :: diagfile
    type( t_diag_dataset ) :: dset
    character(len=12), dimension(3) :: xnames, xlabels, xunits
    character(len=8) :: label_
    integer, parameter :: izero = ichar('0')

    xnames = (/'x1', 'x2', 'x3'/)
    xlabels = (/'x_1', 'x_2', 'x_3'/)
    xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    ! number of slice diagfile: nt2
    do i2 = 0, this%grid_mgr%nt2 - 1
        ! get the time and spatial bounds for this slice
        t2 = this%grid_mgr%tmin2 + i2 * this%grid_mgr%dt2

        if( this%grid_mgr%if_move2 ) then
            xmin_current2 = this%grid_mgr%xmin2 + t2
        else
            xmin_current2 = this%grid_mgr%xmin2
        endif

        xmax_current2 = xmin_current2 + ( this%grid_mgr%xmax2 - this%grid_mgr%xmin2 )

        ! this allocates the diagfile -- must be deallocated later
        call create_diag_file( diagfile )

        call this % get_dir_fname_label_slice( component_idx, i2, dir, filename, label )

        diagfile % ftype = p_diag_grid

        diagfile % filepath = dir
        diagfile % filename = filename
        diagfile % name = label

        ! Iteration info
        diagfile%iter%n = i2
        diagfile%iter%t = t2

        ! Grid info
        diagfile%grid%ndims = 2 ! slice ndim=2

        select case( this%diag_type )

        case( p_boosted_4vector )

            diagfile%grid%units = 'n_0'

            select case( component )

            case(1)
                label_ = '\rho'

                if ( p_x_dim == 1 ) then
                    diagfile%grid%units = 'e \omega_p / c'
                else
                    diagfile%grid%units = 'e \omega_p^'//char(izero+p_x_dim)// &
                                           '/ c^'//char(izero+p_x_dim)
                endif

            case(2)
                label_ = 'j_1'
            case(3)
                label_ = 'j_2'
            case(4)
                label_ = 'j_3'
            end select

            ! diagfile%units = 'm_e c \omega_p e^{-1}'

        case( p_boosted_fields )

            diagfile%grid%units = 'm_e c \omega_p e^{-1}'

            select case( component )
            case(1)
                label_ = 'E_1'
            case(2)
                label_ = 'E_2'
            case(3)
                label_ = 'E_3'
            case(4)
                label_ = 'B_1'
            case(5)
                label_ = 'B_2'
            case(6)
                label_ = 'B_3'
            end select

        end select

        diagFile%grid%label = trim( label_ )
        diagfile%grid%name = diagfile%name

        ! handle i=1 separately
        diagfile%grid%count(1) = this%grid_mgr%nx2
        diagfile%grid%axis(1)%type = diag_axis_linear
        diagfile%grid%axis(1)%min = xmin_current2
        diagfile%grid%axis(1)%max = xmax_current2

        diagfile%grid%axis(1)%name = xnames(1)
        diagfile%grid%axis(1)%label = xlabels(1)
        diagfile%grid%axis(1)%units = xunits(1)

        select case(this%direction(component_idx))
        case(2)
            diagfile%grid%count(2) = this%grid_mgr%gnxperp(2)
            diagfile%grid%axis(2)%type = diag_axis_linear
            diagfile%grid%axis(2)%min = this%grid_mgr%gxperp_min(2)
            diagfile%grid%axis(2)%max = this%grid_mgr%gxperp_max(2)

            diagfile%grid%axis(2)%name = xnames(3)
            diagfile%grid%axis(2)%label = xlabels(3)
            diagfile%grid%axis(2)%units = xunits(3)
        case(3)
            diagfile%grid%count(2) = this%grid_mgr%gnxperp(1)
            diagfile%grid%axis(2)%type = diag_axis_linear
            diagfile%grid%axis(2)%min = this%grid_mgr%gxperp_min(1)
            diagfile%grid%axis(2)%max = this%grid_mgr%gxperp_max(1)

            diagfile%grid%axis(2)%name = xnames(2)
            diagfile%grid%axis(2)%label = xlabels(2)
            diagfile%grid%axis(2)%units = xunits(2)
        end select

        call diagFile % open( p_diag_create, grid%io%comm )

        select case(this%direction(component_idx))
        case(2)
            dataset_dims_slice(1) = this%grid_mgr%nx2
            dataset_dims_slice(2) = this%grid_mgr%gnxperp(2)
        case(3)
            dataset_dims_slice(1) = this%grid_mgr%nx2
            dataset_dims_slice(2) = this%grid_mgr%gnxperp(1)
        end select

        ! don't chunk right now, use contiguous dataset
        call diagFile % start_cdset( label, 2, dataset_dims_slice, &
                                     diag_type_real( p_single ), dset )

        call diagFile % close()

        deallocate( diagFile )
    enddo

end subroutine create_files_component_boosted_diag_slice
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
subroutine create_files_raw_boosted_diag( this, grid )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    class( t_grid ), intent(in) :: grid

    class( t_diag_file ), allocatable :: diagfile
    type( t_diag_dataset ), dimension(this%spec_dim + p_p_dim + 3) :: datasets

    integer :: i, j, n_x_dim, i2, data_type

    if (mpi_node()==0) print *,'Creating RAW boosted files...'

    n_x_dim = this%spec_dim

    !number of diagfile: nt2
    do i2 = 0, this%grid_mgr%nt2 - 1

        ! Prepare path and file name
        call create_diag_file( diagFile )

        call this % populate_raw_diagFile( diagFile, i2 )

        ! create the file
        call diagFile % open( p_diag_create, grid%io%comm )

        ! create datasets
        select case (p_diag_prec)
        case ( p_single )
            data_type = p_diag_float32
        case ( p_double )
            data_type = p_diag_float64
        end select

        if ( diagFile%particles%has_tags ) then
            j = diagFile % particles % nquants - 1
        else
            j = diagFile % particles % nquants
        endif

        ! particle quantities
        do i = 1, j
            call diagFile % start_cdset( diagFile % particles % quants(i), 1, [0], &
                data_type, datasets(i), chunk_size=[p_raw_chunk_size] )
        enddo

        ! particle tags
        if ( diagFile%particles%has_tags ) then
            i = diagFile % particles % nquants
            call diagFile % start_cdset( 'tag', 2, [2,0], p_diag_int32, &
                datasets(i), chunk_size=[2,p_raw_chunk_size] )
        endif

        do i = 1, diagFile % particles % nquants
            call diagFile % end_cdset( datasets(i) )
        enddo

        call diagFile % close()

        deallocate( diagFile )

    enddo

    if (mpi_node()==0) print *,'Done.'

end subroutine create_files_raw_boosted_diag
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
subroutine populate_raw_diagFile_boosted_diag( this, diagFile, i2 )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    class(t_diag_file), intent(inout) :: diagFile
    integer, intent(in) :: i2 ! Must be 0-indexed

    integer :: i, j, n_x_dim
    real(p_double) :: t2

    n_x_dim = this%spec_dim

    ! get the time and spatial bounds for this slice
    t2 = this%grid_mgr%tmin2 + i2 * this%grid_mgr%dt2

    ! Prepare path and file name
    diagFile % ftype = p_diag_particles

    diagFile%filepath = trim(path_mass)//'RAW_BOOST'//p_dir_sep//replace_blanks(trim(this%spec_name))//p_dir_sep
    diagFile%filename = get_filename( i2, '' // 'RAW_BOOST'  // '-' // replace_blanks(trim(this%spec_name)) )

    diagFile%name = this%spec_name

    diagFile%iter%n = i2
    diagFile%iter%t = t2
    diagFile%iter%time_units = '1 / \omega_p'

    diagFile%particles%name = this%spec_name

    ! Species object doesn't have a label (yet) so just set it to the name
    diagFile%particles%label = this%spec_name

    diagFile%particles%np = 0
    diagFile%particles%has_tags = this%add_tag

    j = 0

    ! positions
    do i = 1, n_x_dim
        j = j+1
        diagFile%particles%quants(j) = 'x'//(char(iachar('0')+i))
        diagFile%particles%qlabels(j) = 'x_'//(char(iachar('0')+i))
        diagFile%particles%qunits(j) = 'c/\omega_p'
        diagFile%particles%offset_t(j) = 0.0_p_double
    enddo

    ! momenta
    do i = 1, p_p_dim
        j = j+1
        diagFile%particles%quants(j) = 'p'//(char(iachar('0')+i))
        diagFile%particles%qlabels(j) = 'p_'//(char(iachar('0')+i))
        diagFile%particles%qunits(j) = 'm_e c'
        diagFile%particles%offset_t(j) = -0.5_p_double
    enddo

    ! charge
    j = j+1
    diagFile%particles%quants(j) = 'q'
    diagFile%particles%qlabels(j) = 'q'
    diagFile%particles%qunits(j) = 'e'
    diagFile%particles%offset_t(j) = 0.0_p_double

    ! energy
    j = j+1
    diagFile%particles%quants(j) = 'ene'
    diagFile%particles%qlabels(j) = 'Ene'
    diagFile%particles%qunits(j) = 'm_e c^2'
    diagFile%particles%offset_t(j) = -0.5_p_double

    if ( diagFile%particles%has_tags ) then
        j = j+1
        diagFile%particles%quants(j) = 'tag'
        diagFile%particles%qlabels(j) = 'Tag'
        diagFile%particles%qunits(j) = ''
        diagFile%particles%offset_t(j) = 0.0_p_double
    endif

    diagFile%particles%nquants = j

end subroutine populate_raw_diagFile_boosted_diag
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
subroutine init_raw_func_boosted_diag( this )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    integer :: ierr

    ! set variable list for raw_math_expr
    if (this%raw_math_expr /= '') then

        select case (p_x_dim)
        case (1)
            call setup(this%raw_func, trim(this%raw_math_expr), &
                (/'x1', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)

        case (2)
            call setup(this%raw_func, trim(this%raw_math_expr), &
                (/'x1', 'x2', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)

        case (3)
            call setup(this%raw_func, trim(this%raw_math_expr), &
                (/'x1', 'x2', 'x3', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)

        end select

        ! check if function compiled ok
        if (ierr /= 0) then
            if (mpi_node() == 0) then
                write(0,*) "(*error*) Unable to compile 'raw_math_expr' function:"
                write(0,*) "(*error*) '", trim(this%raw_math_expr), "'"
                write(0,*) "(*error*) aborting..."
            endif
            call abort_program(-1)
        endif
    endif

end subroutine init_raw_func_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine report_boosted_diag( this, tstep, grid )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    type( t_time_step ), intent(in) :: tstep
    class( t_grid ), intent(in) :: grid

    integer :: i, n1

    n1 = n( tstep )

    if( n1 .eq. 0 ) then

        if (mpi_node()==0) print *,'Creating boosted files...'
        call this % create_files( grid )
        if (mpi_node()==0) print *,'Done.'

    else

        do i = 1, this%num_diag_components
            !whether to report, this ndump corresponds to simulation time t1, not t2
            if( .not. test_if_report( tstep, this%ndump_fac ) ) return
            call this % report_component( this%diag_components(i), grid, i )
        enddo

        do i = 1, this%num_diag_components_slice
            !whether to report, this ndump corresponds to simulation time t1, not t2
            if( .not. test_if_report( tstep, this%ndump_fac ) ) return
            call this % report_component_slice( this%diag_components_slice(i), grid, i )
        enddo

        ! adjust position of all buffers
        this%buf_start_idx = this%buf_start_idx + this%num_data_in_buf
        this%num_data_in_buf = 0

    endif

end subroutine report_boosted_diag
!-------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function SliceHasLocalData( gipos, direction, nodeIdx )

  implicit none

  integer, intent(in) :: gipos
  integer, intent(in) :: direction
  integer, dimension(:, :), intent(in) :: nodeIdx

  logical :: SliceHasLocalData

  SliceHasLocalData = ( gipos >= nodeIdx( p_lower, direction ) ) .and. &
                      ( gipos <= nodeIdx( p_upper, direction ) )

end function SliceHasLocalData
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine SliceBounds( gipos, direction, nodeIdx, lperp )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  integer, intent(in) :: gipos
  integer, intent(in) :: direction

  integer, dimension(:,:), intent(in) :: nodeIdx
  integer, intent(out) :: lperp
! integer, dimension(:,:), intent(out) :: gbounds ! we dont have gbounds now


  select case ( direction )
    ! we only have 2 and 3 now
    ! case(1)
    ! lperp = gipos - nodeIdx( p_lower, 1 ) + 1

    ! gbounds(1,p_lower) = nodeIdx( p_lower, 2 )
    ! gbounds(1,p_upper) = nodeIdx( p_upper, 2 )

    ! gbounds(2,p_lower) = nodeIdx( p_lower, 3 )
    ! gbounds(2,p_upper) = nodeIdx( p_upper, 3 )

    case(2)
    ! gbounds(1,p_lower) = nodeIdx( p_lower, 1 )
    ! gbounds(1,p_upper) = nodeIdx( p_upper, 1 )

      lperp = gipos - nodeIdx( p_lower, 2 ) + 1

    ! gbounds(2,p_lower) = nodeIdx( p_lower, 3 )
    ! gbounds(2,p_upper) = nodeIdx( p_upper, 3 )

    case(3)
    ! gbounds(1,p_lower) = nodeIdx( p_lower, 1 )
    ! gbounds(1,p_upper) = nodeIdx( p_upper, 1 )

    ! gbounds(2,p_lower) = nodeIdx( p_lower, 2 )
    ! gbounds(2,p_upper) = nodeIdx( p_upper, 2 )

      lperp = gipos - nodeIdx( p_lower, 3 ) + 1

  end select

end subroutine SliceBounds
!---------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! name: prefix for the data
! component: component of the vdf where the data is stored
! i indexes time, j indexes space
!-------------------------------------------------------------------------------------------------
subroutine report_component_boosted_diag( this, component, grid, component_idx)

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    class( t_grid ), intent(in) :: grid
    integer, intent(in) :: component, component_idx

    character(len=1024) :: dir
    character(len=128 ) :: filename
    character(len=128 ) :: label

    class( t_diag_file ), allocatable :: diagfile
    type( t_diag_dataset ) :: dset
    type( t_diag_chunk ) :: chunk

    ! type(t_vdf) :: write_buf
    type(t_vdf), pointer :: data_buf

    real(p_single), dimension(:), pointer :: write_buf_1d
    real(p_single), dimension(:,:), pointer :: write_buf_2d
    real(p_single), dimension(:,:,:), pointer :: write_buf_3d

    integer :: i2, j, idx, cur_i2, cur_j2, data_idx
    integer :: start_idx, num_data, cur_idx, min_j2, buff_idx
    integer :: cur_idx_within_databuffer

    ! i2 values obtained in the unsorted buffer
    integer, dimension(:), pointer :: i2_vals, num_data_at_i2

    ! integer, dimension(:,:), pointer :: buffer_indices_for_each_i2

    ! sorted indices of those i2 indices
    integer, dimension(:), pointer :: sorted_i2_indices

    ! set pointers to null
    write_buf_1d => null(); write_buf_2d => null(); write_buf_3d => null()
    i2_vals => null(); num_data_at_i2 => null(); sorted_i2_indices => null()

    select case( this%diag_type )

    case( p_boosted_fields )

        ! components 1 to 3 store e
        if( component <= 3 ) then
            data_buf => this%buffered_e
            data_idx = component

        ! components 4 to 6 store b
        else
           data_buf => this%buffered_b
           data_idx = component - 3
       endif

    case( p_boosted_4vector )

       data_buf => this%buffered_4vector
       data_idx = component

    end select

        start_idx = 0

        call alloc(num_data_at_i2, (/ this%grid_mgr%nt2 /),"boost/os-boosted-diag.f03",1796)
        num_data_at_i2(:) = 0

        ! first compute sorted_i2_indices
        if( this%num_data_in_buf > 0 ) then

            call alloc(i2_vals, (/ this%num_data_in_buf /),"boost/os-boosted-diag.f03",1802)
            call alloc(sorted_i2_indices, (/ this%num_data_in_buf /),"boost/os-boosted-diag.f03",1803)

            do j=1, this%num_data_in_buf
                idx = this%buf_start_idx + (j-1)
                i2 = this%grid_mgr%boosted_cell_coords(idx)%t2_idx
                i2_vals(j) = i2
                num_data_at_i2( i2 ) = num_data_at_i2( i2 ) + 1
            enddo

            call heapsort( i2_vals, sorted_i2_indices )

        endif

        call create_diag_file( diagfile )

        ! loop through all lab time slices to reconstruct
        do i2 = 1, this%grid_mgr%nt2

            num_data = num_data_at_i2(i2)
            min_j2 = 1 ! default value in case num_data is 0

            ! we still need to proceed if num_data == 0 because other processes
            ! may have data; in this case, we need to tell the mpi process
            ! doing the writing that we don't have any data, otherwise it
            ! will get stuck waiting for us. however, if we don't have any data
            ! then we don't need to do anything other than give this info to
            ! the writing process.

            if( num_data > 0 ) then

                ! allocate the appropriate dimension buffer. this will hold
                ! all the sorted data for a given i2 until we write it.
                ! the dimensions have to exactly match the ones we eventually
                ! write, which is why we need to keep reallocating this array.
                select case( this%x_dim_ )
                case(1)
                    call alloc(write_buf_1d, (/ num_data /),"boost/os-boosted-diag.f03",1839)
                    chunk % data = c_loc( write_buf_1d )
                case(2)
                    call alloc(write_buf_2d, (/ num_data, data_buf%nx_(2) /),"boost/os-boosted-diag.f03",1842)
                    chunk % data = c_loc( write_buf_2d )
                case(3)
                    call alloc(write_buf_3d, (/ num_data, data_buf%nx_(2), data_buf%nx_(3) /),"boost/os-boosted-diag.f03",1845)
                    chunk % data = c_loc( write_buf_3d )
                end select

            endif

            ! find the minimum j2 attained for all the available data (the heapsort
            ! doesn't respect the j2 ordering posessed by the data prior to sorting )
            do j = 1, num_data
                ! cur_idx_within_databuffer = sorted_i2_indices( start_idx + j )
                cur_idx = this%buf_start_idx + sorted_i2_indices( start_idx + j ) - 1
                cur_j2 = this%grid_mgr%boosted_cell_coords(cur_idx)%x2_idx

                if( j .eq. 1 ) then
                    min_j2 = cur_j2
                else
                    min_j2 = min( min_j2, cur_j2 )
                endif

            enddo

            ! count all the data in the buffer which matches this i2
            do j = 1, num_data

                cur_idx_within_databuffer = sorted_i2_indices( start_idx + j )
                cur_idx = this%buf_start_idx + sorted_i2_indices( start_idx + j ) - 1

                cur_i2 = this%grid_mgr%boosted_cell_coords(cur_idx)%t2_idx

                if( i2 .ne. cur_i2 ) then
                    print *, "Fatal error: wrong i2 in os-boosted-diag.f03"
                    stop
                endif

                cur_j2 = this%grid_mgr%boosted_cell_coords(cur_idx)%x2_idx

                buff_idx = cur_j2 - min_j2 + 1

                ! finally load the write buffer
                select case( this%x_dim_ )
                case(1)
                    write_buf_1d(buff_idx) = data_buf%f1( data_idx, cur_idx_within_databuffer )
                case(2)
                    write_buf_2d(buff_idx, 1:data_buf%nx_(2)) = data_buf%f2( data_idx, cur_idx_within_databuffer, 1:data_buf%nx_(2) )
                case(3)
                    write_buf_3d(buff_idx, 1:data_buf%nx_(2), 1:data_buf%nx_(3)) = data_buf%f3( &
                                    data_idx, cur_idx_within_databuffer, 1:data_buf%nx_(2), 1:data_buf%nx_(3) )
                end select

                ! ! useful debug info -- keeping this here for now
                ! if( num_data > 0 .and. component .eq. 1 .and. i2 .eq. 5 .and. this%diag_type .eq. p_boosted_4vector ) then
                ! print *, 'debug ', mpi_node(), cur_j2, min_j2, num_data, write_buf_2d(buff_idx, 128)
                ! endif

            enddo

            ! update start index to account for progress through the
            ! cell coordinates we've looped through
            start_idx = start_idx + num_data

            call this % get_dir_fname_label( component_idx, i2 - 1, dir, filename, label)

            diagfile % filename = filename
            diagfile % filepath = dir
            diagfile % name = label
            dset % name = label
            dset % ndims = this%x_dim_
            chunk % count(1) = num_data
            chunk % count(2:3) = this%grid_mgr%nxperp(1:2)

            chunk % start(1) = min_j2 - 1
            chunk % start(2:3) = this%grid_mgr%nxperp_min(1:2) - 1
            chunk % stride(1:3) = 1

            ! does not currently work
            ! call diagfile % open( p_diag_update, grid%io%comm, DIAG_MPIIO_INDEPENDENT )

            ! call diagfile % open( p_diag_update, grid%io%comm, DIAG_MPI )
            call diagfile % open( p_diag_update, grid%io%comm )

            call diagfile % open_cdset( dset )

            ! this is necessary because this metadata is only read by
            ! the one process that opens it
            dset % count(1) = this%grid_mgr%nx2
            dset % count(2:3) = this%grid_mgr%gnxperp(1:2)
            dset % data_type = p_diag_float32

            call diagfile % write_par_cdset( dset, chunk )

            call diagfile % end_cdset( dset )

            call diagfile % close()

            if( num_data > 0 ) then
                call freemem(write_buf_1d,"boost/os-boosted-diag.f03",1940)
                call freemem(write_buf_2d,"boost/os-boosted-diag.f03",1941)
                call freemem(write_buf_3d,"boost/os-boosted-diag.f03",1942)
            endif

        enddo ! i2 loop

        deallocate( diagFile )

        call freemem(num_data_at_i2,"boost/os-boosted-diag.f03",1949)

        if( this%num_data_in_buf > 0 ) then

            call freemem(i2_vals,"boost/os-boosted-diag.f03",1953)
            call freemem(sorted_i2_indices,"boost/os-boosted-diag.f03",1954)

        endif

end subroutine report_component_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! name: prefix for the data
! component: component of the vdf where the data is stored
! i indexes time, j indexes space
!-------------------------------------------------------------------------------------------------
subroutine report_component_boosted_diag_slice( this, component, grid, component_idx)

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    class( t_grid ), intent(in) :: grid
    integer, intent(in) :: component, component_idx

    integer, dimension(3,3) :: nodeIdx
    character(len=1024) :: dir
    character(len=128 ) :: filename
    character(len=128 ) :: label

    class( t_diag_file ), allocatable :: diagfile
    type( t_diag_dataset ) :: dset
    type( t_diag_chunk ) :: chunk

    ! type(t_vdf) :: write_buf
    type(t_vdf), pointer :: data_buf

    real(p_single), dimension(:,:), pointer :: write_buf_3d_slice


    integer :: i2, j, idx, cur_i2, cur_j2, data_idx
    integer :: start_idx, num_data, cur_idx, min_j2, buff_idx
    integer :: cur_idx_within_databuffer
    integer :: direction, gipos, lperp !slice position

    ! i2 values obtained in the unsorted buffer
    integer, dimension(:), pointer :: i2_vals, num_data_at_i2

    ! integer, dimension(:,:), pointer :: buffer_indices_for_each_i2

    ! sorted indices of those i2 indices
    integer, dimension(:), pointer :: sorted_i2_indices

    nodeIdx = my_nx_p( grid )
    direction = this%direction(component_idx)
    gipos = this%gipos(component_idx)

    ! set pointers to null
    write_buf_3d_slice => null()
    i2_vals => null(); num_data_at_i2 => null(); sorted_i2_indices => null()

    select case( this%diag_type )

    case( p_boosted_fields )

        ! components 1 to 3 store e
        if( component <= 3 ) then
            data_buf => this%buffered_e
            data_idx = component

        ! components 4 to 6 store b
        else
           data_buf => this%buffered_b
           data_idx = component - 3
       endif

    case( p_boosted_4vector )

       data_buf => this%buffered_4vector
       data_idx = component

    end select

    start_idx = 0

    call alloc(num_data_at_i2, (/ this%grid_mgr%nt2 /),"boost/os-boosted-diag.f03",2036)
    num_data_at_i2(:) = 0

    ! first compute sorted_i2_indices
    if( this%num_data_in_buf > 0 ) then

        call alloc(i2_vals, (/ this%num_data_in_buf /),"boost/os-boosted-diag.f03",2042)
        call alloc(sorted_i2_indices, (/ this%num_data_in_buf /),"boost/os-boosted-diag.f03",2043)

        do j=1, this%num_data_in_buf
            idx = this%buf_start_idx + (j-1)
            i2 = this%grid_mgr%boosted_cell_coords(idx)%t2_idx
            i2_vals(j) = i2
            num_data_at_i2( i2 ) = num_data_at_i2( i2 ) + 1
        enddo

        call heapsort( i2_vals, sorted_i2_indices )

    endif

    call create_diag_file( diagfile )
    ! loop through all lab time slices to reconstruct
    do i2 = 1, this%grid_mgr%nt2

        if (SliceHasLocalData(gipos, direction, nodeIdx)) then

            num_data = num_data_at_i2(i2)
            min_j2 = 1 ! default value in case num_data is 0

            ! we still need to proceed if num_data == 0 because other processes
            ! may have data; in this case, we need to tell the mpi process
            ! doing the writing that we don't have any data, otherwise it
            ! will get stuck waiting for us. however, if we don't have any data
            ! then we don't need to do anything other than give this info to
            ! the writing process.
            if( num_data > 0 ) then

                ! the process has the slice data
                select case(direction)
                    case(2)
                    call alloc(write_buf_3d_slice, (/ num_data, data_buf%nx_(3) /),"boost/os-boosted-diag.f03",2076)
                    chunk % data = c_loc( write_buf_3d_slice )
                    case(3)
                    call alloc(write_buf_3d_slice, (/ num_data, data_buf%nx_(2) /),"boost/os-boosted-diag.f03",2079)
                    chunk % data = c_loc( write_buf_3d_slice )
                end select

            endif

            ! find the minimum j2 attained for all the available data (the heapsort
            ! doesn't respect the j2 ordering posessed by the data prior to sorting )
            do j = 1, num_data
                ! cur_idx_within_databuffer = sorted_i2_indices( start_idx + j )
                cur_idx = this%buf_start_idx + sorted_i2_indices( start_idx + j ) - 1
                cur_j2 = this%grid_mgr%boosted_cell_coords(cur_idx)%x2_idx

                if( j .eq. 1 ) then
                    min_j2 = cur_j2
                else
                    min_j2 = min( min_j2, cur_j2 )
                endif

            enddo

            ! count all the data in the buffer which matches this i2
            do j = 1, num_data

                cur_idx_within_databuffer = sorted_i2_indices( start_idx + j )
                cur_idx = this%buf_start_idx + sorted_i2_indices( start_idx + j ) - 1

                cur_i2 = this%grid_mgr%boosted_cell_coords(cur_idx)%t2_idx

                if( i2 .ne. cur_i2 ) then
                    print *, "Fatal error: wrong i2 in os-boosted-diag.f03"
                    stop
                endif

                cur_j2 = this%grid_mgr%boosted_cell_coords(cur_idx)%x2_idx

                buff_idx = cur_j2 - min_j2 + 1

                ! finally load the write buffer
                call SliceBounds(gipos, direction, nodeIdx, lperp)
                select case(direction)
                case(2)
                    write_buf_3d_slice(buff_idx, 1:data_buf%nx_(3)) = data_buf%f3( &
                                    data_idx, cur_idx_within_databuffer, lperp, 1:data_buf%nx_(3) )
                case(3)
                    write_buf_3d_slice(buff_idx, 1:data_buf%nx_(2)) = data_buf%f3( &
                                    data_idx, cur_idx_within_databuffer, 1:data_buf%nx_(2), lperp )
                end select

                ! ! useful debug info -- keeping this here for now
                ! if( num_data > 0 .and. component .eq. 2 .and. i2 .eq. 5 ) then
                ! print *, 'place 2', mpi_node(), num_data, write_buf_3d_slice(buff_idx, 250)
                ! endif
            enddo

            ! update start index to account for progress through the
            ! cell coordinates we've looped through
            start_idx = start_idx + num_data

            call this % get_dir_fname_label_slice( component_idx, i2 - 1, dir, filename, label)

            diagfile % filename = filename
            diagfile % filepath = dir
            diagfile % name = label
            dset % name = label
            dset % ndims = 2
            select case(direction)
            case(2)
                chunk % count(1) = num_data
                chunk % count(2) = this%grid_mgr%nxperp(2)

                chunk % start(1) = min_j2 - 1
                chunk % start(2) = this%grid_mgr%nxperp_min(2) - 1
            case(3)
                chunk % count(1) = num_data
                chunk % count(2) = this%grid_mgr%nxperp(1)

                chunk % start(1) = min_j2 - 1
                chunk % start(2) = this%grid_mgr%nxperp_min(1) - 1
            end select
            chunk % stride(1:2) = 1

            call diagfile % open( p_diag_update, grid%io%comm )

            call diagfile % open_cdset( dset )
            ! this is necessary because this metadata is only read by
            ! the one process that opens it
            dset % count(1) = this%grid_mgr%nx2
            select case(direction)
            case(2)
                dset % count(2) = this%grid_mgr%gnxperp(2)
            case(3)
                dset % count(2) = this%grid_mgr%gnxperp(1)
            end select
            dset % data_type = p_diag_float32

            call diagfile % write_par_cdset( dset, chunk )
            call diagfile % end_cdset( dset )

            call diagfile % close()

            if( num_data > 0 ) then
                call freemem(write_buf_3d_slice,"boost/os-boosted-diag.f03",2181)
            endif
        else
            num_data = 0 ! no data writing
            min_j2 = 1 ! default value in case num_data is 0

            ! we still need to proceed if num_data == 0 because other processes
            ! may have data; in this case, we need to tell the mpi process
            ! doing the writing that we don't have any data, otherwise it
            ! will get stuck waiting for us. however, if we don't have any data
            ! then we don't need to do anything other than give this info to
            ! the writing process.

            call this % get_dir_fname_label_slice( component_idx, i2 - 1, dir, filename, label)

            diagfile % filename = filename
            diagfile % filepath = dir
            diagfile % name = label
            dset % name = label
            dset % ndims = 2
            select case(direction)
            case(2)
                chunk % count(1) = num_data
                chunk % count(2) = this%grid_mgr%nxperp(2)

                chunk % start(1) = min_j2 - 1
                chunk % start(2) = this%grid_mgr%nxperp_min(2) - 1
            case(3)
                chunk % count(1) = num_data
                chunk % count(2) = this%grid_mgr%nxperp(1)

                chunk % start(1) = min_j2 - 1
                chunk % start(2) = this%grid_mgr%nxperp_min(1) - 1
            end select
            chunk % stride(1:2) = 1

            call diagfile % open( p_diag_update, grid%io%comm )

            call diagfile % open_cdset( dset )
            ! this is necessary because this metadata is only read by
            ! the one process that opens it
            dset % count(1) = this%grid_mgr%nx2
            select case(direction)
            case(2)
                dset % count(2) = this%grid_mgr%gnxperp(2)
            case(3)
                dset % count(2) = this%grid_mgr%gnxperp(1)
            end select
            dset % data_type = p_diag_float32

            call diagfile % write_par_cdset( dset, chunk )
            call diagfile % end_cdset( dset )

            call diagfile % close()

            if( num_data > 0 ) then
                call freemem(write_buf_3d_slice,"boost/os-boosted-diag.f03",2237)
            endif
        endif

    enddo ! i2 loop

    deallocate( diagFile )

    call freemem(num_data_at_i2,"boost/os-boosted-diag.f03",2245)

    if( this%num_data_in_buf > 0 ) then

        call freemem(i2_vals,"boost/os-boosted-diag.f03",2249)
        call freemem(sorted_i2_indices,"boost/os-boosted-diag.f03",2250)

    endif


end subroutine report_component_boosted_diag_slice
!-------------------------------------------------------------------------------------------------




!-------------------------------------------------------------------------------------------------
! this is called at every timestep to extract all data out of the current and
! previous timestep that corresponds to a lab frame timeslice
! in 3D, diag components and components_slice are looped through respectively
!-------------------------------------------------------------------------------------------------
subroutine extract_boosted_diag( this, tstep, t1, grid, g_space, vdf1, vdf2 )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    type(t_time_step), intent(in) :: tstep
    real(p_double), intent(in) :: t1
    class( t_grid ), intent(in) :: grid
    type( t_space ), intent(in) :: g_space
    class(t_vdf), intent(in) :: vdf1
    class(t_vdf), intent(in), optional :: vdf2

    integer :: i1, j, buff_idx, cell_idx !, n
    integer :: x1_idx, x1_idx_prev
    integer :: component_idx, component, direction, gipos
    real(p_double) :: x1_cell, t1_cell, xmin1
    real(p_double) :: dx1, dt1
    real(p_double) :: t1_p, t1_m, x1_p, x1_m
    real(p_double) :: coef_tp, coef_tm, coef_xp, coef_xm
    real(p_double) :: coef_pp, coef_pm, coef_mp, coef_mm
    real(p_double) :: gamma, beta

    real(p_k_fld), dimension(3) :: e_pp, e_pm, e_mp, e_mm
    real(p_k_fld), dimension(3) :: b_pp, b_pm, b_mp, b_mm
    real(p_k_fld), dimension(3) :: buffered_e, buffered_b
    real(p_k_fld) :: buffered_field, buffered_4vector_component !calculate each component
    real(p_k_fld), dimension(4) :: x_pp, x_pm, x_mp, x_mm
    real(p_k_fld), dimension(4) :: buffered_4vector
    integer, dimension(3,3) :: nodeIdx
    integer :: i2, i3

    if (this%ndump_fac .le. 0) return

    nodeIdx(:,1:grid%x_dim) = my_nx_p( grid )
    gamma = this%grid_mgr%gamma
    beta = this%grid_mgr%beta

    ! i1 has the property t = (i1 - 1) * dt1 + tmin
    i1 = n( tstep ) + 1

    buff_idx = this%num_data_in_buf + 1

    cell_idx = this%grid_mgr%starting_cell_idx_at_t1(i1)

    xmin1 = xmin( g_space, 1 ) + this%grid_mgr%dx1 * (grid%my_nx(p_lower,1) - 1)

    dx1 = this%grid_mgr%dx1
    dt1 = this%grid_mgr%dt1

    t1_p = t1
    t1_m = t1 - dt1

    do j = 1, this % grid_mgr % num_t2_points_at_t1(i1)

        t1_cell = this%grid_mgr%boosted_cell_coords(cell_idx)%t1
        x1_cell = this%grid_mgr%boosted_cell_coords(cell_idx)%x1

        ! find the spatial cells to use for this interpolation
        x1_idx = floor( ( x1_cell - xmin1 ) / dx1 )
        x1_idx_prev = floor( ( x1_cell - this%xmin1_prev ) / dx1 )

        ! gridpoints of the interpolation
        x1_m = xmin1 + x1_idx * dx1
        x1_p = x1_m + dx1

! for serious debugging
# 2353 "boost/os-boosted-diag.f03"
        ! m, p --> minus, plus
        coef_tp = (t1_p - t1_cell) / dt1
        coef_xp = (x1_p - x1_cell) / dx1
        coef_tm = (t1_cell - t1_m) / dt1
        coef_xm = (x1_cell - x1_m) / dx1

        coef_pp = coef_tp * coef_xp
        coef_pm = coef_tp * coef_xm
        coef_mp = coef_tm * coef_xp
        coef_mm = coef_tm * coef_xm

        select case( this%diag_type )

        case( p_boosted_fields )
            select case(this%x_dim_)

            case(1)

                e_pp = vdf1%f1(:,x1_idx + 1)
                e_pm = vdf1%f1(:,x1_idx )
                e_mp = this%prev_e%f1(:,x1_idx_prev + 1)
                e_mm = this%prev_e%f1(:,x1_idx_prev )

                b_pp = vdf2%f1(:,x1_idx + 1)
                b_pm = vdf2%f1(:,x1_idx )
                b_mp = this%prev_b%f1(:,x1_idx_prev + 1)
                b_mm = this%prev_b%f1(:,x1_idx_prev )

                call extract_fields_helper( e_pp, e_pm, e_mp, e_mm, &
                                            b_pp, b_pm, b_mp, b_mm, &
                                            coef_pp, coef_pm, coef_mp, coef_mm, &
                                            gamma, beta, &
                                            buffered_e, buffered_b )

                this%buffered_e%f1(:,buff_idx) = buffered_e
                this%buffered_b%f1(:,buff_idx) = buffered_b

            case(2)

                do i2 = 1, this%grid_mgr%nxperp(1)
                    ! if (mpi_node()==0) print *,'i2', i2
                    e_pp = vdf1%f2(:, x1_idx + 1, i2 )
                    e_pm = vdf1%f2(:, x1_idx, i2 )
                    e_mp = this%prev_e%f2(:, x1_idx_prev + 1, i2 )
                    e_mm = this%prev_e%f2(:, x1_idx_prev, i2 )

                    b_pp = vdf2%f2(:, x1_idx + 1, i2 )
                    b_pm = vdf2%f2(:, x1_idx, i2 )
                    b_mp = this%prev_b%f2(:, x1_idx_prev + 1, i2)
                    b_mm = this%prev_b%f2(:, x1_idx_prev, i2 )

                    call extract_fields_helper( e_pp, e_pm, e_mp, e_mm, &
                                                b_pp, b_pm, b_mp, b_mm, &
                                                coef_pp, coef_pm, coef_mp, coef_mm, &
                                                gamma, beta, &
                                                buffered_e, buffered_b )

                    this%buffered_e%f2(:, buff_idx, i2) = buffered_e
                    this%buffered_b%f2(:, buff_idx, i2) = buffered_b
                enddo

            case(3)
                do component_idx = 1, this%num_diag_components
                    do i2 = 1, this%grid_mgr%nxperp(1)
                        do i3 = 1, this%grid_mgr%nxperp(2)

                            e_pp = vdf1%f3(:, x1_idx + 1, i2, i3)
                            e_pm = vdf1%f3(:, x1_idx, i2, i3)
                            e_mp = this%prev_e%f3(:, x1_idx_prev + 1, i2, i3)
                            e_mm = this%prev_e%f3(:, x1_idx_prev, i2, i3)

                            b_pp = vdf2%f3(:, x1_idx + 1, i2, i3)
                            b_pm = vdf2%f3(:, x1_idx, i2, i3)
                            b_mp = this%prev_b%f3(:, x1_idx_prev + 1, i2, i3)
                            b_mm = this%prev_b%f3(:, x1_idx_prev, i2, i3)

                            component = this%diag_components(component_idx)

                            call extract_field_component_helper( e_pp, e_pm, e_mp, e_mm, &
                                                        b_pp, b_pm, b_mp, b_mm, &
                                                        coef_pp, coef_pm, coef_mp, coef_mm, &
                                                        gamma, beta, &
                                                        buffered_field, component )
                            if (component <=3) then
                                this%buffered_e%f3(component, buff_idx, i2, i3) = buffered_field
                            else
                                this%buffered_b%f3(component-3, buff_idx, i2, i3) = buffered_field
                            endif
                        enddo
                    enddo
                enddo

                !loop through slice diag components
                do component_idx = 1, this%num_diag_components_slice
                    direction = this%direction(component_idx)
                    gipos = this%gipos(component_idx)
                    !if this node has slice data, record it
                    if (SliceHasLocalData(gipos, direction, nodeIdx)) then
                        select case(direction)
                        case(2)
                            call SliceBounds(gipos, direction, nodeIdx, i2)
                            do i3 = 1, this%grid_mgr%nxperp(2)
                                e_pp = vdf1%f3(:, x1_idx + 1, i2, i3)
                                e_pm = vdf1%f3(:, x1_idx, i2, i3)
                                e_mp = this%prev_e%f3(:, x1_idx_prev + 1, i2, i3)
                                e_mm = this%prev_e%f3(:, x1_idx_prev, i2, i3)

                                b_pp = vdf2%f3(:, x1_idx + 1, i2, i3)
                                b_pm = vdf2%f3(:, x1_idx, i2, i3)
                                b_mp = this%prev_b%f3(:, x1_idx_prev + 1, i2, i3)
                                b_mm = this%prev_b%f3(:, x1_idx_prev, i2, i3)

                                component = this%diag_components_slice(component_idx)

                                call extract_field_component_helper( e_pp, e_pm, e_mp, e_mm, &
                                                            b_pp, b_pm, b_mp, b_mm, &
                                                            coef_pp, coef_pm, coef_mp, coef_mm, &
                                                            gamma, beta, &
                                                            buffered_field, component )
                                if (component <=3) then
                                    this%buffered_e%f3(component, buff_idx, i2, i3) = buffered_field
                                else
                                    this%buffered_b%f3(component-3, buff_idx, i2, i3) = buffered_field
                                endif
                            enddo
                        case(3)
                            call SliceBounds(gipos, direction, nodeIdx, i3)
                            do i2 = 1, this%grid_mgr%nxperp(1)
                                e_pp = vdf1%f3(:, x1_idx + 1, i2, i3)
                                e_pm = vdf1%f3(:, x1_idx, i2, i3)
                                e_mp = this%prev_e%f3(:, x1_idx_prev + 1, i2, i3)
                                e_mm = this%prev_e%f3(:, x1_idx_prev, i2, i3)

                                b_pp = vdf2%f3(:, x1_idx + 1, i2, i3)
                                b_pm = vdf2%f3(:, x1_idx, i2, i3)
                                b_mp = this%prev_b%f3(:, x1_idx_prev + 1, i2, i3)
                                b_mm = this%prev_b%f3(:, x1_idx_prev, i2, i3)

                                component = this%diag_components_slice(component_idx)

                                call extract_field_component_helper( e_pp, e_pm, e_mp, e_mm, &
                                                            b_pp, b_pm, b_mp, b_mm, &
                                                            coef_pp, coef_pm, coef_mp, coef_mm, &
                                                            gamma, beta, &
                                                            buffered_field, component )
                                if (component <=3) then
                                    this%buffered_e%f3(component, buff_idx, i2, i3) = buffered_field
                                else
                                    this%buffered_b%f3(component-3, buff_idx, i2, i3) = buffered_field
                                endif
                            enddo
                        end select
                    endif
                enddo
            end select


        case( p_boosted_4vector )

            select case( this%x_dim_ )

            case(1)

                x_pp = vdf1%f1(:,x1_idx + 1)
                x_pm = vdf1%f1(:,x1_idx )
                x_mp = this%prev_4vector%f1(:,x1_idx_prev + 1)
                x_mm = this%prev_4vector%f1(:,x1_idx_prev )

                call extract_4vector_helper( x_pp, x_pm, x_mp, x_mm, &
                                             coef_pp, coef_pm, coef_mp, coef_mm, &
                                             gamma, beta, &
                                             buffered_4vector )

                this%buffered_4vector%f1(:,buff_idx) = buffered_4vector

            case(2)

                do i2 = 1, this%grid_mgr%nxperp(1)

                    x_pp = vdf1%f2(:,x1_idx + 1, i2)
                    x_pm = vdf1%f2(:,x1_idx, i2 )
                    x_mp = this%prev_4vector%f2(:,x1_idx_prev + 1, i2)
                    x_mm = this%prev_4vector%f2(:,x1_idx_prev, i2 )

                    call extract_4vector_helper( x_pp, x_pm, x_mp, x_mm, &
                                                 coef_pp, coef_pm, coef_mp, coef_mm, &
                                                 gamma, beta, &
                                                 buffered_4vector )

                    this%buffered_4vector%f2(:,buff_idx,i2) = buffered_4vector

                enddo

            case(3)
                do component_idx = 1, this%num_diag_components
                    do i2 = 1, this%grid_mgr%nxperp(1)
                        do i3 = 1, this%grid_mgr%nxperp(2)

                            x_pp = vdf1%f3(:,x1_idx + 1, i2, i3)
                            x_pm = vdf1%f3(:,x1_idx, i2, i3 )
                            x_mp = this%prev_4vector%f3(:,x1_idx_prev + 1, i2, i3)
                            x_mm = this%prev_4vector%f3(:,x1_idx_prev, i2, i3 )

                            component = this%diag_components(component_idx)
                            call extract_4vector_component_helper( x_pp, x_pm, x_mp, x_mm, &
                                                        coef_pp, coef_pm, coef_mp, coef_mm, &
                                                        gamma, beta, &
                                                        buffered_4vector_component, component )

                            this%buffered_4vector%f3(component,buff_idx,i2, i3) = buffered_4vector_component

                        enddo
                    enddo
                enddo

                !loop through slice diag components
                do component_idx = 1, this%num_diag_components_slice
                    direction = this%direction(component_idx)
                    gipos = this%gipos(component_idx)
                    !if this node has slice data, record it
                    if (SliceHasLocalData(gipos, direction, nodeIdx)) then
                        select case(direction)
                        case(2)
                            call SliceBounds(gipos, direction, nodeIdx, i2)
                            do i3 = 1, this%grid_mgr%nxperp(2)

                                x_pp = vdf1%f3(:,x1_idx + 1, i2, i3)
                                x_pm = vdf1%f3(:,x1_idx, i2, i3 )
                                x_mp = this%prev_4vector%f3(:,x1_idx_prev + 1, i2, i3)
                                x_mm = this%prev_4vector%f3(:,x1_idx_prev, i2, i3 )

                                component = this%diag_components_slice(component_idx)
                                call extract_4vector_component_helper( x_pp, x_pm, x_mp, x_mm, &
                                                            coef_pp, coef_pm, coef_mp, coef_mm, &
                                                            gamma, beta, &
                                                            buffered_4vector_component, component )

                                this%buffered_4vector%f3(component,buff_idx,i2, i3) = buffered_4vector_component
                            enddo
                        case(3)
                            call SliceBounds(gipos, direction, nodeIdx, i3)
                            do i2 = 1, this%grid_mgr%nxperp(1)

                                x_pp = vdf1%f3(:,x1_idx + 1, i2, i3)
                                x_pm = vdf1%f3(:,x1_idx, i2, i3 )
                                x_mp = this%prev_4vector%f3(:,x1_idx_prev + 1, i2, i3)
                                x_mm = this%prev_4vector%f3(:,x1_idx_prev, i2, i3 )

                                component = this%diag_components_slice(component_idx)
                                call extract_4vector_component_helper( x_pp, x_pm, x_mp, x_mm, &
                                                            coef_pp, coef_pm, coef_mp, coef_mm, &
                                                            gamma, beta, &
                                                            buffered_4vector_component, component )

                                this%buffered_4vector%f3(component,buff_idx,i2, i3) = buffered_4vector_component
                            enddo
                        end select
                    endif
                enddo
            end select

        end select

        cell_idx = cell_idx + 1
        buff_idx = buff_idx + 1

    enddo

    ! account for new stored data
    this%num_data_in_buf = this%num_data_in_buf + this%grid_mgr%num_t2_points_at_t1(i1)

    this%xmin1_prev = xmin1

    select case( this%diag_type )

    case( p_boosted_fields )

        call check_vdf_compatability( this%prev_e, vdf1 )
        call check_vdf_compatability( this%prev_b, vdf2 )

        call this%prev_e%copy( vdf1 )
        call this%prev_b%copy( vdf2 )

    case( p_boosted_4vector )

        call check_vdf_compatability( this%prev_4vector, vdf1 )

        call this%prev_4vector%copy( vdf1 )

    end select

end subroutine extract_boosted_diag
!-------------------------------------------------------------------------------------------------

!calculate each component
!-------------------------------------------------------------------------------------------------
subroutine extract_field_component_helper( e_pp, e_pm, e_mp, e_mm, &
                                     b_pp, b_pm, b_mp, b_mm, &
                                     coef_pp, coef_pm, coef_mp, coef_mm, &
                                     gamma, beta, &
                                     buffered_field, component)

    implicit none

    real(p_k_fld), dimension(3), intent(in) :: e_pp, e_pm, e_mp, e_mm
    real(p_k_fld), dimension(3), intent(in) :: b_pp, b_pm, b_mp, b_mm
    real(p_double), intent(in) :: coef_pp, coef_pm, coef_mp, coef_mm
    real(p_double), intent(in) :: gamma, beta
    real(p_k_fld), intent(out) :: buffered_field
    integer, intent(in) :: component

    real(p_k_fld), dimension(3) :: e_interp, b_interp

    ! ! first index is time, second is space
    ! e_interp = real( coef_pp * e_mm + coef_pm * e_mp &
    ! + coef_mp * e_pm + coef_mm * e_pp, p_k_fld )
    ! b_interp = real( coef_pp * b_mm + coef_pm * b_mp &
    ! + coef_mp * b_pm + coef_mm * b_pp, p_k_fld )

    select case(component)
    case(1)
        e_interp(1) = real( coef_pp * e_mm(1) + coef_pm * e_mp(1) &
                + coef_mp * e_pm(1) + coef_mm * e_pp(1), p_k_fld )
        buffered_field = e_interp(1)
    case(2)
        e_interp(2) = real( coef_pp * e_mm(2) + coef_pm * e_mp(2) &
                + coef_mp * e_pm(2) + coef_mm * e_pp(2), p_k_fld )
        b_interp(3) = real( coef_pp * b_mm(3) + coef_pm * b_mp(3) &
                + coef_mp * b_pm(3) + coef_mm * b_pp(3), p_k_fld )
        buffered_field = real( gamma * ( e_interp(2) - beta * b_interp(3) ), p_k_fld )
    case(3)
        e_interp(3) = real( coef_pp * e_mm(3) + coef_pm * e_mp(3) &
                + coef_mp * e_pm(3) + coef_mm * e_pp(3), p_k_fld )
        b_interp(2) = real( coef_pp * b_mm(2) + coef_pm * b_mp(2) &
                + coef_mp * b_pm(2) + coef_mm * b_pp(2), p_k_fld )
        buffered_field = real( gamma * ( e_interp(3) + beta * b_interp(2) ), p_k_fld )
    case(4)
        b_interp(1) = real( coef_pp * b_mm(1) + coef_pm * b_mp(1) &
                   + coef_mp * b_pm(1) + coef_mm * b_pp(1), p_k_fld )
        buffered_field = b_interp(1)
    case(5)
        e_interp(3) = real( coef_pp * e_mm(3) + coef_pm * e_mp(3) &
                + coef_mp * e_pm(3) + coef_mm * e_pp(3), p_k_fld )
        b_interp(2) = real( coef_pp * b_mm(2) + coef_pm * b_mp(2) &
                + coef_mp * b_pm(2) + coef_mm * b_pp(2), p_k_fld )
        buffered_field = real( gamma * ( b_interp(2) + beta * e_interp(3) ), p_k_fld )
    case(6)
        e_interp(2) = real( coef_pp * e_mm(2) + coef_pm * e_mp(2) &
                + coef_mp * e_pm(2) + coef_mm * e_pp(2), p_k_fld )
        b_interp(3) = real( coef_pp * b_mm(3) + coef_pm * b_mp(3) &
                + coef_mp * b_pm(3) + coef_mm * b_pp(3), p_k_fld )
        buffered_field = real( gamma * ( b_interp(3) - beta * e_interp(2) ), p_k_fld )
    end select
    ! now populate the buffers
    ! Eperp ' = gamma[ E_perp - v \times B ] for boost into lab

end subroutine extract_field_component_helper
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! don't mind me, just a nice little helper function at your service
! this function / usage of it is not optimized
!-------------------------------------------------------------------------------------------------
subroutine extract_fields_helper( e_pp, e_pm, e_mp, e_mm, &
                                     b_pp, b_pm, b_mp, b_mm, &
                                     coef_pp, coef_pm, coef_mp, coef_mm, &
                                     gamma, beta, &
                                     buffered_e, buffered_b )

    implicit none

    real(p_k_fld), dimension(3), intent(in) :: e_pp, e_pm, e_mp, e_mm
    real(p_k_fld), dimension(3), intent(in) :: b_pp, b_pm, b_mp, b_mm
    real(p_double), intent(in) :: coef_pp, coef_pm, coef_mp, coef_mm
    real(p_double), intent(in) :: gamma, beta
    real(p_k_fld), dimension(3), intent(out) :: buffered_e, buffered_b

    real(p_k_fld), dimension(3) :: e_interp, b_interp

    ! first index is time, second is space
    e_interp = real( coef_pp * e_mm + coef_pm * e_mp &
                   + coef_mp * e_pm + coef_mm * e_pp, p_k_fld )

    b_interp = real( coef_pp * b_mm + coef_pm * b_mp &
                   + coef_mp * b_pm + coef_mm * b_pp, p_k_fld )

    ! now populate the buffers
    ! Eperp ' = gamma[ E_perp - v \times B ] for boost into lab
    buffered_e(1) = e_interp(1)
    buffered_e(2) = real( gamma * ( e_interp(2) - beta * b_interp(3) ), p_k_fld )
    buffered_e(3) = real( gamma * ( e_interp(3) + beta * b_interp(2) ), p_k_fld )

    buffered_b(1) = b_interp(1)
    buffered_b(2) = real( gamma * ( b_interp(2) + beta * e_interp(3) ), p_k_fld )
    buffered_b(3) = real( gamma * ( b_interp(3) - beta * e_interp(2) ), p_k_fld )

end subroutine extract_fields_helper
!-------------------------------------------------------------------------------------------------


!calculate each component
!-------------------------------------------------------------------------------------------------
subroutine extract_4vector_component_helper( x_pp, x_pm, x_mp, x_mm, &
                                     coef_pp, coef_pm, coef_mp, coef_mm, &
                                     gamma, beta, &
                                     buffered_4vector_component, component )

    implicit none

    real(p_k_fld), dimension(4), intent(in) :: x_pp, x_pm, x_mp, x_mm
    real(p_double), intent(in) :: coef_pp, coef_pm, coef_mp, coef_mm
    real(p_double), intent(in) :: gamma, beta

    integer, intent(in) :: component
    real(p_k_fld), intent(out) :: buffered_4vector_component

    real(p_k_fld), dimension(4) :: x_interp

    ! first index is time, second is space
    select case(component)
    case(1)
        x_interp(1) = real( coef_pp * x_mm(1) + coef_pm * x_mp(1) &
                    + coef_mp * x_pm(1) + coef_mm * x_pp(1), p_k_fld )
        x_interp(2) = real( coef_pp * x_mm(2) + coef_pm * x_mp(2) &
                    + coef_mp * x_pm(2) + coef_mm * x_pp(2), p_k_fld )
        buffered_4vector_component = real( gamma * ( x_interp(1) - beta * x_interp(2) ), p_k_fld )
    case(2)
        x_interp(1) = real( coef_pp * x_mm(1) + coef_pm * x_mp(1) &
                    + coef_mp * x_pm(1) + coef_mm * x_pp(1), p_k_fld )
        x_interp(2) = real( coef_pp * x_mm(2) + coef_pm * x_mp(2) &
                    + coef_mp * x_pm(2) + coef_mm * x_pp(2), p_k_fld )
        buffered_4vector_component = real( gamma * ( x_interp(2) - beta * x_interp(1) ), p_k_fld )
    case(3)
        x_interp(3) = real( coef_pp * x_mm(3) + coef_pm * x_mp(3) &
                    + coef_mp * x_pm(3) + coef_mm * x_pp(3), p_k_fld )
        buffered_4vector_component = x_interp(3)
    case(4)
        x_interp(4) = real( coef_pp * x_mm(4) + coef_pm * x_mp(4) &
                    + coef_mp * x_pm(4) + coef_mm * x_pp(4), p_k_fld )
        buffered_4vector_component = x_interp(4)
    end select

end subroutine extract_4vector_component_helper
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine extract_4vector_helper( x_pp, x_pm, x_mp, x_mm, &
                                     coef_pp, coef_pm, coef_mp, coef_mm, &
                                     gamma, beta, &
                                     buffered_4vector )

    implicit none

    real(p_k_fld), dimension(4), intent(in) :: x_pp, x_pm, x_mp, x_mm
    real(p_double), intent(in) :: coef_pp, coef_pm, coef_mp, coef_mm
    real(p_double), intent(in) :: gamma, beta

    real(p_k_fld), dimension(4), intent(out) :: buffered_4vector

    real(p_k_fld), dimension(4) :: x_interp

    ! first index is time, second is space
    x_interp = real( coef_pp * x_mm + coef_pm * x_mp &
                   + coef_mp * x_pm + coef_mm * x_pp, p_k_fld )

    buffered_4vector(1) = real( gamma * ( x_interp(1) - beta * x_interp(2) ), p_k_fld )
    buffered_4vector(2) = real( gamma * ( x_interp(2) - beta * x_interp(1) ), p_k_fld )
    buffered_4vector(3:4) = x_interp(3:4)

end subroutine extract_4vector_helper
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
! grow the raw particle buffers
!-------------------------------------------------------------------------------------------------
subroutine grow_buffer_boosted_diag( this )

    implicit none

    class(t_boosted_diag), intent(inout) :: this

    real(p_k_part), dimension(:,:), pointer :: temp1
    integer, dimension(:), pointer :: temp2
    integer, dimension(:,:), pointer :: temp3
    integer :: spec_dim, num_par_new

    spec_dim = size( this%buffered_part, 1 )
    num_par_new = this%num_par + p_raw_buf_block

    call alloc(temp1, (/ spec_dim, num_par_new /),"boost/os-boosted-diag.f03",2845)
    call memcpy( temp1, this%buffered_part, spec_dim * this%num_par )
    call freemem(this%buffered_part,"boost/os-boosted-diag.f03",2847)
    this%buffered_part => temp1

    call alloc(temp2, (/ num_par_new /),"boost/os-boosted-diag.f03",2850)
    call memcpy( temp2, this%buffered_diag_idx, this%num_par )
    call freemem(this%buffered_diag_idx,"boost/os-boosted-diag.f03",2852)
    this%buffered_diag_idx => temp2

    if ( associated(this%buffered_tags) ) then
        call alloc(temp3, (/ 2, num_par_new /),"boost/os-boosted-diag.f03",2856)
        call memcpy( temp3, this%buffered_tags, 2 * this%num_par )
        call freemem(this%buffered_tags,"boost/os-boosted-diag.f03",2858)
        this%buffered_tags => temp3
    endif

    this%num_par_max = num_par_new

end subroutine grow_buffer_boosted_diag
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
! push particles back in all dimensions except for x1
!-------------------------------------------------------------------------------------------------
subroutine push_back_boosted_diag( this, x, p, rgamma, t1_m, t1 )

    implicit none

    class(t_boosted_diag), intent(in) :: this
    real(p_k_part), dimension(:), intent(inout) :: x
    real(p_k_part), dimension(:), intent(in) :: p
    real(p_k_part), intent(in) :: rgamma
    real(p_double), intent(in) :: t1_m, t1

    integer :: k
    real(p_k_part) :: vk

    do k = 2, p_x_dim
        vk = p(k) * rgamma
        x(k) = real( x(k) + vk * ( t1_m - t1 ), p_k_part )
    enddo

end subroutine push_back_boosted_diag
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine check_vdf_compatability( vdf1, vdf2 )

    implicit none

    type(t_vdf), intent(in) :: vdf1
    type(t_vdf), intent(in) :: vdf2

    logical :: success
    integer :: x_dim

    x_dim = vdf1%x_dim_

    success = ( (vdf1%x_dim_ .eq. vdf2%x_dim_) &
            .and. all( vdf1%nx_(1:x_dim) .eq. vdf2%nx_(1:x_dim) ) &
            .and. all( vdf1%gc_num_(:,1:x_dim) .eq. vdf2%gc_num_(:,1:x_dim) ) )

    if( .not. success ) then
        write(err_buf__,*) 'ERROR: incompatible vdf shapes in os-boosted-diag';call err__("boost/os-boosted-diag.f03",2912)
        write(err_buf__,*) 'x_dim_: ', vdf1%x_dim_, vdf2%x_dim_;call err__("boost/os-boosted-diag.f03",2913)
        write(err_buf__,*) 'nx_:', vdf1%nx_(1:x_dim), vdf2%nx_(1:x_dim);call err__("boost/os-boosted-diag.f03",2914)
        write(err_buf__,*) 'gc_num_:', vdf1%gc_num_(:,1:x_dim), vdf2%gc_num_(:,1:x_dim);call err__("boost/os-boosted-diag.f03",2915)
        call abort_program(p_err_invalid)
    endif

end subroutine check_vdf_compatability
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine reshape_array_vdf( vdf, nx_new )

    implicit none

    type(t_vdf), intent(inout), pointer :: vdf
    integer, dimension(3), intent(in) :: nx_new

    type(t_vdf) :: vdf_tmp
    integer, dimension(3) :: nx_old

    nx_old(:) = vdf%nx_(:)

    call vdf_tmp % new( vdf, copy = .true. )

    call vdf % cleanup()
    call vdf % new( vdf_tmp%x_dim_, vdf_tmp%f_dim_, nx_new, vdf_tmp%gc_num_, vdf_tmp%dx_, .true. )

    ! copy old data into the fresh vdf
    select case( vdf%x_dim_ )
    case(1)
        vdf%f1( :, 1:nx_old(1) ) = vdf_tmp%f1( :, 1:nx_old(1) )
    case(2)
        vdf%f2( :, 1:nx_old(1), 1:nx_old(2) ) = vdf_tmp%f2( :, 1:nx_old(1), 1:nx_old(2) )
    case(3)
        vdf%f3( :, 1:nx_old(1), 1:nx_old(2), 1:nx_old(3) ) = vdf_tmp%f3( :, 1:nx_old(1), 1:nx_old(2), 1:nx_old(3) )
    end select

    call vdf_tmp % cleanup()

end subroutine reshape_array_vdf
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! Generate a list of indices that result in a sorted array.
! The array to be sorted is this%buffered_diag_idx, which contains the indices of the
! boosted diagnostics that each particle will be written to.
! The function returns idx, n_file, and idx_sorted. The idx buffer is this%num_par long
! and provides an index for each particle that when ordered in this way results in a
! sorted array.
! The n_file buffer contains the index of the last particle put into each diagnostic file.
! The idx_sorted buffer contains a sorted list of particle indices.
!
! For example, if this%num_par is 5 and this%buffered_diag_idx is (/ 2, 1, 2, 3, 1 /),
! the results will be idx = (/ 3, 1, 4, 5, 2 /), n_file = (/ 2, 4, 5 /), and
! idx_sorted = (/ 2, 5, 1, 3, 4 /), where it is assumed that there are three boosted
! diagnostic files.
!-------------------------------------------------------------------------------------------------
subroutine generate_sort_idx_( this, idx, n_file, idx_sorted )

    implicit none

    class(t_boosted_diag), intent(inout) :: this
    integer, dimension(:), intent(out) :: idx, n_file, idx_sorted

    integer :: i, index, isum, ist

    n_file = 0

    ! sort by interpolation cell
    do i = 1, this%num_par
        index = this%buffered_diag_idx(i)
        n_file(index) = n_file(index) + 1
        idx(i) = index
    end do

    isum = 0
    do i = 1, size(n_file)
        ist = n_file(i)
        n_file(i) = isum
        isum = isum + ist
    end do

    do i = 1, this%num_par
        index = idx(i)
        n_file(index) = n_file(index) + 1
        idx(i) = n_file(index)
        idx_sorted(idx(i)) = i
    end do

end subroutine generate_sort_idx_
!-------------------------------------------------------------------------------------------------


end module m_boosted_diag
