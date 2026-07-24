#include "os-config.h"
#include "os-preprocess.fpp"

module m_boosted_grid_mgr

#include "memory/memory.h"


use m_system
use m_parameters
use m_space
use m_grid_define

use m_vdf_define
use m_vdf_interpolate
use m_input_file


implicit none

! note: in all the code below, the 1 and 2 indexing makes life a lot easier.
! just remember that the simulation frame is frame 1, and the frame
! which we're transforming to is frame 2.


!-------------------------------------------------------------------------------------------------
! exact values of t2, x2, are not needed so don't store them
! instead store the indices in the array that the data corresponds
! to since this is what is needed in practice.
! exact values of t1, x1 are needed for interpolation
!-------------------------------------------------------------------------------------------------
type :: t_boosted_cell_coords

    integer :: t2_idx, x2_idx
    real(p_double) :: t1, x1

end type
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! exact values of t2, x2 are required here to process particle trajectories
! indices for t1 are needed for file storage
! indices for x1 boundaries are needed for particle testing
!-------------------------------------------------------------------------------------------------
type :: t_boosted_raw_coords

    real(p_double) :: t2, xmin2, xmax2
    integer :: t2_idx, x1_idx_min, x1_idx_max

end type
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! this object is meant to be instantiated once and shared by
! multiple boosted diagnostics
!-------------------------------------------------------------------------------------------------
type :: t_boosted_grid_mgr

    ! boosted params
    real(p_double) :: signed_gamma, gamma, gamma_times_beta, beta

    ! LOCAL bounds for the boosted simulation
    real(p_double) :: tmin1, tmax1, dt1, xmin1, xmax1, dx1
    integer :: nt1, nx1

    ! moving simulation coordinate system
    logical :: if_move1
    real(p_double) :: move_vel1

    ! bounds for the box that will define the lab data
    ! xmin and xmax are given at t=0
    real(p_double) :: tmin2, tmax2, dt2, xmin2, xmax2, dx2
    integer :: nt2, nx2

    ! for using a moving coordinate lab coordinate system (default)
    logical :: if_move2
    real(p_double) :: move_vel2

    ! transverse data
    ! note: i've commented out data that is not currently used but might
    ! be useful in the future

    ! global transverse spatial params
    real(p_double), dimension(2) :: gxperp_min
    real(p_double), dimension(2) :: gxperp_max
    ! real(p_double), dimension(2) :: dxperp
    integer, dimension(2) :: gnxperp

    ! local transverse spatial parameters (same for lab and boost)
    ! real(p_double), dimension(2) :: xperp_min
    ! real(p_double), dimension(2) :: xperp_max

    ! total number of local transverse cell indices
    integer, dimension(2)        :: nxperp

    ! start position of transverse cell indices w/in global box
    integer, dimension(2)        :: nxperp_min


    ! the indices of the data stored in the buffer,
    ! corresponding to indices of this % boosted_cell_coords
    integer :: buffered_idx_min, buffered_idx_max

    ! array of cell coordinates for each coordinate in frame 2
    ! that belongs to computations in frame 1
    type( t_boosted_cell_coords ), dimension(:), pointer :: boosted_cell_coords => null()
    type( t_boosted_raw_coords ), dimension(:), pointer :: boosted_raw_coords => null()

    ! total number of cells stored
    integer :: total_data, total_data_raw

    ! store the number of points that lie between t2 and t2 - dt2
    ! for each t2 achieved by the simulation
    integer, dimension(:), pointer :: num_t2_points_at_t1 => null()

    ! starting index in the array of boosted_cell_coords that will be
    ! processed at timestep t1
    integer, dimension(:), pointer :: starting_cell_idx_at_t1 => null()

    ! store the number of points that lie between t2 and t2 - dt2
    ! for each t2 achieved by the simulation
    integer, dimension(:), pointer :: num_t2_points_at_t1_raw => null()

    ! starting index in the array of boosted_raw_coords that will be
    ! processed at timestep t1
    integer, dimension(:), pointer :: starting_raw_idx_at_t1 => null()

    contains

    procedure :: init => init_boosted_grid_mgr
    procedure :: init_helper => init_helper_boosted_grid_mgr
    procedure :: init_helper_raw => init_helper_raw_boosted_grid_mgr
    procedure :: cleanup => cleanup_boosted_grid_mgr

    procedure :: set_space1_params => set_space1_params_boosted_grid_mgr
    procedure :: set_space2_params => set_space2_params_boosted_grid_mgr
    ! procedure :: read_input => read_input_boosted_grid_mgr

end type
!-------------------------------------------------------------------------------------------------




public :: t_boosted_grid_mgr, t_boosted_cell_coords


private

integer, parameter :: init_buffer_size = 128


contains


!-------------------------------------------------------------------------------------------------
subroutine init_boosted_grid_mgr( this, ndump_fac, ndump_fac_raw )

    implicit none

    class(t_boosted_grid_mgr), intent(inout) :: this
    integer, intent(in) :: ndump_fac, ndump_fac_raw

    this%gamma = abs( this%signed_gamma )
    this%gamma_times_beta = sqrt(this%gamma**2-1.0_p_double) * sign( 1.0_p_double, this%signed_gamma )
    this%beta = sqrt( 1 - 1 / this%gamma**2) * sign( 1.0_p_double, this%signed_gamma )

    if ( ndump_fac > 0 ) then

        call this % init_helper( .true. )

        allocate( this%boosted_cell_coords( this % total_data ) )

        call this % init_helper( .false. )

    endif

    if ( ndump_fac_raw > 0 ) then

        call this % init_helper_raw( .true. )

        allocate( this%boosted_raw_coords( this % total_data_raw ) )

        call this % init_helper_raw( .false. )

    endif

end subroutine init_boosted_grid_mgr
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
subroutine set_space1_params_boosted_grid_mgr( this, grid, g_space, dt1, tmin1, tmax1 )

    implicit none

    class(t_boosted_grid_mgr), intent(inout) :: this
    class( t_grid ), intent(in) :: grid
    type( t_space ),     intent(in) :: g_space
    real(p_double), intent(in) :: dt1, tmin1, tmax1

    real(p_double) :: gxmin1, gxmax1

    this%tmin1 = tmin1
    this%nt1 = ceiling( (tmax1 - this%tmin1 ) / dt1 ) + 1
    this%tmax1 = this%tmin1 + dt1 * this%nt1
    this%dt1 = dt1

    gxmin1 = g_space%x_bnd_initial(p_lower,1)
    gxmax1 = g_space%x_bnd_initial(p_upper,1)

    this%nx1 = grid%my_nx(3,1)
    this%dx1 = ( gxmax1 - gxmin1 ) / grid%g_nx(1)
    this%xmin1 = gxmin1 + this%dx1 * (grid%my_nx(p_lower,1) - 1)
    this%xmax1 = this%xmin1 + this%dx1 * grid%my_nx(3,1)
    this%if_move1 = if_move(g_space, 1)
    this%move_vel1 = 1.0    ! currently only luminal moving window is supported in osiris

    ! SCR_ROOT( 'grid mgr init', this%xmin1, this%xmin2 )

    ! handle transverse dimensions
    ! global:
    this%gxperp_min = g_space%x_bnd(p_lower,2:3)
    this%gxperp_max = g_space%x_bnd(p_upper,2:3) ! xmax( g_space )(2:)
    this%gnxperp(1:2) = grid%g_nx(2:3)
    ! dxperp- no reason to track for now

    ! local:
    this%nxperp(1:2) = grid%my_nx(3, 2:3)
    this%nxperp_min(1:2) = grid%my_nx(p_lower, 2:3)

end subroutine set_space1_params_boosted_grid_mgr
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! complete initialization of this class is spread out across several functions
! because the data required to populate the class is available at different places in
! the code
!-------------------------------------------------------------------------------------------------
subroutine set_space2_params_boosted_grid_mgr( this, signed_gamma, &
                              tmin2, tmax2, dt2, xmin2, xmax2, nx2, &
                              if_move2, move_vel2 )

    implicit none

    class(t_boosted_grid_mgr), intent(inout) :: this

    real(p_double), intent(in) :: signed_gamma, tmin2, tmax2, dt2, xmin2, xmax2
    integer, intent(in) :: nx2
    logical, intent(in) :: if_move2
    real(p_double), intent(in) :: move_vel2

    integer :: nt2
    real(p_double) :: tmax2_

    if( signed_gamma < 1 .and. signed_gamma > -1 ) then
        SCR_ROOT( 'ERROR in os-boosted-grid-mgr: invalid signed_gamma')
        stop
    endif

    if( tmin2 .ge. tmax2 ) then
        SCR_ROOT( 'ERROR in os-boosted-grid-mgr: tmax2 > tmin1 is required')
        stop
    endif

    if( xmin2 .ge. xmax2 ) then
        SCR_ROOT( 'ERROR in os-boosted-grid-mgr: xmax2 > xmin1 is required')
        stop
    endif

    if( dt2 .le. 0 ) then
        SCR_ROOT( 'ERROR in os-boosted-grid-mgr: dt2 > 0 is required')
        stop
    endif

    if( nx2 .le. 0 ) then
        SCR_ROOT( 'ERROR in os-boosted-grid-mgr: nx2 > 0 is required')
        stop
    endif

    ! correct tmax2 so that dt2 gives the exact spacing
    nt2 = floor( (tmax2 - tmin2 ) / dt2 ) + 1

    ! max time that will be output
    tmax2_ = tmin2 + dt2 * (nt2 - 1)

    this%signed_gamma = signed_gamma
    this%tmin2 = tmin2
    this%tmax2 = tmax2_
    this%dt2 = dt2
    this%xmin2 = xmin2
    this%xmax2 = xmax2
    this%dx2 = (xmax2 - xmin2) / nx2
    this%if_move2 = if_move2
    this%move_vel2 = move_vel2
    this%nx2 = nx2
    this%nt2 = nt2

end subroutine set_space2_params_boosted_grid_mgr
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! Count data: true -> don't set any buffers (they aren't allocated yet),
! just loop through the data and determine how big all the buffers need
! to be
! Count data false: set all buffers
! This function is meant to be called twice, first with count_data = true
! to get the buffer sizes, and then again to actually set the buffers
!-------------------------------------------------------------------------------------------------
subroutine init_helper_boosted_grid_mgr( this, count_data )

    implicit none

    class(t_boosted_grid_mgr), intent(inout) :: this
    logical, intent(in) :: count_data

    real(p_double) :: t1, x2, t2
    real(p_double) :: t1_next, t1_prev
    real(p_double) :: xmin_current2, xmin_current1, xmax_current1
    real(p_double) :: t1_slice, x1_slice
    real(p_double) :: gamma, gamma_times_beta, beta, offset
    integer :: total_data_added, x2_idx
    integer, dimension(this%nt2) :: next_x2_idx
    integer :: next_x2_direction
    integer :: i1, i2

    total_data_added = 0

    gamma = this%gamma
    gamma_times_beta = this%gamma_times_beta
    beta = this%beta

    ! store the number of points that lie between t2 and t2 - dt2
    ! for each t2 achieved by the simulation
    if( .not. associated( this%num_t2_points_at_t1 ) ) then

        call alloc( this%num_t2_points_at_t1, (/this%nt1/) )
        this % num_t2_points_at_t1(:) = 0

        call alloc( this%starting_cell_idx_at_t1, (/this%nt1/) )
        this % starting_cell_idx_at_t1(:) = 0

    endif

    ! next_x2_idx stores the next x2 of each cell to be processed
    if( beta < 0 ) then
        next_x2_idx(:) = this%nx2
        next_x2_direction = -1
    else
        next_x2_idx(:) = 1
        next_x2_direction = 1
    endif

    ! loop through each simulation timestep
    ! start at 2 so that we have the previous data from timestep
    ! 1 ready for interpolation by the time we get to timestep 2
    ! in the simulation
    do i1 = 2, this % nt1

        t1 = this%tmin1 + (i1-1) * this%dt1
        t1_prev = this%tmin1 + (i1-2) * this%dt1
        t1_next = this%tmin1 + (i1-0) * this%dt1

        ! spatial bounds of the simulation at this t1_slice
        if( this%if_move1 ) then

            ! snap to sim grid (does not apply for lab data slices)
            ! this is an exact formula for xmin1 using the current
            ! moving window rounding scheme. if the moving window algorithm
            ! is modified then this code would need to be adjusted accordingly.
            ! the extra dt comes from the fact that the window is moved by an extra dt
            ! starting at t=0; a better moving window protocol probably would be to
            ! skip the window move at t=0. if this were implemented then the extra
            ! term of dt in the following line would be removed.

            offset = this%move_vel1 * ( t1_next - this%tmin1 )
            offset = floor( offset / this%dx1 ) * this%dx1

            xmin_current1 = this%xmin1 + offset
        else
            xmin_current1 = this%xmin1
        endif

        xmax_current1 = xmin_current1 + ( this%xmax1 - this%xmin1 )

        this%starting_cell_idx_at_t1(i1) = total_data_added + 1

        ! loop through each lab frame timeslice
        do i2 = 1, this % nt2

            ! time of this lab timeslice
            t2 = this%tmin2 + (i2-1) * this%dt2

            ! boosted_xmin is xmin at t=tmin. if using the moving window,
            ! the min x value increases with time.
            if( this%if_move2 ) then
                xmin_current2 = this%xmin2 + this%move_vel2 * t2
            else
                xmin_current2 = this%xmin2
            endif

            do while( .true. )

                ! z2 that gives the earliest time when converted to t1
                x2_idx = next_x2_idx( i2 )

                ! all grid points for this time slice have been accounted for
                if( x2_idx < 1 .or. x2_idx > this%nx2 ) exit

                x2 = xmin_current2 + (x2_idx - 1) * this%dx2

                ! coordinates of (t2,x2) transformed into the simulation frame
                ! these don't line up with either the spatial or temporal
                ! coordinates of the simulation, we will eventually need
                ! to interp onto these values
                t1_slice = gamma * t2 + gamma_times_beta * x2
                x1_slice = gamma * x2 + gamma_times_beta * t2

                ! ! useful debug comment -- keeping this here for now
                ! if( i2 .eq. 5 .and. (x2_idx .eq. 199 .or. x2_idx .eq. 200 )) then
                !     print *, 'grid mgr debug', mpi_node(), i2, x2_idx,  't:', t1 - this%dt1, t1, &
                !     t1_slice, 'x: ', xmin_current1, xmax_current1, x1_slice
                ! endif

                ! slice time is less than current slice
                ! keep cycling until we reach the first time that
                ! is within the interpolation region for this schlice
                if( t1_slice <= t1_prev ) then
                    continue

                ! all data has been added;
                ! break out of loop and wait until we get to a larger
                ! t1 to check for more data
                else if( t1_slice > t1 ) then
                    exit

                ! in this case, we found a point for lab slice t2
                ! which lies in the interpolation region ( t1 - dt1, t1 ]
                ! add the data
                else

                    ! the previous condition isn't enough to append the data, we
                    ! also need to make sure that it lies within this domain
                    ! the index still needs to be increased though to account for the
                    ! case that the slice starts outside this domain.
                    if( x1_slice >= xmin_current1 &
                            .and. x1_slice <= xmax_current1 ) then

                        total_data_added = total_data_added + 1

                        if( .not. count_data ) then

                            this%num_t2_points_at_t1(i1) = this%num_t2_points_at_t1(i1) + 1

                            ! load the cell coordinates
                            this%boosted_cell_coords( total_data_added )%t2_idx = i2
                            this%boosted_cell_coords( total_data_added )%x2_idx = x2_idx
                            this%boosted_cell_coords( total_data_added )%t1 = t1_slice
                            this%boosted_cell_coords( total_data_added )%x1 = x1_slice

                        endif
                    endif

                endif

                ! advance to the next coordinate which has a larger time
                ! if beta > 0, then we need to reduce z2_slice
                ! (i.e. lower the next index) to get a larger t1 for
                ! future iterations. opposite of beta < 0.
                next_x2_idx( i2 ) = next_x2_idx( i2 ) + next_x2_direction

            enddo

        enddo

    enddo

    this%total_data = total_data_added

end subroutine init_helper_boosted_grid_mgr
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! Count data: true -> don't set any buffers (they aren't allocated yet),
! just loop through the data and determine how big all the buffers need
! to be
! Count data false: set all buffers
! This function is meant to be called twice, first with count_data = true
! to get the buffer sizes, and then again to actually set the buffers
!-------------------------------------------------------------------------------------------------
subroutine init_helper_raw_boosted_grid_mgr( this, count_data )

    implicit none

    class(t_boosted_grid_mgr), intent(inout) :: this
    logical, intent(in) :: count_data

    real(p_double) :: t1, t2
    real(p_double) :: t1_next, t1_prev
    real(p_double) :: xmin_current2, xmax_current2, xmin_current1, xmax_current1
    real(p_double) :: tmin1, tmax1, xmin1, xmax1, xmin2, xmax2
    real(p_double) :: gamma, gamma_times_beta, beta, offset
    integer :: total_data_added, x1_idx_min, x1_idx_max
    integer :: i1, i2

    total_data_added = 0

    gamma = this%gamma
    gamma_times_beta = this%gamma_times_beta
    beta = this%beta

    ! store the number of points that lie between t2 and t2 - dt2
    ! for each t2 achieved by the simulation
    if( .not. associated( this%num_t2_points_at_t1_raw ) ) then

        call alloc( this%num_t2_points_at_t1_raw, (/this%nt1/) )
        this % num_t2_points_at_t1_raw(:) = 0

        call alloc( this%starting_raw_idx_at_t1, (/this%nt1/) )
        this % starting_raw_idx_at_t1(:) = 0

    endif

    ! loop through each simulation timestep
    ! start at 2 so that we have the previous data from timestep
    ! 1 ready for interpolation by the time we get to timestep 2
    ! in the simulation
    do i1 = 2, this % nt1

        t1 = this%tmin1 + (i1-1) * this%dt1
        t1_prev = this%tmin1 + (i1-2) * this%dt1
        t1_next = this%tmin1 + (i1-0) * this%dt1

        ! spatial bounds of the simulation at this t1_slice
        if( this%if_move1 ) then

            ! snap to sim grid (does not apply for lab data slices)
            ! this is an exact formula for xmin1 using the current
            ! moving window rounding scheme. if the moving window algorithm
            ! is modified then this code would need to be adjusted accordingly.
            ! the extra dt comes from the fact that the window is moved by an extra dt
            ! starting at t=0; a better moving window protocol probably would be to
            ! skip the window move at t=0. if this were implemented then the extra
            ! term of dt in the following line would be removed.

            offset = this%move_vel1 * ( t1_next - this%tmin1 )
            offset = floor( offset / this%dx1 ) * this%dx1

            xmin_current1 = this%xmin1 + offset
        else
            xmin_current1 = this%xmin1
        endif

        xmax_current1 = xmin_current1 + ( this%xmax1 - this%xmin1 )

        this%starting_raw_idx_at_t1(i1) = total_data_added + 1

        ! loop through each lab frame timeslice
        do i2 = 1, this % nt2

            ! time of this lab timeslice
            t2 = this%tmin2 + (i2-1) * this%dt2

            ! boosted_xmin is xmin at t=tmin. if using the moving window,
            ! the min x value increases with time.
            if( this%if_move2 ) then
                xmin_current2 = this%xmin2 + this%move_vel2 * t2
            else
                xmin_current2 = this%xmin2
            endif

            xmax_current2 = xmin_current2 + this%nx2 * this%dx2

            ! compute the min/max simulation times corresponding to this diagnostic
            if ( beta < 0 ) then
                tmin1 = gamma * t2 + gamma_times_beta * xmax_current2
                tmax1 = gamma * t2 + gamma_times_beta * xmin_current2
            else
                tmin1 = gamma * t2 + gamma_times_beta * xmin_current2
                tmax1 = gamma * t2 + gamma_times_beta * xmax_current2
            endif

            ! check if the time is in bounds
            if ( t1 >= tmin1 .and. t1_prev < tmax1 ) then

                ! check if the position is in bounds
                ! find the xmin2 and xmax2 that correspond to these t1 times
                if ( beta < 0 ) then
                    xmin2 = ( t1      / gamma - t2 ) / beta
                    xmax2 = ( t1_prev / gamma - t2 ) / beta
                else
                    xmin2 = ( t1_prev / gamma - t2 ) / beta
                    xmax2 = ( t1      / gamma - t2 ) / beta
                endif

                ! find the xmin1 and xmax1 that correspond to this range
                xmin1 = gamma * xmin2 + gamma_times_beta * t2
                xmax1 = gamma * xmax2 + gamma_times_beta * t2

                ! compute the index for the x1 bounds
                ! we add dt1 here because particles can move that far in a time step
                x1_idx_min = floor( (xmin1 - this%dt1 - xmin_current1) / this%dx1 ) + 1
                x1_idx_max = floor( (xmax1 + this%dt1 - xmin_current1) / this%dx1 ) + 1

                if ( x1_idx_min <= this%nx1 .and. x1_idx_max >= 1 ) then

                    total_data_added = total_data_added + 1

                    if ( .not. count_data ) then
                        this%num_t2_points_at_t1_raw(i1) = this%num_t2_points_at_t1_raw(i1) + 1

                        this%boosted_raw_coords( total_data_added )%t2_idx = i2
                        this%boosted_raw_coords( total_data_added )%t2 = t2
                        this%boosted_raw_coords( total_data_added )%xmin2 = xmin2
                        this%boosted_raw_coords( total_data_added )%xmax2 = xmax2
                        this%boosted_raw_coords( total_data_added )%x1_idx_min = x1_idx_min
                        this%boosted_raw_coords( total_data_added )%x1_idx_max = x1_idx_max
                    endif

                endif

            endif

        enddo

    enddo

    this%total_data_raw = total_data_added

end subroutine init_helper_raw_boosted_grid_mgr
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine cleanup_boosted_grid_mgr( this )

    implicit none

    class(t_boosted_grid_mgr), intent(inout) :: this

    call freemem( this%num_t2_points_at_t1 )
    call freemem( this%starting_cell_idx_at_t1 )
    call freemem( this%num_t2_points_at_t1_raw )
    call freemem( this%starting_raw_idx_at_t1 )
    if (associated(this%boosted_cell_coords)) deallocate( this%boosted_cell_coords )
    if (associated(this%boosted_raw_coords)) deallocate( this%boosted_raw_coords )

end subroutine cleanup_boosted_grid_mgr
!-------------------------------------------------------------------------------------------------



end module m_boosted_grid_mgr