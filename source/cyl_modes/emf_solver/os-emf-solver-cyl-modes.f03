module m_emf_solver_cyl_modes

private

integer, parameter, public :: p_emf_cyl_modes_yee = 0
integer, parameter, public :: p_emf_cyl_modes_fei = 6

interface new_solver_emf_cyl_modes
    module procedure new_solver_emf_cyl_modes
end interface

public :: new_solver_emf_cyl_modes

contains

subroutine new_solver_emf_cyl_modes( this, name )

    use m_emf_define, only : t_emf
    use m_emf_solver_cyl_modes_yee
    use m_emf_solver_cyl_modes_fei
    use m_system

    implicit none

    class( t_emf ), intent(inout) :: this
    character( len = * ), intent(in) :: name

    select case( trim(name) )
    case( "yee" )
        allocate( t_emf_solver_cyl_modes_yee :: this%solver )
    case( "custom", "fei" )
        allocate( t_emf_solver_cyl_modes_fei :: this%solver )
    case default
        if ( mpi_node() == 0 ) then
            write(0,*) "Invalid EMF solver type: '", trim(name), "'"
        endif
        this%solver => null()
    end select

end subroutine new_solver_emf_cyl_modes

end module m_emf_solver_cyl_modes
