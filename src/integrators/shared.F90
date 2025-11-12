!> Module with shared variables and routines for integrators
module shared
    implicit none

    ! Constants
    real(kind=8), parameter :: ZERO = 0.0
    real(kind=8), parameter :: ONE = 1.d0
    real(kind=8), parameter :: TWO = 2.d0
    real(kind=8), parameter :: C12 = 0.5d0
    real(kind=8), parameter :: C13 = ONE/3.d0
    real(kind=8), parameter :: C14 = 0.25d0
    real(kind=8), parameter :: SAFE_LOW = 1d-30

    ! BASIC CONFIGURATION
    integer(kind=4) :: NDIM = 1  ! Number of dimensions
    integer(kind=4) :: NDIM2

    ! Global
    real(kind=8), allocatable :: der(:)
    real(kind=8) :: E_TOL = 1.d-14
    real(kind=8) :: BETA = 0.15d0
    real(kind=8) :: DT_MIN = 1.d-6
    integer(kind=4) :: EXTRA = 0  ! Amount of extra variables (without der) that are not POS
    integer(kind=4) :: EXTRA2

    abstract interface

        !!!! TEMPLATES

        ! Keep going template
        function function_check_keep_tem (y) result(keep_going)
            implicit none
            real(kind=8), dimension(:), intent(in) :: y
            logical :: keep_going
        end function function_check_keep_tem

        ! f__ (t, y__) = (f_i (t, y__), ..., f_n (t, y__)) = der__
        function dydt_tem (t, y) result (derivate)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size (y))      :: derivate
        end function dydt_tem

    end interface
    
    contains

            subroutine init_global (sizey, ndimensions, n_extra)
                implicit none
                integer(kind=4), intent(in) :: sizey, ndimensions
                integer(kind=4), intent(in), optional :: n_extra

                NDIM = ndimensions
                NDIM2 = NDIM * 2

                allocate(der(sizey))

                if (present(n_extra)) EXTRA = n_extra
                EXTRA2 = EXTRA * 2

            end subroutine init_global

            pure function get_index(i) result(idx)
                implicit none
                integer(kind=4), intent(in) :: i
                integer(kind=4) :: idx
                idx = EXTRA2 + NDIM2 * (i - 1) + 1
            end function get_index

    
end module shared