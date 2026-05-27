!> Module with fixed RK integrators
module runge_kutta
    use shared
    
    implicit none
    private
    public :: init_runge_kutta, free_runge_kutta, runge_kutta_caller, runge_kutta_fixed_caller

    ! Workspace arrays
    real(wp), allocatable :: rk(:, :)  ! integrator
    real(wp), allocatable :: rk_imp1(:)  ! implicit 1k integrators
    real(wp), allocatable :: rk_imp_solv(:, :)  ! implicit solver
    real(wp), allocatable :: yscal(:)  ! solver
    real(wp), allocatable :: yaux(:)  ! solver
    real(wp), allocatable :: yhalf(:), derhalf(:)  ! solver
    real(wp), allocatable :: ycaller(:)  ! caller

    ! Constants
    integer(kind=4), parameter :: MAX_N_ITER = 350
    real(wp), parameter :: MAX_DT_FACTOR = 5.0e0_wp

    ! Variables to set
    integer(kind=4) :: OSOL
    integer(kind=4) :: N_STAGES
    integer(kind=4) :: N_IMPLICIT

    ! Derived
    real(wp) :: ONE_OSOL
    real(wp) :: ONE__OSOL_PLUS_ONE
    real(wp) :: ONE__MINUS_ONE_PLUS_TWO_TO_ORD

    ! Pointer to runge_kutta used
    procedure(runge_kutta_tem), pointer :: runge_kutta_ptr => null()

    abstract interface

        subroutine runge_kutta_tem(sizey, y, dydt, t, dt, deri, ynew)
            import :: wp
            import :: dydt_tem
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(wp), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: t, dt
            real(wp), dimension(sizey), intent(in) :: deri
            real(wp), dimension(sizey), intent(out) :: ynew
        end subroutine runge_kutta_tem

        subroutine implicit_funK_tem(sizey, y, dydt, t, dt, sizerk, kin, kout)
            import :: wp
            import :: dydt_tem
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(wp), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: t, dt
            integer(kind=4), intent(in) :: sizerk
            real(wp), dimension(sizey, sizerk), intent(in) :: kin
            real(wp), dimension(sizey, sizerk), intent(out) :: kout
        end subroutine implicit_funK_tem

    end interface

contains

    !!!! HANDLER

    subroutine init_runge_kutta(sizey, which)
        implicit none
        integer(kind=4), intent(in) :: sizey
        integer(kind=4), intent(in) :: which

        N_IMPLICIT = 0

        if (which == 1) then
            runge_kutta_ptr => Euler1
            OSOL = 1
            N_STAGES = 1

        else if (which == 2) then
            runge_kutta_ptr => Euler_back1 ! Implicit
            OSOL = 1
            N_STAGES = 1
            N_IMPLICIT = 1

        else if (which == 3) then
            runge_kutta_ptr => Euler_center2 ! Implicit
            OSOL = 2
            N_STAGES = 1
            N_IMPLICIT = 1

        else if (which == 4) then
            runge_kutta_ptr => Crank_Nicolson2 ! Implicit
            OSOL = 2
            N_STAGES = 2
            N_IMPLICIT = 1

        else if (which == 5) then
            runge_kutta_ptr => Heun2
            OSOL = 2
            N_STAGES = 2

        else if (which == 6) then
            runge_kutta_ptr => midpoint2
            OSOL = 2
            N_STAGES = 2

        else if (which == 7) then
            runge_kutta_ptr => strange2
            OSOL = 2
            N_STAGES = 2

        else if (which == 8) then
            runge_kutta_ptr => Ralston2
            OSOL = 2
            N_STAGES = 2

        else if (which == 9) then
            runge_kutta_ptr => Hammer_Hollingsworth2 ! Implicit
            OSOL = 2
            N_STAGES = 2
            N_IMPLICIT = 1

        else if (which == 10) then
            runge_kutta_ptr => Kraaijevanger_Spijker2 ! Implicit
            OSOL = 2
            N_STAGES = 2
            N_IMPLICIT = 1

        else if (which == 11) then
            runge_kutta_ptr => Qin_Zhang2 ! Implicit
            OSOL = 2
            N_STAGES = 2
            N_IMPLICIT = 1

        else if (which == 12) then
            runge_kutta_ptr => Runge_Kutta3
            OSOL = 3
            N_STAGES = 3

        else if (which == 13) then
            runge_kutta_ptr => Heun3
            OSOL = 3
            N_STAGES = 3

        else if (which == 14) then
            runge_kutta_ptr => Ralston3
            OSOL = 3
            N_STAGES = 3

        else if (which == 15) then
            runge_kutta_ptr => SSPRrk3
            OSOL = 3
            N_STAGES = 3

        else if (which == 16) then
            runge_kutta_ptr => Crouzeix3 ! Implicit
            OSOL = 3
            N_STAGES = 2
            N_IMPLICIT = 1

        else if (which == 17) then
            runge_kutta_ptr => Runge_Kutta_implicit3 ! Implicit
            OSOL = 3
            N_STAGES = 4
            N_IMPLICIT = 1

        else if (which == 18) then
            runge_kutta_ptr => Ralston4 ! Implicit
            OSOL = 4
            N_STAGES = 4
            N_IMPLICIT = 1

        else if (which == 19) then
            runge_kutta_ptr => Lobatto4 ! Implicit
            OSOL = 4
            N_STAGES = 3
            N_IMPLICIT = 1

        else if (which == 20) then
            runge_kutta_ptr => Runge_Kutta4
            OSOL = 4
            N_STAGES = 4

        else if (which == 21) then
            runge_kutta_ptr => Gauss_Legendre4 ! Implicit
            OSOL = 4
            N_STAGES = 2
            N_IMPLICIT = 2

        else if (which == 22) then
            runge_kutta_ptr => Runge_Kutta_four_oct4
            OSOL = 4
            N_STAGES = 4

        else if (which == 23) then
            runge_kutta_ptr => Nystrom5
            OSOL = 5
            N_STAGES = 6

        else if (which == 24) then
            runge_kutta_ptr => Gauss_Legendre6 ! Implicit
            OSOL = 6
            N_STAGES = 3
            N_IMPLICIT = 3

        else if (which == 25) then
            runge_kutta_ptr => Runge_Kutta6
            OSOL = 6
            N_STAGES = 7

        else if (which == 26) then
            runge_kutta_ptr => Abbas6
            OSOL = 6
            N_STAGES = 7

        else
            print *, "Method not available:", which
            stop 1
        end if

        ! Constants
        ONE_OSOL = ONE/real(OSOL, kind=wp)
        ONE__OSOL_PLUS_ONE = ONE/(OSOL + ONE)
        ONE__MINUS_ONE_PLUS_TWO_TO_ORD = ONE/(TWO**OSOL - ONE)

        ! integrator
        allocate (rk(sizey, N_STAGES))

        ! implicit integrators
        if (N_IMPLICIT == 1) allocate (rk_imp1(sizey))
        if (N_IMPLICIT > 0) allocate (rk_imp_solv(sizey, N_IMPLICIT))

        ! solver
        allocate (yscal(sizey))
        allocate (yaux(sizey))
        allocate (yhalf(sizey), derhalf(sizey))

        ! caller
        allocate (ycaller(sizey))
    end subroutine init_runge_kutta

    subroutine free_runge_kutta()
        implicit none
        if (allocated(rk)) deallocate (rk)
        if (allocated(rk_imp1)) deallocate (rk_imp1)
        if (allocated(rk_imp_solv)) deallocate (rk_imp_solv)
        if (allocated(yscal)) deallocate (yscal)
        if (allocated(yaux)) deallocate (yaux)
        if (allocated(yhalf)) deallocate (yhalf)
        if (allocated(derhalf)) deallocate (derhalf)
        if (allocated(ycaller)) deallocate (ycaller)
        nullify (runge_kutta_ptr)
    end subroutine free_runge_kutta

    !!!! Runge Kutta subroutines

    subroutine Euler1(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            ynew(i) = y(i) + dt*deri(i)
        end do

    end subroutine Euler1

    subroutine Euler_back1(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew

        rk_imp1 = ZERO
        call solve_1k_implicit(sizey, y, dydt, t + dt, dt, rk_imp1, ONE, rk(:, 1))

        ynew = y + dt*rk(:, 1)

    end subroutine Euler_back1

    subroutine Euler_center2(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew

        rk_imp1 = ZERO
        call solve_1k_implicit(sizey, y, dydt, t + dt*C1_2, dt, rk_imp1, C1_2, rk(:, 1))

        ynew = y + dt*rk(:, 1)

    end subroutine Euler_center2

    subroutine Crank_Nicolson2(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            rk_imp1(i) = deri(i)*C1_2
        end do
        call solve_1k_implicit(sizey, y, dydt, t + dt, dt, rk_imp1, C1_2, rk(:, 2))

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + rk(i, 2))*C1_2
        end do

    end subroutine Crank_Nicolson2

    subroutine Heun2(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)
        end do
        rk(:, 2) = dydt(t + dt, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + rk(i, 2))*C1_2
        end do

    end subroutine Heun2

    subroutine midpoint2(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)*C1_2
        end do
        rk(:, 2) = dydt(t + dt*C1_2, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*rk(i, 2)
        end do

    end subroutine midpoint2

    subroutine strange2(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)*C3_4
        end do
        rk(:, 2) = dydt(t + dt*C3_4, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + rk(i, 2)*TWO)*C1_3
        end do

    end subroutine strange2

    subroutine Ralston2(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)*C2_3
        end do
        rk(:, 2) = dydt(t + dt*C2_3, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + rk(i, 2)*THREE)*C1_4
        end do

    end subroutine Ralston2

    subroutine Hammer_Hollingsworth2(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            rk_imp1(i) = deri(i)*C1_3
        end do
        call solve_1k_implicit(sizey, y, dydt, t + dt*C2_3, dt, rk_imp1, C1_3, rk(:, 2))

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + rk(i, 2)*THREE)*C1_4
        end do

    end subroutine Hammer_Hollingsworth2

    subroutine Kraaijevanger_Spijker2(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        rk_imp1 = ZERO
        call solve_1k_implicit(sizey, y, dydt, t + dt*C1_2, dt, rk_imp1, C1_2, rk(:, 1))

        do i = 1, sizey
            rk_imp1(i) = -rk(i, 1)*C1_2
        end do
        call solve_1k_implicit(sizey, y, dydt, t + dt*C3_2, dt, rk_imp1, TWO, rk(:, 2))

        do i = 1, sizey
            ynew(i) = y(i) + dt*(-rk(i, 1) + rk(i, 2)*THREE)*C1_2
        end do

    end subroutine Kraaijevanger_Spijker2

    subroutine Qin_Zhang2(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        rk_imp1 = ZERO
        call solve_1k_implicit(sizey, y, dydt, t + dt*C1_4, dt, rk_imp1, C1_4, rk(:, 1))

        do i = 1, sizey
            rk_imp1(i) = rk(i, 1)*C1_2
        end do
        call solve_1k_implicit(sizey, y, dydt, t + dt*C3_4, dt, rk_imp1, C1_4, rk(:, 2))

        do i = 1, sizey
            ynew(i) = y(i) + dt*(rk(i, 1) + rk(i, 2))*C1_2
        end do

    end subroutine Qin_Zhang2

    subroutine Runge_Kutta3(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)*C1_2
        end do
        rk(:, 2) = dydt(t + dt*C1_2, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(-deri(i) + rk(i, 2)*TWO)
        end do
        rk(:, 3) = dydt(t + dt, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + rk(i, 2)*FOUR + rk(i, 3))*C1_6
        end do

    end subroutine Runge_Kutta3

    subroutine Heun3(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)*C1_3
        end do
        rk(:, 2) = dydt(t + dt*C1_3, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*rk(i, 2)*C2_3
        end do
        rk(:, 3) = dydt(t + dt*C2_3, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + rk(i, 3)*THREE)*C1_4
        end do

    end subroutine Heun3

    subroutine Ralston3(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)*C1_2
        end do
        rk(:, 2) = dydt(t + dt*C1_2, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*rk(i, 2)*C3_4
        end do
        rk(:, 3) = dydt(t + dt*C3_4, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i)*TWO + rk(i, 2)*THREE + rk(i, 3)*FOUR)*C1_9
        end do

    end subroutine Ralston3

    subroutine SSPRrk3(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)
        end do
        rk(:, 2) = dydt(t + dt, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(deri(i) + rk(i, 2))*C1_4
        end do
        rk(:, 3) = dydt(t + dt*C1_2, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + rk(i, 2) + rk(i, 3)*FOUR)*C1_6
        end do

    end subroutine SSPRrk3

    subroutine Crouzeix3(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        real(wp), parameter :: aux = C1_2 + SQ3_6
        integer(kind=4) :: i

        rk_imp1 = ZERO
        call solve_1k_implicit(sizey, y, dydt, t + dt*aux, dt, rk_imp1, aux, rk(:, 1))

        do i = 1, sizey
            rk_imp1(i) = -rk(i, 1)*SQ3_6*2
        end do
        call solve_1k_implicit(sizey, y, dydt, t + dt*(C1_2 - SQ3_6), dt, rk_imp1, aux, rk(:, 2))

        do i = 1, sizey
            ynew(i) = y(i) + dt*(rk(i, 1) + rk(i, 2))*C1_2
        end do

    end subroutine Crouzeix3

    subroutine Runge_Kutta_implicit3(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        rk_imp1 = ZERO
        call solve_1k_implicit(sizey, y, dydt, t + dt*C1_2, dt, rk_imp1, C1_2, rk(:, 1))

        do i = 1, sizey
            rk_imp1(i) = rk(i, 1)*C1_6
        end do
        call solve_1k_implicit(sizey, y, dydt, t + dt*C2_3, dt, rk_imp1, C1_2, rk(:, 2))

        do i = 1, sizey
            rk_imp1(i) = (-rk(i, 1) + rk(i, 2))*C1_2
        end do
        call solve_1k_implicit(sizey, y, dydt, t + dt*C1_2, dt, rk_imp1, C1_2, rk(:, 3))

        do i = 1, sizey
            rk_imp1(i) = ((rk(i, 1) - rk(i, 2))*THREE + rk(i, 3))*C1_2
        end do
        call solve_1k_implicit(sizey, y, dydt, t + dt, dt, rk_imp1, C1_2, rk(:, 4))

        do i = 1, sizey
            ynew(i) = y(i) + dt*(rk_imp1(i) + rk(i, 4)*C1_2)
        end do

    end subroutine Runge_Kutta_implicit3

    subroutine Ralston4(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)*C2_5
        end do
        rk(:, 2) = dydt(t + dt*C2_5, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(deri(i)*0.29697761e0_wp + rk(i, 2)*0.15875964e0_wp)
        end do
        rk(:, 3) = dydt(t + dt*0.45573725e0_wp, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(deri(i)*0.21810040e0_wp - rk(i, 2)*3.05096516e0_wp + &
                                            & rk(i, 3)*3.83286476e0_wp)
        end do
        rk(:, 4) = dydt(t + dt, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i)*0.17476028e0_wp - rk(i, 2)*0.55148066e0_wp + rk(i, 3)*1.20553560e0_wp + &
                                            & rk(i, 4)*0.17118478e0_wp)
        end do

    end subroutine Ralston4

    subroutine Lobatto4(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            rk_imp1(i) = deri(i)*C1_4
        end do
        call solve_1k_implicit(sizey, y, dydt, t + dt*C1_2, dt, rk_imp1, C1_4, rk(:, 2))

        do i = 1, sizey
            yaux(i) = y(i) + dt*rk(i, 2)
        end do
        rk(:, 3) = dydt(t + dt, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + rk(i, 2)*FOUR + rk(i, 3))*C1_6
        end do

    end subroutine Lobatto4

    subroutine Runge_Kutta4(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)*C1_2
        end do
        rk(:, 2) = dydt(t + dt*C1_2, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*rk(i, 2)*C1_2
        end do
        rk(:, 3) = dydt(t + dt*C1_2, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*rk(i, 3)
        end do
        rk(:, 4) = dydt(t + dt, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + (rk(i, 2) + rk(i, 3))*TWO + rk(i, 4))*C1_6
        end do

    end subroutine Runge_Kutta4

    subroutine Gauss_Legendre4(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            rk(i, 1) = ZERO
            rk(i, 2) = ZERO
        end do
        call solve_rk_implicit(sizey, N_STAGES, y, dydt, t, dt, FunK_GL4, rk)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(rk(i, 1) + rk(i, 2))*C1_2
        end do

    contains
        subroutine FunK_GL4(sizey, y, dydt, t, dt, sizerk, kin, kout) !! Funk for Gauss_Legendre4
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(wp), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: t, dt
            integer(kind=4), intent(in) :: sizerk
            real(wp), dimension(sizey, sizerk), intent(in) :: kin
            real(wp), dimension(sizey, sizerk), intent(out) :: kout
            integer(kind=4) :: i

            do i = 1, sizey
                yaux(i) = y(i) + dt*(kin(i, 1)*C1_4 + kin(i, 2)*(C1_4 - SQ3_6))
            end do
            kout(:, 1) = dydt(t + dt*(C1_2 - SQ3_6), yaux)

            do i = 1, sizey
                yaux(i) = y(i) + dt*(kin(i, 1)*(C1_4 + SQ3_6) + kin(i, 2)*C1_4)
            end do
            kout(:, 2) = dydt(t + dt*(C1_2 + SQ3_6), yaux)
        end subroutine FunK_GL4

    end subroutine Gauss_Legendre4

    subroutine Runge_Kutta_four_oct4(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            rk(i, 1) = deri(i)*C1_3
            yaux(i) = y(i) + dt*rk(i, 1)
        end do
        rk(:, 2) = dydt(t + dt*C1_3, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(-rk(i, 1) + rk(i, 2))
        end do
        rk(:, 3) = dydt(t + dt*C2_3, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(deri(i) - rk(i, 2) + rk(i, 3))
        end do
        rk(:, 4) = dydt(t + dt, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(deri(i) + (rk(i, 2) + rk(i, 3))*THREE + rk(i, 4))*C1_8
        end do

    end subroutine Runge_Kutta_four_oct4

    subroutine Nystrom5(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*C1_3*deri(i)
        end do
        rk(:, 2) = dydt(t + dt*C1_3, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*C1_25*(FOUR*deri(i) + 6.e0_wp*rk(i, 2))
        end do
        rk(:, 3) = dydt(t + dt*C2_5, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(C1_4*deri(i) - THREE*rk(i, 2) + C15_4*rk(i, 3))
        end do
        rk(:, 4) = dydt(t + dt, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*C1_9*(C2_3*deri(i) + 10.e0_wp*rk(i, 2) + &
                                            & C1_9*(-50.e0_wp*rk(i, 3) + 8.e0_wp*rk(i, 4)))
        end do
        rk(:, 5) = dydt(t + dt*C2_3, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(C1_25*(TWO*deri(i) + 12.e0_wp*rk(i, 2)) + &
                                            & C2_15*rk(i, 3) + C8_75*rk(i, 4))
        end do
        rk(:, 6) = dydt(t + dt*C4_5, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*C1_192*(23.e0_wp*deri(i) + 125.e0_wp*(rk(i, 3) + rk(i, 6)) - 81.e0_wp*rk(i, 5))
        end do

    end subroutine Nystrom5

    subroutine Gauss_Legendre6(sizey, y, dydt, t, dt, deri, ynew) ! Implicit
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            rk(i, 1) = ZERO
            rk(i, 2) = ZERO
            rk(i, 3) = ZERO
        end do
        call solve_rk_implicit(sizey, N_STAGES, y, dydt, t, dt, FunK_GL6, rk)

        do i = 1, sizey
            ynew(i) = y(i) + dt*((rk(i, 1) + rk(i, 3))*FIVE + rk(i, 2)*8.e0_wp)*C1_18
        end do

    contains
        subroutine FunK_GL6(sizey, y, dydt, t, dt, sizerk, kin, kout) !! Funk for Gauss_Legendre6
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(wp), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: t, dt
            integer(kind=4), intent(in) :: sizerk
            real(wp), dimension(sizey, sizerk), intent(in) :: kin
            real(wp), dimension(sizey, sizerk), intent(out) :: kout
            integer(kind=4) :: i

            do i = 1, sizey
                yaux(i) = y(i) + dt*(kin(i, 1)*C5_36 + kin(i, 2)*(C2_3 - SQ15_5)*C1_3 + kin(i, 3)*(C5_36 - SQ15_30))
            end do
            kout(:, 1) = dydt(t + dt*(ONE - SQ15_5)*C1_2, yaux)

            do i = 1, sizey
                yaux(i) = y(i) + dt*(kin(i, 1)*(C5_36 + SQ15_24) + kin(i, 2)*C2_9 + kin(i, 3)*(C5_36 - SQ15_24))
            end do
            kout(:, 2) = dydt(t + dt*C1_2, yaux)

            do i = 1, sizey
                yaux(i) = y(i) + dt*(kin(i, 1)*(C5_36 + SQ15_30) + kin(i, 2)*(C2_9 + SQ15_15) + kin(i, 3)*C5_36)
            end do
            kout(:, 3) = dydt(t + dt*(ONE + SQ15_5)*C1_2, yaux)
        end subroutine FunK_GL6

    end subroutine Gauss_Legendre6

    subroutine Runge_Kutta6(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*C1_4*deri(i)
        end do
        rk(:, 2) = dydt(t + dt*C1_4, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*C1_8*(deri(i) + rk(i, 2))
        end do
        rk(:, 3) = dydt(t + dt*C1_4, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(-C5_6*rk(i, 2) + C4_3*rk(i, 3))
        end do
        rk(:, 4) = dydt(t + dt*C1_2, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(C1_8*(deri(i) + rk(i, 2)) + C1_2*rk(i, 4))
        end do
        rk(:, 5) = dydt(t + dt*C3_4, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(C3_8*rk(i, 2) + C1_4*(rk(i, 3) + rk(i, 5)) - &
                                                & C1_8*rk(i, 4))
        end do
        rk(:, 6) = dydt(t + dt*C3_4, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(C1_7*deri(i) - C2_7*rk(i, 2) + C4_7*(rk(i, 3) + rk(i, 6)))
        end do
        rk(:, 7) = dydt(t + dt, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*(C7_90*(deri(i) + rk(i, 7)) + C32_90*rk(i, 3) + C12_90*rk(i, 4) + &
                                                & C16_90*(rk(i, 5) + rk(i, 6)))
        end do

    end subroutine Runge_Kutta6

    subroutine Abbas6(sizey, y, dydt, t, dt, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew
        integer(kind=4) :: i

        do i = 1, sizey
            yaux(i) = y(i) + dt*deri(i)*C1_3
        end do
        rk(:, 2) = dydt(t + dt*C1_3, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*rk(i, 2)*C2_3
        end do
        rk(:, 3) = dydt(t + dt*C2_3, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(deri(i) + rk(i, 2)*FOUR - rk(i, 3))*C1_12
        end do
        rk(:, 4) = dydt(t + dt*C1_3, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(deri(i)*25.e0_wp - rk(i, 2)*110.e0_wp + &
                                                    & rk(i, 3)*35.e0_wp + rk(i, 4)*90.e0_wp)/48.e0_wp
        end do
        rk(:, 5) = dydt(t + dt*C5_6, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(deri(i)*0.15e0_wp - rk(i, 2)*0.55e0_wp - rk(i, 3)*C1_8 + &
                                                    & rk(i, 4)*C1_2 + rk(i, 5)*0.1e0_wp)
        end do
        rk(:, 6) = dydt(t + dt*C1_6, yaux)

        do i = 1, sizey
            yaux(i) = y(i) + dt*(-deri(i)*195.75e0_wp + rk(i, 2)*495.e0_wp + &
                                                    & rk(i, 3)*53.75e0_wp - rk(i, 4)*590.e0_wp + &
                                                    & rk(i, 5)*32.e0_wp + &
        & rk(i, 6)*400.e0_wp)/195.e0_wp
        end do
        rk(:, 7) = dydt(t + dt, yaux)

        do i = 1, sizey
            ynew(i) = y(i) + dt*((deri(i) + rk(i, 7))*13.e0_wp + (rk(i, 3) + rk(i, 4))*55.e0_wp + (rk(i, 5) + &
                                                        & rk(i, 6))*32.e0_wp)*5.e-3_wp
        end do

    end subroutine Abbas6

    !------------------------------------------------
    !  Solvers Runge Kutta (Implicit and Half Step)
    !------------------------------------------------

    !!!! Implicit: 1D SOLVER
    subroutine solve_1k_implicit(sizey, y, dydt, t, dt, kprev, cte, kout)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: kprev
        real(wp), intent(in) :: cte
        real(wp), dimension(sizey), intent(out) :: kout
        integer(kind=4) :: i, iter

        kout = dydt(t, y + dt*kprev)
        do iter = 1, MAX_N_ITER
            rk_imp_solv(:, 1) = kout
            kout = dydt(t, y + dt*(kprev + cte*kout))
            if (maxval(abs((rk_imp_solv(:, 1) - kout)/(rk_imp_solv(:, 1) + SAFE_LOW))) .le. E_TOL) then
                exit
            end if
        end do
    end subroutine solve_1k_implicit

    !!!! Implicit: ND SOLVER
    subroutine solve_rk_implicit(sizey, nstages, y, dydt, t, dt, impl_funK, rkout)
        implicit none
        integer(kind=4), intent(in) :: sizey, nstages
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        procedure(implicit_funK_tem) :: impl_funK
        real(wp), dimension(sizey, nstages), intent(inout) :: rkout
        integer(kind=4) :: iter

        do iter = 1, MAX_N_ITER
            rk_imp_solv(:, 1:nstages) = rkout
            call impl_funK(sizey, y, dydt, t, dt, nstages, rk_imp_solv(:, 1:nstages), rkout)
            if (maxval(abs((rk_imp_solv(:, 1:nstages) - rkout)/ &
                         & (rk_imp_solv(:, 1:nstages) + SAFE_LOW))) .le. E_TOL) then
                exit
            end if
        end do
    end subroutine solve_rk_implicit

    !!!! Adaptive Timestep: Half step
    subroutine solve_rk_half_step(sizey, y, dydt, t, dt_adap, dt_used, deri, runge_kutta, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t
        real(wp), intent(inout) :: dt_adap
        real(wp), intent(out) :: dt_used
        real(wp), dimension(sizey), intent(in) :: deri
        procedure(runge_kutta_tem), pointer :: runge_kutta
        real(wp), dimension(sizey), intent(out) :: ynew

        real(wp) :: e_calc, ratio, hdt_adap, dt_try
        integer(kind=4) :: iter, i

        iter = 0
        dt_try = dt_adap

        do
            iter = iter + 1
            dt_try = max(dt_try, DT_MIN_NOW)
            hdt_adap = C1_2*dt_try

            ! yscal
            do i = 1, sizey
                yscal(i) = abs(y(i)) + abs(dt_try*deri(i)) + SAFE_LOW
            end do

            ! y(t, dt; der) -> ynew
            call runge_kutta(sizey, y, dydt, t, dt_try, deri, ynew)

            ! y(t, dt/2; der) -> yhalf
            call runge_kutta(sizey, y, dydt, t, hdt_adap, deri, yhalf)

            ! d[yhalf (t + dt/2)]/dt -> derhalf
            derhalf = dydt(t + hdt_adap, yhalf)

            ! yhalf(t + dt/2, dt/2; derhalf) -> yaux
            call runge_kutta(sizey, yhalf, dydt, t + hdt_adap, hdt_adap, derhalf, yaux)

            ! Error
            e_calc = ZERO
            do i = 1, sizey
                e_calc = max(e_calc, abs((ynew(i) - yaux(i))/yscal(i))*ONE__MINUS_ONE_PLUS_TWO_TO_ORD)
            end do
            e_calc = max(e_calc, SAFE_LOW)
            ratio = E_TOL/e_calc

            if (ratio > ONE) then
                dt_used = dt_try
                dt_adap = dt_try*min(BETA*ratio**ONE__OSOL_PLUS_ONE, MAX_DT_FACTOR)
                return

            else
                if (abs(dt_try - DT_MIN_NOW) .le. E_TOL) then
                    dt_used = DT_MIN_NOW
                    dt_adap = dt_try
                    return

                else
                    dt_try = dt_try*min(BETA*ratio**ONE_OSOL, MAX_DT_FACTOR)

                    if ((dt_try /= dt_try) .or. (dt_try .le. DT_MIN_NOW) .or. (iter == MAX_N_ITER)) then
                        dt_used = DT_MIN_NOW
                        dt_adap = DT_MIN_NOW
                        call runge_kutta(sizey, y, dydt, t, dt_try, deri, ynew)
                        return
                    end if
                end if
            end if
        end do

    end subroutine solve_rk_half_step

    ! ------------------------------------------------
    !  Call Runge Kutta integrator (adaptive timestep)
    ! ------------------------------------------------

    subroutine runge_kutta_caller(t, y, dt_adap, dydt, dt, ynew, check_fun)
        implicit none
        real(wp), intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), intent(inout) :: dt_adap  ! This try and also next try
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: dt
        real(wp), dimension(size(y)), intent(out) :: ynew
        procedure(function_check_keep_tem), optional :: check_fun

        integer(kind=4) :: sizey
        real(wp) :: time, t_end, dt_used
        logical :: keep = .True.
        logical :: has_check = .False.

        sizey = size(y)
        has_check = present(check_fun)

        ynew = y
        time = t
        t_end = time + dt

        do while (time < t_end)

            if (has_check) then ! If Check Continue function present
                keep = check_fun(ynew)
                if (.not. keep) then ! If Hard Exit is True
                    dt_adap = time - t ! Replace dt_adap with actual dt used
                    return ! Exit subroutine
                end if
            end if

            ycaller = ynew

            dt_adap = min(dt_adap, t_end - time)
            DT_MIN_NOW = min(DT_MIN, dt_adap)

            der = dydt(time, ycaller)

            call solve_rk_half_step(sizey, ycaller, dydt, time, dt_adap, dt_used, der, runge_kutta_ptr, ynew)

            time = time + dt_used
        end do

    end subroutine runge_kutta_caller

    ! ------------------------------------------------
    !  Call Runge Kutta integrator (fixed timestep)
    ! ------------------------------------------------

    subroutine runge_kutta_fixed_caller(t, y, dt_adap, dydt, dt, ynew, check_fun)
        implicit none
        real(wp), intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), intent(inout) :: dt_adap  ! This is each sub-step
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: dt ! This is full step
        real(wp), dimension(size(y)), intent(out) :: ynew
        procedure(function_check_keep_tem), optional :: check_fun

        integer(kind=4) :: sizey
        real(wp) :: time, t_end
        logical :: keep = .True.
        logical :: has_check = .False.

        sizey = size(y)
        has_check = present(check_fun)

        ynew = y
        time = t
        t_end = time + dt

        do while (time < t_end)

            if (has_check) then ! If Check Continue function present
                keep = check_fun(ynew)
                if (.not. keep) then ! If Hard Exit is True
                    dt_adap = time - t ! Replace dt_adap with actual dt used
                    return ! Exit subroutine
                end if
            end if

            ycaller = ynew

            dt_adap = min(dt_adap, t_end - time)
            DT_MIN_NOW = min(DT_MIN, dt_adap)

            der = dydt(time, ycaller)

            call runge_kutta_ptr(sizey, ycaller, dydt, time, dt_adap, der, ynew)

            time = time + dt_adap
        end do

    end subroutine runge_kutta_fixed_caller

end module runge_kutta
