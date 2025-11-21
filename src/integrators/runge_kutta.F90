!> Module with fixed RK integrators
module runge_kutta
    use shared
    implicit none
    private
    public :: init_runge_kutta, free_runge_kutta, runge_kutta_caller, runge_kutta_fixed_caller

    ! Workspace arrays
    real(kind=8), allocatable :: rk(:,:)  ! integrator
    real(kind=8), allocatable :: rk_imp1(:)  ! implicit 1k integrators
    real(kind=8), allocatable :: rk_imp_solv(:,:)  ! implicit solver
    real(kind=8), allocatable :: yscal(:)  ! solver
    real(kind=8), allocatable :: yaux(:)  ! solver
    real(kind=8), allocatable :: yhalf(:), derhalf(:)  ! solver
    real(kind=8), allocatable :: ycaller(:)  ! caller

    ! Constants
    integer(kind=4), parameter :: MAX_N_ITER = 350
    real(kind=8), parameter :: MAX_DT_FACTOR = 5.0d0

    ! Variables to set
    integer(kind=4) :: OSOL
    integer(kind=4) :: N_STAGES
    integer(kind=4) :: N_IMPLICIT

    ! Derived
    real(kind=8) :: ONE_OSOL
    real(kind=8) :: ONE__OSOL_PLUS_ONE
    real(kind=8) :: ONE__MINUS_ONE_PLUS_TWO_TO_ORD

    ! Pointer to runge_kutta used
    procedure (runge_kutta_tem), pointer :: runge_kutta_ptr => null ()

    abstract interface

        subroutine runge_kutta_tem (sizey, y, dydt, t, dt, deri, ynew)
            import :: dydt_tem
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
        end subroutine runge_kutta_tem

        subroutine implicit_funK_tem (sizey, y, dydt, t, dt, sizerk, kin, kout)
            import :: dydt_tem
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            integer(kind=4), intent(in) :: sizerk
            real(kind=8), dimension(sizey, sizerk), intent(in) :: kin
            real(kind=8), dimension(sizey, sizerk), intent(out) :: kout
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
                runge_kutta_ptr => Runge_Kutta5
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
                print*, "Method not available:", which
                stop 1
            end if

            ! Constants
            ONE_OSOL = ONE / dble(OSOL)
            ONE__OSOL_PLUS_ONE = ONE / (OSOL + ONE)
            ONE__MINUS_ONE_PLUS_TWO_TO_ORD = ONE / (TWO**OSOL - ONE)

            ! integrator
            allocate(rk(sizey, N_STAGES))

            ! implicit integrators
            if (N_IMPLICIT == 1) allocate(rk_imp1(sizey))
            if (N_IMPLICIT > 0) allocate(rk_imp_solv(sizey, N_IMPLICIT))

            ! solver
            allocate(yscal(sizey))
            allocate(yaux(sizey))
            allocate(yhalf(sizey), derhalf(sizey))

            ! caller
            allocate(ycaller(sizey))
        end subroutine init_runge_kutta

        subroutine free_runge_kutta()
            implicit none
            if (allocated(rk)) deallocate(rk)
            if (allocated(rk_imp1)) deallocate(rk_imp1)
            if (allocated(rk_imp_solv)) deallocate(rk_imp_solv)
            if (allocated(yscal)) deallocate(yscal)
            if (allocated(yaux)) deallocate(yaux)
            if (allocated(yhalf)) deallocate(yhalf)
            if (allocated(derhalf)) deallocate(derhalf)
            if (allocated(ycaller)) deallocate(ycaller)
            nullify(runge_kutta_ptr)
        end subroutine free_runge_kutta

        
        !!!! Runge Kutta subroutines

        subroutine Euler1 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            ynew = y + dt * deri

        end subroutine Euler1

        subroutine Euler_back1 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            rk_imp1(1:sizey) = ZERO
            call solve_1k_implicit (sizey, y, dydt, t + dt, dt, rk_imp1(1:sizey), ONE, rk(1:sizey,1))

            ynew = y + dt * rk(1:sizey,1)

        end subroutine Euler_back1

        subroutine Euler_center2 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk_imp1(1:sizey) = ZERO
            call solve_1k_implicit (sizey, y, dydt, t + dt * C1_2, dt, rk_imp1(1:sizey), C1_2, rk(1:sizey,1))

            ynew = y + dt * rk(1:sizey,1)

        end subroutine Euler_center2

        subroutine Crank_Nicolson2 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            rk_imp1(1:sizey) = deri * C1_2
            call solve_1k_implicit (sizey, y, dydt, t + dt, dt, rk_imp1(1:sizey), C1_2, rk(1:sizey,2))

            ynew = y + dt * (deri + rk(1:sizey,2)) * C1_2

        end subroutine Crank_Nicolson2

        subroutine Heun2 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt, y + dt * deri)
            
            ynew = y + dt * (deri + rk(1:sizey,2)) * C1_2

        end subroutine Heun2

        subroutine midpoint2 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_2, y + dt * deri * C1_2)
            
            ynew = y + dt * rk(1:sizey,2)

        end subroutine midpoint2

        subroutine strange2 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C3_4, y + dt * deri * C3_4)
            
            ynew = y + dt * (deri + rk(1:sizey,2) * TWO) * C1_3

        end subroutine strange2

        subroutine Ralston2 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C2_3, y + dt * deri * C2_3)
            
            ynew = y + dt * (deri + rk(1:sizey,2) * THREE) * C1_4

        end subroutine Ralston2

        subroutine Hammer_Hollingsworth2 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk_imp1(1:sizey) = deri * C1_3
            call solve_1k_implicit (sizey, y, dydt, t + dt * C2_3, dt, rk_imp1(1:sizey), C1_3, rk(1:sizey,2))

            ynew = y + dt * (deri + rk(1:sizey,2) * THREE) * C1_4

        end subroutine Hammer_Hollingsworth2
        
        subroutine Kraaijevanger_Spijker2 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            rk_imp1(1:sizey) = ZERO
            call solve_1k_implicit (sizey, y, dydt, t + dt * C1_2, dt, rk_imp1(1:sizey), C1_2, rk(1:sizey,1))

            rk_imp1(1:sizey) = - rk(1:sizey,1) * C1_2
            call solve_1k_implicit (sizey, y, dydt, t + dt * C3_2, dt, rk_imp1(1:sizey), TWO, rk(1:sizey,2))

            ynew = y + dt * (- rk(1:sizey,1) + rk(1:sizey,2) * THREE) * C1_2

        end subroutine Kraaijevanger_Spijker2
        
        subroutine Qin_Zhang2 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            rk_imp1(1:sizey) = ZERO
            call solve_1k_implicit (sizey, y, dydt, t + dt * C1_4, dt, rk_imp1(1:sizey), C1_4, rk(1:sizey,1))

            rk_imp1(1:sizey) = rk(1:sizey,1) * C1_2
            call solve_1k_implicit (sizey, y, dydt, t + dt * C3_4, dt, rk_imp1(1:sizey), C1_4, rk(1:sizey,2))

            ynew = y + dt * (rk(1:sizey,1) + rk(1:sizey,2)) * C1_2

        end subroutine Qin_Zhang2

        subroutine Runge_Kutta3 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_2, y + dt * deri * C1_2)
            rk(1:sizey,3) = dydt (t + dt, y + dt * (- deri + rk(1:sizey,2) * TWO))
            
            ynew = y + dt * (deri + rk(1:sizey,2) * FOUR + rk(1:sizey,3)) * C1_6

        end subroutine Runge_Kutta3

        subroutine Heun3 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_3, y + dt * deri * C1_3)
            rk(1:sizey,3) = dydt (t + dt * C2_3, y + dt * rk(1:sizey,2) * C2_3)
            
            ynew = y + dt * (deri + rk(1:sizey,3) * THREE) * C1_4

        end subroutine Heun3

        subroutine Ralston3 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_2, y + dt * deri * C1_2)
            rk(1:sizey,3) = dydt (t + dt * C3_4, y + dt * rk(1:sizey,2) * C3_4)

            ynew = y + dt * (deri * TWO + rk(1:sizey,2) * THREE + rk(1:sizey,3) * FOUR) * C1_9

        end subroutine Ralston3

        subroutine SSPRrk3 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt, y + dt * der)
            rk(1:sizey,3) = dydt (t + dt * C1_2, y + dt * (deri + rk(1:sizey,2)) * C1_4)

            ynew = y + dt * (deri + rk(1:sizey,2) + rk(1:sizey,3) * FOUR) * C1_6

        end subroutine SSPRrk3

        subroutine Crouzeix3 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            real(kind=8), parameter :: aux = C1_2 + SQ3_6

            rk_imp1(1:sizey) = ZERO
            call solve_1k_implicit (sizey, y, dydt, t + dt * aux, dt, rk_imp1(1:sizey), aux, rk(1:sizey,1))

            rk_imp1(1:sizey) = - rk(1:sizey,1) * SQ3_6 * 2
            call solve_1k_implicit (sizey, y, dydt, t + dt * (C1_2 - SQ3_6), dt, rk_imp1(1:sizey), aux, rk(1:sizey,2))

            ynew = y + dt * (rk(1:sizey,1) + rk(1:sizey,2)) * C1_2
        end subroutine Crouzeix3

        subroutine Runge_Kutta_implicit3 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            rk_imp1(1:sizey) = ZERO
            call solve_1k_implicit (sizey, y, dydt, t + dt * C1_2, dt, rk_imp1(1:sizey), C1_2, rk(1:sizey,1))

            rk_imp1(1:sizey) = rk(1:sizey,1) * C1_6
            call solve_1k_implicit (sizey, y, dydt, t + dt * C2_3, dt, rk_imp1(1:sizey), C1_2, rk(1:sizey,2))

            rk_imp1(1:sizey) = (- rk(1:sizey,1) + rk(1:sizey,2)) * C1_2
            call solve_1k_implicit (sizey, y, dydt, t + dt * C1_2, dt, rk_imp1(1:sizey), C1_2, rk(1:sizey,3))

            rk_imp1(1:sizey) = ((rk(1:sizey,1) - rk(1:sizey,2)) * 3 + rk(1:sizey,3)) * C1_2
            call solve_1k_implicit (sizey, y, dydt, t + dt, dt, rk_imp1(1:sizey), C1_2, rk(1:sizey,4))

            ynew = y + dt * (rk_imp1(1:sizey) + rk(1:sizey,4) * C1_2)
        end subroutine Runge_Kutta_implicit3

        subroutine Ralston4 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C2_5, y + dt * deri * C2_5)
            rk(1:sizey,3) = dydt (t + dt * 0.45573725d0, y + dt * (deri * 0.29697761d0 + rk(1:sizey,2) * 0.15875964d0))
            rk(1:sizey,4) = dydt (t + dt, y + dt * (deri * 0.21810040d0 - rk(1:sizey,2) * 3.05096516d0 + &
                                                & rk(1:sizey,3) * 3.83286476d0))

            ynew = y + dt * (deri * 0.17476028d0 - rk(1:sizey,2) * 0.55148066d0 + rk(1:sizey,3) * 1.20553560d0 + &
                                                & rk(1:sizey,4) * 0.17118478d0)

        end subroutine Ralston4

        subroutine Lobatto4 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
           implicit none
           integer(kind=4), intent(in) :: sizey
           real(kind=8), dimension(sizey),  intent(in) :: y
           procedure(dydt_tem) :: dydt
           real(kind=8), intent(in) :: t, dt
           real(kind=8), dimension(sizey), intent(in) :: deri
           real(kind=8), dimension(sizey), intent(out) :: ynew
        
            rk_imp1(1:sizey) = deri * C1_4
            call solve_1k_implicit (sizey, y, dydt, t + dt * C1_2, dt, rk_imp1(1:sizey), C1_4, rk(1:sizey,2))

            rk(1:sizey,3) = dydt (t + dt, y + dt * rk(1:sizey,2))

            ynew = y + dt * (deri + rk(1:sizey,2) * FOUR + rk(1:sizey,3)) * C1_6

        end subroutine Lobatto4

        subroutine Runge_Kutta4 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_2, y + dt * deri * C1_2)
            rk(1:sizey,3) = dydt (t + dt * C1_2, y + dt * rk(1:sizey,2) * C1_2)
            rk(1:sizey,4) = dydt (t + dt,  y + dt * rk(1:sizey,3))
        
            ynew = y + dt * (deri + (rk(1:sizey,2) + rk(1:sizey,3)) * TWO + rk(1:sizey,4)) * C1_6

        end subroutine Runge_Kutta4

        subroutine Gauss_Legendre4 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            rk(1:sizey,:) = ZERO
            call solve_rk_implicit (sizey, N_STAGES, y, dydt, t, dt, FunK_GL4, rk(1:sizey,:))

            ynew = y + dt * (rk(1:sizey,1) + rk(1:sizey,2)) * C1_2

            contains
                subroutine FunK_GL4 (sizey, y, dydt, t, dt, sizerk, kin, kout) !! Funk for Gauss_Legendre4
                    implicit none
                    integer(kind=4), intent(in) :: sizey
                    real(kind=8), dimension(sizey), intent(in) :: y
                    procedure(dydt_tem) :: dydt
                    real(kind=8), intent(in) :: t, dt
                    integer(kind=4), intent(in) :: sizerk
                    real(kind=8), dimension(sizey, sizerk), intent(in) :: kin
                    real(kind=8), dimension(sizey, sizerk), intent(out) :: kout
        
                    kout(:,1) = dydt (t + dt * (C1_2 - SQ3_6), y + dt * ( &
                        & kin(:,1) * C1_4 + &
                        & kin(:,2) * (C1_4 - SQ3_6)))
                    kout(:,2) = dydt (t + dt * (C1_2 + SQ3_6), y + dt * ( &
                        & kin(:,1) * (C1_4 + SQ3_6) + &
                        & kin(:,2) * C1_4))
                end subroutine FunK_GL4

        end subroutine Gauss_Legendre4 

        subroutine Runge_Kutta_four_oct4 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,1) = deri * C1_3
            rk(1:sizey,2) = dydt (t + dt * C1_3, y + dt * rk(1:sizey,1))
            rk(1:sizey,3) = dydt (t + dt * C2_3, y + dt * (- rk(1:sizey,1) + rk(1:sizey,2)))
            rk(1:sizey,4) = dydt (t + dt, y + dt * (deri - rk(1:sizey,2) + rk(1:sizey,3)))
        
            ynew = y + dt * (deri + (rk(1:sizey,2) + rk(1:sizey,3)) * THREE + rk(1:sizey,4)) * C1_8  

        end subroutine Runge_Kutta_four_oct4

        subroutine Runge_Kutta5 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_4, y + dt * deri * C1_4)
            rk(1:sizey,3) = dydt (t + dt * C1_4, y + dt * (deri + rk(1:sizey,2)) * C1_8)
            rk(1:sizey,4) = dydt (t + dt * C1_2, y + dt * (- rk(1:sizey,2) * C1_2 + rk(1:sizey,3)))
            rk(1:sizey,5) = dydt (t + dt * C3_4, y + dt * (deri + rk(1:sizey,4) * THREE) * C3_16)
            rk(1:sizey,6) = dydt (t + dt, y + dt * (- deri * THREE + rk(1:sizey,2) * TWO + (rk(1:sizey,3) - &
                                                    & rk(1:sizey,4)) * 12.d0 + rk(1:sizey,5) * 8.d0) * C1_7)
        
            ynew = y + dt * ((deri + rk(1:sizey,6)) * 7.d0 + (rk(1:sizey,3) + rk(1:sizey,5)) * 32.d0 + rk(1:sizey,4) * 12.d0)/90.d0

        end subroutine Runge_Kutta5

        subroutine Gauss_Legendre6 (sizey, y, dydt, t, dt, deri, ynew) ! Implicit
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            rk(1:sizey,:) = ZERO
            call solve_rk_implicit (sizey, N_STAGES, y, dydt, t, dt, FunK_GL6, rk(1:sizey,:))

            ynew = y + dt * ((rk(1:sizey,1) + rk(1:sizey,3)) * FIVE + rk(1:sizey,2) * 8.d0) * C1_18

            contains
                subroutine FunK_GL6 (sizey, y, dydt, t, dt, sizerk, kin, kout) !! Funk for Gauss_Legendre6
                    implicit none
                    integer(kind=4), intent(in) :: sizey
                    real(kind=8), dimension(sizey), intent(in) :: y
                    procedure(dydt_tem) :: dydt
                    real(kind=8), intent(in) :: t, dt
                    integer(kind=4), intent(in) :: sizerk
                    real(kind=8), dimension(sizey, sizerk), intent(in) :: kin
                    real(kind=8), dimension(sizey, sizerk), intent(out) :: kout
        
                    kout(:,1) = dydt (t + dt * (ONE - SQ15_5) * C1_2,  y + dt * (&
                            & kin(:,1) * C5_36 + &
                            & kin(:,2) * (C2_3 - SQ15_5) * C1_3 + &
                            & kin(:,3) * (C5_36 - SQ15_30)))
                    kout(:,2) = dydt (t + dt * C1_2, y + dt * (&
                            & kin(:,1) * (C5_36 + SQ15_24)+ &
                            & kin(:,2) * C2_9 + &
                            & kin(:,3) * (C5_36 - SQ15_24)))
                    kout(:,3) = dydt (t + dt * (ONE + SQ15_5) * C1_2, y + dt * (&
                            & kin(:,1) * (C5_36 + SQ15_30) + &
                            & kin(:,2) * (C2_9 + SQ15_15) + &
                            & kin(:,3) * C5_36))
                end subroutine FunK_GL6

        end subroutine Gauss_Legendre6 

        subroutine Runge_Kutta6 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt(t + dt * C1_3, y + dt * C1_3 * deri)
            rk(1:sizey,3) = dydt(t + dt * C1_3, y + dt * (deri + rk(1:sizey,2)) * C1_6)
            rk(1:sizey,4) = dydt(t + dt * C1_2, y + dt * (deri * C1_8 + rk(1:sizey,3) * C3_8))
            rk(1:sizey,5) = dydt(t + dt, y + dt * (deri * C1_2 - rk(1:sizey,3) * C3_2 + rk(1:sizey,4) * TWO))
            rk(1:sizey,6) = dydt(t + dt * C2_3, y + dt * (- deri * C3_7 + rk(1:sizey,3) * C8_7 + rk(1:sizey,4) * C6_7 - &
                                                        & rk(1:sizey,5) * C12_7))
            rk(1:sizey,7) = dydt(t + dt, y + dt * (deri * C7_90 + rk(1:sizey,3) * C16_45 + rk(1:sizey,4) * C2_15 + &
                                                        & rk(1:sizey,5) * (C16_45) + rk(1:sizey,6) * (C7_90)))

            ynew = y + dt * (deri * C7_90 + rk(1:sizey,3) * C16_45 + rk(1:sizey,4) * C2_15 + &
                                                        & rk(1:sizey,5) * C16_45 + rk(1:sizey,6) * C7_90)
                            
        end subroutine Runge_Kutta6

        subroutine Abbas6 (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_3, y + dt * deri * C1_3)
            rk(1:sizey,3) = dydt (t + dt * C2_3, y + dt * rk(1:sizey,2) * C2_3)
            rk(1:sizey,4) = dydt (t + dt * C1_3, y + dt * (deri + rk(1:sizey,2) * FOUR - rk(1:sizey,3)) * C1_12)
            rk(1:sizey,5) = dydt (t + dt * C5_6, y + dt * (deri * 25.d0 - rk(1:sizey,2) * 110.d0 + rk(1:sizey,3) * 35.d0 + &
                                                            & rk(1:sizey,4) * 90.d0)/48.d0)
            rk(1:sizey,6) = dydt (t + dt * C1_6, y + dt * (deri * 0.15d0 - rk(1:sizey,2) * 0.55d0 - rk(1:sizey,3) * C1_8 + &
                                                            & rk(1:sizey,4) * C1_2 + rk(1:sizey,5) * 0.1d0))
            rk(1:sizey,7) = dydt (t + dt, y + dt * (- deri * 195.75d0 + rk(1:sizey,2) * 495.d0 + rk(1:sizey,3) * 53.75d0 - &
                                                            & rk(1:sizey,4) * 590.d0 + rk(1:sizey,5) * 32.d0 + &
            & rk(1:sizey,6) * 400.d0)/195.d0)
        
            ynew = y + dt * ((deri + rk(1:sizey,7)) * 13.d0 + (rk(1:sizey,3) + rk(1:sizey,4)) * 55.d0 + (rk(1:sizey,5) + &
                                                            & rk(1:sizey,6)) * 32.d0) * 5.d-3

        end subroutine Abbas6


        !------------------------------------------------
        !  Solvers Runge Kutta (Implicit and Half Step)
        !------------------------------------------------

        !!!! Implicit: 1D SOLVER
        subroutine solve_1k_implicit (sizey, y, dydt, t, dt, kprev, cte, kout)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: kprev
            real(kind=8), intent(in) :: cte
            real(kind=8), dimension(sizey), intent(out) :: kout
            integer(kind=4) :: i
            
            kout = dydt (t, y + dt * kprev)
            do i = 1, MAX_N_ITER
                rk_imp_solv(1:sizey,1) = kout
                kout = dydt (t, y + dt * (kprev + cte * kout))
                if (maxval(abs((rk_imp_solv(1:sizey,1) - kout) / (rk_imp_solv(1:sizey,1) + SAFE_LOW))) .le. E_TOL) then                    
                    exit
                end if
            end do
        end subroutine solve_1k_implicit

        !!!! Implicit: ND SOLVER
        subroutine solve_rk_implicit (sizey, nstages, y, dydt, t, dt, impl_funK, rkout)
            implicit none
            integer(kind=4), intent(in) :: sizey, nstages
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            procedure(implicit_funK_tem) :: impl_funK
            real(kind=8), dimension(sizey, nstages), intent(inout) :: rkout
            integer(kind=4) :: i

            do i = 1, MAX_N_ITER
                rk_imp_solv(1:sizey, 1:nstages) = rkout
                call impl_funK (sizey, y, dydt, t, dt, nstages, rk_imp_solv(1:sizey, 1:nstages), rkout)
                if (maxval(abs((rk_imp_solv(1:sizey, 1:nstages) - rkout) / &
                             & (rk_imp_solv(1:sizey, 1:nstages) + SAFE_LOW))) .le. E_TOL) then
                    exit
                end if
            end do
        end subroutine solve_rk_implicit
        
        !!!! Adaptive Timestep: Half step
        recursive subroutine solve_rk_half_step (sizey, y, dydt, t, dt_adap, dt_used, deri, runge_kutta, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t
            real(kind=8), intent(inout) :: dt_adap
            real(kind=8), intent(out) :: dt_used
            real(kind=8), dimension(sizey), intent(in) :: deri
            procedure (runge_kutta_tem), pointer :: runge_kutta
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            real(kind=8)                                   :: e_calc, ratio, hdt_adap
            integer(kind=4)                                :: iter = 0
            
            iter = iter + 1
            dt_adap = max (dt_adap, DT_MIN_NOW)
            hdt_adap = C1_2 * dt_adap

            ! yscal
            yscal(:sizey) = abs (y) + abs (dt_adap * deri) + SAFE_LOW

            ! y(t, dt; der) -> ynew
            call runge_kutta (sizey, y, dydt, t, dt_adap, deri, ynew)

            ! y(t, dt/2; der) -> yhalf
            call runge_kutta (sizey, y, dydt, t, hdt_adap, deri, yhalf(:sizey))

            ! d[yhalf (t + dt/2)]/dt -> derhalf
            derhalf(:sizey) = dydt (t + hdt_adap, yhalf(:sizey))

            ! yhalf(t + dt/2, dt/2; derhalf) -> yaux
            call runge_kutta (sizey, yhalf(:sizey), dydt, t + hdt_adap, hdt_adap, derhalf(:sizey), yaux(:sizey))
            
            ! Error
            e_calc = max(maxval(abs((ynew - yaux(:sizey)) / yscal(:sizey)) * ONE__MINUS_ONE_PLUS_TWO_TO_ORD), SAFE_LOW)
            ratio = E_TOL / e_calc

            if (ratio > ONE) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (BETA * ratio**ONE__OSOL_PLUS_ONE, MAX_DT_FACTOR)
                iter = 0

            else
                if (abs(dt_adap - DT_MIN_NOW) .le. E_TOL) then !E_TOL?
                    dt_used = DT_MIN_NOW
                    iter = 0

                else
                    dt_adap = dt_adap * min (BETA * ratio**ONE_OSOL, MAX_DT_FACTOR)

                    if ((dt_adap /= dt_adap) .or. (dt_adap .le. DT_MIN_NOW) .or. (iter == MAX_N_ITER)) then
                        dt_used = DT_MIN_NOW
                        dt_adap = DT_MIN_NOW

                        call runge_kutta (sizey, y, dydt, t, dt_adap, deri, ynew)
                        iter = 0

                    else
                        call solve_rk_half_step (sizey, y, dydt, t, dt_adap, dt_used, deri, runge_kutta, ynew)

                    end if

                end if

            end if

        end subroutine solve_rk_half_step


        ! ------------------------------------------------
        !  Call Runge Kutta integrator (adaptive timestep)
        ! ------------------------------------------------

        subroutine runge_kutta_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), intent(inout) :: dt_adap  ! This try and also next try
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: dt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            procedure(function_check_keep_tem), optional :: check_fun
            
            integer(kind=4) :: sizey
            real(kind=8) :: time, t_end, dt_used
            logical :: keep = .True.
            logical, save :: has_check = .False.

            sizey = size (y)
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

                ycaller(:sizey) = ynew

                dt_adap = min (dt_adap, t_end - time)
                DT_MIN_NOW = min (DT_MIN, dt_adap)

                der(:sizey) = dydt (time, ycaller(:sizey))
                                    
                call solve_rk_half_step (sizey, ycaller(:sizey), dydt, time, dt_adap, dt_used, &
                                    & der(:sizey), runge_kutta_ptr, ynew)

                time = time + dt_used
            end do

        end subroutine runge_kutta_caller


        ! ------------------------------------------------
        !  Call Runge Kutta integrator (fixed timestep)
        ! ------------------------------------------------

        subroutine runge_kutta_fixed_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), intent(inout) :: dt_adap  ! This is each sub-step
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: dt ! This is full step
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            procedure(function_check_keep_tem), optional :: check_fun
            
            integer(kind=4) :: sizey
            real(kind=8) :: time, t_end
            logical :: keep = .True.
            logical, save :: has_check = .False.

            sizey = size (y)
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

                ycaller(:sizey) = ynew

                dt_adap = min (dt_adap, t_end - time)
                DT_MIN_NOW = min (DT_MIN, dt_adap)

                der(:sizey) = dydt (time, ycaller(:sizey))

                call runge_kutta_ptr (sizey, ycaller(:sizey), dydt, time, dt_adap, der(:sizey), ynew)

                time = time + dt_adap
            end do

        end subroutine runge_kutta_fixed_caller
    
end module runge_kutta