!> Module with embedded RK integrators
module embedded
    use shared
    implicit none
    private
    public :: init_embedded, free_embedded, embedded_caller, embedded_fixed_caller

    ! Workspace arrays
    real(kind=8), allocatable :: rk(:,:)  ! integrator
    real(kind=8), allocatable :: yscal(:)  ! solver
    real(kind=8), allocatable :: yaux(:)  ! solver
    real(kind=8), allocatable :: ycaller(:)  ! caller

    ! Constants
    integer(kind=4), parameter :: MAX_N_ITER = 350
    real(kind=8), parameter :: MAX_DT_FACTOR = 5.0d0

    ! Variables to set
    integer(kind=4) :: OSOL
    real(kind=8) :: ONE_OSOL
    integer(kind=4) :: OAUX
    real(kind=8) :: ONE_OAUX

    ! Pointer to embedded used
    procedure (embedded_tem), pointer :: embedded_ptr => null ()

    abstract interface

        subroutine embedded_tem (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            import :: dydt_tem
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
        end subroutine embedded_tem
    
    end interface
    
    contains
    
        !!!! HANDLER

        subroutine init_embedded(sizey, which)
            implicit none
            integer(kind=4), intent(in) :: sizey
            integer(kind=4), intent(in) :: which
            integer(kind=4) :: n_stages

            if (which == 1) then
                embedded_ptr => Fehlberg1_2
                OSOL = 2
                OAUX = 1
                n_stages = 2
            
            else if (which == 2) then
                embedded_ptr => Heun_Euler2_1
                OSOL = 1
                OAUX = 2
                n_stages = 2
            
            else if (which == 3) then
                embedded_ptr => Fehlberg2_1
                OSOL = 2
                OAUX = 1
                n_stages = 2
            
            else if (which == 4) then
                embedded_ptr => Bogacki_Shampine3_2
                OSOL = 3
                OAUX = 2
                n_stages = 4
            
            else if (which == 5) then
                embedded_ptr => Zonneveld4_3
                OSOL = 4
                OAUX = 3
                n_stages = 5
            
            else if (which == 6) then
                embedded_ptr => Merson4_3
                OSOL = 4
                OAUX = 3
                n_stages = 5
            
            else if (which == 7) then
                embedded_ptr => Fehlberg4_5
                OSOL = 4
                OAUX = 5
                n_stages = 6
            
            else if (which == 8) then
                embedded_ptr => Cash_Karp5_4
                OSOL = 5
                OAUX = 4
                n_stages = 6
            
            else if (which == 9) then
                embedded_ptr => Dormand_Prince5_4
                OSOL = 5
                OAUX = 4
                n_stages = 7
            
            else if (which == 10) then
                embedded_ptr => Verner6_5
                OSOL = 6
                OAUX = 5
                n_stages = 8
            
            else if (which == 11) then
                embedded_ptr => Fehlberg7_8
                OSOL = 7
                OAUX = 8
                n_stages = 13
            
            else if (which == 12) then
                embedded_ptr => Dormand_Prince8_7
                OSOL = 8
                OAUX = 7
                n_stages = 13
            
            else
                print*, "Method not available:", which
                stop 1
            end if

            ! Constants
            ONE_OSOL = ONE / dble(OSOL)
            ONE_OAUX = ONE / dble(OAUX)

            ! integrator
            allocate(rk(sizey, n_stages))

            ! solver
            allocate(yscal(sizey))
            allocate(yaux(sizey))

            ! caller
            allocate(ycaller(sizey))
        end subroutine init_embedded

        subroutine free_embedded()
            implicit none
            if (allocated(rk)) deallocate(rk)
            if (allocated(yscal)) deallocate(yscal)
            if (allocated(yaux)) deallocate(yaux)
            if (allocated(ycaller)) deallocate(ycaller)
            nullify(embedded_ptr)
        end subroutine free_embedded

        
        !!!! Embedded subroutines

        subroutine Fehlberg1_2 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_2, y + dt * deri * C1_2)

            ynew = y + dt * rk(1:sizey,2)
            
            yauxi = y + dt * C1_2 * (deri + rk(1:sizey,2))

        end subroutine Fehlberg1_2

        subroutine Heun_Euler2_1 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt, y + dt * deri)
            
            ynew = y + dt * (deri + rk(1:sizey,2)) * C1_2
            yauxi = y + dt * deri

        end subroutine Heun_Euler2_1

        subroutine Fehlberg2_1 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C3_4, y + dt * deri * C3_4)

            ynew = y + dt * (deri * C1_3 + rk(1:sizey,2) * C2_3)
            
            yauxi = y + dt * rk(1:sizey,2)

        end subroutine Fehlberg2_1

        subroutine Bogacki_Shampine3_2 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew

            rk(1:sizey,2) = dydt (t + dt * C1_2, y + dt * deri * C1_2)
            rk(1:sizey,3) = dydt (t + dt * C3_4, y + dt * rk(1:sizey,2) * C3_4)

            ynew = y + dt * (deri * TWO + rk(1:sizey,2) * THREE + rk(1:sizey,3) * FOUR) * C1_9

            rk(1:sizey,4) = dydt (t + dt, ynew)           
            
            yauxi = y + dt * (deri * 7.d0/24.d0 + rk(1:sizey,2) * C1_4 + rk(1:sizey,3) * C1_3 + rk(1:sizey,4) * C1_8)

        end subroutine Bogacki_Shampine3_2

        subroutine Zonneveld4_3 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            real(kind=8) :: dt_2

            dt_2 = dt * C1_2
            
            rk(1:sizey,2) = dydt (t + dt_2, y + dt * deri * C1_2)
            rk(1:sizey,3) = dydt (t + dt_2, y + dt * rk(1:sizey,2) * C1_2)
            rk(1:sizey,4) = dydt (t + dt, y + dt * rk(1:sizey,3))
            rk(1:sizey,5) = dydt (t + dt * C3_4, y + dt * ( &
                & deri * FIVE + rk(1:sizey,2) * 7.d0 + rk(1:sizey,3) * 13.d0 - rk(1:sizey,4)) * 0.03125d0)
            
            rk(1:sizey,1) = rk(1:sizey,2) + rk(1:sizey,3)  ! Used as kaux

            ynew = y + dt * (deri + rk(1:sizey,1) * TWO + rk(1:sizey,4)) * C1_6
            yauxi = y + dt * (- deri * THREE + rk(1:sizey,1) * 14.d0 + rk(1:sizey,4) * 13.d0 - rk(1:sizey,5) * 32.d0) * C1_6

        end subroutine Zonneveld4_3

        subroutine Merson4_3 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_3, y + dt * deri * C1_3)
            rk(1:sizey,3) = dydt (t + dt * C1_3, y + dt * (deri + rk(1:sizey,2)) * C1_6)
            rk(1:sizey,4) = dydt (t + dt * C1_2, y + dt * (deri + rk(1:sizey,3) * THREE) * C1_8)
            rk(1:sizey,5) = dydt (t + dt, y + dt * (deri - rk(1:sizey,3) * THREE + rk(1:sizey,4) * FOUR) * C1_2)
            
            ynew = y + dt * (deri + rk(1:sizey,4) * FOUR + rk(1:sizey,5)) * C1_6
            yauxi = y + dt * (deri + rk(1:sizey,3) * THREE + rk(1:sizey,4) * FOUR + rk(1:sizey,5) * TWO) * 0.1d0

        end subroutine Merson4_3

        subroutine Fehlberg4_5 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_4, y + dt * deri * C1_4)
            rk(1:sizey,3) = dydt (t + dt * C3_8, y + dt * (deri * THREE + rk(1:sizey,2) * 9.d0) * 0.03125d0)
            rk(1:sizey,4) = dydt (t + dt * 12.d0/13.d0, y + dt * (deri * 1932.d0 - rk(1:sizey,2) * 7200.d0 + &
                                                & rk(1:sizey,3) * 7296.d0)/2197.d0)
            rk(1:sizey,5) = dydt (t + dt, y + dt * ((deri * 8341.d0 + rk(1:sizey,3) * 29440.d0 - &
                                                & rk(1:sizey,4) * 845.d0)/4104.d0 - rk(1:sizey,2) * 8.d0))
            rk(1:sizey,6) = dydt (t + dt * C1_2, y + dt * ((- deri * 1216.d0 + rk(1:sizey,4) * 1859.d0)/4104.d0 + &
                                                & rk(1:sizey,2) * TWO - rk(1:sizey,3) * 3544.d0/2565.d0 - rk(1:sizey,5) * 0.275d0))
            
            ynew = y + dt * ((deri * 475.d0 + rk(1:sizey,4) * 2197.d0)/4104.d0 + rk(1:sizey,3) * 1408.d0/2565.d0 - &
                                                & rk(1:sizey,5) * C1_5)
            yauxi = y + dt * ((deri * 6688.d0 + rk(1:sizey,4) * 28561.d0 + rk(1:sizey,6) * 2052.d0)/56430.d0 + & 
                                                & rk(1:sizey,3) * 6656.d0/12825.d0 - rk(1:sizey,5) * 0.18d0)

        end subroutine Fehlberg4_5

        subroutine Cash_Karp5_4 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_5, y + dt * deri * C1_5)
            rk(1:sizey,3) = dydt (t + dt * 0.3d0, y + dt * (deri * THREE + rk(1:sizey,2) * 9.d0) * 0.025d0)
            rk(1:sizey,4) = dydt (t + dt * 0.6d0, y + dt * (deri * THREE - rk(1:sizey,2) * 9.d0 + rk(1:sizey,3) * 12.d0) * 0.1d0)
            rk(1:sizey,5) = dydt (t + dt, y + dt * (- deri * 11.d0 + rk(1:sizey,2) * 135.d0 - rk(1:sizey,3) * 140.d0 + &
                                                & rk(1:sizey,4) * 70.d0)/54.d0)
            rk(1:sizey,6) = dydt (t + dt * 0.875d0, y + dt * (deri * 3262.d0 + rk(1:sizey,2) * 37800.d0 + &
                                                & rk(1:sizey,3) * 4600.d0 + rk(1:sizey,4) * 44275.d0 + &
                                                & rk(1:sizey,5) * 6831.d0)/110592.d0)
            
            ynew = y + dt * (deri * 37.d0/378.d0 + rk(1:sizey,3) * 250.d0/621.d0 + rk(1:sizey,4) * 125.d0/594.d0 + &
                                                & rk(1:sizey,6) * 512.d0/1771.d0)
            yauxi = y + dt * ((deri * 5650.d0 + rk(1:sizey,4) * 13525.d0)/55296.d0 + rk(1:sizey,3) * 18575.d0/48384.d0 + &
                                                & rk(1:sizey,5) * 277.d0/14336.d0 + rk(1:sizey,6) * C1_4)

        end subroutine Cash_Karp5_4

        subroutine Dormand_Prince5_4 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_5, y + dt * deri * C1_5)
            rk(1:sizey,3) = dydt (t + dt * 0.3d0, y + dt * (deri * 0.075d0 + rk(1:sizey,2) * 0.225d0))
            rk(1:sizey,4) = dydt (t + dt * 0.8d0, y + dt * (deri * 44.d0 - rk(1:sizey,2) * 168.d0 + rk(1:sizey,3) * 160.d0)/45.d0)
            rk(1:sizey,5) = dydt (t + dt * C8_9, y + dt * (deri * 19372.d0 - rk(1:sizey,2) * 76080.d0 + &
                                                & rk(1:sizey,3) * 64448.d0 - rk(1:sizey,4) * 1908.d0)/6561.d0)
            rk(1:sizey,6) = dydt (t + dt, y + dt * ((deri * 9017.d0 - rk(1:sizey,2) * 34080.d0 + &
                                                & rk(1:sizey,4) * 882.d0)/3168.d0 + rk(1:sizey,3) * 46732.d0/5247.d0 - &
                                                & rk(1:sizey,5) * 5103.d0/18656.d0))
            
            ynew = y + dt * ((deri * 35.d0 + rk(1:sizey,4) * 250.d0)/384.d0 + rk(1:sizey,3) * 500.d0/1113.d0 - &
                                                & rk(1:sizey,5) * 2187.d0/6784.d0 + rk(1:sizey,6) * 11.d0/84.d0)
            
            rk(1:sizey,7) = dydt (t + dt, ynew)

            yauxi = y + dt * ((deri * 5179.d0 + rk(1:sizey,4) * 35370.d0)/57600.d0 + rk(1:sizey,3) * 7571.d0/16695.d0 - &
                                                & rk(1:sizey,5) * 92097.d0/339200.d0 + rk(1:sizey,6) * 187.d0/2100.d0 + &
                                                & rk(1:sizey,7) * 0.025d0)

        end subroutine Dormand_Prince5_4

        subroutine Verner6_5 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2) = dydt (t + dt * C1_6, y + dt * deri * C1_6)
            rk(1:sizey,3) = dydt (t + dt * C4_15, y + dt * (deri * FOUR + rk(1:sizey,2) * 16.d0)/75.d0)
            rk(1:sizey,4) = dydt (t + dt * C2_3, y + dt * (deri * FIVE - rk(1:sizey,2) * 16.d0 + rk(1:sizey,3) * 15.d0) * C1_6)
            rk(1:sizey,5) = dydt (t + dt * C5_6, y + dt * ((- deri * 165.d0 - rk(1:sizey,3) * 425.d0)/64.d0 + &
                                                & (rk(1:sizey,2) * 880.d0 + rk(1:sizey,4) * 85.d0)/96.d0))
            rk(1:sizey,6) = dydt (t + dt, y + dt * ((deri * 612.d0 + rk(1:sizey,5) * 88.d0)/255.d0 - rk(1:sizey,2) * 8.d0 + &
                                                &  (rk(1:sizey,3) * 4015.d0 - rk(1:sizey,4) * 187.d0)/612.d0))
            rk(1:sizey,7) = dydt (t + dt * C1_15, y + dt * ((- deri * 8263.d0 + rk(1:sizey,2) * 24800.d0)/15000.d0 - &
                                                & rk(1:sizey,3) * 643.d0/680.d0 - rk(1:sizey,4) * 0.324d0 + &
                                                & rk(1:sizey,5) * 2484.d0/10625.d0))
            rk(1:sizey,8) = dydt (t + dt, y + dt * (deri * 3501.d0/1720.d0 + (297275.d0 * rk(1:sizey,3) - &
                                                & 367200 * rk(1:sizey,2))/52632.d0 - rk(1:sizey,4) * 319.d0/2322.d0 + &
                                                & rk(1:sizey,5) * 24068.d0/84065.d0 + rk(1:sizey,7) * 3850.d0/26703.d0))
            
            ynew = y + dt * (deri * 0.075d0 + rk(1:sizey,3) * 875.d0/2244.d0 + (rk(1:sizey,4) * 3703.d0 + &
                                                & rk(1:sizey,7) * 125.d0)/11592.d0 + rk(1:sizey,5) * 264.d0/1955.d0 + &
                                                & rk(1:sizey,8) * 43.d0/616.d0)
            yauxi = y + dt * ((deri * 13.d0 + rk(1:sizey,4) * 50.d0)/160.d0 + (rk(1:sizey,3) * 2375.d0 + &
                                                & rk(1:sizey,6) * 408.d0)/5984.d0 + rk(1:sizey,5) * 12.d0/85.d0)

        end subroutine Verner6_5

        subroutine Fehlberg7_8 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2)   = dydt (t + dt * C2_27, y + dt * deri * C2_27)
            rk(1:sizey,3)   = dydt (t + dt * C1_9, y + dt * (deri + rk(1:sizey,2) * THREE)/36.d0)
            rk(1:sizey,4)   = dydt (t + dt * C1_6, y + dt * (deri + rk(1:sizey,3) * THREE)/24.d0)
            rk(1:sizey,5)   = dydt (t + dt * C5_12, y + dt * (deri * C5_12 + (rk(1:sizey,4) - rk(1:sizey,3)) * 1.5625d0))
            rk(1:sizey,6)   = dydt (t + dt * C1_2, y + dt * (deri + rk(1:sizey,4) * FIVE + rk(1:sizey,5) * FOUR) * 0.05d0)
            rk(1:sizey,7)   = dydt (t + dt * C5_6, y + dt * (- deri * 25.d0 + rk(1:sizey,4) * 125.d0 - rk(1:sizey,5) * 260.d0 + &
                                                & rk(1:sizey,6) * 250.d0)/108.d0)
            rk(1:sizey,8)   = dydt (t + dt * C1_6, y + dt * (deri * 93.d0 + rk(1:sizey,5) * 244.d0 - rk(1:sizey,6) * 200.d0 + &
                                                & rk(1:sizey,7) * 13.d0)/900.d0)
            rk(1:sizey,9)   = dydt (t + dt * C2_3, y + dt * (deri * TWO + (- rk(1:sizey,4) * 795.d0 + rk(1:sizey,5) * 1408.d0 - &
                                                & rk(1:sizey,6) * 1070.d0 + rk(1:sizey,7) * 67.d0)/90.d0 + rk(1:sizey,8) * THREE))
            rk(1:sizey,10)  = dydt (t + dt * C1_3,  y + dt * ((- deri * 91.d0 + rk(1:sizey,4) * 23.d0  + &
                                                & rk(1:sizey,6) * 622.d0)/108.d0 - rk(1:sizey,5) * 976.d0/135.d0 + &
                                                & (- rk(1:sizey,7) * 19.d0 + rk(1:sizey,8) * 170.d0 - rk(1:sizey,9) * 5.d0)/60.d0))

            rk(1:sizey,1)  = - rk(1:sizey,4) * 8525.d0 + rk(1:sizey,5) * 17984.d0  !! Used as kaux

            rk(1:sizey,11)  = dydt (t + dt, y + dt * (deri * 2383.d0 + rk(1:sizey,1) - rk(1:sizey,6) * 15050.d0 + &
                                                & rk(1:sizey,7) * 2133.d0 + rk(1:sizey,8) * 2250.d0 + rk(1:sizey,9) * 1125.d0 + &
                                                & rk(1:sizey,10) * 1800.d0)/4100.d0)
            rk(1:sizey,12)  = dydt (t, y + dt * (((deri - rk(1:sizey,7))/205.d0 + ((- rk(1:sizey,6) + rk(1:sizey,10)) * TWO + &
                                                & (- rk(1:sizey,8) + rk(1:sizey,9)))/41.d0) * THREE))
            rk(1:sizey,13)  = dydt (t + dt, y + dt * ((- deri * 1777.d0 + rk(1:sizey,1) - rk(1:sizey,6) * 14450.d0 + & 
                                                & rk(1:sizey,7) * 2193.d0 + rk(1:sizey,8) * 2550.d0 + rk(1:sizey,9) * 825.d0 + &
                                                & rk(1:sizey,10) * 1900.d0)/4100.d0 + rk(1:sizey,12)))

            rk(1:sizey,1) = rk(1:sizey,6) * 272.d0 + (rk(1:sizey,7) + rk(1:sizey,8)) * 216.d0 + (rk(1:sizey,9) + &
                                                & rk(1:sizey,10)) * 27.d0

            ynew = y + dt * ((deri + rk(1:sizey,11)) * 41.d0 + rk(1:sizey,1)) * C1_840
            yauxi = y + dt * (rk(1:sizey,1) + (rk(1:sizey,12) + rk(1:sizey,13)) * 41.d0) * C1_840

        end subroutine Fehlberg7_8
        
        subroutine Dormand_Prince8_7 (sizey, y, dydt, t, dt, deri, yauxi, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey),  intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: yauxi
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            rk(1:sizey,2)  = dydt (t + dt * C1_18, y + dt * (C1_18 * deri))
            rk(1:sizey,3)  = dydt (t + dt * C1_12, y + dt * (C1_48 * deri + 0.0625d0 * rk(1:sizey,2)))
            rk(1:sizey,4)  = dydt (t + dt * C1_8, y + dt * (0.03125d0 * deri + 0.09375d0 * rk(1:sizey,3)))
            rk(1:sizey,5)  = dydt (t + dt * 0.3125d0, y + dt * (0.3125d0 * deri - 1.171875d0 * rk(1:sizey,3) + &
                                                & 1.171875d0 * rk(1:sizey,4)))
            rk(1:sizey,6)  = dydt (t + dt * C3_8, y + dt * (0.0375d0 * deri + 0.1875d0 * rk(1:sizey,4) + 0.15d0 * rk(1:sizey,5)))
            rk(1:sizey,7)  = dydt (t + dt * 0.1475d0, y + dt * (29443841.d0/614563906.d0 * deri + &
                                                & 77736538.d0/692538347.d0 * rk(1:sizey,4) - &
                                                & 28693883.d0/1125000000.d0 * rk(1:sizey,5) + &
                                                & 23124283.d0/1800000000.d0 * rk(1:sizey,6)))
            rk(1:sizey,8)  = dydt (t + dt * 0.465d0, y + dt * (16016141.d0/946692911.d0 * deri + &
                                                & 61564180.d0/158732637.d0 * rk(1:sizey,4) + &
                                                & 22789713.d0/633445777.d0 * rk(1:sizey,5) + &
                                                & 545815736.d0/2771057229.d0 * rk(1:sizey,6) - &
                                                & 180193667.d0/1043307555.d0 * rk(1:sizey,7)))
            rk(1:sizey,9)  = dydt (t + dt * 5490023248.d0/9719169821.d0, y + dt * (39632708.d0/573591083.d0 * deri - &
                                                & 433636366.d0/683701615.d0 * rk(1:sizey,4) - &
                                                & 421739975.d0/2616292301.d0 * rk(1:sizey,5) + &
                                                & 100302831.d0/723423059.d0 * rk(1:sizey,6) + &
                                                & 790204164.d0/839813087.d0 * rk(1:sizey,7) + &
                                                & 800635310.d0/3783071287.d0 * rk(1:sizey,8)))
            rk(1:sizey,10) = dydt (t + dt * 0.65d0, y + dt * (246121993.d0/1340847787.d0 * deri - &
                                                & 37695042795.d0/15268766246.d0 * rk(1:sizey,4) - &
                                                & 309121744.d0/1061227803.d0 * rk(1:sizey,5) - &
                                                & 12992083.d0/490766935.d0 * rk(1:sizey,6) + &
                                                & 6005943493.d0/2108947869.d0 * rk(1:sizey,7) + &
                                                & 393006217.d0/1396673457.d0 * rk(1:sizey,8) + &
                                                & 123872331.d0/1001029789.d0 * rk(1:sizey,9)))
            rk(1:sizey,11) = dydt (t + dt * 1201146811.d0/1299019798.d0, y + dt * (-1028468189.d0/846180014.d0 * deri + &
                                                & 8478235783.d0/508512852.d0 * rk(1:sizey,4) + &
                                                & 1311729495.d0/1432422823.d0 * rk(1:sizey,5) - &
                                                & 10304129995.d0/1701304382.d0 * rk(1:sizey,6) - &
                                                & 48777925059.d0/3047939560.d0 * rk(1:sizey,7) + &
                                                & 15336726248.d0/1032824649.d0 * rk(1:sizey,8) - &
                                                & 45442868181.d0/3398467696.d0 * rk(1:sizey,9) + &
                                                & 3065993473.d0/597172653.d0 * rk(1:sizey,10)))
            rk(1:sizey,12) = dydt (t + dt, y + dt * (185892177.d0/718116043.d0 * deri - &
                                                & 3185094517.d0/667107341.d0 * rk(1:sizey,4) - &
                                                & 477755414.d0/1098053517.d0 * rk(1:sizey,5) - &
                                                & 703635378.d0/230739211.d0 * rk(1:sizey,6) + &
                                                & 5731566787.d0/1027545527.d0 * rk(1:sizey,7) + &
                                                & 5232866602.d0/850066563.d0 * rk(1:sizey,8) - &
                                                & 4093664535.d0/808688257.d0 * rk(1:sizey,9) + &
                                                & 3962137247.d0/1805957418.d0 * rk(1:sizey,10) + &
                                                & 65686358.d0/487910083.d0 * rk(1:sizey,11)))
            rk(1:sizey,13) = dydt (t + dt, y + dt * (403863854.d0/491063109.d0 * deri - &
                                                & 5068492393.d0/434740067.d0 * rk(1:sizey,4) - &
                                                & 411421997.d0/543043805.d0 * rk(1:sizey,5) + &
                                                & 652783627.d0/914296604.d0 * rk(1:sizey,6) + &
                                                & 11173962825.d0/925320556.d0 * rk(1:sizey,7) - &
                                                & 13158990841.d0/6184727034.d0 * rk(1:sizey,8) + &
                                                & 3936647629.d0/1978049680.d0 * rk(1:sizey,9) - &
                                                & 160528059.d0/685178525.d0 * rk(1:sizey,10) + &
                                                & 248638103.d0/1413531060.d0 * rk(1:sizey,11)))

            ynew = y + dt * (13451932.d0/455176623.d0 * deri - 808719846.d0/976000145.d0 * rk(1:sizey,6) + &
                                                & 1757004468.d0/5645159321.d0 * rk(1:sizey,7) + &
                                                & 656045339.d0/265891186.d0 * rk(1:sizey,8) - &
                                                & 3867574721.d0/1518517206.d0 * rk(1:sizey,9) + &
                                                & 465885868.d0/322736535.d0 * rk(1:sizey,10) + &
                                                & 53011238.d0/667516719.d0 * rk(1:sizey,11) + C2_45 * rk(1:sizey,12))
            yauxi = y + dt * (14005451.d0/335480064.d0 * deri - 59238493.d0/1068277825.d0 * rk(1:sizey,6) + &
                                                & 181606767.d0/758867731.d0 * rk(1:sizey,7) + &
                                                & 561292985.d0/797845732.d0 * rk(1:sizey,8) - &
                                                & 1041891430.d0/1371343529.d0 * rk(1:sizey,9) + &
                                                & 760417239.d0/1151165299.d0 * rk(1:sizey,10) + &
                                                & 118820643.d0/751138087.d0 * rk(1:sizey,11) - &
                                                & 528747749.d0/2220607170.d0 * rk(1:sizey,12) + C1_4 * rk(1:sizey,13))

        end subroutine Dormand_Prince8_7


        !------------------------------------------------
        !  Solver Embedded (adaptive timestep)
        !------------------------------------------------

        recursive subroutine solve_embedded (sizey, y, dydt, t, dt_adap, dt_used, deri, embedded, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t
            real(kind=8), intent(inout) :: dt_adap
            real(kind=8), intent(out) :: dt_used
            real(kind=8), dimension(sizey), intent(in) :: deri
            procedure (embedded_tem), pointer :: embedded
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            integer(kind=4), save :: iter = 0
            real(kind=8) :: e_calc, ratio

            iter = iter + 1
            dt_adap = max (dt_adap, DT_MIN_NOW)  

            ! yscal
            yscal(:sizey) = abs (y) + abs (dt_adap * deri(:sizey)) + SAFE_LOW          
            
            ! y(t, dt) -> yaux | ynew
            call embedded (sizey, y, dydt, t, dt_adap, deri, yaux(:sizey), ynew)
            
            ! Error
            e_calc = max (maxval (abs ((ynew - yaux(:sizey)) / yscal(:sizey))), SAFE_LOW)
            ratio = E_TOL / e_calc

            if (ratio > ONE) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (BETA * ratio**ONE_OSOL, MAX_DT_FACTOR)
                iter = 0

            else
                if (abs(dt_adap - DT_MIN_NOW) .le. E_TOL) then !E_TOL?
                    dt_used = DT_MIN_NOW
                    iter = 0

                else
                    dt_adap = dt_adap * min (BETA * ratio**ONE_OAUX, MAX_DT_FACTOR)

                    if ((dt_adap /= dt_adap) .or. (dt_adap < DT_MIN_NOW) .or. (iter == MAX_N_ITER)) then
                        dt_adap = DT_MIN_NOW
                        dt_used = DT_MIN_NOW
                        
                        call embedded (sizey, y, dydt, t, dt_adap, deri, yaux(:sizey), ynew)
                        iter = 0

                    else
                        call solve_embedded (sizey, y, dydt, t, dt_adap, dt_used, deri, embedded, ynew)

                    end if

                end if

            end if

        end subroutine solve_embedded
                                         
        
        ! ------------------------------------------------
        !  Call embedded integrator (adaptive timestep)
        ! ------------------------------------------------

        subroutine embedded_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
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
            logical :: has_check = .False.

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
                                    
                call solve_embedded (sizey, ycaller(:sizey), dydt, time, dt_adap, dt_used, &
                                    & der(:sizey), embedded_ptr, ynew)

                time = time + dt_used
            end do

        end subroutine embedded_caller


        ! ------------------------------------------------
        !  Call embedded integrator (fixed timestep)
        ! ------------------------------------------------

        subroutine embedded_fixed_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
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
            logical :: has_check = .False.

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

                call embedded_ptr (sizey, ycaller(:sizey), dydt, time, dt_adap, der(:sizey), yaux(:sizey), ynew)  ! yaux is dummy

                time = time + dt_adap
            end do

        end subroutine embedded_fixed_caller
    
end module embedded