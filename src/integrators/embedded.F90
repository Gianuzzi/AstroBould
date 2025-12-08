!> Module with embedded RK integrators
module embedded
    use shared
    implicit none
    private
    public :: init_embedded, free_embedded, embedded_caller, embedded_fixed_caller

    ! Workspace arrays
    real(wp), allocatable :: rk(:, :)  ! integrator
    real(wp), allocatable :: yscal(:)  ! solver
    real(wp), allocatable :: yaux(:)  ! solver
    real(wp), allocatable :: ycaller(:)  ! caller

    ! Constants
    integer(kind=4), parameter :: MAX_N_ITER = 350
    real(wp), parameter :: MAX_DT_FACTOR = 5.0e0_wp

    ! Variables to set
    integer(kind=4) :: OSOL
    real(wp) :: ONE_OSOL
    integer(kind=4) :: OAUX
    real(wp) :: ONE_OAUX

    ! Pointer to embedded used
    procedure(embedded_tem), pointer :: embedded_ptr => null()

    abstract interface

        subroutine embedded_tem(sizey, y, dydt, t, dt, deri, yauxi, ynew)
            import :: wp
            import :: dydt_tem
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(wp), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: t, dt
            real(wp), dimension(sizey), intent(in) :: deri
            real(wp), dimension(sizey), intent(out) :: yauxi
            real(wp), dimension(sizey), intent(out) :: ynew
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
            print *, "Method not available:", which
            stop 1
        end if

        ! Constants
        ONE_OSOL = ONE/real(OSOL, kind=wp)
        ONE_OAUX = ONE/real(OAUX, kind=wp)

        ! integrator
        allocate (rk(sizey, n_stages))

        ! solver
        allocate (yscal(sizey))
        allocate (yaux(sizey))

        ! caller
        allocate (ycaller(sizey))
    end subroutine init_embedded

    subroutine free_embedded()
        implicit none
        if (allocated(rk)) deallocate (rk)
        if (allocated(yscal)) deallocate (yscal)
        if (allocated(yaux)) deallocate (yaux)
        if (allocated(ycaller)) deallocate (ycaller)
        nullify (embedded_ptr)
    end subroutine free_embedded

        !!!! Embedded subroutines

    subroutine Fehlberg1_2(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C1_2, y + dt*deri*C1_2)

        ynew = y + dt*rk(1:sizey, 2)

        yauxi = y + dt*C1_2*(deri + rk(1:sizey, 2))

    end subroutine Fehlberg1_2

    subroutine Heun_Euler2_1(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt, y + dt*deri)

        ynew = y + dt*(deri + rk(1:sizey, 2))*C1_2
        yauxi = y + dt*deri

    end subroutine Heun_Euler2_1

    subroutine Fehlberg2_1(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C3_4, y + dt*deri*C3_4)

        ynew = y + dt*(deri*C1_3 + rk(1:sizey, 2)*C2_3)

        yauxi = y + dt*rk(1:sizey, 2)

    end subroutine Fehlberg2_1

    subroutine Bogacki_Shampine3_2(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C1_2, y + dt*deri*C1_2)
        rk(1:sizey, 3) = dydt(t + dt*C3_4, y + dt*rk(1:sizey, 2)*C3_4)

        ynew = y + dt*(deri*TWO + rk(1:sizey, 2)*THREE + rk(1:sizey, 3)*FOUR)*C1_9

        rk(1:sizey, 4) = dydt(t + dt, ynew)

        yauxi = y + dt*(deri*7.e0_wp/24.e0_wp + rk(1:sizey, 2)*C1_4 + rk(1:sizey, 3)*C1_3 + rk(1:sizey, 4)*C1_8)

    end subroutine Bogacki_Shampine3_2

    subroutine Zonneveld4_3(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew
        real(wp) :: dt_2

        dt_2 = dt*C1_2

        rk(1:sizey, 2) = dydt(t + dt_2, y + dt*deri*C1_2)
        rk(1:sizey, 3) = dydt(t + dt_2, y + dt*rk(1:sizey, 2)*C1_2)
        rk(1:sizey, 4) = dydt(t + dt, y + dt*rk(1:sizey, 3))
        rk(1:sizey, 5) = dydt(t + dt*C3_4, y + dt*( &
            & deri*FIVE + rk(1:sizey, 2)*7.e0_wp + rk(1:sizey, 3)*13.e0_wp - rk(1:sizey, 4))*0.03125e0_wp)

        rk(1:sizey, 1) = rk(1:sizey, 2) + rk(1:sizey, 3)  ! Used as kaux

        ynew = y + dt*(deri + rk(1:sizey, 1)*TWO + rk(1:sizey, 4))*C1_6
        yauxi = y + dt*(-deri*THREE + rk(1:sizey, 1)*14.e0_wp + rk(1:sizey, 4)*13.e0_wp - &
                                        & rk(1:sizey, 5)*32.e0_wp)*C1_6

    end subroutine Zonneveld4_3

    subroutine Merson4_3(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C1_3, y + dt*deri*C1_3)
        rk(1:sizey, 3) = dydt(t + dt*C1_3, y + dt*(deri + rk(1:sizey, 2))*C1_6)
        rk(1:sizey, 4) = dydt(t + dt*C1_2, y + dt*(deri + rk(1:sizey, 3)*THREE)*C1_8)
        rk(1:sizey, 5) = dydt(t + dt, y + dt*(deri - rk(1:sizey, 3)*THREE + rk(1:sizey, 4)*FOUR)*C1_2)

        ynew = y + dt*(deri + rk(1:sizey, 4)*FOUR + rk(1:sizey, 5))*C1_6
        yauxi = y + dt*(deri + rk(1:sizey, 3)*THREE + rk(1:sizey, 4)*FOUR + rk(1:sizey, 5)*TWO)*0.1e0_wp

    end subroutine Merson4_3

    subroutine Fehlberg4_5(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C1_4, y + dt*deri*C1_4)
        rk(1:sizey, 3) = dydt(t + dt*C3_8, y + dt*(deri*THREE + rk(1:sizey, 2)*9.e0_wp)*0.03125e0_wp)
        rk(1:sizey, 4) = dydt(t + dt*12.e0_wp/13.e0_wp, y + dt*(deri*1932.e0_wp - rk(1:sizey, 2)*7200.e0_wp + &
                                            & rk(1:sizey, 3)*7296.e0_wp)/2197.e0_wp)
        rk(1:sizey, 5) = dydt(t + dt, y + dt*((deri*8341.e0_wp + rk(1:sizey, 3)*29440.e0_wp - &
                                            & rk(1:sizey, 4)*845.e0_wp)/4104.e0_wp - rk(1:sizey, 2)*8.e0_wp))
        rk(1:sizey, 6) = dydt(t + dt*C1_2, y + dt*((-deri*1216.e0_wp + rk(1:sizey, 4)*1859.e0_wp)/4104.e0_wp + &
                                            & rk(1:sizey, 2)*TWO - rk(1:sizey, 3)*3544.e0_wp/2565.e0_wp - &
                                            & rk(1:sizey, 5)*0.275e0_wp))

        ynew = y + dt*((deri*475.e0_wp + rk(1:sizey, 4)*2197.e0_wp)/4104.e0_wp + rk(1:sizey, 3)*1408.e0_wp/2565.e0_wp - &
                                            & rk(1:sizey, 5)*C1_5)
        yauxi = y + dt*((deri*6688.e0_wp + rk(1:sizey, 4)*28561.e0_wp + rk(1:sizey, 6)*2052.e0_wp)/56430.e0_wp + &
                                            & rk(1:sizey, 3)*6656.e0_wp/12825.e0_wp - rk(1:sizey, 5)*0.18e0_wp)

    end subroutine Fehlberg4_5

    subroutine Cash_Karp5_4(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C1_5, y + dt*deri*C1_5)
        rk(1:sizey, 3) = dydt(t + dt*0.3e0_wp, y + dt*(deri*THREE + rk(1:sizey, 2)*9.e0_wp)*0.025e0_wp)
        rk(1:sizey, 4) = dydt(t + dt*0.6e0_wp, y + dt*(deri*THREE - rk(1:sizey, 2)*9.e0_wp + &
                                            & rk(1:sizey, 3)*12.e0_wp)*0.1e0_wp)
        rk(1:sizey, 5) = dydt(t + dt, y + dt*(-deri*11.e0_wp + rk(1:sizey, 2)*135.e0_wp - rk(1:sizey, 3)*140.e0_wp + &
                                            & rk(1:sizey, 4)*70.e0_wp)/54.e0_wp)
        rk(1:sizey, 6) = dydt(t + dt*0.875e0_wp, y + dt*(deri*3262.e0_wp + rk(1:sizey, 2)*37800.e0_wp + &
                                            & rk(1:sizey, 3)*4600.e0_wp + rk(1:sizey, 4)*44275.e0_wp + &
                                            & rk(1:sizey, 5)*6831.e0_wp)/110592.e0_wp)

        ynew = y + dt*(deri*37.e0_wp/378.e0_wp + rk(1:sizey, 3)*250.e0_wp/621.e0_wp + &
                                            & rk(1:sizey, 4)*125.e0_wp/594.e0_wp + rk(1:sizey, 6)*512.e0_wp/1771.e0_wp)
        yauxi = y + dt*((deri*5650.e0_wp + rk(1:sizey, 4)*13525.e0_wp)/55296.e0_wp + &
                                            & rk(1:sizey, 3)*18575.e0_wp/48384.e0_wp + &
                                            & rk(1:sizey, 5)*277.e0_wp/14336.e0_wp + rk(1:sizey, 6)*C1_4)

    end subroutine Cash_Karp5_4

    subroutine Dormand_Prince5_4(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C1_5, y + dt*deri*C1_5)
        rk(1:sizey, 3) = dydt(t + dt*0.3e0_wp, y + dt*(deri*0.075e0_wp + rk(1:sizey, 2)*0.225e0_wp))
        rk(1:sizey, 4) = dydt(t + dt*0.8e0_wp, y + dt*(deri*44.e0_wp - rk(1:sizey, 2)*168.e0_wp + &
                                            & rk(1:sizey, 3)*160.e0_wp)/45.e0_wp)
        rk(1:sizey, 5) = dydt(t + dt*C8_9, y + dt*(deri*19372.e0_wp - rk(1:sizey, 2)*76080.e0_wp + &
                                            & rk(1:sizey, 3)*64448.e0_wp - rk(1:sizey, 4)*1908.e0_wp)/6561.e0_wp)
        rk(1:sizey, 6) = dydt(t + dt, y + dt*((deri*9017.e0_wp - rk(1:sizey, 2)*34080.e0_wp + &
                                            & rk(1:sizey, 4)*882.e0_wp)/3168.e0_wp + &
                                            & rk(1:sizey, 3)*46732.e0_wp/5247.e0_wp - rk(1:sizey, 5)*5103.e0_wp/18656.e0_wp))

        ynew = y + dt*((deri*35.e0_wp + rk(1:sizey, 4)*250.e0_wp)/384.e0_wp + rk(1:sizey, 3)*500.e0_wp/1113.e0_wp - &
                                            & rk(1:sizey, 5)*2187.e0_wp/6784.e0_wp + rk(1:sizey, 6)*11.e0_wp/84.e0_wp)

        rk(1:sizey, 7) = dydt(t + dt, ynew)

        yauxi = y + dt*((deri*5179.e0_wp + rk(1:sizey, 4)*35370.e0_wp)/57600.e0_wp + &
                                            & rk(1:sizey, 3)*7571.e0_wp/16695.e0_wp - &
                                            & rk(1:sizey, 5)*92097.e0_wp/339200.e0_wp + &
                                            & rk(1:sizey, 6)*187.e0_wp/2100.e0_wp + rk(1:sizey, 7)*0.025e0_wp)

    end subroutine Dormand_Prince5_4

    subroutine Verner6_5(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C1_6, y + dt*deri*C1_6)
        rk(1:sizey, 3) = dydt(t + dt*C4_15, y + dt*(deri*FOUR + rk(1:sizey, 2)*16.e0_wp)/75.e0_wp)
        rk(1:sizey, 4) = dydt(t + dt*C2_3, y + dt*(deri*FIVE - rk(1:sizey, 2)*16.e0_wp + &
                                            & rk(1:sizey, 3)*15.e0_wp)*C1_6)
        rk(1:sizey, 5) = dydt(t + dt*C5_6, y + dt*((-deri*165.e0_wp - rk(1:sizey, 3)*425.e0_wp)/64.e0_wp + &
                                            & (rk(1:sizey, 2)*880.e0_wp + rk(1:sizey, 4)*85.e0_wp)/96.e0_wp))
        rk(1:sizey, 6) = dydt(t + dt, y + dt*((deri*612.e0_wp + rk(1:sizey, 5)*88.e0_wp)/255.e0_wp - &
                                            & rk(1:sizey, 2)*8.e0_wp + &
                                            &  (rk(1:sizey, 3)*4015.e0_wp - rk(1:sizey, 4)*187.e0_wp)/612.e0_wp))
        rk(1:sizey, 7) = dydt(t + dt*C1_15, y + dt*((-deri*8263.e0_wp + rk(1:sizey, 2)*24800.e0_wp)/15000.e0_wp - &
                                            & rk(1:sizey, 3)*643.e0_wp/680.e0_wp - rk(1:sizey, 4)*0.324e0_wp + &
                                            & rk(1:sizey, 5)*2484.e0_wp/10625.e0_wp))
        rk(1:sizey, 8) = dydt(t + dt, y + dt*(deri*3501.e0_wp/1720.e0_wp + (297275.e0_wp*rk(1:sizey, 3) - &
                                            & 367200*rk(1:sizey, 2))/52632.e0_wp - rk(1:sizey, 4)*319.e0_wp/2322.e0_wp + &
                                            & rk(1:sizey, 5)*24068.e0_wp/84065.e0_wp + &
                                            & rk(1:sizey, 7)*3850.e0_wp/26703.e0_wp))

        ynew = y + dt*(deri*0.075e0_wp + rk(1:sizey, 3)*875.e0_wp/2244.e0_wp + (rk(1:sizey, 4)*3703.e0_wp + &
                                            & rk(1:sizey, 7)*125.e0_wp)/11592.e0_wp + rk(1:sizey, 5)*264.e0_wp/1955.e0_wp + &
                                            & rk(1:sizey, 8)*43.e0_wp/616.e0_wp)
        yauxi = y + dt*((deri*13.e0_wp + rk(1:sizey, 4)*50.e0_wp)/160.e0_wp + (rk(1:sizey, 3)*2375.e0_wp + &
                                            & rk(1:sizey, 6)*408.e0_wp)/5984.e0_wp + rk(1:sizey, 5)*12.e0_wp/85.e0_wp)

    end subroutine Verner6_5

    subroutine Fehlberg7_8(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C2_27, y + dt*deri*C2_27)
        rk(1:sizey, 3) = dydt(t + dt*C1_9, y + dt*(deri + rk(1:sizey, 2)*THREE)/36.e0_wp)
        rk(1:sizey, 4) = dydt(t + dt*C1_6, y + dt*(deri + rk(1:sizey, 3)*THREE)/24.e0_wp)
        rk(1:sizey, 5) = dydt(t + dt*C5_12, y + dt*(deri*C5_12 + (rk(1:sizey, 4) - rk(1:sizey, 3))*1.5625e0_wp))
        rk(1:sizey, 6) = dydt(t + dt*C1_2, y + dt*(deri + rk(1:sizey, 4)*FIVE + rk(1:sizey, 5)*FOUR)*0.05e0_wp)
        rk(1:sizey, 7) = dydt(t + dt*C5_6, y + dt*(-deri*25.e0_wp + rk(1:sizey, 4)*125.e0_wp - &
                                            & rk(1:sizey, 5)*260.e0_wp + rk(1:sizey, 6)*250.e0_wp)/108.e0_wp)
        rk(1:sizey, 8) = dydt(t + dt*C1_6, y + dt*(deri*93.e0_wp + rk(1:sizey, 5)*244.e0_wp - &
                                            & rk(1:sizey, 6)*200.e0_wp + rk(1:sizey, 7)*13.e0_wp)/900.e0_wp)
        rk(1:sizey, 9) = dydt(t + dt*C2_3, y + dt*(deri*TWO + (-rk(1:sizey, 4)*795.e0_wp + &
                                            & rk(1:sizey, 5)*1408.e0_wp - rk(1:sizey, 6)*1070.e0_wp + &
                                            & rk(1:sizey, 7)*67.e0_wp)/90.e0_wp + &
                                            & rk(1:sizey, 8)*THREE))
        rk(1:sizey, 10) = dydt(t + dt*C1_3, y + dt*((-deri*91.e0_wp + rk(1:sizey, 4)*23.e0_wp + &
                                            & rk(1:sizey, 6)*622.e0_wp)/108.e0_wp - rk(1:sizey, 5)*976.e0_wp/135.e0_wp + &
                                            & (-rk(1:sizey, 7)*19.e0_wp + rk(1:sizey, 8)*170.e0_wp - &
                                            & rk(1:sizey, 9)*5.e0_wp)/60.e0_wp))

        rk(1:sizey, 1) = -rk(1:sizey, 4)*8525.e0_wp + rk(1:sizey, 5)*17984.e0_wp  !! Used as kaux

        rk(1:sizey, 11) = dydt(t + dt, y + dt*(deri*2383.e0_wp + rk(1:sizey, 1) - rk(1:sizey, 6)*15050.e0_wp + &
                                            & rk(1:sizey, 7)*2133.e0_wp + rk(1:sizey, 8)*2250.e0_wp + &
                                            & rk(1:sizey, 9)*1125.e0_wp + &
                                            & rk(1:sizey, 10)*1800.e0_wp)/4100.e0_wp)
        rk(1:sizey, 12) = dydt(t, y + dt*(((deri - rk(1:sizey, 7))/205.e0_wp + ((-rk(1:sizey, 6) + rk(1:sizey, 10))*TWO + &
                                            & (-rk(1:sizey, 8) + rk(1:sizey, 9)))/41.e0_wp)*THREE))
        rk(1:sizey, 13) = dydt(t + dt, y + dt*((-deri*1777.e0_wp + rk(1:sizey, 1) - rk(1:sizey, 6)*14450.e0_wp + &
                                            & rk(1:sizey, 7)*2193.e0_wp + rk(1:sizey, 8)*2550.e0_wp + &
                                            & rk(1:sizey, 9)*825.e0_wp + &
                                            & rk(1:sizey, 10)*1900.e0_wp)/4100.e0_wp + rk(1:sizey, 12)))

        rk(1:sizey, 1) = rk(1:sizey, 6)*272.e0_wp + (rk(1:sizey, 7) + rk(1:sizey, 8))*216.e0_wp + (rk(1:sizey, 9) + &
                                            & rk(1:sizey, 10))*27.e0_wp

        ynew = y + dt*((deri + rk(1:sizey, 11))*41.e0_wp + rk(1:sizey, 1))*C1_840
        yauxi = y + dt*(rk(1:sizey, 1) + (rk(1:sizey, 12) + rk(1:sizey, 13))*41.e0_wp)*C1_840

    end subroutine Fehlberg7_8

    subroutine Dormand_Prince8_7(sizey, y, dydt, t, dt, deri, yauxi, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t, dt
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: yauxi
        real(wp), dimension(sizey), intent(out) :: ynew

        rk(1:sizey, 2) = dydt(t + dt*C1_18, y + dt*(C1_18*deri))
        rk(1:sizey, 3) = dydt(t + dt*C1_12, y + dt*(C1_48*deri + 0.0625e0_wp*rk(1:sizey, 2)))
        rk(1:sizey, 4) = dydt(t + dt*C1_8, y + dt*(0.03125e0_wp*deri + 0.09375e0_wp*rk(1:sizey, 3)))
        rk(1:sizey, 5) = dydt(t + dt*0.3125e0_wp, y + dt*(0.3125e0_wp*deri - 1.171875e0_wp*rk(1:sizey, 3) + &
                                            & 1.171875e0_wp*rk(1:sizey, 4)))
        rk(1:sizey, 6) = dydt(t + dt*C3_8, y + dt*(0.0375e0_wp*deri + 0.1875e0_wp*rk(1:sizey, 4) + &
                                            & 0.15e0_wp*rk(1:sizey, 5)))
        rk(1:sizey, 7) = dydt(t + dt*0.1475e0_wp, y + dt*(29443841.e0_wp/614563906.e0_wp*deri + &
                                            & 77736538.e0_wp/692538347.e0_wp*rk(1:sizey, 4) - &
                                            & 28693883.e0_wp/1125000000.e0_wp*rk(1:sizey, 5) + &
                                            & 23124283.e0_wp/1800000000.e0_wp*rk(1:sizey, 6)))
        rk(1:sizey, 8) = dydt(t + dt*0.465e0_wp, y + dt*(16016141.e0_wp/946692911.e0_wp*deri + &
                                            & 61564180.e0_wp/158732637.e0_wp*rk(1:sizey, 4) + &
                                            & 22789713.e0_wp/633445777.e0_wp*rk(1:sizey, 5) + &
                                            & 545815736.e0_wp/2771057229.e0_wp*rk(1:sizey, 6) - &
                                            & 180193667.e0_wp/1043307555.e0_wp*rk(1:sizey, 7)))
        rk(1:sizey, 9) = dydt(t + dt*5490023248.e0_wp/9719169821.e0_wp, y + dt*(39632708.e0_wp/573591083.e0_wp*deri - &
                                            & 433636366.e0_wp/683701615.e0_wp*rk(1:sizey, 4) - &
                                            & 421739975.e0_wp/2616292301.e0_wp*rk(1:sizey, 5) + &
                                            & 100302831.e0_wp/723423059.e0_wp*rk(1:sizey, 6) + &
                                            & 790204164.e0_wp/839813087.e0_wp*rk(1:sizey, 7) + &
                                            & 800635310.e0_wp/3783071287.e0_wp*rk(1:sizey, 8)))
        rk(1:sizey, 10) = dydt(t + dt*0.65e0_wp, y + dt*(246121993.e0_wp/1340847787.e0_wp*deri - &
                                            & 37695042795.e0_wp/15268766246.e0_wp*rk(1:sizey, 4) - &
                                            & 309121744.e0_wp/1061227803.e0_wp*rk(1:sizey, 5) - &
                                            & 12992083.e0_wp/490766935.e0_wp*rk(1:sizey, 6) + &
                                            & 6005943493.e0_wp/2108947869.e0_wp*rk(1:sizey, 7) + &
                                            & 393006217.e0_wp/1396673457.e0_wp*rk(1:sizey, 8) + &
                                            & 123872331.e0_wp/1001029789.e0_wp*rk(1:sizey, 9)))
        rk(1:sizey, 11) = dydt(t + dt*1201146811.e0_wp/1299019798.e0_wp, y + dt* &
                                            & (-1028468189.e0_wp/846180014.e0_wp*deri + &
                                            & 8478235783.e0_wp/508512852.e0_wp*rk(1:sizey, 4) + &
                                            & 1311729495.e0_wp/1432422823.e0_wp*rk(1:sizey, 5) - &
                                            & 10304129995.e0_wp/1701304382.e0_wp*rk(1:sizey, 6) - &
                                            & 48777925059.e0_wp/3047939560.e0_wp*rk(1:sizey, 7) + &
                                            & 15336726248.e0_wp/1032824649.e0_wp*rk(1:sizey, 8) - &
                                            & 45442868181.e0_wp/3398467696.e0_wp*rk(1:sizey, 9) + &
                                            & 3065993473.e0_wp/597172653.e0_wp*rk(1:sizey, 10)))
        rk(1:sizey, 12) = dydt(t + dt, y + dt*(185892177.e0_wp/718116043.e0_wp*deri - &
                                            & 3185094517.e0_wp/667107341.e0_wp*rk(1:sizey, 4) - &
                                            & 477755414.e0_wp/1098053517.e0_wp*rk(1:sizey, 5) - &
                                            & 703635378.e0_wp/230739211.e0_wp*rk(1:sizey, 6) + &
                                            & 5731566787.e0_wp/1027545527.e0_wp*rk(1:sizey, 7) + &
                                            & 5232866602.e0_wp/850066563.e0_wp*rk(1:sizey, 8) - &
                                            & 4093664535.e0_wp/808688257.e0_wp*rk(1:sizey, 9) + &
                                            & 3962137247.e0_wp/1805957418.e0_wp*rk(1:sizey, 10) + &
                                            & 65686358.e0_wp/487910083.e0_wp*rk(1:sizey, 11)))
        rk(1:sizey, 13) = dydt(t + dt, y + dt*(403863854.e0_wp/491063109.e0_wp*deri - &
                                            & 5068492393.e0_wp/434740067.e0_wp*rk(1:sizey, 4) - &
                                            & 411421997.e0_wp/543043805.e0_wp*rk(1:sizey, 5) + &
                                            & 652783627.e0_wp/914296604.e0_wp*rk(1:sizey, 6) + &
                                            & 11173962825.e0_wp/925320556.e0_wp*rk(1:sizey, 7) - &
                                            & 13158990841.e0_wp/6184727034.e0_wp*rk(1:sizey, 8) + &
                                            & 3936647629.e0_wp/1978049680.e0_wp*rk(1:sizey, 9) - &
                                            & 160528059.e0_wp/685178525.e0_wp*rk(1:sizey, 10) + &
                                            & 248638103.e0_wp/1413531060.e0_wp*rk(1:sizey, 11)))

        ynew = y + dt*(13451932.e0_wp/455176623.e0_wp*deri - 808719846.e0_wp/976000145.e0_wp*rk(1:sizey, 6) + &
                                            & 1757004468.e0_wp/5645159321.e0_wp*rk(1:sizey, 7) + &
                                            & 656045339.e0_wp/265891186.e0_wp*rk(1:sizey, 8) - &
                                            & 3867574721.e0_wp/1518517206.e0_wp*rk(1:sizey, 9) + &
                                            & 465885868.e0_wp/322736535.e0_wp*rk(1:sizey, 10) + &
                                            & 53011238.e0_wp/667516719.e0_wp*rk(1:sizey, 11) + C2_45*rk(1:sizey, 12))
        yauxi = y + dt*(14005451.e0_wp/335480064.e0_wp*deri - 59238493.e0_wp/1068277825.e0_wp*rk(1:sizey, 6) + &
                                            & 181606767.e0_wp/758867731.e0_wp*rk(1:sizey, 7) + &
                                            & 561292985.e0_wp/797845732.e0_wp*rk(1:sizey, 8) - &
                                            & 1041891430.e0_wp/1371343529.e0_wp*rk(1:sizey, 9) + &
                                            & 760417239.e0_wp/1151165299.e0_wp*rk(1:sizey, 10) + &
                                            & 118820643.e0_wp/751138087.e0_wp*rk(1:sizey, 11) - &
                                            & 528747749.e0_wp/2220607170.e0_wp*rk(1:sizey, 12) + C1_4*rk(1:sizey, 13))

    end subroutine Dormand_Prince8_7

    !------------------------------------------------
    !  Solver Embedded (adaptive timestep)
    !------------------------------------------------

    recursive subroutine solve_embedded(sizey, y, dydt, t, dt_adap, dt_used, deri, embedded, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt
        real(wp), intent(in) :: t
        real(wp), intent(inout) :: dt_adap
        real(wp), intent(out) :: dt_used
        real(wp), dimension(sizey), intent(in) :: deri
        procedure(embedded_tem), pointer :: embedded
        real(wp), dimension(sizey), intent(out) :: ynew

        integer(kind=4), save :: iter = 0
        real(wp) :: e_calc, ratio

        iter = iter + 1
        dt_adap = max(dt_adap, DT_MIN_NOW)

        ! yscal
        yscal(:sizey) = abs(y) + abs(dt_adap*deri(:sizey)) + SAFE_LOW

        ! y(t, dt) -> yaux | ynew
        call embedded(sizey, y, dydt, t, dt_adap, deri, yaux(:sizey), ynew)

        ! Error
        e_calc = max(maxval(abs((ynew - yaux(:sizey))/yscal(:sizey))), SAFE_LOW)
        ratio = E_TOL/e_calc

        if (ratio > ONE) then
            dt_used = dt_adap
            dt_adap = dt_adap*min(BETA*ratio**ONE_OSOL, MAX_DT_FACTOR)
            iter = 0

        else
            if (abs(dt_adap - DT_MIN_NOW) .le. E_TOL) then !E_TOL?
                dt_used = DT_MIN_NOW
                iter = 0

            else
                dt_adap = dt_adap*min(BETA*ratio**ONE_OAUX, MAX_DT_FACTOR)

                if ((dt_adap /= dt_adap) .or. (dt_adap < DT_MIN_NOW) .or. (iter == MAX_N_ITER)) then
                    dt_adap = DT_MIN_NOW
                    dt_used = DT_MIN_NOW

                    call embedded(sizey, y, dydt, t, dt_adap, deri, yaux(:sizey), ynew)
                    iter = 0

                else
                    call solve_embedded(sizey, y, dydt, t, dt_adap, dt_used, deri, embedded, ynew)

                end if

            end if

        end if

    end subroutine solve_embedded

    ! ------------------------------------------------
    !  Call embedded integrator (adaptive timestep)
    ! ------------------------------------------------

    subroutine embedded_caller(t, y, dt_adap, dydt, dt, ynew, check_fun)
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

            ycaller(:sizey) = ynew

            dt_adap = min(dt_adap, t_end - time)
            DT_MIN_NOW = min(DT_MIN, dt_adap)

            der(:sizey) = dydt(time, ycaller(:sizey))

            call solve_embedded(sizey, ycaller(:sizey), dydt, time, dt_adap, dt_used, der(:sizey), embedded_ptr, ynew)

            time = time + dt_used
        end do

    end subroutine embedded_caller

    ! ------------------------------------------------
    !  Call embedded integrator (fixed timestep)
    ! ------------------------------------------------

    subroutine embedded_fixed_caller(t, y, dt_adap, dydt, dt, ynew, check_fun)
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

            ycaller(:sizey) = ynew

            dt_adap = min(dt_adap, t_end - time)
            DT_MIN_NOW = min(DT_MIN, dt_adap)

            der(:sizey) = dydt(time, ycaller(:sizey))

            call embedded_ptr(sizey, ycaller(:sizey), dydt, time, dt_adap, der(:sizey), yaux(:sizey), ynew)  ! yaux is dummy

            time = time + dt_adap
        end do

    end subroutine embedded_fixed_caller

end module embedded
