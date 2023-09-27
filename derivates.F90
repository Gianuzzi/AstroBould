module derivates
    use const
    use parameters
    use gravity
    use stokes
    implicit none
    !y = /omega, mi, Ri, xi, yi, vxi, vyi, .../
    real(kind=8) :: yb(1 + 6 * (Ntot+1)), ybnew(1 + 6 * (Ntot+1))
    real(kind=8) :: ya(1 + 6 * (Ntot+1)), yanew(1 + 6 * (Ntot+1))
    real(kind=8) :: yr(1 + 6 * (Ntot+1)), yrnew(1 + 6 * (Ntot+1))

    !! Dampings
    real(kind=8), parameter :: tau_o = inf ! Omega
    real(kind=8), parameter :: tau_m = inf ! M0

    procedure (dydt_temp), pointer :: dydt => null ()
    type dydt_ptr
        procedure(dydt_temp), pointer, nopass :: f_i => null ()
    end type dydt_ptr
    abstract interface
        ! Here must be every f_i defined explicitly
        function dydt_temp (t, y) result (der) 
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size (y))      :: der
        end function dydt_temp
    end interface

    contains

        function domega_dt(t, omega) result(domega)
            implicit none
            real(kind=8), intent(in) :: t, omega
            real(kind=8) :: domega
            if (tau_o >= inf) then
                domega = cero
            else
                domega = -exp(-t / tau_o) * omega / tau_o
            end if
        end function domega_dt

        function dm0_dt(t, m0) result(dm0)
            implicit none
            real(kind=8), intent(in) :: t, m0
            real(kind=8) :: dm0
            if (tau_m >= inf) then
                dm0 = cero
            else
                dm0 = -exp(-t / tau_m) * m0 / tau_m
            end if
        end function dm0_dt

        function dydt_bar_ex(t, y) result(dydt)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: rib(0:Nboul,2), rb(2), vb(2)
            integer(kind=4) :: i, ineqs

            !y = /omega, mi, Ri, xi, yi, vxi, vyi, .../
            dydt = cero

            do i = 0, Nboul
                rib(i,1) = cos(omega * t + theta_b(i)) * r_b(i)
                rib(i,2) = sin(omega * t + theta_b(i)) * r_b(i)
                ! dydt(ineqs+4) = - omega * rib(i,2)   ! Velocidad x del boulder i
                ! dydt(ineqs+5) = omega * rib(i,1)   ! Velocidad y del boulder i
                ! dydt(ineqs+6 : ineqs+7) = - omega2 * rib(i,:)
            end do

            do i = Nboul+1, Ntot
                ineqs = i*neqs
                rb = y(ineqs+4 : ineqs+5)
                vb = y(ineqs+6 : ineqs+7)
                dydt(ineqs+4 : ineqs+5) = vb
                call accbar(m(0:Nboul), rb, rib, dydt(ineqs+6 : ineqs+7))
                call accsto(t, rb, vb, dydt(ineqs+6 : ineqs+7))
            end do
        end function dydt_bar_ex

        function dydt_bar_im(t, y) result(dydt)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: rib(0:Nboul,2), vib(0:Nboul,2), rb(2), vb(2)
            real(kind=8) :: m(0:Ntot), omega2, radius(0:Nboul)
            integer(kind=4) :: i, ineqs

            !y = /omega, mi, Ri, xi, yi, vxi, vyi, .../
            dydt = cero
            dydt(1) = domega_dt(t, y(1))
            dydt(2) = dm0_dt(t, y(2))

            omega2 = y(1) * y(1)
            do i = 0, Nboul ! Integradores usan desde i=1
                ineqs = i * neqs
                m(i)      = y(ineqs+2)
                radius(i) = y(ineqs+3)
                rib(i,:)  = y(ineqs+4 : ineqs+5)
                vib(i,:)  = y(ineqs+6 : ineqs+7)
                dydt(ineqs+4 : ineqs+5) = y(ineqs+6 : ineqs+7)
                dydt(ineqs+6 : ineqs+7) = -omega2 * rib(i,:) + dydt(1) * (/-rib(i,2), rib(i,1)/) ! a_i = -w^2 * r_i + dw/dt * (-ry,rx)  
            end do

            do i = Nboul+1, Ntot
                ineqs = i * neqs
                m(i) = y(i*neqs+2)
                rb   = y(ineqs+4 : ineqs+5)
                vb   = y(ineqs+6 : ineqs+7)
                dydt(ineqs+4 : ineqs+5) = vb
                call accbar(m(0:Nboul), rb, rib, dydt(ineqs+6 : ineqs+7))
                call accsto(t, rb, vb, dydt(ineqs+6 : ineqs+7))
            end do
        end function dydt_bar_im

        function dydt_ast_ex(t, y) result(dydt)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: ra(2), va(2)
            real(kind=8) :: ria(1:Nboul,2), via(1:Nboul,2), aia(1:Nboul,2)
            real(kind=8) :: rcm(2), vcm(2)
            integer(kind=4) :: i, ineqs
            
            !y = /omega, mi, Ri, xi, yi, vxi, vyi, .../
            dydt = cero

            do i = 1, Nboul
                ria(i,1) = cos(omega * t + theta_a(i))
                ria(i,2) = sin(omega * t + theta_a(i))
            end do
            ria = ria * R0
            via(:,1) = -omega * ria(:,2)
            via(:,2) =  omega * ria(:,1)
            aia = -omega2 * ria

           rcm(1) = cos(omega * t + theta_acm)
           rcm(2) = sin(omega * t + theta_acm)
           vcm(1) = -omega * rcm(2)
           vcm(2) =  omega * rcm(1)

            do i = Nboul+1, Ntot
                ineqs = i * neqs
                ra = y(ineqs+4 : ineqs+5)
                va = y(ineqs+6 : ineqs+7)
                dydt(ineqs+4 : ineqs+5) = va
                call accast(omega, m(0:Nboul), ra, ria, dydt(ineqs+6 : ineqs+7))
                call accsto(t, ra-rcm, va-vcm, dydt(ineqs+6 : ineqs+7))
            end do
        end function dydt_ast_ex

        function dydt_ast_im(t, y) result(dydt)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8)    :: ria(1:Nboul,2), via(1:Nboul,2), aia(1:Nboul,2)
            real(kind=8)    :: ra(2), va(2)
            real(kind=8)    :: rcm(2), vcm(2), acm(2)
            real(kind=8)    :: m(0:Ntot), omega2, radius(0:Nboul)
            integer(kind=4) :: i, ineqs

            !y = /omega, mi, Ri, xi, yi, vxi, vyi, .../
            dydt = cero
            dydt(1) = domega_dt(t, y(1))
            dydt(2) = dm0_dt(t, y(2))

            omega2 = y(1) * y(1)
            m(0)      = y(2)
            radius(0) = y(3)
            do i = 1, Nboul
                ineqs = i * neqs
                m(i)      = y(ineqs+2)
                radius(i) = y(ineqs+3)
                ria(i,:)  = y(ineqs+4 : ineqs+5)
                via(i,:)  = y(ineqs+6 : ineqs+7)
                aia(i,:)  = -omega2 * ria(i,:) + dydt(1) * (/-ria(i,2), ria(i,1)/) ! a_i = -w^2 * r_i + dw/dt * (-ry,rx)
                dydt(ineqs+4 : ineqs+5) = via(i,:)
                dydt(ineqs+6 : ineqs+7) = aia(i,:)
            end do

            do i = Nboul+1, Ntot
                ineqs = i * neqs
                m(i) = y(ineqs+2)
                ra   = y(ineqs+4 : ineqs+5)
                va   = y(ineqs+6 : ineqs+7)
                dydt(ineqs+4 : ineqs+5) = va
                call accast(y(1), m(0:Nboul), ra, ria, dydt(ineqs+6 : ineqs+7))
                call rcmfromast(ria, via, aia, m(0:Nboul), rcm, vcm, acm)
                call accsto(t, ra-rcm, va-vcm, dydt(ineqs+6 : ineqs+7))
            end do
        end function dydt_ast_im

        function dydt_rot(t, y) result(dydt)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8), dimension(2)             :: rr, vr, rb, vb
            real(kind=8)                           :: m(0:Ntot), radius(0:Nboul)
            integer(kind=4)                        :: i, ineqs

            !y = /omega, mi, Ri, xi, yi, vxi, vyi, .../
            dydt = cero
            dydt(1) = domega_dt(t, y(1))
            dydt(2) = dm0_dt(t, y(2))

            do i = 0, Nboul
                ineqs = i * neqs
                m(i)      = y(ineqs+2)
                radius(i) = y(ineqs+3)
            end do

            do i = Nboul+1, Ntot
                ineqs = i * neqs
                m(i) = y(ineqs+2)
                rr   = y(ineqs+4 : ineqs+5)
                vr   = y(ineqs+6 : ineqs+7)
                dydt(ineqs+4 : ineqs+5) = vr
                call accrot(y(1), m(0:Nboul), rr, vr, ria, dydt(ineqs+6 : ineqs+7)) ! ria es invariante
                ! rcm es invariante ? (no, pero no importa)
                rb    = rr - rcm
                vb(1) = vr(1) - y(1) * rr(2) - vcm(1)
                vb(2) = vr(2) + y(1) * rr(1) - vcm(2)
                call accsto(t, rb, vb, dydt(ineqs+6 : ineqs+7))
            end do
        end function dydt_rot
end module derivates