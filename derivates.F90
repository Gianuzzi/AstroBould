module derivates
    use const
    use parameters
    use forces
    implicit none
    !y = /omega, mi, xi, yi, vxi, vyi, .../

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
                domega = -exp(-(t-t0) / tau_o) * omega / tau_o
            end if
        end function domega_dt

        function dm0_dt(t, m0) result(dm0)
            implicit none
            real(kind=8), intent(in) :: t, m0
            real(kind=8) :: dm0
            if (tau_m >= inf) then
                dm0 = cero
            else
                dm0 = -exp(-(t-t0) / tau_m) * m0 / tau_m
            end if
        end function dm0_dt

        function dydt_bar_ex(t, y) result(dydt)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: rib(0:Nboul,2), rb(2), vb(2), ab(2)
            integer(kind=4) :: i, ineqs

            !y = /omega, mi, xi, yi, vxi, vyi, .../
            dydt = cero
            ! print*, 't = ', t

            do i = 0, Nboul
                rib(i,1) = cos(omega * (t-t0) + theta_b(i)) * r_b(i)
                rib(i,2) = sin(omega * (t-t0) + theta_b(i)) * r_b(i)
            end do

            do i = Nboul+1, Ntot
                ineqs = i*neqs
                rb = y(ineqs+3 : ineqs+4)
                vb = y(ineqs+5 : ineqs+6)
                call apply_force(t-t0, omega, m, rb, vb, rib, ab)
                dydt(ineqs+3 : ineqs+4) = vb
                dydt(ineqs+5 : ineqs+6) = ab
            end do
        end function dydt_bar_ex

        function dydt_bar_im(t, y) result(dydt)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: rib(0:Nboul,2), vib(0:Nboul,2), rb(2), vb(2), ab(2)
            real(kind=8) :: m(0:Nboul), omega2!, radius(0:)
            integer(kind=4) :: i, ineqs

            !y = /omega, mi, xi, yi, vxi, vyi, .../
            dydt = cero
            ! dydt(1) = domega_dt(t, y(1)) !not used
            ! dydt(2) = dm0_dt(t, y(2)) !not used

            omega2 = y(1) * y(1)
            do i = 0, Nboul ! Integradores usan desde i=1
                ineqs = i * neqs
                m(i)      = y(ineqs+2)
                rib(i,:)  = y(ineqs+3 : ineqs+4)
                vib(i,:)  = y(ineqs+5 : ineqs+6)
                dydt(ineqs+3 : ineqs+4) = y(ineqs+5 : ineqs+6)
                dydt(ineqs+5 : ineqs+6) = -omega2 * rib(i,:) !+ dydt(1) * (/-rib(i,2), rib(i,1)/) ! a_i = -w^2 * r_i + dw/dt * (-ry,rx)
            end do

            do i = Nboul+1, Ntot
                ineqs = i * neqs
                rb = y(ineqs+3 : ineqs+4)
                vb = y(ineqs+5 : ineqs+6)
                call apply_force(t-t0, y(1), m, rb, vb, rib, ab)
                dydt(ineqs+3 : ineqs+4) = vb
                dydt(ineqs+5 : ineqs+6) = ab
            end do
        end function dydt_bar_im
end module derivates