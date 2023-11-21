module forces
    use const
    use parameters
    implicit none

    contains

        !! Get distances
        subroutine set_db0(rb, db0)
            real(kind=8), intent(in) :: rb(2)
            real(kind=8), intent(out) :: db0

            db0 = sqrt(rb(1)*rb(1) + rb(2)*rb(2))
        end subroutine set_db0

        subroutine set_da0(rb, rb0, raux, da0)
            real(kind=8), intent(in) :: rb(2), rb0(2)
            real(kind=8), intent(out) :: raux(2), da0

            raux = rb - rb0
            da0 = sqrt(raux(1)**2 + raux(2)**2)
        end subroutine set_da0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! GRAVITY !!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! Potencial
            
        real(kind=8) function potbar(rb,rib) result(p)
            implicit none
            real(kind=8), intent(in) :: rb(2), rib(0:Nboul,2)
            real(kind=8) :: d
            integer(kind=4) :: i

            d = sqrt((rb(1)-rib(0,1))**2 + (rb(2)-rib(0,2))**2)
            if (d < R0) then
                p = inf
                return
            end if
            p = - Gmi(0) / d
            ! p = cero
            do i = 1, Nboul
                d = sqrt((rb(1)-rib(i,1))**2 + (rb(2)-rib(i,2))**2)
                p = p - Gmi(i)/d
            end do
        end function potbar

        real(kind=8) function potast(ra,ria) result(p)
            implicit none
            real(kind=8), intent(in) :: ra(2), ria(1:Nboul,2)
            real(kind=8) :: d
            integer(kind=4) :: i
            
            d = sqrt(ra(1)*ra(1) + ra(2)*ra(2))
            if (d < R0) then
                p = inf
                return
            end if
            p = - Gmi(0) / d
            ! p = cero
            do i = 1, Nboul
                d = sqrt((ra(1)-ria(i,1))**2 + (ra(2)-ria(i,2))**2)
                p = p - Gmi(i) / d + omega2 * mucm(i) * (ria(i,1) * ra(1) + ria(i,2) * ra(2)) ! + w^2 * dot(rcm,r)
            end do
        end function potast

        real(kind=8) function potrot(rr,ria) result(p)
            implicit none
            real(kind=8), intent(in) :: rr(2), ria(1:Nboul,2)
            
            p = potast(rr, ria)
            p = p - omega2 * (rr(1)*rr(1) + rr(2)*rr(2)) * uno2 ! - w^2 * r^2 / 2
        end function potrot


        !!! Acceleration

        subroutine accbar(m,rb,rib,ab)
            implicit none
            real(kind=8), intent(in)  :: m(0:Nboul), rb(2), rib(0:Nboul,2)
            real(kind=8), intent(out) :: ab(2)
            real(kind=8) :: d 
            integer(kind=4) :: i
            call set_da0(rb, rib(0,:), raux, da0)
            if (da0 < rmin) hexit = .True.
            ab = - m(0) * raux / (da0*da0*da0)
            do i = 1, Nboul
                d = sqrt((rb(1)-rib(i,1))**2 + (rb(2)-rib(i,2))**2)
                ab = ab - m(i) * (rb - rib(i,:)) / (d*d*d)
            end do
            ab = ab * G
        end subroutine accbar

        subroutine accast(omega,m,ra,ria,aa)
            implicit none
            real(kind=8), intent(in)  :: m(0:Nboul), ra(2), ria(1:Nboul,2), omega
            real(kind=8), intent(out) :: aa(2)
            real(kind=8) :: d, omega2, mucm(1:Nboul)
            integer(kind=4) :: i
            mucm = m(1:Nboul) / m(0)
            omega2 = omega * omega
            d = sqrt(ra(1)*ra(1) + ra(2)*ra(2))
            if (d < rmin) hexit = .True.
            aa = - m(0) * ra / (d*d*d)
            ! aa = cero
            do i = 1, Nboul
                d = sqrt((ra(1)-ria(i,1))**2 + (ra(2)-ria(i,2))**2)
                aa = aa - m(i) * (ra - ria(i,:)) / (d*d*d) - omega2 * mucm(i) * ria(i,:) ! - w^2 * rcm
            end do
            aa = aa * G
        end subroutine accast

        subroutine accrot(omega,m,rr,vr,ria,ar)
            implicit none
            real(kind=8), intent(in)  :: m(0:Nboul), rr(2), vr(2), ria(1:Nboul,2), omega
            real(kind=8), intent(out) :: ar(2)   
            call accast(omega, m, rr, ria, ar)
            ar = ar +  omega * omega * rr             ! Centrifugal (+ w^2 * r)
            ar = ar - dos * omega * (/-vr(2), vr(1)/) ! Coriolis    (- 2w * (-vy, vx))
        end subroutine accrot


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!! STOKES !!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Parameters

        subroutine set_C_and_Alpha(tau_a,tau_e,C,alpha)
            implicit none
            real(kind=8), intent(in)  :: tau_a, tau_e
            real(kind=8), intent(out) :: C, alpha
            C = uno / (dos * tau_a) + uno / tau_e
            alpha = (dos * tau_a) / ((dos * tau_a) + tau_e)
        end subroutine set_C_and_Alpha

        subroutine getfact_stok(t,f_stk)
            implicit none
            real(kind=8), intent(in)    :: t
            real(kind=8), intent(out) :: f_stk
            f_stk = uno2 * (uno + tanh(1.d1 * (uno - t / t_stokes)))
        end subroutine getfact_stok

        !!! Acceleration

        subroutine accsto(t,rb,vb,rdd, GM)
            implicit none
            real(kind=8), intent(in) :: t, rb(2), vb(2)
            real(kind=8), intent(in) :: GM
            real(kind=8), intent(inout) :: rdd(2)
            real(kind=8) :: vc(2)

            if (.not. lostokes) return
            call getfact_stok(t,f_stk)
            call set_db0(rb, db0) ! Así ya queda db0 calculado
            vc  = sqrt(GM / (db0*db0*db0)) * (/-rb(2),rb(1)/)
            rdd = rdd - C_stk * (vb - a_stk * vc) * f_stk ! -C(v * - alpha*vc)
        end subroutine accsto


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! GEO-POTENTIAL !!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Acceleration

        subroutine accJ2(J2, rdd)
            implicit none
            real(kind=8), intent(in) :: J2
            real(kind=8), intent(inout) :: rdd(2)

            ! call set_da0(rb, rb0, raux, da0) !! Se supone que ya está calculado
            if (.not. loJ2) return
            rdd = rdd - GM0 * J2 * raux * 1.5d0 / (da0*da0*da0*da0*da0)
        end subroutine accJ2

end module forces
