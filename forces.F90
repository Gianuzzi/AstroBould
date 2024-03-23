module forces
    use const
    use parameters
    implicit none

    private :: set_db0, set_da0, accbar, accsto, accJ2, accnaisto !, acc_triax

    real(kind=8) :: J2coef = cero
    real(kind=8) :: ReC20coe = cero, ReC22coe = cero
    real(kind=8) :: da03 = cero, da05 = cero, da07 = cero
    real(kind=8) :: db03 = cero

    contains

        !! Get distances
        subroutine set_db0(rb, db0)
            implicit none
            real(kind=8), intent(in) :: rb(2)
            real(kind=8), intent(out) :: db0

            db0 = sqrt(rb(1)*rb(1) + rb(2)*rb(2))
            db03 = db0 * db0 * db0
        end subroutine set_db0

        subroutine set_da0(rb, rb0, raux, da0)
            implicit none
            real(kind=8), intent(in) :: rb(2), rb0(2)
            real(kind=8), intent(out) :: raux(2), da0

            raux = rb - rb0
            da0 = sqrt(raux(1)**2 + raux(2)**2)
            if (da0 < rmin) hexit = 1
            if (da0 > rmax) hexit = 2
            da03 = da0 * da0 * da0
            da05 = da03 * da0 * da0
            da07 = da05 * da0 * da0
        end subroutine set_da0

        subroutine apply_force(t, omega, m, rb, vb, rib, rdd)
            implicit none
            real(kind=8), intent(in) :: t, omega, m(0:), rb(2), vb(2), rib(0:,:)
            real(kind=8), intent(inout) :: rdd(2)
            
            call set_da0(rb, rib(0,:2), raux, da0)
            call accbar(m, rb, rib, rdd)
            if (lostokes .or. lostokes_naive) then
                call set_db0(rb, db0)
                if (lostokes) call accsto(GM, t, rb, vb, rdd)
                if (lostokes_naive) call accnaisto(omega, rb, vb, rdd)
            end if
            if (loJ2) call accJ2(GM0, rdd)
            ! if (loTriAx) call acc_triax(GM0, t, rdd)
        end subroutine apply_force

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! GRAVITY !!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! Potencial
            
        real(kind=8) function potbar(rb, rib) result(p)
            implicit none
            real(kind=8), intent(in) :: rb(2), rib(0:,:)
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

        real(kind=8) function potast(ra, ria) result(p)
            implicit none
            real(kind=8), intent(in) :: ra(2), ria(Nboul,2)
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

        real(kind=8) function potrot(rr, ria) result(p)
            implicit none
            real(kind=8), intent(in) :: rr(2), ria(Nboul,2)
            
            p = potast(rr, ria)
            p = p - omega2 * (rr(1)*rr(1) + rr(2)*rr(2)) * uno2 ! - w^2 * r^2 / 2
        end function potrot


        !!! Acceleration

        subroutine accbar(m, rb, rib, ab)
            implicit none
            real(kind=8), intent(in) :: m(0:), rb(2), rib(0:,:)
            real(kind=8), intent(out) :: ab(2)
            real(kind=8) :: d 
            integer(kind=4) :: i
            ab = - m(0) * raux / da03
            do i = 1, Nboul
                d = sqrt((rb(1)-rib(i,1))**2 + (rb(2)-rib(i,2))**2)
                ab = ab - m(i) * (rb - rib(i,:)) / (d*d*d)
            end do
            ab = ab * G
        end subroutine accbar

        subroutine accast(omega, m, ra, ria, aa)
            implicit none
            real(kind=8), intent(in) :: m(0:), ra(2), ria(Nboul,2), omega
            real(kind=8), intent(out) :: aa(2)
            real(kind=8) :: d, omega2, mucm(Nboul)
            integer(kind=4) :: i
            mucm = m(1:Nboul) / m(0)
            omega2 = omega * omega
            d = sqrt(ra(1)*ra(1) + ra(2)*ra(2))
            if (d < rmin) hexit = 1
            if (d > rmax) hexit = 2
            aa = - m(0) * ra / (d*d*d)
            ! aa = cero
            do i = 1, Nboul
                d = sqrt((ra(1)-ria(i,1))**2 + (ra(2)-ria(i,2))**2)
                aa = aa - m(i) * (ra - ria(i,:)) / (d*d*d) - omega2 * mucm(i) * ria(i,:) ! - w^2 * rcm
            end do
            aa = aa * G
        end subroutine accast

        subroutine accrot(omega, m, rr, vr, ria, ar)
            implicit none
            real(kind=8), intent(in) :: m(0:), rr(2), vr(2), ria(Nboul,2), omega
            real(kind=8), intent(out) :: ar(2)   
            call accast(omega, m, rr, ria, ar)
            ar = ar +  omega * omega * rr             ! Centrifugal (+ w^2 * r)
            ar = ar - dos * omega * (/-vr(2), vr(1)/) ! Coriolis    (- 2w * (-vy, vx))
        end subroutine accrot


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!! STOKES !!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Parameters

        subroutine set_C_and_Alpha(tau_a, tau_e, C, alpha)
            implicit none
            real(kind=8), intent(in) :: tau_a, tau_e
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

        subroutine accsto(GM, t, rb, vb, rdd)
            implicit none
            real(kind=8), intent(in) :: GM
            real(kind=8), intent(in) :: t, rb(2), vb(2)
            real(kind=8), intent(inout) :: rdd(2)
            real(kind=8) :: vc(2)

            call getfact_stok(t,f_stk)
            vc  = sqrt(GM / db03) * (/-rb(2),rb(1)/)
            rdd = rdd - C_stk * (vb - a_stk * vc) * f_stk ! -C(v * - alpha*vc)
        end subroutine accsto


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!! NAIVE STOKES !!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine accnaisto(omega, rb, vb, rdd)
            implicit none
            real(kind=8), intent(in) :: omega
            real(kind=8), intent(in) :: rb(2), vb(2)
            real(kind=8), intent(inout) :: rdd(2)
            real(kind=8) :: gamma, vrad

            vrad = dot_product(vb, rb) / db0 ! This is v_radial
            gamma = - eta * omega * vrad ! This is a_radial
            rdd = rdd + gamma * rb / db0
        end subroutine accnaisto

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! GEO-POTENTIAL !!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Parameters

        subroutine set_J2(J2)
            real(kind=8), intent(in) :: J2
            J2coef = 1.5d0 * J2
        end subroutine set_J2

        !!! Acceleration

        subroutine accJ2(GM0, rdd)
            implicit none
            real(kind=8), intent(in) :: GM0
            real(kind=8), intent(inout) :: rdd(2)
            
            rdd = rdd - GM0 * J2coef * raux / da05
        end subroutine accJ2

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!! TRI-AXIAL !!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! !!! Parameters

        ! subroutine set_tri_coefs(a, b, c, C20, C22, Re)
        !     implicit none
        !     real(kind=8), intent(in) :: a, b, c
        !     real(kind=8), intent(out) :: C20, C22, Re
        !     real(kind=8) :: Re2

        !     if ((a < b) .or. (b < c)) then
        !         write(*,*) 'ERROR: a < b < c'
        !         stop
        !     end if

        !     Re = (a * b * c)**(uno3)
        !     Re2 = Re * Re
        !     C20 = (dos * c*c - a*a - b*b) / (10.d0 * Re2)
        !     C22 = (a*a - b*b) / (20.d0 * Re2)
        
        !     ReC20coe = Re2 * C20 * 1.5d0
        !     ReC22coe = Re2 * C22 * 3.0d0
        ! end subroutine set_tri_coefs

        ! !!! Acceleration

        ! subroutine acc_triax(GM0, t, rdd)
        !     implicit none
        !     real(kind=8), intent(in) :: GM0, t
        !     real(kind=8), intent(inout) :: rdd(2)
        !     real(kind=8) :: ac_rot(2) = cero
        !     real(kind=8) :: fac1, fac2, xx, yy, co, se
            
        !     xx = raux(1) * raux(1)
        !     yy = raux(2) * raux(2)
        !     fac1 = -3.d0 * xx + 7.d0 * yy
        !     fac2 = -7.d0 * xx + 3.d0 * yy

        !     ac_rot(1) = raux(1) * ((ReC20coe / da05) + (ReC22coe * fac1 / da07))
        !     ac_rot(2) = raux(2) * ((ReC20coe / da05) + (ReC22coe * fac2 / da07))
        !     ac_rot = ac_rot * GM0
            
        !     !! Ahora hay que rotarlo
        !     co = cos(omega * t)
        !     se = sin(omega * t)
        !     rdd = rdd + (/ac_rot(1) * co - ac_rot(2) * se, &
        !                 & ac_rot(1) * se + ac_rot(2) * co/)
        ! end subroutine acc_triax

end module forces
