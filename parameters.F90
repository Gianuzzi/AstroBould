module parameters
    use const
    implicit none
    !!! Parámetros y auxiliares
    integer(kind=4), parameter :: Nboul = 1, Npart = 1, Ntot = Npart + Nboul !!! Boulders, Particles
    integer(kind=4), parameter :: neqs = 6, NP = (Nboul+1) * neqs + 1 !! Np es posicion de primer partícula
    integer(kind=4) :: Nactive
    real(kind=8) :: m(0:Ntot), mu(Ntot), GM, mcm, Gmi(0:Ntot), mucm(0:Ntot)
    real(kind=8) :: radius(0:Ntot)
    real(kind=8) :: rib(0:Ntot,2), ria(1:Ntot,2), rcm(2)
    real(kind=8) :: vib(0:Ntot,2), via(1:Ntot,2), vcm(2)
    real(kind=8) :: aib(0:Ntot,2), aia(1:Ntot,2), acm(2)
    real(kind=8) :: rb(2), ra(2), rr(2)
    real(kind=8) :: vb(2), va(2), vr(2)
    real(kind=8) :: ab(2), aa(2), ar(2)
    real(kind=8) :: da0
    real(kind=8) :: da, de, amax, amin, emax, emin
    real(kind=8) :: xc(6), xe(6), xc0(6), xe0(6), Res0
    real(kind=8) :: ea, ee, ei, eM, ew, eO, eR
    real(kind=8) :: theta_a(Nboul), theta_acm
    real(kind=8) :: r_b(0:Nboul), theta_b(0:Nboul)
    real(kind=8) :: omega, wk, lambda, R0, Prot, a_corot
    real(kind=8) :: omega2, lambda2, fac_omega=uno, rmax=inf


    contains

        subroutine rcmfromast(ria,via,aia,m,rcm,vcm,acm)
            implicit none
            real(kind=8), intent(in)  :: ria(1:Nboul,2), via(1:Nboul,2), aia(1:Nboul,2)
            real(kind=8), intent(in)  :: m(0:Nboul)
            real(kind=8), intent(out) :: rcm(2), vcm(2), acm(2)
            real(kind=8) :: mucm(1:Nboul)
            integer(kind=4) :: i

            rcm = cero ; vcm = cero ; acm = cero
            mucm = m(1:Nboul) / sum(m)
            do i = 1, Nboul!size(ria,1)
                rcm = rcm + mucm(i) * ria(i,:) ! Posición del centro de masas
                vcm = vcm + mucm(i) * via(i,:) ! Velocidad del centro de masas
                acm = acm + mucm(i) * aia(i,:) ! Aceleración del centro de masas
            end do
        end subroutine rcmfromast

        subroutine ast2bar(ra,va,aa,rcm,vcm,acm,rb,vb,ab)
            implicit none
            real(kind=8), intent(in)  :: ra(2), va(2), aa(2)
            real(kind=8), intent(in)  :: rcm(2), vcm(2), acm(2)
            real(kind=8), intent(out) :: rb(2), vb(2), ab(2)
            rb = ra - rcm ! Posición
            vb = va - vcm ! Velocidad
            ab = aa - acm ! Aceleración
        end subroutine ast2bar

        subroutine bar2ast(rb,vb,ab,rcm,vcm,acm,ra,va,aa)
            implicit none
            real(kind=8), intent(in)  :: rb(2), vb(2), ab(2)
            real(kind=8), intent(in)  :: rcm(2), vcm(2), acm(2)
            real(kind=8), intent(out) :: ra(2), va(2), aa(2)
            ra = rb + rcm ! Posición
            va = vb + vcm ! Velocidad
            aa = ab + acm ! Aceleración
        end subroutine bar2ast

        subroutine rot2ast(rr,vr,ar,omega,ra,va,aa)
            implicit none
            real(kind=8), intent(in)  :: rr(2), vr(2), ar(2), omega
            real(kind=8), intent(out) :: ra(2), va(2), aa(2)
            real(kind=8) :: omega2
            omega2 = omega * omega
            ra(1) = rr(1)
            ra(2) = rr(2)
            va(1) = vr(1) - omega * rr(2)
            va(2) = vr(2) + omega * rr(1)
            aa(1) = ar(1) - omega * dos * vr(2) - omega2 * rr(1)
            aa(2) = ar(2) + omega * dos * vr(1) - omega2 * rr(2)
        end subroutine rot2ast

        subroutine ast2rot(ra,va,aa,omega,rr,vr,ar)
            implicit none
            real(kind=8), intent(in)  :: ra(2), va(2), aa(2), omega
            real(kind=8), intent(out) :: rr(2), vr(2), ar(2)
            real(kind=8) :: omega2
            omega2 = omega * omega
            rr(1) = ra(1)
            rr(2) = ra(2)
            vr(1) = va(1) + omega * rr(2)
            vr(2) = va(2) - omega * rr(1)
            ar(1) = aa(1) + dos * omega * vr(2) + omega2 * rr(1)
            ar(2) = aa(2) - dos * omega * vr(1) + omega2 * rr(2)
        end subroutine ast2rot

        subroutine ang_mom_bar(rib, vib, mi, radius, omega, ang_mom)
            implicit none
            real(kind=8), intent(in)  :: rib(0:Ntot,2), vib(0:Ntot,2), mi(0:Ntot), radius(0:Ntot)
            real(kind=8), intent(in)  :: omega
            real(kind=8), intent(out) :: ang_mom
            integer(kind=4) :: i
            ang_mom = cero
            do i = 0, Nboul
                ang_mom = ang_mom + mi(i) * (rib(i,1) * vib(i,2) - rib(i,2) * vib(i,1)) ! Traslacional
                ang_mom = ang_mom + 0.4d0 * mi(i) * radius(i) * radius(i) * omega ! Rotacional
            end do
            do i = Nboul+1, Ntot
                ang_mom = ang_mom + mi(i) * (rib(i,1) * vib(i,2) - rib(i,2) * vib(i,1)) ! Traslacional
            end do
        end subroutine ang_mom_bar

        subroutine ang_mom_ast(m, radius, omega, ang_mom)
            implicit none
            real(kind=8), intent(in)  :: m(0:Ntot), radius(0:Ntot)  
            real(kind=8), intent(in)  :: omega
            real(kind=8), intent(out) :: ang_mom
            real(kind=8) :: omegar2
            integer(kind=4) :: i
            omegar2 = omega * radius(0) * radius(0)
            ang_mom = cero
            ! ang_mom = ang_mom + sum(m(0:Nboul)) * dot_product(rcm, rcm) * omega ! Rotacional
            ang_mom = ang_mom + 0.4d0 * m(0) * omegar2 ! Rotacional
            do i = 1, Nboul
                ang_mom = ang_mom + m(i) * omegar2 ! Traslacional
                ang_mom = ang_mom + 0.4d0 * m(i) * radius(i) * radius(i) * omega ! Rotacional
            end do
            do i = Nboul+1, Ntot
                ang_mom = ang_mom + m(i) * (ria(i,1) * via(i,2) - ria(i,2) * via(i,1)) ! Traslacional
            end do
        end subroutine ang_mom_ast
end module parameters