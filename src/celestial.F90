!> Module with coordinates/elements, and some extra cel-mech routines.
module celestial
    use constants, only: wp, cero, uno, uno2, uno3, dos, G, pi, twopi, tini, epsilon, sqepsilon
    
    implicit none
    
    contains

        ! Get period
        pure function get_Period(mass, a) result(per)
            implicit none
            real(wp), intent(in) :: mass, a
            real(wp) :: per

            per = cero
            if (mass > cero) per = twopi * sqrt(a * a * a / (G * mass))
        end function get_Period

        ! Get corotation a from rotating body
        pure function get_a_corot(mass, omega) result(acorot)
            implicit none
            real(wp), intent(in) :: mass, omega
            real(wp) :: acorot

            if (omega < tini) then
                acorot = cero
            else 
                acorot = (G * mass / (omega * omega))**(1/3.)
            end if
        end function get_a_corot

        ! Get center of mass from masses, positions and velocities
        pure subroutine get_center_of_mass(mass, rib, vib, mcm, rcm, vcm)
            implicit none
            real(wp), intent(in) :: mass(:), rib(:,:), vib(:,:)
            real(wp), intent(out) :: mcm, rcm(2), vcm(2)
            integer(kind=4) :: i

            rcm = cero
            vcm = cero
            mcm = cero
            do i = 1, size(mass)
                rcm = rcm + mass(i) * rib(i,:)
                vcm = vcm + mass(i) * vib(i,:)
                mcm = mcm + mass(i)
            end do
            rcm = rcm / mcm ! rcm = sum_i m_i * r_i / M
            vcm = vcm / mcm ! vcm = sum_i m_i * v_i / M
        end subroutine get_center_of_mass

        ! Get acceleration and potential energy from single mass
        pure subroutine get_acc_and_pot_single(mass, rib, xy_target, dr_max, acc, pot, inside)
            implicit none
            real(wp), intent(in) :: mass, rib(2), xy_target(2), dr_max
            real(wp), intent(inout) :: acc(2), pot
            logical, intent(inout), optional :: inside
            real(wp) :: dx, dy, dr, dr2

            dx = xy_target(1) - rib(1)
            dy = xy_target(2) - rib(2)
            dr2 = dx * dx + dy * dy
            dr = sqrt(dr2)
            
            if (dr > tini) then
                acc = acc - G * (mass / (dr2 * dr)) * (/dx, dy/)
                pot = pot - G * mass / dr
            end if

            if (present(inside)) inside = dr < max(dr_max, tini)

        end subroutine get_acc_and_pot_single

        ! Get coordinates from a body
        subroutine coord(msum, a, e, inc, capm, omega, capom, xc)
            implicit none
            real(wp), intent(in)  :: msum, a, e, inc, capm, omega, capom
            real(wp), intent(out) :: xc(6)
            real(wp) :: sp, cp, so, co, si, ci
            real(wp) :: d11, d12, d13, d21, d22, d23
            real(wp) :: cape, dummy, scap, ccap, sqe, sqgma
            real(wp) :: ri, xfac1, xfac2, vfac1, vfac2
            
            ! Generate rotation matrices (on p. 42 of Fitzpatrick)
            sp = sin(omega)
            cp = cos(omega)
            so = sin(capom)
            co = cos(capom)
            si = sin(inc)
            ci = cos(inc)
            d11 = cp * co - sp * so * ci
            d12 = cp * so + sp * co * ci
            d13 = sp * si
            d21 = -sp * co - cp * so * ci
            d22 = -sp * so + cp * co * ci
            d23 = cp * si
            
            ! Get the other quantities depending on orbit type (i.e. ialpha)
            call aver(capm, e, cape, dummy)
            scap = sin(cape)
            ccap = cos(cape)
            sqe = sqrt(uno - e * e)
            sqgma = sqrt(G * msum * a)
            xfac1 = a * (ccap - e)
            xfac2 = a * sqe * scap
            ri = uno / (a * (uno - e * ccap))
            vfac1 = -ri * sqgma * scap
            vfac2 = ri * sqgma * sqe * ccap
            
            xc(1) = d11 * xfac1 + d21 * xfac2
            xc(2) = d12 * xfac1 + d22 * xfac2
            xc(3) = d13 * xfac1 + d23 * xfac2
            xc(4) = d11 * vfac1 + d21 * vfac2
            xc(5) = d12 * vfac1 + d22 * vfac2
            xc(6) = d13 * vfac1 + d23 * vfac2
        end subroutine coord

        ! Get true f and g from a body
        subroutine aver(dm, e, u, f)
            implicit none
            real(wp), intent(in)  :: dm, e
            real(wp), intent(out) :: u, f
            real(wp) :: u0, dif, seno, cose
            integer(kind=4), parameter :: MAX_ITER = 100
            integer(kind=4) :: i
            
            u0 = dm
            dif = uno
            i = 0
            do while (dif > epsilon)
                u = dm + e * sin(u0)
                dif = abs(u - u0)
                u0 = u
                i = i + 1
                if (i > MAX_ITER) then
                    print*, 'aver: too many iterations. Error:', dif
                    exit
                end if
            end do
            seno = sqrt(uno + e) * sin(uno2 * u)
            cose = sqrt(uno - e) * cos(uno2 * u)
            f = dos * atan2(seno, cose)
        end subroutine aver

        ! Get elements from a body
        pure subroutine elem(msum, xc, a, e, inc, capm, omega, capom)
            implicit none
            real(wp), intent(in)  :: msum, xc(6)
            real(wp), intent(out) :: a, e, inc, capm, omega, capom
            real(wp) :: gmsum
            real(wp) :: x, y, z, vx, vy, vz
            real(wp) :: hx, hy, hz, h2, h, fac, u
            real(wp) :: r, v, v2, vdotr, energy
            real(wp) :: cape, cw, sw, w, face, capf, tmpf
            integer(kind=4) :: ialpha
                
            gmsum = G * msum
            x = xc(1)
            y = xc(2)
            z = xc(3)
            vx = xc(4)
            vy = xc(5)
            vz = xc(6)
            
            hx = y * vz - z * vy
            hy = z * vx - x * vz
            hz = x * vy - y * vx
            h2 = hx * hx + hy * hy + hz * hz
            h = sqrt(h2)
            inc = acos(hz / h)
            
            fac = sqrt(hx*hx + hy*hy) / h
            if (fac < epsilon) then
                capom = cero
                u = atan2(y, x)
                if (abs(inc - pi) < 10.e0_wp * epsilon) then
                    u = -u
                end if
            else
                capom = atan2(hx, -hy)
                u = atan2(z / sin(inc), x * cos(capom) + y * sin(capom))
            end if
            
            if (capom < cero) then
                capom = capom + twopi
            end if
            if (u < cero) then
                u = u + twopi
            end if
            
            r = sqrt(x * x + y * y + z * z)
            v2 = vx * vx + vy * vy + vz * vz
            v = sqrt(v2)
            vdotr = x * vx + y * vy + z * vz
            energy = uno2 * v2 - gmsum / r
            
            if (abs(energy * r / gmsum) < sqepsilon) then
                ialpha = 0
            else
                if (energy < cero) then
                    ialpha = -1
                else if (energy > cero) then
                    ialpha = 1
                end if
            end if
            
            !! Ellipse
            if (ialpha == -1) then
                a = -uno2 * gmsum / energy
                fac = uno - h2 / (gmsum * a)
                if (fac > epsilon) then
                    e = sqrt(fac)
                    face = (a - r) / (a * e)
                    if (face > uno) then
                        cape = cero
                    else
                        if (face > -uno) then
                            cape = acos(face)
                        else
                            cape = pi
                        end if
                    end if
                    if (vdotr < cero) then
                        cape = twopi - cape
                    end if
                    cw = (cos(cape) - e) / (uno - e * cos(cape))
                    sw = sqrt(uno - e * e) * sin(cape) / (uno - e * cos(cape))
                    w = atan2(sw, cw)
                    if (w < cero) then
                        w = w + twopi
                    end if
                else
                    e = cero
                    w = u
                    cape = u
                end if
                capm = cape - e * sin(cape)
                omega = u - w
                if (omega < cero) then
                    omega = omega + twopi
                end if
                omega = omega - int(omega / twopi) * twopi
            end if
            
            !! Hypérbola
            if (ialpha == 1) then
                a = uno2 * gmsum / energy
                fac = h2 / (gmsum * a)
                if (fac > epsilon) then
                    e = sqrt(uno + fac)
                    tmpf = (a + r) / (a * e)
                    if (tmpf < uno) then
                        tmpf = uno
                    end if
                    capf = log(tmpf + sqrt(tmpf * tmpf - uno))
                    if (vdotr < cero) then
                        capf = -capf
                    end if
                    cw = (e - cosh(capf)) / (e * cosh(capf) - uno)
                    sw = sqrt(e * e - uno) * sinh(capf) / (e * cosh(capf) - uno)
                    w = atan2(sw, cw)
                    if (w < cero) then
                        w = w + twopi
                    end if
                else
                    e = uno
                    tmpf = uno2 * h2 / gmsum
                    w = acos(dos * tmpf / r - uno)
                    if (vdotr < cero) then
                        w = twopi - w
                    end if
                    tmpf = (a + r) / (a * e)
                    capf = log(tmpf + sqrt(tmpf * tmpf - uno))
                end if
                capm = e * sinh(capf) - capf
                omega = u - w
                if (omega < cero) then
                    omega = omega + twopi
                end if
                omega = omega - int(omega / twopi) * twopi
            end if

            !! Parábola    
            if (ialpha == 0) then
                a = uno2 * h2 / gmsum
                e = uno
                w = acos(dos * a / r - uno)
                if (vdotr < cero) then
                    w = twopi - w
                end if
                tmpf = tan(uno2 * w)
                capm = tmpf * (uno + tmpf * tmpf * uno3)
                omega = u - w
                if (omega < cero) then 
                    omega = omega + twopi
                end if
                omega = omega - int(omega / twopi) * twopi
            end if
        end subroutine elem

        ! Get astrocentric from jacobi coordinates
        pure subroutine jacobi_to_astroc (mass_array, coords_j, coords_a)
            implicit none
            real(wp), dimension(:), intent(in) :: mass_array
            real(wp), dimension(:,:), intent(in) :: coords_j
            real(wp), dimension(:,:), intent(out) :: coords_a
            real(wp) :: mass_acum, mass_acum_prev
            integer(kind=4) :: i

            coords_a = coords_j
            mass_acum_prev = mass_array(1)  ! Arrancamos desde 1
            do i = 2, size(coords_j,1)
                mass_acum = mass_acum_prev + mass_array(i)
                coords_a(i,1:4) = coords_j(i,1:4) + &
                                & (mass_acum_prev * (coords_a(i-1,1:4) - coords_j(i-1,1:4)) + &
                                & mass_array(i) * coords_a(i-1,1:4)) / mass_acum
                mass_acum_prev = mass_acum
            end do
        end subroutine jacobi_to_astroc

        ! Get jacobi from astrocentric coordinates
        pure subroutine astroc_to_jacob (mass_array, coords_a, coords_j)
            implicit none
            real(wp), dimension(:), intent(in) :: mass_array
            real(wp), dimension(:,:), intent(in) :: coords_a
            real(wp), dimension(:,:), intent(out) :: coords_j
            real(wp) :: mass_acum, mass_acum_prev
            integer(kind=4) :: i

            coords_j = coords_a
            mass_acum_prev = mass_array(1)  ! Arrancamos desde 1
            do i = 2, size(coords_a,1)
                mass_acum = mass_acum_prev + mass_array(i)
                coords_j(i,1:4) = coords_a(i,1:4) - &
                                & (mass_acum_prev * (coords_a(i-1,1:4) - coords_j(i-1,1:4)) + &
                                & mass_array(i) * coords_a(i-1,1:4)) / mass_acum
                mass_acum_prev = mass_acum
            end do
        end subroutine astroc_to_jacob

        ! Get barycentric from astrocentric coordinates
        pure subroutine baric_to_astroc (mass_array, coords_b, coords_a)
            implicit none
            real(wp), dimension(:), intent(in) :: mass_array
            real(wp), dimension(:,:), intent(in) :: coords_b
            real(wp), dimension(:,:), intent(out) :: coords_a
            integer(kind=4) :: i
            real(wp), dimension(4) :: xyvxvy

            xyvxvy = cero
            
            do i = 2, size(coords_b,1)
                xyvxvy = xyvxvy - coords_b(i,1:4) * mass_array(i) / mass_array(1)
            end do
            do i = 2, size(coords_b,1)
                coords_a(i,1:4) = coords_b(i,1:4) - xyvxvy
            end do
            coords_a(1, 1:4) = cero
        end subroutine baric_to_astroc

        ! Get astrocentric from barycentric coordinates
        pure subroutine astroc_to_baric (mass_array, coords_a, coords_b)
            implicit none
            real(wp), dimension(:), intent(in) :: mass_array
            real(wp), dimension(:,:), intent(in) :: coords_a
            real(wp), dimension(:,:), intent(out) :: coords_b
            integer(kind=4) :: i
            real(wp) :: mass_tot
            real(wp), dimension(4) :: xyvxvy

            mass_tot = mass_array(1)
            xyvxvy = cero

            do i = 2, size(coords_a,1)
                mass_tot = mass_tot + mass_array(i)
                xyvxvy = xyvxvy + coords_a(i,1:4) * mass_array(i)
            end do
            xyvxvy = - xyvxvy / mass_tot
            do i = 1, size(coords_a,1)
                coords_b(i,1:4) = coords_a(i,1:4) + xyvxvy
            end do
        end subroutine astroc_to_baric

        ! Osculating to geometrical
        pure subroutine coord2geom(mass, radius, J2, xc, a, e, inc, capm, omega, capom, n)
            implicit none
            real(wp), intent(in) :: mass, radius, J2
            real(wp), dimension(6), intent(in) :: xc
            real(wp), intent(out) :: a, e, inc, capm, omega, capom
            real(wp), intent(out) :: n
            ! Coordinates
            real(wp):: x, y, z
            real(wp):: vx, vy, vz
            !! Polar
            real(wp) :: r, L
            real(wp) :: vr, vL
            logical :: coplanar
            ! Corrections
            real(wp) :: rc, Lc
            real(wp) :: vrc, vLc
            real(wp) :: zc, vzc
            ! Extra Geometric elements
            real(wp) :: lambda, varpi
            ! For extra a corrections
            logical, parameter :: extra_cor = .True.
            real(wp) :: aux, det, aux_a, aux_b, aux_c
            real(wp) :: Hz
            real(wp) :: r0
            ! Iterations
            integer(kind=4) :: i
            integer(kind=4), parameter :: MAX_ITER = 100
            real(wp), parameter :: err_rel = 1.e-9_wp
            real(wp) :: a_prev

            ! Check
            if (abs(J2) < tini) then
                call elem(mass, xc, a, e, inc, capm, omega, capom)
                return
            end if

            ! Get coordinates
            x = xc(1)
            y = xc(2)
            z = xc(3)
            vx = xc(4)
            vy = xc(5)
            vz = xc(6)

            coplanar = (abs(z) < tini) .and. (abs(vz) < tini)

            r = sqrt(x * x + y * y)
            L = atan2(y, x)
            vr = vx * cos(L) + vy * sin(L)
            vL = (-vx * sin(L) + vy * cos(L)) / r

            Hz = x * vy - y * vx

            ! Parameters or extra 'a' corrections
            aux_a = uno
            aux_b = - Hz * Hz / (G * mass)  ! Negative
            aux_c = 1.5e0_wp * radius * radius * J2
            det = aux_b**2 - 4.e0_wp * aux_a * aux_c
            if (det < cero) then
                r0 = r
            else 
                aux = sqrt(det)
                if (aux > aux_b) then
                    r0 = (-aux_b + aux) / (dos * aux_a)
                else
                    r0 = (-aux_b - aux) / (dos * aux_a)
                end if
            end if

            ! Initial condition
            a = r
            e = cero
            inc = cero

            rc = cero
            Lc = cero
            zc = cero
            vrc = cero
            vLc = cero
            vzc = cero

            ! Loop
            do i = 1, MAX_ITER
                a_prev = a
                call get_geom(a, e, inc, lambda, varpi, capom, rc, Lc, zc, vrc, vLc, vzc, n)
                if ((abs(a_prev - a) / a_prev) < err_rel) exit
            end do

            ! Get output elements
            capm = lambda - varpi

            if (inc < (uno2 * pi)) then
                omega = varpi - capom
            else
                omega = capom - varpi
            end if            

            contains

                pure subroutine get_geom(a, e, inc, lambda, varpi, capom, rc, Lc, zc, vrc, vLc, vzc, n)
                    implicit none
                    real(wp), intent(inout) :: a, e, inc, lambda, varpi, capom
                    real(wp), intent(inout) :: rc, Lc, zc, vrc, vLc, vzc
                    real(wp), intent(out) :: n  ! The only important

                    real(wp) :: kappa, nu, eta, chi
                    real(wp) :: n2, kappa2, nu2, eta2, chi2
                    real(wp) :: alpha_1, alpha_2, alpha2

                    real(wp), parameter :: cn1 = 3.e0_wp/4.e0_wp
                    real(wp), parameter :: cn2 = -9.e0_wp/32.e0_wp
                    real(wp), parameter :: cn3 = 27.e0_wp/128.e0_wp
                    real(wp), parameter :: cn4 = 3.e0_wp
                    real(wp), parameter :: cn5 = -12.e0_wp

                    real(wp), parameter :: ck1 = -cn1
                    real(wp), parameter :: ck2 = cn2
                    real(wp), parameter :: ck3 = -cn3
                    real(wp), parameter :: ck4 = -9.e0_wp

                    real(wp), parameter :: cnu1 = 9.e0_wp/4.e0_wp
                    real(wp), parameter :: cnu2 = -81.e0_wp/32.e0_wp
                    real(wp), parameter :: cnu3 = 729.e0_wp/128.e0_wp
                    real(wp), parameter :: cnu4 = 6.e0_wp
                    real(wp), parameter :: cnu5 = -51.e0_wp/4.e0_wp

                    real(wp) :: nKep, R2_a2

                    real(wp) :: aux_x, aux_y

                    ! Compute frequencies
                    
                    nKep = sqrt(G * mass / (a * a * a))
                    R2_a2 = radius * radius / (a * a)

                    n = nKep * (uno + J2 * R2_a2 * (cn1 &
                                            & + cn2 * R2_a2 * J2 &
                                            & + cn3 * R2_a2**2 * J2**2 &
                                            & + cn4 * e**2 &
                                            & + cn5 * inc**2))
                    n2 = n * n

                    kappa = nKep * (uno + J2 * R2_a2 * (ck1 &
                                                & + ck2 * R2_a2 * J2 &
                                                & + ck3 * R2_a2**2 * J2**2 &
                                                & + ck4 * inc**2))
                    kappa2 = kappa * kappa

                    nu = nKep * (uno + J2 * R2_a2 * (cnu1 &
                                            & + cnu2 * R2_a2 * J2 &
                                            & + cnu3 * R2_a2**2 * J2**2 &
                                            & + cnu4 * e**2 &
                                            & + cnu5 * inc**2))
                    nu2 = nu * nu
                    
                    eta2 = nKep * nKep * (uno - dos * J2 * R2_a2)
                    eta = sqrt(eta2)

                    chi2 = nKep * nKep * (uno + 7.5e0_wp * J2 * R2_a2)
                    chi = sqrt(chi2)

                    alpha_1 = (dos * nu + kappa) * uno3
                    alpha_2 = (dos * nu - kappa)
                    alpha2 = alpha_1 * alpha_2

                    ! Compute geometric elements

                    a = (r - rc) / (uno - (vL - vLc - n) / (dos * n))

                    e = sqrt( ((vL - vLc - n) / (dos * n))**2 + ((vr - vrc) / (a * kappa))**2 )

                    inc = sqrt( ((z - zc) / a)**2 + ((vz - vzc) / (a * nu))**2 )

                    lambda = L - Lc - dos * (n * (vr - vrc) / (a * kappa2))

                    aux_y = vr - vrc
                    aux_x = a * kappa * (uno - (r - rc) / a)
                    varpi = lambda - atan2(aux_y, aux_x)

                    if (.not. coplanar) then
                        aux_y = nu * (z - zc)
                        aux_x = vz - vzc
                        capom = lambda - atan2(aux_y, aux_x)
                    else
                        capom = cero
                    end if


                    ! Compute corrections

                    rc = a * (e**2 * (1.5e0_wp * (eta2 / kappa2) &
                                    & - uno &
                                    & - uno2 * eta2 / kappa2 * cos(dos * (lambda - varpi))) &
                    & + inc**2 * (0.75e0_wp * chi2 / kappa2 &
                                & - uno &
                                & + 0.25e0_wp * chi2 / alpha2 * cos(dos * (lambda - capom))) &
                    &)
                    
                    Lc = n * (e**2 * (0.75e0_wp + uno2 * eta2 / kappa2) / kappa * sin(dos * (lambda - varpi)) &
                            & - inc**2 * 0.25e0_wp * chi2 / (alpha2 * nu) * sin(dos * (lambda - capom)) &
                        &)

                    zc = a * inc * e * chi2 / kappa * (uno2 / alpha_1 * sin(dos * lambda - varpi - capom) &
                                                    & - 1.5e0_wp / alpha_2 * sin(varpi - capom))
                    
                    vrc = a * (e**2 * eta2 / kappa * sin(dos * (lambda - varpi)) &
                            & - inc**2 * uno2 * chi2 / alpha2 * nu * sin(dos * (lambda - capom)))
                    
                    vLc = n * (e**2 * (3.5e0_wp &
                                    & - 3.e0_wp * eta2 / kappa2 &
                                    & - uno2 * kappa2 / n**2 &
                                    & + (1.5e0_wp + eta2 / kappa2) * cos(dos * (lambda - varpi))) &
                            & + inc**2 * uno2 * (4.e0_wp &
                                                & - kappa2 / n2 &
                                                & - 3.e0_wp * chi2 / kappa2 &
                                                & - chi2 / alpha2 * cos(dos * (lambda - capom))) &
                            &)
                    
                    vzc = a * e * inc * uno2 * chi2 / kappa * (&
                                    & (kappa + nu) / alpha_1 * cos(dos * lambda - varpi - capom) &
                                    & + 3.e0_wp * (kappa - nu) / alpha_2 * cos(varpi - capom) &
                                                            &)

                    ! Extra a corrections
                    if (extra_cor) a = r0 * (uno + e**2 + inc**2)
                    
                end subroutine get_geom
        
        end subroutine coord2geom
    
end module celestial