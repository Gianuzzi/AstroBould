!> Module with coordinates, and some cel-mech routines
module celestial
    use constants, only: cero, uno, uno2, dos, G, pi, twopi, tini, epsilon
    implicit none
    
    contains

        ! Get period
        function get_Period(mass, a) result(per)
            implicit none
            real(kind=8), intent(in) :: mass, a
            real(kind=8) :: per

            per = cero
            if (mass > cero) per = twopi * sqrt(a * a * a / (G * mass))
        end function get_Period

        ! Get corotation a from rotating body
        function get_a_corot(mass, omega) result(acorot)
            implicit none
            real(kind=8), intent(in) :: mass, omega
            real(kind=8) :: acorot

            if (omega < tini) then
                acorot = cero
            else 
                acorot = (G * mass / (omega * omega))**(1/3.)
            end if
        end function get_a_corot

        ! Get center of mass from masses, positions and velocities
        subroutine get_center_of_mass(mass, rib, vib, mcm, rcm, vcm)
            implicit none
            real(kind=8), intent(in) :: mass(:), rib(:,:), vib(:,:)
            real(kind=8), intent(out) :: mcm, rcm(2), vcm(2)
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

        ! Get acceleration and potential energy from single mass (MISSING G)
        subroutine get_acc_and_pot_single(mass, rib, xy_target, dr_max, acc, pot, inside)
            implicit none
            real(kind=8), intent(in) :: mass, rib(2), xy_target(2), dr_max
            real(kind=8), intent(inout) :: acc(2), pot
            logical, intent(inout), optional :: inside
            real(kind=8) :: dx, dy, dr, dr2

            dx = xy_target(1) - rib(1)
            dy = xy_target(2) - rib(2)
            dr2 = dx * dx + dy * dy
            dr = sqrt(dr2)
            if (dr < max(dr_max, tini)) then
                if (present(inside)) inside = .True.
                return
            end if
            acc = acc - (mass / (dr2 * dr)) * (/dx, dy/)
            pot = pot - mass / dr
            if (present(inside)) inside = .False.
        end subroutine get_acc_and_pot_single

        ! Get coordinates from a body
        subroutine coord(msum, a, e, inc, capm, omega, capom, xc)
            implicit none
            real(kind=8), intent(in)  :: msum, a, e, inc, capm, omega, capom
            real(kind=8), intent(out) :: xc(6)
            real(kind=8) :: sp, cp, so, co, si, ci
            real(kind=8) :: d11, d12, d13, d21, d22, d23
            real(kind=8) :: cape, dummy, scap, ccap, sqe, sqgma
            real(kind=8) :: ri, xfac1, xfac2, vfac1, vfac2
            
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
            real(kind=8), intent(in)  :: dm, e
            real(kind=8), intent(out) :: u, f
            real(kind=8) :: u0, dif, seno, cose
            integer(kind=4) :: i, MAX_ITER = 50
            
            u0 = dm
            dif = uno
            i = 0
            do while (dif > epsilon)
                u = dm + e * sin(u0)
                dif = abs(u - u0)
                u0 = u
                i = i + 1
                if (i > MAX_ITER) then
                    print*, 'aver: too many iterations'
                    exit
                end if
            end do
            seno = sqrt(uno + e) * sin(uno2 * u)
            cose = sqrt(uno - e) * cos(uno2 * u)
            f = dos * atan2(seno, cose)
        end subroutine aver

        ! Get elements from a body
        subroutine elem(msum, xc, a, e, inc, capm, omega, capom)
            implicit none
            real(kind=8), intent(in)  :: msum, xc(6)
            real(kind=8), intent(out) :: a, e, inc, capm, omega, capom
            real(kind=8) :: gmsum
            real(kind=8) :: x, y, z, vx, vy, vz
            real(kind=8) :: hx, hy, hz, h2, h, fac, u
            real(kind=8) :: r, v, v2, vdotr, energy
            real(kind=8) :: cape, cw, sw, w, face, capf, tmpf
            integer(kind=4) :: ialpha = 0 ! Just initialization
                
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
                if (abs(inc - pi) < 10.d0 * epsilon) then
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
            
            if (abs(energy * r / gmsum) < sqrt(epsilon)) then
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
                capm = tmpf * (uno + tmpf * tmpf / 3.d0)
                omega = u - w
                if (omega < cero) then 
                    omega = omega + twopi
                end if
                omega = omega - int(omega / twopi) * twopi
            end if
        end subroutine elem

        ! Get astrocentric from jacobi coordinates
        subroutine jacobi_to_astroc (mass_array, coords_j, coords_a)
            implicit none
            real(kind=8), dimension(:), intent(in) :: mass_array
            real(kind=8), dimension(:,:), intent(in) :: coords_j
            real(kind=8), dimension(:,:), intent(out) :: coords_a
            real(kind=8) :: mass_acum, mass_acum_prev
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
        subroutine astroc_to_jacob (mass_array, coords_a, coords_j)
            implicit none
            real(kind=8), dimension(:), intent(in) :: mass_array
            real(kind=8), dimension(:,:), intent(in) :: coords_a
            real(kind=8), dimension(:,:), intent(out) :: coords_j
            real(kind=8) :: mass_acum, mass_acum_prev
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
        subroutine baric_to_astroc (mass_array, coords_b, coords_a)
            implicit none
            real(kind=8), dimension(:), intent(in) :: mass_array
            real(kind=8), dimension(:,:), intent(in) :: coords_b
            real(kind=8), dimension(:,:), intent(out) :: coords_a
            integer(kind=4) :: i
            real(kind=8), dimension(4) :: xyvxvy = cero
            
            do i = 2, size(coords_b,1)
                xyvxvy = xyvxvy - coords_b(i,1:4) * mass_array(i) / mass_array(1)
            end do
            do i = 2, size(coords_b,1)
                coords_a(i,1:4) = coords_b(i,1:4) - xyvxvy
            end do
            coords_a(1, 1:4) = cero
        end subroutine baric_to_astroc

        ! Get astrocentric from barycentric coordinates
        subroutine astroc_to_baric (mass_array, coords_a, coords_b)
            implicit none
            real(kind=8), dimension(:), intent(in) :: mass_array
            real(kind=8), dimension(:,:), intent(in) :: coords_a
            real(kind=8), dimension(:,:), intent(out) :: coords_b
            integer(kind=4) :: i
            real(kind=8) :: mass_tot
            real(kind=8), dimension(4) :: xyvxvy = cero

            mass_tot = mass_array(1)
            do i = 2, size(coords_a,1)
                mass_tot = mass_tot + mass_array(i)
                xyvxvy = xyvxvy + coords_a(i,1:4) * mass_array(i)
            end do
            xyvxvy = - xyvxvy / mass_tot
            do i = 1, size(coords_a,1)
                coords_b(i,1:4) = coords_a(i,1:4) + xyvxvy
            end do
        end subroutine astroc_to_baric
    
end module celestial