!> Module with main derivate function.
module derivates
    use, intrinsic :: ieee_arithmetic
    use constants, only: wp, G, cero, uno, uno2, uno3, dos, myepsilon, tini, megno_factor
    use auxiliary, only: cross2D_z, rotate2D
    use parameters, only: sim, &
                          & asteroid_data, &  !! |axis_a, axis_b, inertia|
                          & boulders_coords, boulders_data, & !! (Nb, 4) |mass,radius,theta_Ast0,dist_Ast|
                          & m_arr, R_arr, &
                          & hard_exit
    use accelerations, only: use_damp, damp_time, damp_coef_1, damp_coef_2, damp_model, &
                            & use_drag, use_drag_moons, drag_coef, drag_time, &
                            & use_stokes, use_stokes_moons, stokes_C, stokes_alpha, stokes_time, &
                            & use_ellipsoid, K_coef, L_coef, &
                            & use_manual_J2_from_cm, use_manual_J2_from_primary, J2K_coef, & 
                            & use_boulder_z, Gmboulder_z_coef, dz2_boulder_z_coef

    implicit none

    abstract interface
        ! Here must be every f_i defined explicitly
        function dydt_template(t, y) result(der)
            import :: wp
            implicit none
            real(wp), intent(in)               :: t
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(size(y))       :: der
        end function dydt_template

        function dydt_single_template(t, y) result(der)
            import :: wp
            implicit none
            real(wp), intent(in) :: t
            real(wp), intent(in) :: y
            real(wp) :: der
        end function dydt_single_template

    end interface

contains

    pure function get_index(i) result(idx)
        implicit none
        integer(kind=4), intent(in) :: i
        integer(kind=4) :: idx
        idx = 4*i - 1
    end function get_index

    pure function get_variational_index(i, first_particle, Ntotal) result(idx)
        implicit none
        integer(kind=4), intent(in) :: i, first_particle, Ntotal  ! Without asteroid
        integer(kind=4) :: idx
        idx = 4*(Ntotal + 1) - 1 + 7*(i - first_particle)
    end function get_variational_index

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! DERIVATIVES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function dydt(t, y) result(der)
        !y = /theta, omega, xA, yA, vxA, vyA, Moon, Part, .../
        implicit none
        real(wp), intent(in)               :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y))       :: der
        real(wp) :: theta, omega
        real(wp) :: coords_A(4), coords_M(4), coords_P(4), dr_vec(2), dr, dr2
        real(wp) :: acc_grav(2), acc_grav_m(2), torque
        integer(kind=4) :: i, idx
        integer(kind=4) :: j, jdx
        integer(kind=4) :: vdx  ! For variational
        integer(kind=4) :: N_total, N_particles, last_moon, first_particle 
        real(wp) :: c2th, s2th  ! For triaxial
        real(wp) :: Q_eff, dQdx, dQdy  ! For triaxial
        real(wp) :: inv_dr, inv_dr2, inv_dr3, inv_dr7  ! For triaxial and extra forces
        real(wp) :: theta_moon  ! For triaxial
        real(wp) :: xy_rotated(2)  ! For triaxial
        real(wp) :: Gmast, Gmcomb, Gmj, Gmi  ! For extra/COM forces
        real(wp) :: dr_ver(2), dv_vec(2)  ! For extra forces
        real(wp) :: vel_circ(2), v2  ! For extra forces
        real(wp) :: two_ener  ! For extra forces
        real(wp) :: aux_J2K, aux_inv_dr3_boulder_z  ! For extra forces
        real(wp) :: mean_movement  ! For extra forces
        real(wp) :: vel_radial(2), acc_radial_drag(2)  ! For extra forces
        real(wp) :: damp_f, drag_f, stokes_f  ! For extra forces
        real(wp) :: prod, dist, glob_prod, glob_dist  ! for MEGNO
        real(wp) :: aux_real

        der = cero  ! init der at cero

        last_moon = 1 + sim%Nmoon_active  ! This would be the last moon (+ 1 bc asteroid is 1)
        first_particle = last_moon + 1
        N_particles = sim%Npart_active
        N_total = last_moon + N_particles

        ! Calculate the angle of the asteroid
        theta = y(1)
        omega = y(2)
        der(1) = omega

        ! Set the position derivates: dX/dt = V
        do i = 1, N_total
            idx = get_index(i)
            der(idx:idx + 1) = y(idx + 2:idx + 3)
        end do

        coords_A = y(3:6)  ! Asteroid pos + vel
        torque = cero  ! Init torque at cero
        acc_grav = cero  ! Init acceleration at cero

        ! ---> Omega Damping <---
        if (use_damp) then
            !! Damping
            damp_f = uno2*(uno + tanh(1.e1_wp*(uno - t/damp_time)))
            select case (damp_model)
            case (1) ! domega/dt = tau
                der(2) = der(2) + damp_coef_1*damp_f
            case (2) ! domega/dt = -exp(- (t-t0) / tau) * omega0 / tau = - omega / tau
                der(2) = der(2) - omega/damp_coef_1*damp_f
            case (3) ! domega/dt = A * B * (t-t0)**(B-1) * omega0 * exp (A * (t-t0)**B) = (A * B * (t-t0)**(B-1)) * omega
                der(2) = der(2) +  &
                        & damp_coef_1*(t - cero + tini)**(damp_coef_2 - uno)*omega*damp_f
            end select
        end if

        ! ---> Forces / Escapes/ Collisions acting from COM of asteroid <---
        !! Includes triaxial

        !!! Set-up
        Gmast = G*m_arr(1)
        if (use_stokes) stokes_f = uno2*(uno + tanh(1.e1_wp*(uno - t/stokes_time)))
        if (use_drag) drag_f = uno2*(uno + tanh(1.e1_wp*(uno - t/drag_time)))
        if (use_ellipsoid) then
            c2th = cos(dos*theta)
            s2th = sin(dos*theta)
        end if

        !! Moons (massive)
        do j = 2, last_moon
            jdx = get_index(j)
            coords_M = y(jdx:jdx + 3)  ! Moon

            !! ASTEROID AND MOON
            dr_vec = coords_M(1:2) - coords_A(1:2)  ! From Asteroid to Moon
            dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
            dr = sqrt(dr2)

            ! Check if collision or Escape
            if ((dr < R_arr(1) + R_arr(j)) .or. &
              & (dr < sim%min_distance) .or. &
              & (dr > sim%max_distance .and. sim%max_distance > cero)) then
                hard_exit = .True.
            end if

            if (dr < myepsilon) cycle  ! Skip to avoid NaNs

            ! Extra needed
            Gmj = G*m_arr(j)
            inv_dr3 = uno/(dr2*dr)

            ! INIT Acceleration to 0
            acc_grav = cero

            ! ---> Triaxial <---
            if (use_ellipsoid) then

                ! Anti-rotate target to check if inside
                xy_rotated = rotate2D(dr_vec, -theta)
                if (((xy_rotated(1) + R_arr(j))/asteroid_data(1))**2 &
                & + ((xy_rotated(2) + R_arr(j))/asteroid_data(2))**2 < uno) then
                    hard_exit = .True.
                end if

                ! Extra needed
                inv_dr2 = inv_dr3*dr

                ! Q = (x²-y²) cos(2th) + 2xy sin(2th)
                ! Q_eff = 5 Q / r⁴
                Q_eff = 5*( (dr_vec(1)**2 - dr_vec(2)**2)*c2th + dos*dr_vec(1)*dr_vec(2)*s2th )*inv_dr2*inv_dr2
                dQdx = dos*(dr_vec(1)*c2th + dr_vec(2)*s2th)  ! 2x cos(2th) + 2y sin(2th)
                dQdy = -dos*(dr_vec(2)*c2th - dr_vec(1)*s2th)  ! - 2y cos(2th) + 2x sin(2th)

                ! Accelerations by triaxial
                ! a_x = -G mA / r³ (x - K x / r² - L (dQ/dx / r² - x 5 Q / r⁴))
                ! a_y = -G mA / r³ (y - K y / r² - L (dQ/dy / r² - y 5 Q / r⁴))
                acc_grav(1) = -(Gmast*inv_dr3)*( &
                                &  dr_vec(1) &
                                &  - K_coef*dr_vec(1)*inv_dr2 &
                                &  - L_coef*(dQdx*inv_dr2 - dr_vec(1)*Q_eff) &
                                &)
                acc_grav(2) = -(Gmast*inv_dr3)*( &
                                &  dr_vec(2) &
                                &  - K_coef*dr_vec(2)*inv_dr2 &
                                &  - L_coef*(dQdy*inv_dr2 - dr_vec(2)*Q_eff) &
                                &)

                ! Torque to Asteroid
                theta_moon = atan2(dr_vec(2), dr_vec(1))
                torque = torque - dos*m_arr(j)*L_coef*inv_dr3*sin(dos*(theta_moon - theta))

            ! ---> Manual J2 from asteroid CM <---
            else if (use_manual_J2_from_cm) then

                ! Accelerations by J2: G mA (x, y) (J2K / r²) / r³
                acc_grav = Gmast*dr_vec*J2K_coef/dr2*inv_dr3

            ! ---> Manual boulder_z from asteroid CM <---
            else if (use_boulder_z) then

                aux_inv_dr3_boulder_z = uno/(dr2 + dz2_boulder_z_coef)**(1.5e0_wp)

                ! Acceleration by boulders Z: -2 G (x, y) boulder_z / r³
                acc_grav = -Gmboulder_z_coef*dr_vec*aux_inv_dr3_boulder_z

            end if

            ! Acceleration of moon j from asteroid
            der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav

            !! ASTEROID FROM MOON: This force is "felt" by the asteroid CM
            ! Acceleration of asteroid from moon j [REACTION]
            der(5:6) = der(5:6) - acc_grav*m_arr(j)/m_arr(1)

            
            ! -------- NON CONSERVATIVE -------

            ! ---> Drag and/or Stokes <---
            if (use_drag_moons .or. use_stokes_moons) then
                Gmcomb = Gmast + Gmj
                inv_dr = inv_dr3*dr2
                dr_ver = dr_vec*inv_dr
                dv_vec = coords_M(3:4) - coords_A(3:4)  ! Velocity from Asteroid to Moon
                v2 = dot_product(dv_vec, dv_vec)

                ! Get energy
                two_ener = dos*Gmcomb*inv_dr - v2

                ! Check if unbound
                if (two_ener > cero) then ! Can calculate only in this case
                    mean_movement = abs(two_ener)**(1.5e0_wp)/Gmcomb ! n

                    ! ---> Drag <---
                    if (use_drag_moons) then
                        vel_radial = dot_product(dr_ver, dv_vec)
                        acc_radial_drag = -drag_coef*mean_movement*vel_radial

                        ! Acceleration of moon j by drag
                        !! acc = -a_r_drag * (x, y) / r * factor
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_radial_drag*dr_ver*drag_f

                    end if

                    ! ---> Stokes <---
                    if (use_stokes_moons) then
                        vel_circ = mean_movement*(/-dr_vec(2), dr_vec(1)/)  ! v_circ = n (-y, x)

                        ! Acceleration of moon j by Stokes
                        !! acc = -C * (v - alpha * vc) * factor
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) - stokes_C*(dv_vec - stokes_alpha*vel_circ)*stokes_f

                    end if

                end if

            end if            

        end do

        !! Particles (massless)
        do j = first_particle, N_total
            jdx = get_index(j)
            coords_P = y(jdx:jdx + 3)  ! Particle

            !! ASTEROID AND PARTICLE
            dr_vec = coords_P(1:2) - coords_A(1:2)  ! From Asteroid to Particle
            dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
            dr = sqrt(dr2)

            ! Check if collision or Escape
            if ((dr < sim%min_distance) .or. &
              & (dr > sim%max_distance .and. sim%max_distance > cero)) then
                hard_exit = .True.
            end if

            if (dr < myepsilon) cycle  ! Skip to avoid NaNs

            ! Extra needed
            inv_dr3 = uno/(dr2*dr)

            ! INIT Acceleration to 0
            acc_grav = cero

            ! ---> Triaxial <---
            if (use_ellipsoid) then

                ! Anti-rotate target to check if inside
                xy_rotated = rotate2D(dr_vec, -theta)
                if ((xy_rotated(1)/asteroid_data(1))**2 &
                & + (xy_rotated(2)/asteroid_data(2))**2 < uno) then
                    hard_exit = .True.
                end if

                ! Extra needed
                inv_dr2 = inv_dr3*dr

                ! Q = (x²-y²) cos(2th) + 2xy sin(2th)
                ! Q_eff = 5 Q / r⁴
                Q_eff = 5*( (dr_vec(1)**2 - dr_vec(2)**2)*c2th + dos*dr_vec(1)*dr_vec(2)*s2th )*inv_dr2*inv_dr2
                dQdx = dos*(dr_vec(1)*c2th + dr_vec(2)*s2th)  ! 2x cos(2th) + 2y sin(2th)
                dQdy = -dos*(dr_vec(2)*c2th - dr_vec(1)*s2th)  ! - 2y cos(2th) + 2x sin(2th)

                ! Accelerations by triaxial
                ! a_x = -G mA / r³ (x - K x / r² - L (dQ/dx / r² - x 5 Q / r⁴))
                ! a_y = -G mA / r³ (y - K y / r² - L (dQ/dy / r² - y 5 Q / r⁴))
                acc_grav(1) = -(Gmast*inv_dr3)*( dr_vec(1) &
                                &  - K_coef*dr_vec(1)*inv_dr2 &
                                &  - L_coef*(dQdx*inv_dr2 - dr_vec(1)*Q_eff) )
                acc_grav(2) = -(Gmast*inv_dr3)*( dr_vec(2) &
                                &  - K_coef*dr_vec(2)*inv_dr2 &
                                &  - L_coef*(dQdy*inv_dr2 - dr_vec(2)*Q_eff) )
                
                ! ! Variational [MEGNO]
                ! if (sim%megno_active) then
                !     vdx = get_variational_index(j, first_particle, N_total)
                !     coords_P = y(vdx:vdx + 3)  ! Variational particle
                    

                !     ! dvx = mu / r³ * (r² * (-dy (K - 3 r²) x y + dx (4 x⁴ + 5 x² y² + y⁴ - K (2 x² + y²))) &
                !     !                 & + L (6 dx x⁴ - 17 dy x³ y + 9 dx x² y² - 7 dy x y³ - 7 dx y⁴) cos(2th) &
                !     !                 & + 2 L (4 dy x⁴ + 4 dx x³ y - 3 dy x² y² + 9 dx x y³ - 2 dy y⁴) sin(2th))
                !     der(vdx + 2) = der(vdx + 2) + aux_real*m_arr(1) * (&
                !                 & dr2*( &
                !                         & -coords_P(2)*(K_coef - 3*dr2)*dr_vec(1)*dr_vec(2) &
                !                         & + coords_P(1)*(&
                !                             & 3*dr_vec(1)**4 &
                !                             & + 5*(dr_vec(1)*dr_vec(2))**2 &
                !                             & + dr2*dr2 &
                !                             & - K_coef*(dr_vec(1)**2 + dr2 ) &
                !                             & ) &
                !                         & ) &
                !                 & + L_coef*( &
                !                         & 6*coords_P(1)*dr_vec(1)**4 &
                !                         & -17*coords_P(2)*dr_vec(1)**3*dr_vec(2) &
                !                         & + 9*coords_P(1)*(dr_vec(1)*dr_vec(2))**2 &
                !                         & -7*coords_P(2)*dr_vec(1)*dr_vec(2)**3 &
                !                         & -7*coords_P(1)*dr_vec(2)**4 &
                !                         & )*c2th &
                !                 & + 2*L_coef*( &
                !                           & 4*coords_P(2)*dr_vec(1)**4 &
                !                           & + 4*coords_P(1)*dr_vec(1)**3*dr_vec(2) &
                !                           & -3*coords_P(2)*(dr_vec(1)*dr_vec(2))**2 &
                !                           & +9*coords_P(1)*dr_vec(1)*dr_vec(2)**3 &
                !                           & -2*coords_P(2)*dr_vec(2)**4 &
                !                           & )*s2th)
                    
                !     ! dvy = mu / r³ * (r² * (-dx (K - 3 r²) x y + dy (x⁴ + 5 x² y² + 4 y⁴ - K (x² + 2 y²))) &
                !     !                 & + L (7 dy x⁴ + 7 dx x³ y - 9 dy x² y² + 17 dy x y³ - 6 dy y⁴) cos(2th) &
                !     !                 & - 2 L (2 dx x⁴ - 9 dy x³ y + 3 dx x² y² - 4 dx x y³ - 4 dx y⁴) sin(2th))
                !     der(vdx + 3) = der(vdx + 3) + aux_real*m_arr(1) * (&
                !                 & dr2*( &
                !                         & -coords_P(1)*(K_coef - 3*dr2)*dr_vec(1)*dr_vec(2) &
                !                         & + coords_P(2)*(&
                !                             & dr2*dr2 &
                !                             & + 5*(dr_vec(1)*dr_vec(2))**2 &
                !                             & + 3*dr_vec(2)**4 &
                !                             & - K_coef*(dr2 + dr_vec(2)**2) &
                !                             & ) &
                !                         & ) &
                !                 & + L_coef*( &
                !                         & 7*coords_P(2)*dr_vec(1)**4 &
                !                         & + 7*coords_P(1)*dr_vec(1)**3*dr_vec(2) &
                !                         & -9*coords_P(2)*(dr_vec(1)*dr_vec(2))**2 &
                !                         & + 17*coords_P(1)*dr_vec(1)*dr_vec(2)**3 &
                !                         & -6*coords_P(2)*dr_vec(2)**4 &
                !                         & )*c2th &
                !                 & + 2*L_coef*( &
                !                           & 2*coords_P(1)*dr_vec(1)**4 &
                !                           & -9*coords_P(2)*dr_vec(1)**3*dr_vec(2) &
                !                           & + 3*coords_P(1)*(dr_vec(1)*dr_vec(2))**2 &
                !                           & -4*coords_P(2)*dr_vec(1)*dr_vec(2)**3 &
                !                           & -4*coords_P(1)*dr_vec(2)**4 &
                !                           & )*s2th)

                ! end if

            ! ---> Manual J2 from asteroid CM <---
            else if (use_manual_J2_from_cm) then

                ! Acceleration of particle j by J2
                acc_grav = Gmast*dr_vec*J2K_coef/dr2*inv_dr3  !! Add G (x, y) (J2K / r²) / r³

                ! Variational [MEGNO]
                if (sim%megno_active) then
                    vdx = get_variational_index(j, first_particle, N_total)
                    coords_P = y(vdx:vdx + 3)  ! Variational particle
                    
                    inv_dr7 = inv_dr3 * inv_dr3 / dr
                    aux_real = Gmast*J2K_coef*inv_dr7

                    ! dvx = mu J2k / r⁷ * (-5 dy x y + dx (-4 x² + y²))
                    der(vdx + 2) = der(vdx + 2) + aux_real*(&
                            & -5*coords_P(2)*dr_vec(1)*dr_vec(2) &
                            & + coords_P(1)*(-5*dr_vec(1)*dr_vec(1) + dr2))
                    
                    ! dvy = mu J2k / r⁷ * (-5 dx x y + dy (x² - 4 y²))
                    der(vdx + 3) = der(vdx + 3) + aux_real*(&
                            & -5*coords_P(1)*dr_vec(1)*dr_vec(2) &
                            & + coords_P(2)*(dr2 - 5*dr_vec(2)*dr_vec(2)))

                end if

            ! ---> Manual boulder_z from asteroid CM <---
            else if (use_boulder_z) then

                aux_inv_dr3_boulder_z = uno/(dr2 + dz2_boulder_z_coef)**(1.5e0_wp)

                ! Acceleration of particle j by boulders Z
                acc_grav = - Gmboulder_z_coef*dr_vec*aux_inv_dr3_boulder_z  !! Add - 2 G (x, y) boulder_z / r³

                ! Variational [MEGNO]
                if (sim%megno_active) then
                    vdx = get_variational_index(j, first_particle, N_total)
                    coords_P = y(vdx:vdx + 3)  ! Variational particle
                    
                    aux_real = Gmboulder_z_coef/(dr2 + dz2_boulder_z_coef)**(2.5e0_wp)

                    ! dvx = - mu / (r² + rz²)²·⁵ * (-3 dy x y + dx (-2 x² + y² + rz²))
                    der(vdx + 2) = der(vdx + 2) + aux_real*(&
                            & 3*coords_P(2)*dr_vec(1)*dr_vec(2) &
                            & - coords_P(1)*(-3*dr_vec(1)*dr_vec(1) + dr2 + dz2_boulder_z_coef))
                    
                    ! dvy = mu / (r² + rz²)²·⁵ * (3 dx x y - dy (x² - 2 y² + rz²))
                    der(vdx + 3) = der(vdx + 3) + aux_real*(&
                            & 3*coords_P(1)*dr_vec(1)*dr_vec(2) &
                            & - coords_P(2)*(dr2 - 3*dr_vec(2)*dr_vec(2) + dz2_boulder_z_coef))

                end if

            end if

            ! Acceleration of particle j from asteroid CM
            der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav

            ! -------- NON CONSERVATIVE -------

            ! ---> Drag and/or Stokes <---
            if (use_drag .or. use_stokes) then
                inv_dr = inv_dr3*dr2
                dr_ver = dr_vec*inv_dr
                dv_vec = coords_P(3:4) - coords_A(3:4)  ! Velocity from Asteroid to Particle
                v2 = dot_product(dv_vec, dv_vec)

                ! Get energy
                if (use_manual_J2_from_cm) then
                    aux_J2K = J2K_coef/dr2  ! J2K_coef is negative
                    two_ener = dos*Gmast*inv_dr*(uno - aux_J2K) - v2  ! Check if unbound
                    mean_movement = sqrt(Gmast*inv_dr3)*(uno - aux_J2K*uno3)
                else
                    two_ener = dos*Gmast*inv_dr - v2  ! Check if unbound
                    mean_movement = abs(two_ener)**(1.5e0_wp)/Gmast ! n
                end if

                if (two_ener > cero) then ! Can calculate only in this case

                    ! ---> Drag <---
                    if (use_drag) then
                        vel_radial = dot_product(dr_ver, dv_vec)
                        acc_radial_drag = -drag_coef*mean_movement*vel_radial

                        ! Acceleration of particle j by drag
                        !! acc = -a_r_drag * (x, y) / r * factor
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_radial_drag*dr_ver*drag_f

                    end if
                    ! ---> Stokes <---
                    if (use_stokes) then
                        vel_circ = mean_movement*(/-dr_vec(2), dr_vec(1)/)  ! v_circ = n (-y, x)

                        ! Acceleration of particle j by Stokes
                        !! acc = -C * (v - alpha * vc) * factor
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) - stokes_C*(dv_vec - stokes_alpha*vel_circ)*stokes_f

                    end if

                end if

            end if

        end do

        ! ---> GRAVITY (if not triaxial, only boulders) <---

        ! First, Asteroid to all
        if (.not. use_ellipsoid) then  ! Only if NOT triaxial
            ! First we do only boulder 0 (primary) for possible J2

            !! Get boulder coords
            boulders_coords(0, 1) = boulders_data(0, 4)*cos(theta + boulders_data(0, 3))  ! x
            boulders_coords(0, 2) = boulders_data(0, 4)*sin(theta + boulders_data(0, 3))  ! y
            boulders_coords(0, 3) = -omega*boulders_coords(0, 2)  ! vx
            boulders_coords(0, 4) = omega*boulders_coords(0, 1)  ! vy

            boulders_coords(0, :) = boulders_coords(0, :) + coords_A  ! Move to Asteroid

            ! Aux needed
            Gmi = G*boulders_data(0, 1)

            !! Moons (massive)
            do j = 2, last_moon ! +1 porque j=1 es asteroid
                jdx = get_index(j)
                coords_M = y(jdx:jdx + 3)  ! Moon

                !! BOULDER AND MOON
                dr_vec = coords_M(1:2) - boulders_coords(0, 1:2)  ! From Boulder 0 (primary) to Moon
                dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                dr = sqrt(dr2)

                if (dr < myepsilon) cycle  ! Skip to avoid NaNs

                ! Check if collision
                if (dr < boulders_data(0, 2) + R_arr(j)) hard_exit = .True.

                ! Gravitational acceleration of moon j by boulder 0 (and possible J2)
                if (use_manual_J2_from_primary) then  

                    ! Including extra J2
                    acc_grav = -Gmi*dr_vec/(dr2*dr)*(uno - J2K_coef/dr2)  !! -G m0 (x, y) (1 - K / r²) / r³
                else

                    ! Purely central
                    acc_grav = -Gmi*dr_vec/(dr2*dr)  !! -G m0 (x, y) / r³
                end if

                ! Acceleration of moon j by boulder 0
                der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav

                !! ASTEROID FROM MOON: This force is "felt" by the asteroid CM
                ! Acceleration of asteroid from moon j [REACTION]
                der(5:6) = der(5:6) - acc_grav*m_arr(j)/m_arr(1)  ! Force moved to Asteroid CM

                ! Torque to Asteroid
                torque = torque + cross2D_z(boulders_coords(0, 1:2) - coords_A(1:2), -acc_grav*m_arr(j))  !! r x F
            
            end do

            !! Particles (massless)
            do j = first_particle, N_total
                jdx = get_index(j)
                coords_P = y(jdx:jdx + 3)  ! Particle

                !! BOULDER AND PARTICLE
                dr_vec = coords_P(1:2) - boulders_coords(0, 1:2)  ! From Boulder 0 (primary) to Particle
                dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                dr = sqrt(dr2)

                if (dr < myepsilon) cycle  ! Skip to avoid NaNs

                ! Check if collision
                if (dr < boulders_data(0, 2)) hard_exit = .True.

                ! Gravitational acceleration of particle j by boulder 0 (and possible J2)
                if (use_manual_J2_from_primary) then  

                    ! Including extra J2
                    acc_grav = -Gmi*dr_vec/(dr2*dr)*(uno - J2K_coef/dr2)  !! -G m0 (x, y) (1 - K / r²) / r³

                    ! ! Variational [MEGNO]
                    ! if (sim%megno_active) then
                    !     vdx = get_variational_index(j, first_particle, N_total)

                    !     coords_P = y(vdx:vdx + 3)  ! Variational particle

                    !     inv_dr7 = uno / (dr2*dr2*dr2*dr)
                        
                    !     ! dvx = mu / r⁷ * (dy (-5 K + 3 r²) x y + dx (2 x² (-2 K + x²) + (K + x²) y² - y⁴)
                    !     der(vdx+2) = der(vdx+2) + aux_real*inv_dr7*(&
                    !         & coords_P(2)*(-5*J2K_coef + 3*dr2)*dr_vec(1)*dr_vec(2) &
                    !         & + coords_P(1)*( &
                    !             & 3*dr_vec(1)**4 &
                    !             & + (dr_vec(1)*dr_vec(2))**2 &
                    !             & -dr2*dr2 &
                    !             & + J2K_coef*(-5*dr_vec(1)*dr_vec(1) + dr2) &
                    !             & ) &
                    !         & )

                    !     ! dvy = mu / r⁷ * (dx (-5 K + 3 l^2) x y + dy (-x^4 + x^2 y^2 + 2 y^4 + K (x^2 - 4 y^2)))
                    !     der(vdx+3) = der(vdx+3) + aux_real*inv_dr7*(&
                    !         & coords_P(1)*(-5*J2K_coef + 3*dr2)*dr_vec(1)*dr_vec(2) &
                    !         & + coords_P(2)*( &
                    !             & -dr2*dr2 &
                    !             & + (dr_vec(1)*dr_vec(2))**2 &
                    !             & + 3*dr_vec(2)**4 &
                    !             & + J2K_coef*(-5*dr_vec(2)*dr_vec(2) + dr2) &
                    !             & ) &
                    !         & )

                    ! end if
                else

                    ! Purely central
                    acc_grav = -Gmi*dr_vec/(dr2*dr)  !! -G m0 (x, y) / r³

                    ! Variational [MEGNO]
                    if (sim%megno_active) then
                        vdx = get_variational_index(j, first_particle, N_total)

                        coords_P = y(vdx:vdx + 3)  ! Variational particle
                        
                        der(vdx + 2:vdx + 3) = der(vdx + 2:vdx + 3) - Gmi * &
                            & ( coords_P(1:2)*dr2 - 3 * dot_product(dr_vec,coords_P(1:2))*dr_vec) / (dr2*dr2*dr)

                    end if
                end if

                ! Acceleration of particle j by boulder 0
                der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav

            end do

            ! Now the rest of the boulders
            do i = 1, sim%Nboulders

                !! Get boulder coords
                boulders_coords(i, 1) = boulders_data(i, 4)*cos(theta + boulders_data(i, 3))  ! x
                boulders_coords(i, 2) = boulders_data(i, 4)*sin(theta + boulders_data(i, 3))  ! y
                boulders_coords(i, 3) = -omega*boulders_coords(i, 2)  ! vx
                boulders_coords(i, 4) = omega*boulders_coords(i, 1)  ! vy

                boulders_coords(i, :) = boulders_coords(i, :) + coords_A  ! Move to Asteroid

                ! Aux needed
                Gmi = G*boulders_data(i, 1)

                !! Moons (massive)
                do j = 2, last_moon ! +1 porque j=1 es asteroid
                    jdx = get_index(j)
                    coords_M = y(jdx:jdx + 3)  ! Moon

                    !! BOULDER AND MOON
                    dr_vec = coords_M(1:2) - boulders_coords(i, 1:2)  ! From Boulder i (primary) to Moon
                    dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                    dr = sqrt(dr2)

                    if (dr < myepsilon) cycle  ! Skip to avoid NaNs

                    ! Check if collision
                    if (dr < boulders_data(i, 2) + R_arr(j)) hard_exit = .True.

                    ! Gravitational acceleration of moon j by boulder i
                    acc_grav = -Gmi*dr_vec/(dr2*dr)  !! -G mi (x, y) / r³

                    ! Acceleration of moon j by boulder i
                    der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav

                    !! ASTEROID FROM MOON: This force is "felt" by the asteroid CM
                    ! Acceleration of asteroid from moon j [REACTION]
                    der(5:6) = der(5:6) - acc_grav*m_arr(j)/m_arr(1)  ! Force moved to Asteroid CM

                    ! Torque to Asteroid
                    torque = torque + cross2D_z(boulders_coords(i, 1:2) - coords_A(1:2), -acc_grav*m_arr(j))  !! r x F
                
                end do

                !! Particles (massless)
                do j = first_particle, N_total
                    jdx = get_index(j)
                    coords_P = y(jdx:jdx + 3)  ! Particle

                    !! BOULDER AND PARTICLE
                    dr_vec = coords_P(1:2) - boulders_coords(i, 1:2)  ! From Boulder i (primary) to Particle
                    dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                    dr = sqrt(dr2)

                    if (dr < myepsilon) cycle  ! Skip to avoid NaNs

                    ! Check if collision
                    if (dr < boulders_data(i, 2)) hard_exit = .True.

                    ! Gravitational acceleration of particle j by boulder i
                    acc_grav = -Gmi*dr_vec/(dr2*dr)  !! -G mi (x, y) / r³

                    ! Acceleration of particle j by boulder i
                    der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav

                    ! Variational [MEGNO]
                    if (sim%megno_active) then
                        vdx = get_variational_index(j, first_particle, N_total)

                        coords_P = y(vdx:vdx + 3)  ! Variational particle
                        
                        der(vdx + 2:vdx + 3) = der(vdx + 2:vdx + 3) - Gmi * &
                            & ( coords_P(1:2)*dr2 - 3 * dot_product(dr_vec,coords_P(1:2))*dr_vec) / (dr2*dr2*dr)

                    end if

                end do
            
            end do

        end if

        !! Update Omega with Torque
        der(2) = der(2) + torque/asteroid_data(3) ! Torque/Inertia

        ! Second, Moons to particles
        do i = 2, last_moon
            idx = get_index(i)
            coords_M = y(idx:idx + 3)  ! Moon i (M)

            !! Particles (massless)
            do j = first_particle, N_total
                jdx = get_index(j)
                coords_P = y(jdx:jdx + 3)  ! Particle

                !! MOON AND PARTICLE
                dr_vec = coords_P(1:2) - coords_M(1:2)  ! From Moon i (M) to Particle j (P)
                dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                dr = sqrt(dr2)

                if (dr < myepsilon) cycle  ! Skip to avoid NaNs

                ! Check if collision
                if (dr < R_arr(i)) hard_exit = .True.

                ! Auxiliar
                aux_real = -G*m_arr(i)/(dr2*dr)

                ! Acceleration of particle j by moon i
                der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + aux_real*dr_vec  ! -G mMoon (x, y) / r³

                ! Variational [MEGNO]
                if (sim%megno_active) then
                    vdx = get_variational_index(j, first_particle, N_total)

                    coords_P = y(vdx:vdx + 3)  ! Variational particle

                    der(vdx + 2:vdx + 3) = der(vdx + 2:vdx + 3) + aux_real / dr2 * &
                        & ( coords_P(1:2)*dr2 - 3 * dot_product(dr_vec,coords_P(1:2))*dr_vec)

                end if

            end do

        end do

        ! Third, Moons to Moons (if requested)
        if (sim%use_moon_gravity) then
            do i = 2, last_moon - 1
                idx = get_index(i)
                coords_M = y(idx:idx + 3)  ! Moon i (M)

                !! Other Moons (massive)
                do j = i + 1, last_moon ! Ahora sí +1 porque j=1 es asteroid
                    jdx = get_index(j)
                    coords_P = y(jdx:jdx + 3)  ! Moon 2 (P)

                    !! MOON 1 (M) AND MOON 2 (P)
                    dr_vec = coords_P(1:2) - coords_M(1:2)  ! From Moon i (M) to Moon j (P)
                    dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                    dr = sqrt(dr2)

                    if (dr < myepsilon) cycle  ! Skip to avoid NaNs

                    ! Check if collision
                    if (dr < R_arr(i) + R_arr(j)) hard_exit = .True.

                    ! Moons acceleration per unit mass
                    acc_grav_m = G*dr_vec/(dr2*dr)  !! G (x, y) / r³

                    ! Both accelerations
                    der(idx + 2:idx + 3) = der(idx + 2:idx + 3) + acc_grav_m*m_arr(j)
                    der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) - acc_grav_m*m_arr(i)

                end do

            end do

        end if

        ! Fourth, Extra variational if MEGNO
        if (sim%megno_active) then

            ! Initialize to 0
            glob_prod = cero
            glob_dist = cero

            ! Loop trhough particles
            do i = first_particle, N_total
                vdx = get_variational_index(i, first_particle, N_total)
                ! Update 'd positions' with 'd velocities'
                der(vdx:vdx + 1) = y(vdx + 2:vdx + 3)

                ! Get and update prod  and dist
                prod = y(vdx)*der(vdx) + y(vdx+1)*der(vdx+1) + y(vdx+2)*der(vdx+2) + y(vdx+3)*der(vdx+3)
                if (.not. ieee_is_finite(prod) .or. abs(prod) < tini) prod = cero

                dist = y(vdx)*y(vdx) + y(vdx+1)*y(vdx+1) + y(vdx+2)*y(vdx+2) + y(vdx+3)*y(vdx+3)
                if (.not. ieee_is_finite(dist) .or. dist < tini) dist = tini
                
                glob_prod = glob_prod + prod
                glob_dist = glob_dist + dist

                ! Calculate dot{lambda}
                der(vdx + 4) = prod / dist

                ! Calculate dot{Y}
                der(vdx + 5) = prod / dist * t / megno_factor

                ! Calculate dot{<Y>}
                if (t > 0) der(vdx + 6) = dos * y(vdx + 5) / t
                
                ! print*, t, i-first_particle, y(vdx + 5), y(vdx + 5) / max(t, tini)

            end do

            ! Now, we compute the global 

            ! Calculate dot{lambda}
            der(vdx + 7) = glob_prod / glob_dist

            ! Calculate dot{Y}
            der(vdx + 8) = glob_prod / glob_dist * t / megno_factor

            ! Calculate dot{<Y>}
            if (t > 0) der(vdx + 9) = dos * y(vdx + 8) / t

        endif

    end function dydt

end module derivates
