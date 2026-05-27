!> Module with main derivate function.
module derivates
    use, intrinsic :: ieee_arithmetic
    use iso_fortran_env, only: int64
    use omp_lib
    use constants, only: wp, G, cero, uno, uno2, uno3, dos, twopi, tini, megno_factor, infinito
    use auxiliary, only: cross2D_z, rotate2D
    use parameters, only: sim, &
                          & asteroid_data, &
                          & boulders_coords, boulders_data, &
                          & m_arr, R_arr, &
                          & hard_exit, &
                          & get_index
    use accelerations, only: use_damp, damp_time, damp_coef_1, damp_coef_2, damp_model, &
                            & use_drag, use_drag_moons, drag_coef, drag_time, &
                            & use_stokes, use_stokes_moons, stokes_C, stokes_alpha, stokes_time, &
                            & use_ellipsoid, K_coef, L_coef, &
                            & use_manual_J2_from_cm, use_manual_J2_from_primary, J2K_coef, &
                            & use_boulder_z, Gmboulder_z_coef, dz2_boulder_z_coef
    use collisions, only: init_collisions, collisions_brute, collisions_grid, collisions_verlet

    implicit none
    private
    public :: dydt, set_dydt, dydt_grav_f, dydt_coll_f

    abstract interface
        function dydt_template(t, y) result(der)
            import :: wp
            implicit none
            real(wp), intent(in) :: t
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(size(y)) :: der
        end function dydt_template

        subroutine dydt_grav_template(t, y, der, first_particle, N_total)
            import :: wp
            implicit none
            real(wp), intent(in) :: t
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(:), intent(inout) :: der
            integer(kind=4), intent(in) :: first_particle, N_total
        end subroutine dydt_grav_template

    end interface

    procedure(dydt_grav_template), pointer :: dydt_grav => null()

contains

    subroutine set_dydt(Ntotal, sinodic)
        implicit none
        integer(kind=4), intent(in) :: Ntotal
        logical, intent(in) :: sinodic

        if (sinodic) then
            dydt_grav => dydt_grav_sinodic
        else
            dydt_grav => dydt_grav_inertial
        end if

        call init_collisions(Ntotal, sinodic)
    end subroutine set_dydt

    pure function get_variational_index(i, first_particle, Ntotal) result(idx)
        implicit none
        integer(kind=4), intent(in) :: i, first_particle, Ntotal
        integer(kind=4) :: idx
        idx = 4*(Ntotal + 1) - 1 + 7*(i - first_particle)
    end function get_variational_index

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                        HELPER SUBROUTINES                               !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine set_pos_derivatives(y, der, N_total)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in) :: N_total

        integer(kind=4) :: i, idx

        der(1) = y(2)

        do i = 1, N_total
            idx = get_index(i)
            der(idx:idx + 1) = y(idx + 2:idx + 3)
        end do

    end subroutine set_pos_derivatives

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                       GRAVITATIONAL DERIVATIVES                         !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine dydt_grav_inertial(t, y, der, first_particle, N_total)
        implicit none
        real(wp),               intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in) :: first_particle, N_total

        ! ── Scalar shared state ──────────
        real(wp) :: theta, omega
        real(wp) :: coords_A(4)
        real(wp) :: torque
        real(wp) :: Gmast, Gmi, Gmj, Gmcomb
        real(wp) :: c2th, s2th
        real(wp) :: damp_f, drag_f, stokes_f
        real(wp) :: rcoll, rescape

        ! ── Per-iteration locals (PRIVATE in parallel regions) ──────
        real(wp) :: coords_M(4), coords_P(4)
        real(wp) :: dr_vec(2), dr, dr2
        real(wp) :: acc_grav(2), acc_grav_m(2)
        real(wp) :: inv_dr, inv_dr2, inv_dr3, inv_dr7
        real(wp) :: Q_eff, dQdx, dQdy
        real(wp) :: theta_moon, xy_rotated(2)
        real(wp) :: dr_ver(2), dv_vec(2)
        real(wp) :: vel_circ(2), vel_radial(2), acc_radial_drag(2)
        real(wp) :: v2, two_ener, mean_movement
        real(wp) :: aux_J2K, aux_inv_dr3_boulder_z, aux_real
        integer(kind=4) :: i, idx, j, jdx, vdx, last_moon

        last_moon = first_particle - 1
        theta = y(1)
        omega = y(2)
        coords_A = y(3:6)
        torque = cero
        Gmast = G*m_arr(1)

        if (use_stokes)   stokes_f = uno2*(uno + tanh(1.e1_wp*(uno - t/stokes_time)))
        if (use_drag)     drag_f = uno2*(uno + tanh(1.e1_wp*(uno - t/drag_time)))
        if (use_ellipsoid) then
            c2th = cos(dos*theta)
            s2th = sin(dos*theta)
        end if
        if (sim%max_distance <= cero) then
            rescape = infinito
        else
            rescape = sim%max_distance
        end if

        ! ── Omega damping (serial — single scalar write) ──────
        if (use_damp) then
            damp_f = uno2*(uno + tanh(1.e1_wp*(uno - t/damp_time)))
            select case (damp_model)
                case (1)
                    der(2) = der(2) + damp_coef_1*damp_f
                case (2)
                    der(2) = der(2) - omega/damp_coef_1*damp_f
                case (3)
                    der(2) = der(2) + &
                            & damp_coef_1*(t - cero + tini)**(damp_coef_2 - uno)*omega*damp_f
            end select
        end if

        ! ===================================================================
        ! ── Asteroid COM → moons  (serial: der(5:6) and torque conflict) ───
        ! ===================================================================
        do j = 2, last_moon
            jdx = get_index(j)
            coords_M = y(jdx:jdx + 3)

            dr_vec = coords_M(1:2) - coords_A(1:2)
            dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)

            rcoll = min(sim%min_distance, R_arr(1) + R_arr(j))
            if ((dr2 < rcoll*rcoll) .or. (dr2 > rescape*rescape)) then
                hard_exit = .True.
                cycle
            end if

            dr = sqrt(dr2)
            Gmj = G*m_arr(j)
            inv_dr3 = uno/(dr2*dr)
            acc_grav = cero

            if (use_ellipsoid) then
                xy_rotated = rotate2D(dr_vec, -theta)
                if (((xy_rotated(1) + R_arr(j))/asteroid_data(1))**2 &
                &  + ((xy_rotated(2) + R_arr(j))/asteroid_data(2))**2 < uno) then
                    hard_exit = .True.
                end if
                inv_dr2 = inv_dr3*dr
                Q_eff = 5*( (dr_vec(1)**2 - dr_vec(2)**2)*c2th &
                           &   + dos*dr_vec(1)*dr_vec(2)*s2th )*inv_dr2*inv_dr2
                dQdx = dos*(dr_vec(1)*c2th + dr_vec(2)*s2th)
                dQdy = -dos*(dr_vec(2)*c2th - dr_vec(1)*s2th)
                acc_grav(1) = -(Gmast*inv_dr3)*( dr_vec(1) &
                              &  - K_coef*dr_vec(1)*inv_dr2 &
                              &  - L_coef*(dQdx*inv_dr2 - dr_vec(1)*Q_eff) )
                acc_grav(2) = -(Gmast*inv_dr3)*( dr_vec(2) &
                              &  - K_coef*dr_vec(2)*inv_dr2 &
                              &  - L_coef*(dQdy*inv_dr2 - dr_vec(2)*Q_eff) )
                theta_moon = atan2(dr_vec(2), dr_vec(1))
                torque = torque - dos*m_arr(j)*L_coef*inv_dr3*sin(dos*(theta_moon - theta))

            else if (use_manual_J2_from_cm) then
                acc_grav = Gmast*dr_vec*J2K_coef/dr2*inv_dr3

            else if (use_boulder_z) then
                aux_inv_dr3_boulder_z = uno/(dr2 + dz2_boulder_z_coef)**(1.5e0_wp)
                acc_grav = -Gmboulder_z_coef*dr_vec*aux_inv_dr3_boulder_z
            end if

            der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav
            der(5:6) = der(5:6) - acc_grav*m_arr(j)/m_arr(1)

            if (use_drag_moons .or. use_stokes_moons) then
                Gmcomb = Gmast + Gmj
                inv_dr = inv_dr3*dr2
                dr_ver = dr_vec*inv_dr
                dv_vec = coords_M(3:4) - coords_A(3:4)
                v2 = dot_product(dv_vec, dv_vec)
                two_ener = dos*Gmcomb*inv_dr - v2
                if (two_ener > cero) then
                    mean_movement = abs(two_ener)**(1.5e0_wp)/Gmcomb
                    if (use_drag_moons) then
                        vel_radial = dot_product(dr_ver, dv_vec)
                        acc_radial_drag = -drag_coef*mean_movement*vel_radial
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) &
                                              & + acc_radial_drag*dr_ver*drag_f
                    end if
                    if (use_stokes_moons) then
                        vel_circ = mean_movement*(/-dr_vec(2), dr_vec(1)/)
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) &
                                              & - stokes_C*(dv_vec - stokes_alpha*vel_circ)*stokes_f
                    end if
                end if
            end if

        end do  ! asteroid COM → moons (serial)

        ! ===================================================================
        ! ── Asteroid COM → particles  (PARALLEL: each j owns unique slots) ─
        ! ===================================================================
        !$OMP PARALLEL DO                             &
        !$OMP   DEFAULT(SHARED)                       &
        !$OMP   PRIVATE(j, jdx, vdx,                  &
        !$OMP           coords_P, dr_vec, dr, dr2,    &
        !$OMP           rcoll, inv_dr, inv_dr2,       &
        !$OMP           inv_dr3, inv_dr7,             &
        !$OMP           acc_grav, Q_eff, dQdx, dQdy,  &
        !$OMP           xy_rotated,                   &
        !$OMP           aux_inv_dr3_boulder_z,        &
        !$OMP           aux_J2K, aux_real,            &
        !$OMP           dr_ver, dv_vec,               &
        !$OMP           vel_circ, vel_radial,         &
        !$OMP           acc_radial_drag,              &
        !$OMP           v2, two_ener, mean_movement)  &
        !$OMP   SCHEDULE(DYNAMIC, 256)
        do j = first_particle, N_total
            jdx = get_index(j)
            coords_P = y(jdx:jdx + 3)

            dr_vec = coords_P(1:2) - coords_A(1:2)
            dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)

            rcoll = min(sim%min_distance, R_arr(1) + R_arr(j))
            if ((dr2 < rcoll*rcoll) .or. (dr2 > rescape*rescape)) then
                !$OMP ATOMIC WRITE
                hard_exit = .True.
                cycle
            end if

            dr = sqrt(dr2)
            inv_dr3 = uno/(dr2*dr)
            acc_grav = cero

            if (use_ellipsoid) then
                xy_rotated = rotate2D(dr_vec, -theta)
                if (((xy_rotated(1) + R_arr(j))/asteroid_data(1))**2 &
                &  + ((xy_rotated(2) + R_arr(j))/asteroid_data(2))**2 < uno) then
                    !$OMP ATOMIC WRITE
                    hard_exit = .True.
                    cycle
                end if
                inv_dr2 = inv_dr3*dr
                Q_eff = 5*( (dr_vec(1)**2 - dr_vec(2)**2)*c2th &
                          &   + dos*dr_vec(1)*dr_vec(2)*s2th )*inv_dr2*inv_dr2
                dQdx = dos*(dr_vec(1)*c2th + dr_vec(2)*s2th)
                dQdy = -dos*(dr_vec(2)*c2th - dr_vec(1)*s2th)
                acc_grav(1) = -(Gmast*inv_dr3)*( dr_vec(1) &
                              &  - K_coef*dr_vec(1)*inv_dr2 &
                              &  - L_coef*(dQdx*inv_dr2 - dr_vec(1)*Q_eff) )
                acc_grav(2) = -(Gmast*inv_dr3)*( dr_vec(2) &
                              &  - K_coef*dr_vec(2)*inv_dr2 &
                              &  - L_coef*(dQdy*inv_dr2 - dr_vec(2)*Q_eff) )

            else if (use_manual_J2_from_cm) then
                acc_grav = Gmast*dr_vec*J2K_coef/dr2*inv_dr3
                if (sim%megno_active) then
                    vdx = get_variational_index(j, first_particle, N_total)
                    coords_P = y(vdx:vdx + 3)
                    inv_dr7 = inv_dr3 * inv_dr3 / dr
                    aux_real = Gmast*J2K_coef*inv_dr7
                    der(vdx + 2) = der(vdx + 2) + aux_real*( &
                            & -5*coords_P(2)*dr_vec(1)*dr_vec(2) &
                            & + coords_P(1)*(-5*dr_vec(1)*dr_vec(1) + dr2))
                    der(vdx + 3) = der(vdx + 3) + aux_real*( &
                            & -5*coords_P(1)*dr_vec(1)*dr_vec(2) &
                            & + coords_P(2)*(dr2 - 5*dr_vec(2)*dr_vec(2)))
                end if

            else if (use_boulder_z) then
                aux_inv_dr3_boulder_z = uno/(dr2 + dz2_boulder_z_coef)**(1.5e0_wp)
                acc_grav = -Gmboulder_z_coef*dr_vec*aux_inv_dr3_boulder_z
                if (sim%megno_active) then
                    vdx = get_variational_index(j, first_particle, N_total)
                    coords_P = y(vdx:vdx + 3)
                    aux_real = Gmboulder_z_coef/(dr2 + dz2_boulder_z_coef)**(2.5e0_wp)
                    der(vdx + 2) = der(vdx + 2) + aux_real*( &
                            &  3*coords_P(2)*dr_vec(1)*dr_vec(2) &
                            &  - coords_P(1)*(-3*dr_vec(1)*dr_vec(1) + dr2 + dz2_boulder_z_coef))
                    der(vdx + 3) = der(vdx + 3) + aux_real*( &
                            &  3*coords_P(1)*dr_vec(1)*dr_vec(2) &
                            &  - coords_P(2)*(dr2 - 3*dr_vec(2)*dr_vec(2) + dz2_boulder_z_coef))
                end if
            end if

            der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav

            if (use_drag .or. use_stokes) then
                inv_dr = inv_dr3*dr2
                dr_ver = dr_vec*inv_dr
                dv_vec = coords_P(3:4) - coords_A(3:4)
                v2 = dot_product(dv_vec, dv_vec)
                if (use_manual_J2_from_cm) then
                    aux_J2K = J2K_coef/dr2
                    two_ener = dos*Gmast*inv_dr*(uno - aux_J2K) - v2
                    mean_movement = sqrt(Gmast*inv_dr3)*(uno - aux_J2K*uno3)
                else
                    two_ener = dos*Gmast*inv_dr - v2
                    mean_movement = abs(two_ener)**(1.5e0_wp)/Gmast
                end if
                if (two_ener > cero) then
                    if (use_drag) then
                        vel_radial = dot_product(dr_ver, dv_vec)
                        acc_radial_drag = -drag_coef*mean_movement*vel_radial
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) &
                                              & + acc_radial_drag*dr_ver*drag_f
                    end if
                    if (use_stokes) then
                        vel_circ = mean_movement*(/-dr_vec(2), dr_vec(1)/)
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) &
                                              & - stokes_C*(dv_vec - stokes_alpha*vel_circ)*stokes_f
                    end if
                end if
            end if

        end do  ! asteroid COM → particles
        !$OMP END PARALLEL DO

        ! ================================================================
        ! ── Boulder gravity  (only when NOT triaxial) ───────────────────
        ! ================================================================
        if (.not. use_ellipsoid) then

            ! ── Boulder 0 (serial coordinate setup + serial moon loop) ───
            boulders_coords(0, 1) = boulders_data(0, 4)*cos(theta + boulders_data(0, 3))
            boulders_coords(0, 2) = boulders_data(0, 4)*sin(theta + boulders_data(0, 3))
            boulders_coords(0, 3) = -omega*boulders_coords(0, 2)
            boulders_coords(0, 4) = omega*boulders_coords(0, 1)
            boulders_coords(0, :) = boulders_coords(0, :) + coords_A
            Gmi = G*boulders_data(0, 1)

            do j = 2, last_moon  ! serial: der(5:6) and torque conflict
                jdx = get_index(j)
                coords_M = y(jdx:jdx + 3)
                dr_vec = coords_M(1:2) - boulders_coords(0, 1:2)
                dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                rcoll = boulders_data(0, 2) + R_arr(j)
                if (dr2 < rcoll*rcoll) then
                    hard_exit = .True.
                    cycle
                end if
                dr = sqrt(dr2)
                if (use_manual_J2_from_primary) then
                    acc_grav = -Gmi*dr_vec/(dr2*dr)*(uno - J2K_coef/dr2)
                else
                    acc_grav = -Gmi*dr_vec/(dr2*dr)
                end if
                der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav
                der(5:6) = der(5:6) - acc_grav*m_arr(j)/m_arr(1)
                torque = torque + cross2D_z(boulders_coords(0, 1:2) - coords_A(1:2), -acc_grav*m_arr(j))
            end do

            ! ── Boulder 0 → particles  (PARALLEL) ──────────────────
            !$OMP PARALLEL DO                             &
            !$OMP   DEFAULT(SHARED)                       &
            !$OMP   PRIVATE(j, jdx, vdx,                  &
            !$OMP           coords_P, dr_vec, dr, dr2,    &
            !$OMP           rcoll, acc_grav, aux_real)    &
            !$OMP   SCHEDULE(STATIC)
            do j = first_particle, N_total
                jdx = get_index(j)
                coords_P = y(jdx:jdx + 3)
                dr_vec = coords_P(1:2) - boulders_coords(0, 1:2)
                dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                rcoll = boulders_data(0, 2) + R_arr(j)
                if (dr2 < rcoll*rcoll) then
                    !$OMP ATOMIC WRITE
                    hard_exit = .True.
                    cycle
                end if
                dr = sqrt(dr2)
                if (use_manual_J2_from_primary) then
                    acc_grav = -Gmi*dr_vec/(dr2*dr)*(uno - J2K_coef/dr2)
                else
                    acc_grav = -Gmi*dr_vec/(dr2*dr)
                    if (sim%megno_active) then
                        vdx = get_variational_index(j, first_particle, N_total)
                        coords_P = y(vdx:vdx + 3)
                        der(vdx + 2:vdx + 3) = der(vdx + 2:vdx + 3) - Gmi * &
                            & ( coords_P(1:2)*dr2 - 3*dot_product(dr_vec, coords_P(1:2))*dr_vec) &
                            & / (dr2*dr2*dr)
                    end if
                end if
                der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav
            end do
            !$OMP END PARALLEL DO

            ! ── Remaining boulders (serial outer i, parallel inner j) ─────
            do i = 1, sim%Nboulders
                boulders_coords(i, 1) = boulders_data(i, 4)*cos(theta + boulders_data(i, 3))
                boulders_coords(i, 2) = boulders_data(i, 4)*sin(theta + boulders_data(i, 3))
                boulders_coords(i, 3) = -omega*boulders_coords(i, 2)
                boulders_coords(i, 4) = omega*boulders_coords(i, 1)
                boulders_coords(i, :) = boulders_coords(i, :) + coords_A
                Gmi = G*boulders_data(i, 1)

                do j = 2, last_moon  ! serial: conflicts on der(5:6) and torque
                    jdx = get_index(j)
                    coords_M = y(jdx:jdx + 3)
                    dr_vec = coords_M(1:2) - boulders_coords(i, 1:2)
                    dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                    rcoll = boulders_data(i, 2) + R_arr(j)
                    if (dr2 < rcoll*rcoll) then
                        hard_exit = .True.
                        cycle
                    end if
                    dr = sqrt(dr2)
                    acc_grav = -Gmi*dr_vec/(dr2*dr)
                    der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav
                    der(5:6) = der(5:6) - acc_grav*m_arr(j)/m_arr(1)
                    torque = torque + cross2D_z(boulders_coords(i, 1:2) - coords_A(1:2), -acc_grav*m_arr(j))
                end do

                !$OMP PARALLEL DO                             &
                !$OMP   DEFAULT(SHARED)                       &
                !$OMP   PRIVATE(j, jdx, vdx,                  &
                !$OMP           coords_P, dr_vec, dr, dr2,    &
                !$OMP           rcoll, acc_grav, aux_real)    &
                !$OMP   SCHEDULE(STATIC)
                do j = first_particle, N_total
                    jdx = get_index(j)
                    coords_P = y(jdx:jdx + 3)
                    dr_vec = coords_P(1:2) - boulders_coords(i, 1:2)
                    dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                    rcoll = boulders_data(i, 2) + R_arr(j)
                    if (dr2 < rcoll*rcoll) then
                        !$OMP ATOMIC WRITE
                        hard_exit = .True.
                        cycle
                    end if
                    dr = sqrt(dr2)
                    acc_grav = -Gmi*dr_vec/(dr2*dr)
                    der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav
                    if (sim%megno_active) then
                        vdx = get_variational_index(j, first_particle, N_total)
                        coords_P = y(vdx:vdx + 3)
                        der(vdx + 2:vdx + 3) = der(vdx + 2:vdx + 3) - Gmi * &
                            & ( coords_P(1:2)*dr2 - 3*dot_product(dr_vec, coords_P(1:2))*dr_vec) &
                            & / (dr2*dr2*dr)
                    end if
                end do
                !$OMP END PARALLEL DO

            end do  ! boulder loop

        end if  ! .not. use_ellipsoid

        ! ── Torque → asteroid spin (serial) ───────────────────────────
        der(2) = der(2) + torque/asteroid_data(3)

        ! =================================================================
        ! ── Mutual moon gravity  (serial: both i and j sides written) ────
        ! =================================================================
        if (sim%use_moon_gravity) then
            do i = 2, last_moon - 1
                idx = get_index(i)
                coords_M = y(idx:idx + 3)
                do j = i + 1, last_moon
                    jdx = get_index(j)
                    coords_P = y(jdx:jdx + 3)
                    dr_vec = coords_P(1:2) - coords_M(1:2)
                    dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                    rcoll = R_arr(i) + R_arr(j)
                    if (dr2 < rcoll*rcoll) then
                        hard_exit = .True.
                        cycle
                    end if
                    dr = sqrt(dr2)
                    acc_grav_m = G*dr_vec/(dr2*dr)
                    der(idx + 2:idx + 3) = der(idx + 2:idx + 3) + acc_grav_m*m_arr(j)
                    der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) - acc_grav_m*m_arr(i)
                end do
            end do
        end if

        ! ===================================================================
        ! ── Moon → particles  (serial outer i, PARALLEL inner j) ───────────
        ! ===================================================================
        do i = 2, last_moon
            idx = get_index(i)
            coords_M = y(idx:idx + 3)

            !$OMP PARALLEL DO                             &
            !$OMP   DEFAULT(SHARED)                       &
            !$OMP   PRIVATE(j, jdx, vdx,                  &
            !$OMP           coords_P, dr_vec, dr, dr2,    &
            !$OMP           rcoll, aux_real)              &
            !$OMP   SCHEDULE(STATIC)
            do j = first_particle, N_total
                jdx = get_index(j)
                coords_P = y(jdx:jdx + 3)
                dr_vec = coords_P(1:2) - coords_M(1:2)
                dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                rcoll = R_arr(i) + R_arr(j)
                if (dr2 < rcoll*rcoll) then
                    !$OMP ATOMIC WRITE
                    hard_exit = .True.
                    cycle
                end if
                dr = sqrt(dr2)
                aux_real = -G*m_arr(i)/(dr2*dr)
                der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + aux_real*dr_vec
                if (sim%megno_active) then
                    vdx = get_variational_index(j, first_particle, N_total)
                    coords_P = y(vdx:vdx + 3)
                    der(vdx + 2:vdx + 3) = der(vdx + 2:vdx + 3) + aux_real/dr2 * &
                        & ( coords_P(1:2)*dr2 - 3*dot_product(dr_vec, coords_P(1:2))*dr_vec )
                end if
            end do
            !$OMP END PARALLEL DO

        end do  ! moon → particles

    end subroutine dydt_grav_inertial


    subroutine dydt_grav_sinodic(t, y, der, dummy, N_total)
        implicit none
        real(wp),               intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in) :: dummy, N_total

        real(wp) :: theta, omega
        real(wp) :: coords_P(4), dr_vec(2), dr, dr2
        real(wp) :: acc_grav(2)
        real(wp) :: inv_dr, inv_dr2, inv_dr3, inv_dr7
        real(wp) :: Q_eff
        real(wp) :: Gmast, Gmi
        real(wp) :: dr_ver(2), dv_vec(2)
        real(wp) :: vel_circ(2), vel_radial(2), acc_radial_drag(2)
        real(wp) :: v2, two_ener, mean_movement
        real(wp) :: aux_J2K, aux_inv_dr3_boulder_z, aux_real
        real(wp) :: drag_f, stokes_f
        real(wp) :: rcoll, rescape
        integer(kind=4) :: i, j, jdx, vdx
        integer(kind=4), parameter :: first_particle = 2

        theta = y(1)
        omega = y(2)
        Gmast = G*m_arr(1)

        if (use_stokes) stokes_f = uno2*(uno + tanh(1.e1_wp*(uno - t/stokes_time)))
        if (use_drag)   drag_f = uno2*(uno + tanh(1.e1_wp*(uno - t/drag_time)))
        if (sim%max_distance <= cero) then
            rescape = infinito
        else
            rescape = sim%max_distance
        end if

        ! ====================================================================
        ! ── Asteroid COM → particles  (PARALLEL) ────────────────────────────
        ! ====================================================================
        !$OMP PARALLEL DO                              &
        !$OMP   DEFAULT(SHARED)                        &
        !$OMP   PRIVATE(j, jdx, vdx,                   &
        !$OMP           coords_P, dr_vec, dr, dr2,     &
        !$OMP           rcoll, inv_dr, inv_dr2,        &
        !$OMP           inv_dr3, inv_dr7,              &
        !$OMP           acc_grav, Q_eff,               &
        !$OMP           aux_inv_dr3_boulder_z,         &
        !$OMP           aux_J2K, aux_real,             &
        !$OMP           dr_ver, dv_vec,                &
        !$OMP           vel_circ, vel_radial,          &
        !$OMP           acc_radial_drag,               &
        !$OMP           v2, two_ener, mean_movement)   &
        !$OMP   SCHEDULE(DYNAMIC, 256)
        do j = first_particle, N_total
            jdx = get_index(j)
            coords_P = y(jdx:jdx + 3)

            dr_vec = coords_P(1:2)
            dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)

            rcoll = min(sim%min_distance, R_arr(1) + R_arr(j))
            if ((dr2 < rcoll*rcoll) .or. (dr2 > rescape*rescape)) then
                !$OMP ATOMIC WRITE
                hard_exit = .True.
                cycle
            end if

            dr = sqrt(dr2)
            inv_dr3 = uno/(dr2*dr)
            acc_grav = cero

            if (use_ellipsoid) then
                if ((dr_vec(1)/asteroid_data(1))**2 + (dr_vec(2)/asteroid_data(2))**2 < uno) then
                    !$OMP ATOMIC WRITE
                    hard_exit = .True.
                    cycle
                end if
                inv_dr2 = inv_dr3*dr
                Q_eff = 5*( dr_vec(1)**2 - dr_vec(2)**2 )*inv_dr2*inv_dr2
                acc_grav(1) = -(Gmast*inv_dr3)*dr_vec(1)*( uno &
                              &  - K_coef*inv_dr2 &
                              &  - L_coef*(dos*inv_dr2 - Q_eff) )
                acc_grav(2) = -(Gmast*inv_dr3)*dr_vec(2)*( uno &
                              &  - K_coef*inv_dr2 &
                              &  - L_coef*(-dos*inv_dr2 - Q_eff) )

            else if (use_manual_J2_from_cm) then
                acc_grav = Gmast*dr_vec*J2K_coef/dr2*inv_dr3
                if (sim%megno_active) then
                    vdx = get_variational_index(j, first_particle, N_total)
                    coords_P = y(vdx:vdx + 3)
                    inv_dr7 = inv_dr3 * inv_dr3 / dr
                    aux_real = Gmast*J2K_coef*inv_dr7
                    der(vdx + 2) = der(vdx + 2) + aux_real*( &
                            & -5*coords_P(2)*dr_vec(1)*dr_vec(2) &
                            & + coords_P(1)*(-5*dr_vec(1)*dr_vec(1) + dr2))
                    der(vdx + 3) = der(vdx + 3) + aux_real*( &
                            & -5*coords_P(1)*dr_vec(1)*dr_vec(2) &
                            & + coords_P(2)*(dr2 - 5*dr_vec(2)*dr_vec(2)))
                end if

            else if (use_boulder_z) then
                aux_inv_dr3_boulder_z = uno/(dr2 + dz2_boulder_z_coef)**(1.5e0_wp)
                acc_grav = -Gmboulder_z_coef*dr_vec*aux_inv_dr3_boulder_z
                if (sim%megno_active) then
                    vdx = get_variational_index(j, first_particle, N_total)
                    coords_P = y(vdx:vdx + 3)
                    aux_real = Gmboulder_z_coef/(dr2 + dz2_boulder_z_coef)**(2.5e0_wp)
                    der(vdx + 2) = der(vdx + 2) + aux_real*( &
                            &  3*coords_P(2)*dr_vec(1)*dr_vec(2) &
                            &  - coords_P(1)*(-3*dr_vec(1)*dr_vec(1) + dr2 + dz2_boulder_z_coef))
                    der(vdx + 3) = der(vdx + 3) + aux_real*( &
                            &  3*coords_P(1)*dr_vec(1)*dr_vec(2) &
                            &  - coords_P(2)*(dr2 - 3*dr_vec(2)*dr_vec(2) + dz2_boulder_z_coef))
                end if
            end if

            der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav

            ! Coriolis and centrifugal
            der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) &
                                  & + dos*y(2)*(/ y(jdx + 3), -y(jdx + 2) /) &
                                  & + y(2)**2 * y(jdx:jdx + 1)

            if (sim%megno_active) then
                vdx = get_variational_index(j, first_particle, N_total)
                coords_P = y(vdx:vdx + 3)
                der(vdx + 2) = der(vdx + 2) + dos*omega*coords_P(4) + omega*omega*coords_P(1)
                der(vdx + 3) = der(vdx + 3) - dos*omega*coords_P(3) + omega*omega*coords_P(2)
            end if

            if (use_drag .or. use_stokes) then
                inv_dr = inv_dr3*dr2
                dr_ver = dr_vec*inv_dr
                dv_vec = coords_P(3:4)
                v2 = dot_product(dv_vec, dv_vec)
                if (use_manual_J2_from_cm) then
                    aux_J2K = J2K_coef/dr2
                    two_ener = dos*Gmast*inv_dr*(uno - aux_J2K) - v2
                    mean_movement = sqrt(Gmast*inv_dr3)*(uno - aux_J2K*uno3)
                else
                    two_ener = dos*Gmast*inv_dr - v2
                    mean_movement = abs(two_ener)**(1.5e0_wp)/Gmast
                end if
                if (two_ener > cero) then
                    if (use_drag) then
                        vel_radial = dot_product(dr_ver, dv_vec)
                        acc_radial_drag = -drag_coef*mean_movement*vel_radial
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) &
                                              & + acc_radial_drag*dr_ver*drag_f
                    end if
                    if (use_stokes) then
                        vel_circ = mean_movement*(/-dr_vec(2), dr_vec(1)/)
                        der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) &
                                              & - stokes_C*(dv_vec - stokes_alpha*vel_circ)*stokes_f
                    end if
                end if
            end if

        end do
        !$OMP END PARALLEL DO

        ! ====================================================================
        ! ── Boulder gravity  (only when NOT triaxial) ───────────────────────
        ! ====================================================================
        if (.not. use_ellipsoid) then

            Gmi = G*boulders_data(0, 1)

            !$OMP PARALLEL DO                             &
            !$OMP   DEFAULT(SHARED)                       &
            !$OMP   PRIVATE(j, jdx, vdx,                  &
            !$OMP           coords_P, dr_vec, dr, dr2,    &
            !$OMP           rcoll, acc_grav, aux_real)    &
            !$OMP   SCHEDULE(STATIC)
            do j = first_particle, N_total
                jdx = get_index(j)
                coords_P = y(jdx:jdx + 3)
                dr_vec = coords_P(1:2) - boulders_coords(0, 1:2)
                dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                rcoll = boulders_data(0, 2) + R_arr(j)
                if (dr2 < rcoll*rcoll) then
                    !$OMP ATOMIC WRITE
                    hard_exit = .True.
                    cycle
                end if
                dr = sqrt(dr2)
                if (use_manual_J2_from_primary) then
                    acc_grav = -Gmi*dr_vec/(dr2*dr)*(uno - J2K_coef/dr2)
                else
                    acc_grav = -Gmi*dr_vec/(dr2*dr)
                    if (sim%megno_active) then
                        vdx = get_variational_index(j, first_particle, N_total)
                        coords_P = y(vdx:vdx + 3)
                        der(vdx + 2:vdx + 3) = der(vdx + 2:vdx + 3) - Gmi * &
                            & ( coords_P(1:2)*dr2 - 3*dot_product(dr_vec, coords_P(1:2))*dr_vec ) &
                            & / (dr2*dr2*dr)
                    end if
                end if
                der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav
            end do
            !$OMP END PARALLEL DO

            do i = 1, sim%Nboulders
                Gmi = G*boulders_data(i, 1)

                !$OMP PARALLEL DO                             &
                !$OMP   DEFAULT(SHARED)                       &
                !$OMP   PRIVATE(j, jdx, vdx,                  &
                !$OMP           coords_P, dr_vec, dr, dr2,    &
                !$OMP           rcoll, acc_grav, aux_real)    &
                !$OMP   SCHEDULE(STATIC)
                do j = first_particle, N_total
                    jdx = get_index(j)
                    coords_P = y(jdx:jdx + 3)
                    dr_vec = coords_P(1:2) - boulders_coords(i, 1:2)
                    dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                    rcoll = boulders_data(i, 2) + R_arr(j)
                    if (dr2 < rcoll*rcoll) then
                        !$OMP ATOMIC WRITE
                        hard_exit = .True.
                        cycle
                    end if
                    dr = sqrt(dr2)
                    acc_grav = -Gmi*dr_vec/(dr2*dr)
                    der(jdx + 2:jdx + 3) = der(jdx + 2:jdx + 3) + acc_grav
                    if (sim%megno_active) then
                        vdx = get_variational_index(j, first_particle, N_total)
                        coords_P = y(vdx:vdx + 3)
                        der(vdx + 2:vdx + 3) = der(vdx + 2:vdx + 3) - Gmi * &
                            & ( coords_P(1:2)*dr2 - 3*dot_product(dr_vec, coords_P(1:2))*dr_vec ) &
                            & / (dr2*dr2*dr)
                    end if
                end do
                !$OMP END PARALLEL DO

            end do

        end if  ! .not. use_ellipsoid

    end subroutine dydt_grav_sinodic


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                       COLLISION DERIVATIVES                             !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine dydt_coll(t, y, der, first_particle, N_total)
        implicit none
        real(wp),               intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in) :: first_particle, N_total

        integer(kind=4) :: last_moon, N_particles

        last_moon = first_particle - 1
        N_particles = N_total - last_moon

        if (sim%use_moon_soft_sphere_col .and. (last_moon > 3)) then
            if (last_moon <= sim%grid_col_min_bodies + 1) then
                call collisions_brute(y, der, 2, last_moon, .True.)
            else if (sim%use_verlet_col .and. sim%use_verlet_with_moons) then
                call collisions_verlet(t, y, der, 2, last_moon, .True.)
            else
                call collisions_grid(y, der, 2, last_moon, .True.)
            end if
        end if

        if (sim%use_part_soft_sphere_col) then
            if (N_particles <= sim%grid_col_min_bodies) then
                call collisions_brute(y, der, first_particle, N_total, .False.)
            else if (sim%use_verlet_col .and. .not. sim%use_verlet_with_moons) then
                call collisions_verlet(t, y, der, first_particle, N_total, .False.)
            else
                call collisions_grid(y, der, first_particle, N_total, .False.)
            end if
        end if

    end subroutine dydt_coll


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                          MEGNO DERIVATIVES                              !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine dydt_megno(t, y, der, first_particle, N_total)
        implicit none
        real(wp),               intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in) :: first_particle, N_total

        integer(kind=4) :: i, vdx
        real(wp) :: prod, dist, glob_prod, glob_dist

        if (.not. sim%megno_active) return

        glob_prod = cero
        glob_dist = cero

        do i = first_particle, N_total
            vdx = get_variational_index(i, first_particle, N_total)

            der(vdx:vdx + 1) = y(vdx + 2:vdx + 3)

            prod = y(vdx)*der(vdx) + y(vdx+1)*der(vdx+1) &
                 + y(vdx+2)*der(vdx+2) + y(vdx+3)*der(vdx+3)
            if (.not. ieee_is_finite(prod) .or. abs(prod) < tini) prod = cero

            dist = y(vdx)*y(vdx) + y(vdx+1)*y(vdx+1) &
                 + y(vdx+2)*y(vdx+2) + y(vdx+3)*y(vdx+3)
            if (.not. ieee_is_finite(dist) .or. dist < tini) dist = tini

            glob_prod = glob_prod + prod
            glob_dist = glob_dist + dist

            der(vdx + 4) = prod / dist
            der(vdx + 5) = prod / dist * t / megno_factor
            if (t > 0) der(vdx + 6) = dos * y(vdx + 5) / t
        end do

        der(vdx + 7) = glob_prod / glob_dist
        der(vdx + 8) = glob_prod / glob_dist * t / megno_factor
        if (t > 0) der(vdx + 9) = dos * y(vdx + 8) / t

    end subroutine dydt_megno


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                     TOP-LEVEL DERIVATIVE FUNCTIONS                      !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function dydt(t, y) result(der)
        implicit none
        real(wp),               intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: der

        integer(kind=4) :: last_moon, first_particle, N_total

        der = cero

        last_moon = 1 + sim%Nmoon_active
        first_particle = last_moon + 1
        N_total = last_moon + sim%Npart_active

        call dydt_grav(t, y, der, first_particle, N_total)
        call dydt_coll(t, y, der, first_particle, N_total)
        call set_pos_derivatives(y, der, N_total)
        call dydt_megno(t, y, der, first_particle, N_total)

    end function dydt

    function dydt_grav_f(t, y) result(der)
        implicit none
        real(wp),               intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: der

        integer(kind=4) :: last_moon, first_particle, N_total

        der = cero

        last_moon = 1 + sim%Nmoon_active
        first_particle = last_moon + 1
        N_total = last_moon + sim%Npart_active

        call dydt_grav(t, y, der, first_particle, N_total)

    end function dydt_grav_f

    function dydt_coll_f(t, y) result(der)
        implicit none
        real(wp),               intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: der

        integer(kind=4) :: last_moon, first_particle, N_total

        der = cero

        last_moon = 1 + sim%Nmoon_active
        first_particle = last_moon + 1
        N_total = last_moon + sim%Npart_active

        call dydt_coll(t, y, der, first_particle, N_total)

    end function dydt_coll_f

end module derivates