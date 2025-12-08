!> Module with acceleration initializers, and some usage examples.

module accelerations
    use auxiliary, only: cross2D_z
    use constants, only: wp, cero, uno, uno2, uno3, dos, G, tini

    implicit none

    logical :: use_stokes = .False., use_stokes_moons = .False.
    real(wp) :: stokes_time = cero, stokes_C = cero, stokes_alpha = cero  ! Stokes
    logical :: use_drag = .False., use_drag_moons = .False.
    real(wp) :: drag_coef = cero, drag_time = cero  ! Drag
    logical :: use_ellipsoid = .False.
    real(wp) :: C20_coef = cero, C22_coef = cero, Re_coef = cero ! Ellipsoid basics
    real(wp) :: K_coef = cero, L_coef = cero ! Ellipsoid deep
    logical :: use_manual_J2 = .False.
    real(wp) :: J2K_coef = cero ! Manual J2 basics
    logical :: use_damp = .False.
    real(wp) :: damp_coef_1 = cero, damp_coef_2 = cero, damp_time = cero ! Omega Damping
    integer(kind=4) :: damp_model = -1 ! Omega Damping

contains

    ! xy_vec: Position vector [x, y]
    ! dr_vec: Relative position vector between B and A [xy_vecA - xy_vecB]. From OTHER to ME
    ! dv_vec: Relative velocity vector between B and A [vxvy_vecA - vxvy_vecB]
    ! dr: |dr_vec|
    ! total_mass: mass A + mass B

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!! GRAVITY !!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    pure subroutine gravity(mass, dr_vec, dr, acc)
        implicit none
        real(wp), intent(in) :: mass, dr_vec(2), dr
        real(wp), intent(inout) :: acc(2)

        acc = acc - G*mass*dr_vec/(dr*dr*dr)  ! G m (x, y) / r
    end subroutine gravity

    pure subroutine torque_grav_Z(xy_vec, acc_applied, torque_m) ! torque / mass
        implicit none
        real(wp), intent(in) :: xy_vec(2), acc_applied(2)
        real(wp), intent(inout) :: torque_m

        torque_m = torque_m + cross2D_z(xy_vec, acc_applied) ! r x acc  ! mass to be added later
    end subroutine torque_grav_Z

    ! Both
    pure subroutine gravity_and_torque_Z(mass, xy_vec, dr_vec, dr, torque_m, acc)
        implicit none
        real(wp), intent(in) :: mass, xy_vec(2), dr_vec(2), dr
        real(wp), intent(inout) :: torque_m, acc(2)
        real(wp) :: acc_grav(2)

        acc_grav = cero  ! Init
        call gravity(mass, dr_vec, dr, acc_grav)
        call torque_grav_Z(xy_vec, acc_grav, torque_m)

        acc = acc + acc_grav
    end subroutine gravity_and_torque_Z

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!! TRIAXIAL !!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!! Init Parameters
    subroutine init_ellipsoid(axis_a, axis_b, axis_c)
        implicit none
        real(wp), intent(in) :: axis_a, axis_b, axis_c

        if ((axis_a > cero) .and. (axis_c > cero)) then
            use_ellipsoid = .True.
            Re_coef = (axis_a*axis_b*axis_c)**uno3
            C20_coef = (dos*axis_c**2 - axis_a**2 - axis_b**2)/(10.e0_wp*Re_coef**2)
            C22_coef = (axis_a**2 - axis_b**2)/(20.e0_wp*Re_coef**2)
            K_coef = 1.5e0_wp*Re_coef**2*C20_coef
            L_coef = 3.e0_wp*Re_coef**2*C22_coef
        else
            use_ellipsoid = .False.
        end if
    end subroutine init_ellipsoid

    !!! Acceleration Example
    subroutine ellipsoid_acceleration(mass, theta, dr_vec, dr, acc)
        implicit none
        real(wp), intent(in) :: mass, theta, dr_vec(2), dr
        real(wp), intent(inout) :: acc(2)
        real(wp) :: c2th, s2th
        real(wp) :: Q_param, dQdx, dQdy
        real(wp) :: Q_param_eff
        real(wp) :: inv_dr2, inv_dr3

        c2th = cos(dos*theta)
        s2th = sin(dos*theta)

        Q_param = (dr_vec(1)**2 - dr_vec(2)**2)*c2th + dos*dr_vec(1)*dr_vec(2)*s2th ! (x²-y²) cos(2th) + 2xy sin(2th)
        dQdx = dos*(dr_vec(1)*c2th + dr_vec(2)*s2th)  ! 2x cos(2th) + 2y sin(2th)
        dQdy = -dos*(dr_vec(2)*c2th - dr_vec(1)*s2th)  ! - 2y cos(2th) + 2x sin(2th)

        inv_dr3 = uno/(dr*dr*dr)
        inv_dr2 = dr*inv_dr3

        ! Q_param = (x²-y²) cos(2th) + 2xy sin(2th)
        ! Q_param_eff = 5 * Q_param / r⁴
        Q_param_eff = 5.e0_wp*Q_param*inv_dr2*inv_dr2 ! Q_ef = 5 * Q / r⁴

        ! a_unit_massx = G / r³ (x - K x / r² - L (dQ/dx / r² - x 5 Q / r⁴))
        ! a_unit_massy = G / r³ (y - K y / r² - L (dQ/dy / r² - y 5 Q / r⁴))
        acc(1) = acc(1) - (G*mass*inv_dr3)* &
                    & (dr_vec(1) &
                    &  - K_coef*dr_vec(1)*inv_dr2 &
                    &  - L_coef*(dQdx*inv_dr2 - dr_vec(1)*Q_param_eff))
        acc(2) = acc(2) - (G*mass*inv_dr3)* &
                    & (dr_vec(2) &
                    &  - K_coef*dr_vec(2)*inv_dr2 &
                    &  - L_coef*(dQdy*inv_dr2 - dr_vec(2)*Q_param_eff))
    end subroutine ellipsoid_acceleration

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! MANUAL J2 !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!! Init Parameters
    subroutine init_manual_J2(J2, radius)
        implicit none
        real(wp), intent(in) :: J2, radius

        if ((J2 > cero) .and. (radius > cero)) then
            use_manual_J2 = .True.
            J2K_coef = -1.5e0_wp*radius**2*J2  ! Minus, bc C20 = -J2
        else
            use_manual_J2 = .False.
        end if
    end subroutine init_manual_J2

    !!! Acceleration Example
    subroutine manual_J2_acceleration(mass, dr_vec, dr, acc)
        implicit none
        real(wp), intent(in) :: mass, dr_vec(2), dr
        real(wp), intent(inout) :: acc(2)
        real(wp) :: inv_dr2, inv_dr3

        inv_dr3 = uno/(dr*dr*dr)
        inv_dr2 = dr*inv_dr3

        ! a_unit_massx = G / r³ (x - K x / r²)
        ! a_unit_massy = G / r³ (y - K y / r²)
        acc(1) = acc(1) - (G*mass*inv_dr3)*dr_vec(1)*(uno - J2K_coef*inv_dr2)
        acc(2) = acc(2) - (G*mass*inv_dr3)*dr_vec(2)*(uno - J2K_coef*inv_dr2)
    end subroutine manual_J2_acceleration

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!! STOKES !!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!! Init Parameters
    subroutine init_stokes(tau_a, tau_e, active_timescale, apply_to_moons)
        implicit none
        real(wp), intent(in) :: tau_a, tau_e, active_timescale
        logical, intent(in) :: apply_to_moons

        if (active_timescale > cero .and. abs(tau_a) > cero .and. abs(tau_e) > cero) then
            use_stokes = .True.
            stokes_time = active_timescale
            stokes_C = uno/(dos*tau_a) + uno/tau_e
            stokes_alpha = (dos*tau_a)/((dos*tau_a) + tau_e)
        else
            use_stokes = .False.
        end if
        use_stokes_moons = use_stokes .and. apply_to_moons
    end subroutine init_stokes

    !!! Acceleration Example
    subroutine stokes(time, total_mass, dr_vec, dv_vec, dr, acc)
        implicit none
        real(wp), intent(in) :: total_mass, time, dr_vec(2), dv_vec(2), dr
        real(wp), intent(inout) :: acc(2)
        real(wp) :: stokes_factor, vel_circ(2)

        stokes_factor = uno2*(uno + tanh(1.e1_wp*(uno - time/stokes_time)))
        vel_circ = sqrt(G*total_mass/(dr*dr*dr))*(/-dr_vec(2), dr_vec(1)/)  ! v_circ = n (-y, x)

        acc = acc - stokes_C*(dv_vec - stokes_alpha*vel_circ)*stokes_factor ! = -C * (v - alpha * vc) * factor
    end subroutine stokes

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!! NAIVE-STOKES (DRAG) !!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!! Init Parameters
    subroutine init_drag(drag_coefficient, active_timescale, apply_to_moons)
        implicit none
        real(wp), intent(in) :: drag_coefficient, active_timescale
        logical, intent(in) :: apply_to_moons

        if (active_timescale > cero .and. abs(drag_coefficient) > cero) then
            use_drag = .True.
            drag_coef = drag_coefficient
            drag_time = active_timescale
        else
            use_drag = .False.
        end if
        use_drag_moons = use_drag .and. apply_to_moons
    end subroutine init_drag

    !!! Acceleration Example
    subroutine drag(time, total_mass, dr_vec, dv_vec, dr, acc)
        implicit none
        real(wp), intent(in) :: time, total_mass, dr_vec(2), dv_vec(2), dr
        real(wp), intent(inout) :: acc(2)
        real(wp) :: drag_factor, aux_real
        real(wp) :: Gmcm, v2, acc_radial, vel_radial, mean_movement

        drag_factor = uno2*(uno + tanh(1.e1_wp*(uno - time/drag_time)))
        Gmcm = G*total_mass
        v2 = dot_product(dv_vec, dv_vec)

        ! Debemos chequear que la partícula no esté "desligada"
        aux_real = dos*Gmcm/dr - v2
        if (aux_real < cero) return ! No se puede calcular

        mean_movement = aux_real**(1.5e0_wp)/Gmcm ! n
        vel_radial = dot_product(dr_vec, dv_vec)/dr
        acc_radial = -drag_coef*mean_movement*vel_radial

        acc = acc + acc_radial*dr_vec/dr*drag_factor ! = -a_r * (x, y) / r * factor
    end subroutine drag

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! ROTATION DAMPING !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!! Init Parameters
    subroutine init_damping(coefficient_1, coefficient_2, active_timescale, model)
        implicit none
        real(wp), intent(in) :: coefficient_1, coefficient_2, active_timescale
        integer(kind=4), intent(in) :: model

        if (active_timescale > cero .and. abs(coefficient_1) > cero) then
            use_damp = .True.
            damp_coef_1 = coefficient_1
            damp_coef_2 = coefficient_2
            damp_time = active_timescale
            damp_model = model
        else
            use_damp = .False.
        end if
    end subroutine init_damping

    !!! Acceleration Example (These are like the torque, but not mass needed)

    !!!! omega(t) = omega0 + tau * t
    subroutine damping_linear(time, acc_omega)
        implicit none
        real(wp), intent(in) :: time
        real(wp), intent(inout) :: acc_omega
        real(wp) :: damp_factor

        damp_factor = uno2*(uno + tanh(1.e1_wp*(uno - time/damp_time)))
        !! domega/dt = tau
        acc_omega = acc_omega + damp_coef_1*damp_factor
    end subroutine damping_linear

    !!!! omega(t) = omega0 * exp (- (t - t0) / tau) = omega0 * exp (- (t - t0) / tau)
    subroutine damping_exp(time, omega, acc_omega)
        implicit none
        real(wp), intent(in) :: time, omega
        real(wp), intent(inout) :: acc_omega
        real(wp) :: damp_factor

        damp_factor = uno2*(uno + tanh(1.e1_wp*(uno - time/damp_time)))
        !! domega/dt = -exp(- (t-t0) / tau) * omega0 / tau = - omega / tau
        acc_omega = acc_omega - omega/damp_coef_1*damp_factor
    end subroutine damping_exp

    !!!! omega(t) = omega0 * exp (A * (t-t0)**B)
    subroutine damping_expoly(time, omega, acc_omega, initial_time)
        implicit none
        real(wp), intent(in) :: time, omega
        real(wp), intent(inout) :: acc_omega
        real(wp), intent(in), optional :: initial_time
        real(wp) :: damp_factor, t0 = cero

        if (present(initial_time)) t0 = initial_time
        damp_factor = uno2*(uno + tanh(1.e1_wp*(uno - time/damp_time)))

        !! domega/dt = A * B * (t-t0)**(B-1) * omega0 * exp (A * (t-t0)**B) = (A * B * (t-t0)**(B-1)) * omega
        acc_omega = acc_omega + damp_coef_1*(time - t0 + tini)**(damp_coef_2 - uno)*omega*damp_factor
    end subroutine damping_expoly

end module accelerations
