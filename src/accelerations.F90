module accelerations
    use auxiliary, only: cross2D_z
    use constants, only: cero, uno, uno2, dos, G, tini
    implicit none
    real(kind=8) :: stokes_time = cero, stokes_C = cero, stokes_alpha = cero  ! Stokes
    real(kind=8) :: drag_coef = cero, drag_time = cero  ! Drag
    real(kind=8) :: J2_coef = cero ! J2
    real(kind=8) :: damp_coef_1 = cero, damp_coef_2 = cero, damp_time = cero ! Omega Damping
    integer(kind=8) :: damp_model = -1 ! Omega Damping

    contains

        ! xy_vec: Position vector [x, y]
        ! dr_vec: Relative position vector between B and A [xy_vecA - xy_vecB]. From OTHER to ME
        ! dv_vec: Relative velocity vector between B and A [vxvy_vecA - vxvy_vecB]
        ! dr: |dr_vec|
        ! total_mass: mass A + mass B

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!! GRAVITY !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine gravity(mass, dr_vec, dr, acc)
            implicit none
            real(kind=8), intent(in) :: mass, dr_vec(2), dr
            real(kind=8), intent(inout) :: acc(2)

            acc = acc - G * mass * dr_vec / (dr * dr * dr)  ! G m (x, y) / r
        end subroutine gravity

        subroutine torque_grav_Z(xy_vec, acc_applied, torque_m) ! torque / mass
            implicit none
            real(kind=8), intent(in) :: xy_vec(2), acc_applied(2)
            real(kind=8), intent(inout) :: torque_m

            torque_m = torque_m + cross2D_z(xy_vec, acc_applied) ! r x acc  ! mass to be added later
        end subroutine torque_grav_Z

        ! Both
        subroutine gravity_and_torque_Z(mass, xy_vec, dr_vec, dr, torque_m, acc)
            implicit none
            real(kind=8), intent(in) :: mass, xy_vec(2), dr_vec(2), dr
            real(kind=8), intent(inout) :: torque_m, acc(2)
            real(kind=8) :: acc_grav(2) = cero
            
            call gravity(mass, dr_vec, dr, acc_grav)
            call torque_grav_Z(xy_vec, acc_grav, torque_m)

            acc = acc + acc_grav
        end subroutine gravity_and_torque_Z
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!! STOKES !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Init Parameters
        subroutine init_stokes(tau_a, tau_e, charac_timescale)
            implicit none
            real(kind=8), intent(in) :: tau_a, tau_e, charac_timescale
            
            stokes_time = charac_timescale
            stokes_C = uno / (dos * tau_a) + uno / tau_e
            stokes_alpha = (dos * tau_a) / ((dos * tau_a) + tau_e)
        end subroutine init_stokes

        !!! Acceleration
        subroutine stokes(time, total_mass, dr_vec, dv_vec, dr, acc)
            implicit none
            real(kind=8), intent(in) :: total_mass, time, dr_vec(2), dv_vec(2), dr
            real(kind=8), intent(inout) :: acc(2)
            real(kind=8) :: stokes_factor, vel_circ(2)

            stokes_factor = uno2 * (uno + tanh(1.d1 * (uno - time / stokes_time)))
            vel_circ  = sqrt(G * total_mass / (dr * dr * dr)) * (/-dr_vec(2), dr_vec(1)/)  ! v_circ = n (-y, x)

            acc = acc - stokes_C * (dv_vec - stokes_alpha * vel_circ) * stokes_factor ! = -C * (v - alpha * vc) * factor
        end subroutine stokes

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!! NAIVE-STOKES (DRAG) !!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Init Parameters
        subroutine init_drag(drag_coefficient, characteristic_timescale)
            implicit none
            real(kind=8), intent(in) :: drag_coefficient, characteristic_timescale
            
            drag_coef = drag_coefficient
            drag_time = characteristic_timescale
        end subroutine init_drag


        !!! Acceleration
        subroutine drag(time, total_mass, dr_vec, dv_vec, dr, acc)
            implicit none
            real(kind=8), intent(in) :: time, total_mass, dr_vec(2), dv_vec(2), dr 
            real(kind=8), intent(inout) :: acc(2)
            real(kind=8) :: drag_factor, aux_real
            real(kind=8) :: Gmcm, v2, acc_radial, vel_radial, mean_movement 

            drag_factor = uno2 * (uno + tanh(1.d1 * (uno - time / drag_time)))
            Gmcm = G * total_mass
            v2 = dot_product(dv_vec, dv_vec)

            ! Debemos chequear que la partícula no esté "desligada"
            aux_real = dos * Gmcm / dr - v2
            if (aux_real < cero) return ! No se puede calcular

            mean_movement = aux_real**(1.5d0) / Gmcm ! n
            vel_radial = dot_product(dr_vec, dv_vec) / dr 
            acc_radial = - drag_coef * mean_movement * vel_radial
            
            acc = acc + acc_radial * dr_vec / dr * drag_factor ! = -a_r * (x, y) / r * factor
        end subroutine drag

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! GEO-POTENTIAL !!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Init Parameters
        subroutine init_J2(J2_coefficient)
            implicit none
            real(kind=8), intent(in) :: J2_coefficient
            
            J2_coef = J2_coefficient * 1.5d0
        end subroutine init_J2

        !!! Acceleration
        subroutine J2_acceleration(mass, dr_vec, dr, acc)
            implicit none
            real(kind=8), intent(in) :: mass, dr_vec(2), dr
            real(kind=8), intent(inout) :: acc(2)
            
            acc = acc - G * mass * J2_coef * dr_vec / (dr * dr * dr * dr * dr)
        end subroutine J2_acceleration

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!! ROTATION DAMPING !!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Init Parameters
        subroutine init_damping(coefficient_1, coefficient_2, characteristic_timescale, model)
            implicit none
            real(kind=8), intent(in) :: coefficient_1, coefficient_2, characteristic_timescale
            integer(kind=4), intent(in) :: model
            
            damp_coef_1 = coefficient_1
            damp_coef_2 = coefficient_2
            damp_time = characteristic_timescale
            damp_model = model
        end subroutine init_damping

        !!! Acceleration  (These are like the torque, but not mass needed)

        !!!! omega(t) = omega0 + tau * t
        subroutine damping_linear(time, acc_omega)
            implicit none
            real(kind=8), intent(in) :: time
            real(kind=8), intent(inout) :: acc_omega
            real(kind=8) :: damp_factor

            damp_factor = uno2 * (uno + tanh(1.d1 * (uno - time / damp_time)))
            acc_omega = acc_omega + damp_coef_1 * damp_factor ! domega/dt = tau
        end subroutine damping_linear

        !!!! omega(t) = omega0 * exp (- (t - t0) / tau) = omega0 * exp (- (t - t0) / tau)
        subroutine damping_exp(time, omega, acc_omega)
            implicit none
            real(kind=8), intent(in) :: time, omega
            real(kind=8), intent(inout) :: acc_omega
            real(kind=8) :: damp_factor

            damp_factor = uno2 * (uno + tanh(1.d1 * (uno - time / damp_time)))
            acc_omega = acc_omega - omega / damp_coef_1 * damp_factor ! domega/dt = -exp(- (t-t0) / tau) * omega0 / tau = - omega / tau
        end subroutine damping_exp

        !!!! omega(t) = omega0 * exp (A * (t-t0)**B)
        subroutine damping_expoly(time, omega, acc_omega, initial_time)
            implicit none
            real(kind=8), intent(in) :: time, omega
            real(kind=8), intent(inout) :: acc_omega
            real(kind=8), intent(in), optional :: initial_time
            real(kind=8) :: damp_factor, t0 = cero

            if (present(initial_time)) t0 = initial_time
            damp_factor = uno2 * (uno + tanh(1.d1 * (uno - time / damp_time)))

            ! domega/dt = A * B * (t-t0)**(B-1) * omega0 * exp (A * (t-t0)**B) = (A * B * (t-t0)**(B-1)) * omega
            acc_omega = acc_omega + damp_coef_1 * (time - t0 + tini)**(damp_coef_2 - uno) * omega * damp_factor
        end subroutine damping_expoly
    


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!! ALL !!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! Here we try to repeat calculations as minimum as possible

        subroutine acc_all(time, mass_other, mass_me, xy_vec, dr_vec, dv_vec, dr, omega, torque, acc_omega, acc)
            implicit none
            real(kind=8), intent(in) :: time
            real(kind=8), intent(in) :: mass_other, mass_me
            real(kind=8), intent(in) :: xy_vec(2), dr_vec(2), dv_vec(2), dr
            real(kind=8), intent(in) :: omega
            real(kind=8), intent(inout) :: torque, acc_omega, acc(2)
            real(kind=8) :: Gmtot, Gmass, dr_ver(2), dr3, factor
            real(kind=8) :: acc_grav(2)
            real(kind=8) :: vel_circ(2), v2
            real(kind=8) :: aux_real
            real(kind=8) :: mean_movement
            real(kind=8) :: vel_radial(2), acc_radial(2)

            Gmtot = G * (mass_me + mass_other)
            Gmass = G * mass_other
            dr_ver = dr_vec / dr
            dr3 = dr * dr * dr

            !! Gravity
            acc_grav = - Gmass * dr_vec / dr3  ! G m (x, y) / r
            acc = acc + acc_grav

            !! Torque
            torque = torque + cross2D_z(xy_vec, acc_grav * mass_me)! r x F  ! mass to be added later

            !! Stokes
            factor = uno2 * (uno + tanh(1.d1 * (uno - time / stokes_time)))
            vel_circ  = sqrt(Gmtot / dr3) * (/-dr_vec(2), dr_vec(1)/)  ! v_circ = n (-y, x)

            acc = acc - stokes_C * (dv_vec - stokes_alpha * vel_circ) * factor ! = -C * (v - alpha * vc) * factor
            
            !! Naive-Stokes (Drag)
            factor = uno2 * (uno + tanh(1.d1 * (uno - time / drag_time)))
            v2 = dot_product(dv_vec, dv_vec)
            
            aux_real = dos * Gmtot / dr - v2  ! Check if unbound
            if (aux_real > cero) then ! Can calculate only in this case
                mean_movement = aux_real**(1.5d0) / Gmtot ! n
                vel_radial = dot_product(dr_ver, dv_vec)
                acc_radial = - drag_coef * mean_movement * vel_radial                
            
                acc = acc + acc_radial * dr_ver * factor ! = -a_r * (x, y) / r * factor
            end if

            !! J2
            acc = acc - G * mass_other * J2_coef * dr_ver / (dr3 * dr)

            !! Damping
            factor = uno2 * (uno + tanh(1.d1 * (uno - time / damp_time)))
            select case (damp_model)
                case (1) ! domega/dt = tau
                    acc_omega = acc_omega + damp_coef_1 * factor 
                case (2) ! domega/dt = -exp(- (t-t0) / tau) * omega0 / tau = - omega / tau     
                    acc_omega = acc_omega - omega / damp_coef_1 * factor
                case (3) ! domega/dt = A * B * (t-t0)**(B-1) * omega0 * exp (A * (t-t0)**B) = (A * B * (t-t0)**(B-1)) * omega
                    acc_omega = acc_omega + &
                              & damp_coef_1 * (time - cero + tini)**(damp_coef_2 - uno) * omega * factor
            end select

        end subroutine acc_all

end module accelerations
