module forces
    use parameters
    implicit none
    
    procedure (dydt_template), pointer :: dydt => null ()
    procedure (dydt_single_template), pointer :: domegadt => null ()
    procedure (dydt_single_template), pointer :: dmassdt => null ()
    
    abstract interface
        ! Here must be every f_i defined explicitly
        function dydt_template (t, y) result (der) 
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size (y))      :: der
        end function dydt_template

        function dydt_single_template (t, y) result (der) 
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), intent(in) :: y
            real(kind=8) :: der
        end function dydt_single_template
    end interface
    
    contains


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!! DERIVATIVES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function dydt_single_null (t, variable) result(dydt)
            ! No hay derivadas
            implicit none
            real(kind=8), intent(in) :: t, variable
            real(kind=8) :: dydt
            
            dydt = cero
        end function dydt_single_null
    
        function dydt_null (t, y) result(dydt)
            ! No hay derivadas
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y)) :: dydt
            
            dydt = cero
        end function

        !! domega/dt

        function domega_dt_linear (t, omega) result(domegadt)
            ! omega(t) = omega0 + omega_linear_damping_slope * t
            ! domega/dt = omega_linear_damping_slope
            implicit none
            real(kind=8), intent(in) :: t, omega
            real(kind=8) :: domegadt
            
            domegadt = omega_linear_damping_slope
        end function domega_dt_linear
    
        function domega_dt_exponential (t, omega) result(domegadt)
            ! omega(t) = omega0 * exp (-(t-t0)/omega_exp_damping_time)
            ! domega/dt = -exp(- (t-t0) / omega_exp_damping_time) * omega0 / omega_exp_damping_time = 
            !           = - omega / omega_exp_damping_time
            implicit none
            real(kind=8), intent(in) :: t, omega
            real(kind=8) :: domegadt
            
            ! domegadt = -exp (- (t - initial_time) / omega_exp_damping_time) * &
            !             & asteroid_initial_conditions(8) / omega_exp_damping_time
            domegadt = - omega / omega_exp_damping_time
        end function domega_dt_exponential

        function domega_dt_expoly (t, omega) result(domegadt)
            ! omega(t) = omega0 * exp (A * (t-t0)**B)
            ! domega/dt = A * B * (t-t0)**(B-1) * omega0 * exp (A * (t-t0)**B) = (A * B * (t-t0)**(B-1)) * omega
            implicit none
            real(kind=8), intent(in) :: t, omega
            real(kind=8) :: domegadt
            
            ! domegadt = omega_exp_poly_A * omega_exp_poly_B * (t - initial_time + tini)**(omega_exp_poly_B - uno) * &
            !            & exp(omega_exp_poly_A * (t - initial_time + tini)**omega_exp_poly_B) * asteroid_initial_conditions(8)
            domegadt = omega_exp_poly_AB * (t - initial_time + tini)**(omega_exp_poly_B - uno) * omega
        end function domega_dt_expoly

        !! dmass/dt

        function dmass_dt_exponential(t, mass) result(dmassdt)
            ! Exponential mass damping
            implicit none
            real(kind=8), intent(in) :: t, mass
            real(kind=8) :: dmassdt
            
            dmassdt = -exp (- (t - initial_time) / mass_exp_damping_time) * mass / mass_exp_damping_time
        end function dmass_dt_exponential

        !! dY/dt

        function dydt_explicit_v1 (t, y) result(dydt)
            !y = / x0, y0, vx0, vy0, Bould, Part, .../
            ! Utiliza las posiciones exactas para los boulders, con w constante
            !! En caso de modificar la velocidad angular por fuera de la función,
            !! se debe modificar asteroid_theta_correction para que las posiciones sean correctas
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: rib(0:Nboulders,2), vib(0:Nboulders,2)
            real(kind=8) :: rb(2), vb(2), ab(2)
            integer(kind=4) :: i, ineqs, particle_i
            real(kind=8) :: dummy_real, dummy_real2(2)

            ! Initialize
            dydt = cero            

            ! Explicit positions, velocities and accelerations for boulders
            do i = 0, Nboulders
                ineqs = i * equation_size
                rib(i,1) = cos(asteroid_omega * (t - initial_time) + theta_ast_arr(i) + asteroid_theta_correction) * dist_ast_arr(i)
                rib(i,2) = sin(asteroid_omega * (t - initial_time) + theta_ast_arr(i) + asteroid_theta_correction) * dist_ast_arr(i)
                vib(i,1) = -asteroid_omega * rib(i,2)
                vib(i,2) = asteroid_omega * rib(i,1)
                dydt(ineqs+1) = vib(i,1) ! for the sake of the integrator
                dydt(ineqs+2) = vib(i,2) ! for the sake of the integrator
                dydt(ineqs+3) = -asteroid_omega2 * rib(i,1) ! for the sake of the integrator
                dydt(ineqs+4) = -asteroid_omega2 * rib(i,2) ! for the sake of the integrator
            end do

            ! Particles accelerations
            !$OMP PARALLEL DEFAULT(SHARED) &
            !$OMP PRIVATE(i,rb,vb,ab,particle_i)
            !$OMP DO SCHEDULE (STATIC)
            do i = 1, Nactive
                ab = cero
                particle_i = (i + Nboulders) * equation_size
                rb(1) = y(particle_i+1)
                rb(2) = y(particle_i+2)
                vb(1) = y(particle_i+3)
                vb(2) = y(particle_i+4)
                call apply_force(&
                    & t, &
                    & mass_ast_arr, particles_mass(i), &
                    & particles_hexit(i), &
                    & rb, vb, &
                    & rib, vib, &
                    & asteroid_mass, asteroid_pos, asteroid_vel, & ! No varían en cada step en este modelo
                    & asteroid_inertia, & ! inertia no varía en cada step en este modelo
                    & dummy_real, dummy_real2, ab) ! Solo aceleración a partículas
                dydt(particle_i+1) = vb(1)
                dydt(particle_i+2) = vb(2)
                dydt(particle_i+3) = ab(1)
                dydt(particle_i+4) = ab(2)
            end do
            !$OMP END DO
            !$OMP END PARALLEL

        end function dydt_explicit_v1
        
        function dydt_implicit_v1 (t, y) result(dydt) ! Tienen un error sistemático que no detecto aún
            !y = / x0, y0, vx0, vy0, Bould, Part, .../
            !! This is the hardest one to implement, because omega is not constant,
            !! and cant be calculated properly. Also, because it is not in the parameters array
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: rib(0:Nboulders,2), vib(0:Nboulders,2)
            real(kind=8) :: rb(2), vb(2), ab(2)
            real(kind=8) :: rcm(2), vcm(2), raux(0:Nboulders,2), vaux(0:Nboulders,2)
            real(kind=8) :: omega, omega2, omegadot, inertia, angmom
            integer(kind=4) :: i, ineqs, particle_i
            real(kind=8) :: dummy_real2(2), aux_real

            ! Initialize
            dydt = cero

            ! Calculate the center of mass of the asteroid
            rcm = cero
            vcm = cero
            do i = 0, Nboulders 
                ineqs = i * equation_size
                rib(i,1) = y(ineqs+1)
                rib(i,2) = y(ineqs+2)
                vib(i,1) = y(ineqs+3)
                vib(i,2) = y(ineqs+4)
                rcm = rcm + mass_ast_arr(i) * rib(i,:)
                vcm = vcm + mass_ast_arr(i) * vib(i,:)
            end do
            rcm = rcm / asteroid_mass ! rcm = sum_i m_i * r_i / M
            vcm = vcm / asteroid_mass ! vcm = sum_i m_i * v_i / M

            if (sqrt(sum(rcm * rcm)) < error_tolerance) rcm = cero ! Avoid numerical errors
            if (sqrt(sum(vcm * vcm)) < error_tolerance) vcm = cero ! Avoid numerical errors

            ! Calculate angmom and inertia
            aux_real = cero
            angmom = cero
            inertia = cero
            do i = 0, Nboulders
                raux(i,:) = rib(i,:) - rcm
                vaux(i,:) = vib(i,:) - vcm
                aux_real = 0.4d0 * mass_ast_arr(i) * radius_ast_arr(i)**2 ! Inertia Sphere
                angmom = angmom + mass_ast_arr(i) * cross2D(rib(i,:), vib(i,:)) ! Traslacional (+ Rotational)
                ! angmom = angmom + aux_real * cross2D(raux(i,:), vaux(i,:))/sum(raux(i,:) * raux(i,:)) ! Rotacional (Sphere)
                inertia = inertia + inertia_ast_arr(i) + mass_ast_arr(i) * sum(raux(i,:) * raux(i,:)) ! Sphere + Steiner
            end do
            !! Get Omega from L/I
            omega = angmom / inertia
            omega2 = omega * omega

            ! Particles accelerations (and possible torque for the asteroid)
            omegadot = cero ! Init
            !$OMP PARALLEL DEFAULT(SHARED) &
            !$OMP PRIVATE(i,rb,vb,ab,particle_i)
            !$OMP DO REDUCTION(+:omegadot) SCHEDULE (STATIC)
            do i = 1, Nactive
                ab = cero
                particle_i = (i + Nboulders) * equation_size
                rb(1) = y(particle_i+1)
                rb(2) = y(particle_i+2)
                vb(1) = y(particle_i+3)
                vb(2) = y(particle_i+4)
                call apply_force(&
                    & t, &
                    & mass_ast_arr, particles_mass(i), &
                    & particles_hexit(i), &
                    & rb, vb, &
                    & rib, vib, &
                    & asteroid_mass, rcm, vcm, & ! rcm y vcm pueden haber variado, así que se recalcularon
                    & inertia, & ! inertia puede haber variado, así que se recalculó
                    & omegadot, dummy_real2, ab)
                dydt(particle_i+1) = vb(1)
                dydt(particle_i+2) = vb(2)
                dydt(particle_i+3) = ab(1)
                dydt(particle_i+4) = ab(2)
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            ! Update omega
            omegadot = domegadt(t, omega) + omegadot 

            ! Set bodies accelerations (acá al final porque necesitamos omegadot)
            do i = 0, Nboulders 
                ineqs = i * equation_size
                dydt(ineqs+1) = y(ineqs+3)
                dydt(ineqs+2) = y(ineqs+4)
                dydt(ineqs+3) = -omega2 * raux(i,1) - omegadot * raux(i,2) ! a_i = -w^2 * r_i + dw/dt * (-ry,rx)
                dydt(ineqs+4) = -omega2 * raux(i,2) + omegadot * raux(i,1) ! a_i = -w^2 * r_i + dw/dt * (-ry,rx)
            end do

        end function dydt_implicit_v1
        
        function dydt_explicit_v2 (t, y) result(dydt)
            !y = /theta, omega, xA, yA, vxA, vyA, Part, .../
            ! Utiliza las posiciones exactas para los boulders, con w constante
            !! En caso de modificar la velocidad angular por fuera de la función,
            !! se debe modificar asteroid_theta_correction para que las posiciones sean correctas
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: rib(0:Nboulders,2), vib(0:Nboulders,2)
            real(kind=8) :: rb(2), vb(2), ab(2)
            real(kind=8) :: rcm(2), vcm(2)
            real(kind=8) :: theta
            integer(kind=4) :: i, particle_i
            real(kind=8) :: dummy_real, dummy_real2(2)

            dydt = cero

            ! Calculate the angle of the asteroid
            theta = asteroid_omega * (t - initial_time) + asteroid_theta_correction ! theta = w(t-t0) + theta0
            dydt(1) = asteroid_omega ! for the sake of the integrator
            ! Set the asteroid derivatives
            dydt(3:4) = y(5:6)

            ! Calculate the positions of the boulders
            do i = 0, Nboulders
                rib(i,1) = cos(theta + theta_ast_arr(i)) * dist_ast_arr(i)
                rib(i,2) = sin(theta + theta_ast_arr(i)) * dist_ast_arr(i)
                vib(i,1) = -asteroid_omega * rib(i,2)
                vib(i,2) = asteroid_omega * rib(i,1)
            end do
            rib(0:,1) = y(3) + rib(0:,1) ! xA + rcos(theta + theta_i)
            rib(0:,2) = y(4) + rib(0:,2) ! yA + rsin(theta + theta_i)
            vib(0:,1) = y(5) + vib(0:,1) ! vxA - w r sin(theta + theta_i)
            vib(0:,2) = y(6) + vib(0:,2) ! vyA + w r cos(theta + theta_i)

            ! Define center of mass and velocity
            rcm = y(3:4)
            vcm = y(5:6)

            ! Particles accelerations
            !$OMP PARALLEL DEFAULT(SHARED) &
            !$OMP PRIVATE(i,rb,vb,ab,particle_i)
            !$OMP DO SCHEDULE (STATIC)
            do i = 1, Nactive
                ab = cero
                particle_i = i * equation_size + 2
                rb(1) = y(particle_i+1)
                rb(2) = y(particle_i+2)
                vb(1) = y(particle_i+3)
                vb(2) = y(particle_i+4)
                call apply_force(&
                    & t, &
                    & mass_ast_arr, particles_mass(i), &
                    & particles_hexit(i), &
                    & rb, vb, &
                    & rib, vib, &
                    & asteroid_mass, rcm, vcm, &
                    & asteroid_inertia, & ! inertia no varía
                    & dummy_real, dummy_real2, ab) ! Solo aceleración a partículas
                dydt(particle_i+1) = vb(1)
                dydt(particle_i+2) = vb(2)
                dydt(particle_i+3) = ab(1)
                dydt(particle_i+4) = ab(2)
            end do
            !$OMP END DO
            !$OMP END PARALLEL

        end function dydt_explicit_v2
        
        function dydt_implicit_v2 (t, y) result(dydt)
            !y = /theta, omega, xA, yA, vxA, vyA, Part, .../
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: rib(0:Nboulders,2), vib(0:Nboulders,2)
            real(kind=8) :: rb(2), vb(2), ab(2)
            real(kind=8) :: rcm(2), vcm(2)
            real(kind=8) :: theta, omega, omegadot
            integer(kind=4) :: i, particle_i
            real(kind=8) :: dummy_real2(2)

            dydt = cero

            ! Calculate the angle of the asteroid
            theta = y(1)
            omega = y(2)
            dydt(1) = omega
            ! Set the asteroid derivatives
            dydt(3:4) = y(5:6)
            dydt(5:6) = (/cero, cero/)

            ! Calculate the positions of the boulders
            do i = 0, Nboulders
                rib(i,1) = cos(theta + theta_ast_arr(i)) * dist_ast_arr(i)
                rib(i,2) = sin(theta + theta_ast_arr(i)) * dist_ast_arr(i)
                vib(i,1) = -omega * rib(i,2)
                vib(i,2) = omega * rib(i,1)
            end do
            rib(0:,1) = y(3) + rib(0:,1) ! xA + rcos(theta + theta_i)
            rib(0:,2) = y(4) + rib(0:,2) ! yA + rsin(theta + theta_i)
            vib(0:,1) = y(5) + vib(0:,1) ! vxA - w r sin(theta + theta_i)
            vib(0:,2) = y(6) + vib(0:,2) ! vyA + w r cos(theta + theta_i)

            ! Define center of mass and velocity
            rcm = y(3:4)
            vcm = y(5:6)
            
            ! Particles accelerations and omegadot
            omegadot = cero ! Init
            !$OMP PARALLEL DEFAULT(SHARED) &
            !$OMP PRIVATE(i,rb,vb,ab,particle_i)
            !$OMP DO REDUCTION(+:omegadot) SCHEDULE (STATIC)
            do i = 1, Nactive
                ab = cero
                particle_i = i * equation_size + 2
                rb(1) = y(particle_i+1)
                rb(2) = y(particle_i+2)
                vb(1) = y(particle_i+3)
                vb(2) = y(particle_i+4)
                call apply_force(&
                    & t, &
                    & mass_ast_arr, particles_mass(i), &
                    & particles_hexit(i), &
                    & rb, vb, &
                    & rib, vib, &
                    & asteroid_mass, rcm, vcm, &
                    & asteroid_inertia, & ! inertia no varía
                    & omegadot, dummy_real2, ab)
                dydt(particle_i+1) = vb(1)
                dydt(particle_i+2) = vb(2)
                dydt(particle_i+3) = ab(1)
                dydt(particle_i+4) = ab(2)
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            dydt(2) = domegadt(t, omega) + omegadot ! Update omega

        end function dydt_implicit_v2      
        
        function dydt_no_boulders (t, y) result(dydt)
            !y = /xA, yA, vxA, vyA, Part, .../
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: dydt
            real(kind=8) :: r0b(2), v0b(2)
            real(kind=8) :: rb(2), vb(2), ab(2)
            integer(kind=4) :: i, particle_i

            dydt = cero

            ! Set the asteroid derivatives
            dydt(1:2) = y(3:4)
            dydt(3:4) = (/cero, cero/)
            
            ! Calculate the positions of m0
            r0b = y(1:2)
            v0b = y(3:4)

            !$OMP PARALLEL DEFAULT(SHARED) &
            !$OMP PRIVATE(i,rb,vb,ab,particle_i)
            !$OMP DO SCHEDULE (STATIC)
            do i = 1, Nactive
                ab = cero
                particle_i = i * equation_size
                rb(1) = y(particle_i+1)
                rb(2) = y(particle_i+2)
                vb(1) = y(particle_i+3)
                vb(2) = y(particle_i+4)
                call apply_force_no_boulders(&
                    & t, &
                    & mass_ast_arr(0), &
                    & particles_hexit(i), &
                    & rb, vb, &
                    & r0b, v0b, &
                    & ab)
                dydt(particle_i+1) = vb(1)
                dydt(particle_i+2) = vb(2)
                dydt(particle_i+3) = ab(1)
                dydt(particle_i+4) = ab(2)
            end do
            !$OMP END DO
            !$OMP END PARALLEL
            
        end function dydt_no_boulders
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!! ACCELERATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine apply_force(t, m, mp, hexit_p, rb, vb, rib, vib, mcm, rcm, vcm, inertia, domegadt, acm, ab)
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), intent(in) :: m(0:Nboulders), mp
            integer(kind=4), intent(inout) :: hexit_p
            real(kind=8), intent(in) :: rb(2), vb(2)
            real(kind=8), intent(in) :: rib(0:Nboulders,2), vib(0:Nboulders,2)
            real(kind=8), intent(inout) :: mcm, rcm(2), vcm(2)
            real(kind=8), intent(inout) :: inertia
            real(kind=8), intent(inout) :: domegadt, acm(2), ab(2) ! YA DEBEN VENIR INICIALIZADOS
            real(kind=8) :: r_from_i(0:Nboulders,2), r_from_cm(2), v_from_cm(2)
            real(kind=8) :: dist_from_i(0:Nboulders), dist_from_cm, torque
            integer(kind=4) :: i

            ! Calculate distance and vector to boulders
            do i = 0, Nboulders
                r_from_i(i,:) = rb - rib(i,:)
                dist_from_i(i) = sqrt(r_from_i(i,1)*r_from_i(i,1) + r_from_i(i,2)*r_from_i(i,2))
            end do

            !Lets see what we need...
            if ((use_torque .and. (mp > tini)) .or. use_naive_stokes .or. use_stokes) then
                ! If cm not given, calculate it. mcm is the condition
                if (mcm .le. cero) then
                    mcm = sum(m)
                    rcm = cero
                    vcm = cero
                    do i = 0, Nboulders
                        rcm = rcm + m(i) * rib(i,:)
                        vcm = vcm + m(i) * vib(i,:)
                    end do
                    rcm = rcm / mcm
                    vcm = vcm / mcm
                end if
                r_from_cm = rb - rcm
                v_from_cm = vb - vcm
            else 
                r_from_cm = rb
                v_from_cm = vb
            end if
            dist_from_cm = sqrt(r_from_cm(1)*r_from_cm(1) + r_from_cm(2)*r_from_cm(2))
            
            ! Check if collision or escape
            if (any(dist_from_i < radius_ast_arr)) hexit_p = 1
            if (dist_from_cm < min_distance) hexit_p = 1
            if (dist_from_cm > max_distance) hexit_p = 2
            if (hexit_p > 0) return
        
            ! Calculate gravity (and torque if needed)
            if (use_torque .and. (mp > tini)) then
                torque = cero
                call gravitational_acceleration_and_torque(m, mp, rib, rcm, r_from_i, dist_from_i, torque, ab)
                if (inertia .le. cero) then
                    inertia = cero
                    do i = 0, Nboulders
                        inertia = inertia + m(i) * sum(rib(i,:) * rib(i,:))
                    end do
                end if
                domegadt = domegadt + torque / inertia
            else
                call gravitational_acceleration(m, r_from_i, dist_from_i, ab)
            end if

            ! Calculate stokes
            if (use_stokes) call stokes_acceleration(mcm, t, r_from_cm, v_from_cm, dist_from_cm, ab)

            ! Calculate naive stokes
            if (use_naive_stokes) call naive_stokes_acceleration(mcm, t, r_from_cm, v_from_cm, dist_from_cm, ab)

            ! Calculate J2
            if (use_J2) call J2_acceleration(m(0), r_from_i(0,:), dist_from_i(0), ab)

        end subroutine apply_force

        subroutine apply_force_no_boulders(t, m0, hexit_p, rb, vb, r0b, v0b, ab)
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), intent(in) :: m0
            integer(kind=4), intent(inout) :: hexit_p
            real(kind=8), intent(in) :: rb(2), vb(2)
            real(kind=8), intent(in) :: r0b(2), v0b(2)
            real(kind=8), intent(inout) :: ab(2) ! YA DEBEN VENIR INICIALIZADOS
            real(kind=8) :: r_from_0(2), v_from_0(2)
            real(kind=8) :: dist_from_0

            ! Calculate distance and vector to boulders
            r_from_0 = rb - r0b
            v_from_0 = vb - v0b
            dist_from_0 = sqrt(r_from_0(1)*r_from_0(1) + r_from_0(2)*r_from_0(2))
            
            ! Check if collision or escape
            if (dist_from_0 < min_distance) hexit_p = 1
            if (dist_from_0 > max_distance) hexit_p = 2
            if (hexit_p > 0) return
        
            ! Calculate gravity (and torque if needed)
            call gravitational_acceleration((/m0/), r_from_0, (/dist_from_0/), ab)

            ! Calculate stokes
            if (use_stokes) call stokes_acceleration(m0, t, r_from_0, v_from_0, dist_from_0, ab)

            ! Calculate naive stokes
            if (use_naive_stokes) call naive_stokes_acceleration(m0, t, r_from_0, v_from_0, dist_from_0, ab)

            ! Calculate J2
            if (use_J2) call J2_acceleration(m0, r_from_0, dist_from_0, ab)

        end subroutine apply_force_no_boulders
        
        !!!
        !!! CONFIGURATION
        !!!
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!! GRAVITY !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine gravitational_acceleration (m_bouls, r_from_i, dist_from_i, ab)
            implicit none
            real(kind=8), intent(in) :: m_bouls(0:Nboulders), r_from_i(0:Nboulders,2), dist_from_i(0:Nboulders)
            real(kind=8), intent(inout) :: ab(2)
            real(kind=8) :: acc_grav(2)
            integer(kind=4) :: i

            acc_grav = cero
            do i = 0, Nboulders
                acc_grav = acc_grav - m_bouls(i) * r_from_i(i,:) / dist_from_i(i)**3
            end do
            ab = ab + acc_grav * G
        end subroutine gravitational_acceleration

        !!!! INCLUDING TORQUE
        subroutine gravitational_acceleration_and_torque (m_bouls, mp, rib, rcm, r_from_i, dist_from_i, torque, ab)
            implicit none
            real(kind=8), intent(in) :: m_bouls(0:Nboulders), mp
            real(kind=8), intent(in) :: rib(0:Nboulders,2), rcm(2)
            real(kind=8), intent(in) :: r_from_i(0:Nboulders,2), dist_from_i(0:Nboulders)
            real(kind=8), intent(inout) :: torque, ab(2)
            real(kind=8) :: acc_grav(2), aux(2), torque_p
            integer(kind=4) :: i

            acc_grav = cero
            torque_p = cero
            do i = 0, Nboulders
                aux = - m_bouls(i) * r_from_i(i,:) / dist_from_i(i)**3 ! acc_pB 
                acc_grav = acc_grav + aux
                torque_p = torque_p + cross2D(rib(i,:) - rcm, - aux * mp) ! r x F_Bp = r x ( - acc_pB * m_p)
            end do
            torque = torque + torque_p * G
            ab = ab + acc_grav * G
        end subroutine gravitational_acceleration_and_torque
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!! STOKES !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! Init Parameters

        subroutine set_stokes_C_and_alpha (tau_a, tau_e, C, alpha)
            implicit none
            real(kind=8), intent(in) :: tau_a, tau_e
            real(kind=8), intent(out) :: C, alpha
            
            C = uno / (dos * tau_a) + uno / tau_e
            alpha = (dos * tau_a) / ((dos * tau_a) + tau_e)
        end subroutine set_stokes_C_and_alpha

        !!! Acceleration

        subroutine stokes_acceleration (mcm, t, r_from_cm, v_from_cm, dist_from_cm, ab)
            implicit none
            real(kind=8), intent(in) :: mcm, t, r_from_cm(2), v_from_cm(2), dist_from_cm
            real(kind=8), intent(inout) :: ab(2)
            real(kind=8) :: stokes_factor, vel_circ(2)

            stokes_factor = uno2 * (uno + tanh(1.d1 * (uno - t / stokes_charac_time)))
            vel_circ  = sqrt(G * mcm / dist_from_cm**3) * (/-r_from_cm(2), r_from_cm(1)/)
            ab = ab - stokes_C * (v_from_cm - stokes_alpha * vel_circ) * stokes_factor ! =-C(v * - alpha*vc)
        end subroutine stokes_acceleration

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! NAIVE-STOKES !!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine naive_stokes_acceleration (mcm, t, r_from_cm, v_from_cm, dist_from_cm, ab)
            implicit none
            real(kind=8), intent(in) :: mcm, t, r_from_cm(2), v_from_cm(2), dist_from_cm 
            real(kind=8), intent(inout) :: ab(2)
            real(kind=8) :: drag_factor, aux_real
            real(kind=8) :: Gmcm, v2, acc_radial, vel_radial, mean_movement 

            drag_factor = uno2 * (uno + tanh(1.d1 * (uno - t / drag_charac_time)))
            Gmcm = G * mcm
            v2 = dot_product(v_from_cm, v_from_cm)

            ! Debemos chequear que la partícula no esté "desligada"
            aux_real = dos * Gmcm / dist_from_cm - v2
            if (aux_real < cero) return ! No se puede calcular

            mean_movement = aux_real**(1.5d0) / Gmcm ! n
            vel_radial = dot_product(v_from_cm, r_from_cm) / dist_from_cm 
            acc_radial = - drag_coefficient * mean_movement * vel_radial
            
            ab = ab + acc_radial * r_from_cm / dist_from_cm * drag_factor
            user_real_var2 = acc_radial * r_from_cm / dist_from_cm * drag_factor
        end subroutine naive_stokes_acceleration

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!! GEO-POTENTIAL !!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine J2_acceleration(m0, r_from_primary, distance_from_primary, ab)
            implicit none
            real(kind=8), intent(in) :: m0, r_from_primary(2), distance_from_primary
            real(kind=8), intent(inout) :: ab(2)
            
            ab = ab - G * m0 * J2_effective * r_from_primary / distance_from_primary**5
        end subroutine J2_acceleration       
    
end module forces
