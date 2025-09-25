module derivates
    use constants, only: G, cero, uno, uno2, dos, tini
    use auxiliary, only: cross2D_z
    use parameters, only: sim, system, &
                          & boulders_coords, boulders_data, & !! (Nb, 4) |mass,radius,theta_Ast0,dist_Ast|
                          & m_arr, R_arr, hard_exit
    use accelerations, only: use_damp, damp_time, damp_coef_1, damp_coef_2, damp_model, &
                            & use_drag, drag_coef, drag_time, &
                            & use_stokes, stokes_C, stokes_alpha, stokes_time, &
                            & use_J2, J2_coef


    implicit none

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

        function get_index(i) result(idx)
            implicit none
            integer(kind=4), intent(in) :: i
            integer(kind=4) :: idx
            idx = 4 * i - 1
        end function get_index


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!! DERIVATIVES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        function dydt (t, y) result(der)
            !y = /theta, omega, xA, yA, vxA, vyA, Part, .../
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size(y))       :: der
            real(kind=8) :: theta, omega
            real(kind=8) :: coords_A(4), coords_M(4), coords_P(4), dr_vec(2), dr, dr2
            real(kind=8) :: acc_grav_m(2), torque
            integer(kind=4) :: i, idx
            integer(kind=4) :: j, jdx
            integer(kind=4) :: N_total, last_moon
            real(kind=8) :: Gmast, Gmcomb, dr_ver(2), dv_vec(2), dr3  ! For extra forces
            real(kind=8) :: vel_circ(2), v2  ! For extra forces
            real(kind=8) :: aux_real  ! For extra forces
            real(kind=8) :: mean_movement  ! For extra forces
            real(kind=8) :: vel_radial(2), acc_radial(2)  ! For extra forces
            real(kind=8) :: damp_f, drag_f, stokes_f  ! For extra forces
            logical :: use_extra
            
            der = cero  ! init der at cero

            last_moon = system%Nmoons_active + 1  ! This would be the last moon (+ 1 bc asteroid is 1)
            N_total = last_moon + system%Nparticles_active

            ! Calculate the angle of the asteroid
            theta = y(1)
            omega = y(2)
            der(1) = omega
            

            ! Set the position derivates: dX/dt = V
            do i = 1, N_total
                idx = get_index(i)
                der(idx:idx+1) = y(idx+2:idx+3)
            end do
            
            coords_A = y(3:6)  ! Asteroid pos + vel
            torque = cero  ! Init torque at cero


            ! ---> Omega Damping <---
            if (use_damp) then
                !! Damping
                damp_f = uno2 * (uno + tanh(1.d1 * (uno - t / damp_time)))
                select case (damp_model)
                    case (1) ! domega/dt = tau
                        der(2) = der(2) + damp_coef_1 * damp_f
                    case (2) ! domega/dt = -exp(- (t-t0) / tau) * omega0 / tau = - omega / tau     
                        der(2) = der(2) - omega / damp_coef_1 * damp_f
                    case (3) ! domega/dt = A * B * (t-t0)**(B-1) * omega0 * exp (A * (t-t0)**B) = (A * B * (t-t0)**(B-1)) * omega
                        der(2) = der(2) +  &
                                & damp_coef_1 * (t - cero + tini)**(damp_coef_2 - uno) * omega * damp_f
                end select
            end if

            ! ---> Forces acting from COM of asteroid <---
            if (use_J2 .or. use_drag .or. use_stokes) then
                Gmast = G * m_arr(1)
                if (use_stokes) stokes_f = uno2 * (uno + tanh(1.d1 * (uno - t / stokes_time)))
                if (use_drag) drag_f = uno2 * (uno + tanh(1.d1 * (uno - t / drag_time)))
                use_extra = .True.
            else 
                use_extra = .False.
            end if

            !! Moons (massive)
            do j = 2, last_moon
                jdx = get_index(j)
                coords_P = y(jdx:jdx+3)  ! Moon

                !! ASTEROID AND MOON
                dr_vec = coords_P(1:2) - coords_A(1:2)  ! From Asteroid to Moon
                dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                dr = sqrt(dr2)

                ! Check if collision or Escape
                if ((dr < R_arr(1) + R_arr(j)) .or. &
                  & (dr < sim%min_distance) .or. &
                  & (dr > sim%max_distance .and. sim%max_distance > cero)) hard_exit = .True.

                ! Check if more needed
                if (use_extra) then
                    Gmcomb = Gmast + (G * m_arr(j))
                    dr_ver = dr_vec / dr
                    dr3 = dr2 * dr
                end if

                ! ---> J2 <---
                if (use_J2) der(jdx+2:jdx+3) = der(jdx+2:jdx+3) - Gmcomb * J2_coef * dr_ver / (dr3 * dr)

                ! ---> Drag <---
                if (use_drag) then
                    dv_vec = coords_P(3:4) - coords_A(3:4)  ! Velocity from Asteroid to Moon or Particle 
                    v2 = dot_product(dv_vec, dv_vec)
                    
                    aux_real = dos * Gmcomb / dr - v2  ! Check if unbound
                    if (aux_real > cero) then ! Can calculate only in this case
                        mean_movement = aux_real**(1.5d0) / Gmcomb ! n
                        vel_radial = dot_product(dr_ver, dv_vec)
                        acc_radial = - drag_coef * mean_movement * vel_radial                
                        der(jdx+2:jdx+3) = der(jdx+2:jdx+3) + acc_radial * dr_ver * drag_f ! = -a_r * (x, y) / r * factor
                    end if
                end if

                ! ---> Stokes <---
                if (use_stokes) then
                    vel_circ  = sqrt(Gmcomb / dr3) * (/-dr_vec(2), dr_vec(1)/)  ! v_circ = n (-y, x)
                    der(jdx+2:jdx+3) = der(jdx+2:jdx+3) - stokes_C * (dv_vec - stokes_alpha * vel_circ) * stokes_f ! = -C * (v - alpha * vc) * factor
                end if

            end do

            !! Particles (massless)
            do j = last_moon + 1, N_total
                jdx = get_index(j)
                coords_P = y(jdx:jdx+3)  ! Particle

                !! ASTEROID AND PARTICLE
                dr_vec = coords_P(1:2) - coords_A(1:2)  ! From Asteroid to Particle
                dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                dr = sqrt(dr2)

                ! Check if collision or Escape
                if ((dr < sim%min_distance) .or. &
                  & (dr > sim%max_distance .and. sim%max_distance > cero)) hard_exit = .True.

                ! Check if more needed
                if (use_extra) then
                    dr_ver = dr_vec / dr
                    dr3 = dr2 * dr
                end if

                ! ---> J2 <---
                if (use_J2) der(jdx+2:jdx+3) = der(jdx+2:jdx+3) - Gmast * J2_coef * dr_ver / (dr3 * dr)

                ! ---> Drag <---
                if (use_drag) then
                    dv_vec = coords_P(3:4) - coords_A(3:4)  ! Velocity from Asteroid to Particle 
                    v2 = dot_product(dv_vec, dv_vec)
                    
                    aux_real = dos * Gmast / dr - v2  ! Check if unbound
                    if (aux_real > cero) then ! Can calculate only in this case
                        mean_movement = aux_real**(1.5d0) / Gmast ! n
                        vel_radial = dot_product(dr_ver, dv_vec)
                        acc_radial = - drag_coef * mean_movement * vel_radial                
                        der(jdx+2:jdx+3) = der(jdx+2:jdx+3) + acc_radial * dr_ver * drag_f ! = -a_r * (x, y) / r * factor
                    end if
                end if

                ! ---> Stokes <---
                if (use_stokes) then
                    vel_circ  = sqrt(Gmast / dr3) * (/-dr_vec(2), dr_vec(1)/)  ! v_circ = n (-y, x)
                    der(jdx+2:jdx+3) = der(jdx+2:jdx+3) - stokes_C * (dv_vec - stokes_alpha * vel_circ) * stokes_f ! = -C * (v - alpha * vc) * factor
                end if

            end do

            ! ---> GRAVITY <---

            ! First, Asteroid to all
            do i = 0, system%asteroid%Nboulders

                !! Get boulder coords
                boulders_coords(i,1) = boulders_data(i,4) * cos(theta + boulders_data(i,3))  ! y
                boulders_coords(i,2) = boulders_data(i,4) * sin(theta + boulders_data(i,3))  ! x
                boulders_coords(i,3) = -omega * boulders_coords(i,2)  ! vx
                boulders_coords(i,4) = omega * boulders_coords(i,1)  ! vy

                boulders_coords(i,:) = boulders_coords(i,:) + coords_A  ! Move to Asteroid

                !! Particles (massless)
                do j = last_moon + 1, N_total
                    jdx = get_index(j)
                    coords_P = y(jdx:jdx+3)  ! Particle

                    !! BOULDER AND PARTICLE
                    dr_vec = coords_P(1:2) - boulders_coords(i,1:2)  ! From Boulder to Particle
                    dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                    dr = sqrt(dr2)

                    ! Check if collision
                    if (dr < boulders_data(i,2)) hard_exit = .True.

                    ! Particle acceleration
                    der(jdx+2:jdx+3) = der(jdx+2:jdx+3) - G * boulders_data(i,1) * dr_vec / (dr2 * dr)  ! G mBoul (x, y) / r³

                end do

                !! Moons (massive)
                do j = 2, last_moon ! +1 porque j=1 es asteroid
                    jdx = get_index(j)
                    coords_M = y(jdx:jdx+3)  ! Moon

                    !! BOULDER AND MOON
                    dr_vec = coords_M(1:2) - boulders_coords(i,1:2)  ! From Boulder to Moon
                    dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                    dr = sqrt(dr2)

                    ! Check if collision
                    if (dr < boulders_data(i,2) + R_arr(j)) hard_exit = .True.

                    ! Moon acceleration per unit mass (WITHOUT MASSES)
                    acc_grav_m = - G * dr_vec / (dr2 * dr)  !! G (x, y) / r³
                    der(jdx+2:jdx+3) = der(jdx+2:jdx+3) + acc_grav_m * boulders_data(i,1)  ! a_m * mass

                    !! ASTEROID FROM MOON: This force is "felt" by the asteroid CM
                    ! Asteroid acceleration
                    der(5:6) = der(5:6) - acc_grav_m * boulders_data(i,1) * m_arr(j) / m_arr(1)  ! Force moved to Asteroid CM

                    ! Torque to Asteroid
                    torque = torque - cross2D_z(boulders_coords(i,1:2) - coords_A(1:2), acc_grav_m * m_arr(j) * boulders_data(i,1))  !! r x F

                end do

            end do

            !! Update Omega with Torque
            der(2) = der(2) + torque / system%asteroid%inertia

            ! Second, Moons to particles
            do i = 2, last_moon
                idx = get_index(i)
                coords_M = y(idx:idx+3)  ! Moon 1 (M)

                !! Particles (massless)
                do j = last_moon + 1, N_total
                    jdx = get_index(j)
                    coords_P = y(jdx:jdx+3)  ! Particle

                    !! MOON AND PARTICLE
                    dr_vec = coords_P(1:2) - coords_M(1:2)  ! From Moon 1 (M) to Particle
                    dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                    dr = sqrt(dr2)

                    ! Check if collision
                    if (dr < R_arr(i)) hard_exit = .True.

                    ! Particle acceleration
                    der(jdx+2:jdx+3) = der(jdx+2:jdx+3) - G * m_arr(i) * dr_vec / (dr2 * dr)  ! G mMoon (x, y) / r³

                end do
            
            end do

            ! Third, Moons to Moons (if requested)
            if (sim%use_moon_gravity) then
                do i = 2, last_moon - 1
                    idx = get_index(i)
                    coords_M = y(idx:idx+3)  ! Moon 1 (M)

                    !! Other Moons (massive)
                    do j = i + 1, last_moon ! Ahora sí +1 porque j=1 es asteroid
                        jdx = get_index(j)
                        coords_P = y(jdx:jdx+3)  ! Moon 2 (P)

                        !! MOON 1 (M) AND MOON 2 (P)
                        dr_vec = coords_P(1:2) - coords_M(1:2)  ! From Moon 1 (M) to Moon 2 (P)
                        dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                        dr = sqrt(dr2)

                        ! Check if collision
                        if (dr < R_arr(i) + R_arr(j)) hard_exit = .True.

                        ! Moons acceleration per unit mass
                        acc_grav_m = G * dr_vec / (dr2 * dr)  !! G (x, y) / r³

                        ! Both accelerations
                        der(idx+2:idx+3) = der(idx+2:idx+3) + acc_grav_m * m_arr(j)
                        der(jdx+2:jdx+3) = der(jdx+2:jdx+3) - acc_grav_m * m_arr(i)

                    end do

                end do
            
            end if

        end function dydt  
    
end module derivates
