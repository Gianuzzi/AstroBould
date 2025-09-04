module derivates
    use constants, only: G
    use auxiliary, only: cross2D_z
    use parameters, only: system, &
                          & boulders_coords, boulders_data, & !! (Nb, 4) |mass,radius,theta_Ast0,dist_Ast|
                          & m_arr, R_arr, hexit_arr
    use accelerations
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
            real(kind=8) :: acc_grav(2), torque
            integer(kind=4) :: i, idx
            integer(kind=4) :: j, jdx
            integer(kind=4) :: N_total, last_moon
            
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
            
            coords_A = y(3:6)  ! Asteroid pos
            torque = cero  ! Init torque at cero

            ! ---> GRAVITY <---

            ! First, Asteroid to all
            do i = 0, system%asteroid%Nboulders

                !! Get boulder coords
                boulders_coords(i,1) = boulders_data(i,4) * cos(theta + boulders_data(i,3))
                boulders_coords(i,2) = boulders_data(i,4) * sin(theta + boulders_data(i,3))
                boulders_coords(i,3) = -omega * boulders_coords(i,2)  ! vx
                boulders_coords(i,4) = omega * boulders_coords(i,1)  ! vy

                boulders_coords(i,:) = boulders_coords(i,:) + coords_A

                !! Particles (massless)
                do j = last_moon + 1, N_total
                    jdx = get_index(j)
                    coords_P = y(jdx:jdx+3)  ! Particle

                    !! BOULDER AND PARTICLE
                    dr_vec = coords_P(1:2) - boulders_coords(i,1:2)  ! From Boulder to Particle
                    dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                    dr = sqrt(dr2)

                    ! Check if collision or escape
                    if (dr < boulders_data(i,2)) hexit_arr(j) = 1
                    ! if (dr < min_distance) hexit_arr(j) = 1
                    ! if (dr > max_distance) hexit_arr(j) = 2
                    if (hexit_arr(j) > 0) cycle  ! Collision !!

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

                    ! Check if collision or escape
                    if (dr < boulders_data(i,2) + R_arr(j)) hexit_arr(j) = 1
                    ! if (dr < min_distance) hexit_arr(j) = 1
                    ! if (dr > max_distance) hexit_arr(j) = 2
                    if (hexit_arr(j) > 0) cycle  ! Collision !!

                    ! Moon acceleration
                    acc_grav = - G * boulders_data(i,1) * dr_vec / (dr2 * dr)  !! G mBoul (x, y) / r³
                    der(jdx+2:jdx+3) = der(jdx+2:jdx+3) + acc_grav

                    ! Torque
                    torque = torque - cross2D_z(boulders_coords(i,1:2) - coords_A(1:2), acc_grav * m_arr(j))  !! r x F  !!! m_arr from PARAMS

                    !! ASTEROID FROM MOON ( Just one direction here)
                    dr_vec = coords_M(1:2) - coords_A(1:2)  ! From Asteroid to Moon
                    dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                    dr = sqrt(dr2)
                    ! Asteroid acceleration
                    der(5:6) = der(5:6) + G * m_arr(j) * dr_vec / (dr2 * dr)  ! m_arr from PARAMS

                end do

            end do

            !! Update Omega with Torque
            der(2) = torque / system%asteroid%inertia

            ! Second, Moons to all (but no asteroid)
            do i = 2, last_moon  ! Para no repetir, agregamos un exit debajo
                idx = get_index(i)
                coords_M = y(idx:idx+3)  ! Moon 1 (M)

                !! Particles (massless)
                do j = last_moon + 1, N_total
                    jdx = 4 * j - 1 
                    coords_P = y(jdx:jdx+3)  ! Particle

                    !! BOULDER AND PARTICLE
                    dr_vec = coords_P(1:2) - coords_M(1:2)  ! From Moon 1 (M) to Particle
                    dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                    dr = sqrt(dr2)

                    ! Check if collision or escape
                    if (dr < R_arr(i)) hexit_arr(j) = 1
                    ! if (dr < min_distance) hexit_arr(j) = 1
                    ! if (dr > max_distance) hexit_arr(j) = 2
                    if (hexit_arr(j) > 0) cycle  ! Collision !!

                    ! Particle acceleration
                    der(jdx+2:jdx+3) = der(jdx+2:jdx+3) - G * m_arr(i) * dr_vec / (dr2 * dr)  ! G mMoon (x, y) / r³

                end do

                if (i .eq. last_moon) exit  ! No repitamos la última

                !! Other Moons (massive)
                do j = i + 1, last_moon ! Ahora sí +1 porque j=1 es asteroid
                    jdx = get_index(j)
                    coords_P = y(jdx:jdx+3)  ! Moon 2 (P)

                    !! MOON 1 (M) AND MOON 2 (P)
                    dr_vec = coords_P(1:2) - coords_M(1:2)  ! From Moon 1 (M) to Moon 2 (P)
                    dr2 = dr_vec(1) * dr_vec(1) + dr_vec(2) * dr_vec(2)
                    dr = sqrt(dr2)

                    ! Check if collision or escape
                    if (dr < R_arr(i) + R_arr(j)) hexit_arr(j) = 1
                    ! if (dr < min_distance) hexit_arr(j) = 1
                    ! if (dr > max_distance) hexit_arr(j) = 2
                    if (hexit_arr(j) > 0) cycle  ! Collision !!

                    ! Moons acceleration per unit mass
                    acc_grav = G * dr_vec / (dr2 * dr)  !! G (x, y) / r³

                    ! Both accelerations
                    der(idx+2:idx+3) = der(idx+2:idx+3) + acc_grav * m_arr(j)
                    der(jdx+2:jdx+3) = der(jdx+2:jdx+3) - acc_grav * m_arr(i)

                end do

            end do


            !! Forces acting from COM of asteroid


        ! !!! Init Parameters

        ! subroutine set_stokes_C_and_alpha (tau_a, tau_e, C, alpha)
        !     implicit none
        !     real(kind=8), intent(in) :: tau_a, tau_e
        !     real(kind=8), intent(out) :: C, alpha
            
        !     C = uno / (dos * tau_a) + uno / tau_e
        !     alpha = (dos * tau_a) / ((dos * tau_a) + tau_e)
        ! end subroutine set_stokes_C_and_alpha

        ! !!! Acceleration

        ! subroutine stokes_acceleration (mcm, t, r_from_cm, v_from_cm, dist_from_cm, ab)
        !     implicit none
        !     real(kind=8), intent(in) :: mcm, t, r_from_cm(2), v_from_cm(2), dist_from_cm
        !     real(kind=8), intent(inout) :: ab(2)
        !     real(kind=8) :: stokes_factor, vel_circ(2)

        !     stokes_factor = uno2 * (uno + tanh(1.d1 * (uno - t / stokes_charac_time)))
        !     vel_circ  = sqrt(G * mcm / dist_from_cm**3) * (/-r_from_cm(2), r_from_cm(1)/)
        !     ab = ab - stokes_C * (v_from_cm - stokes_alpha * vel_circ) * stokes_factor ! =-C(v * - alpha*vc)
        ! end subroutine stokes_acceleration

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!! NAIVE-STOKES !!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! subroutine naive_stokes_acceleration (mcm, t, r_from_cm, v_from_cm, dist_from_cm, ab)
        !     implicit none
        !     real(kind=8), intent(in) :: mcm, t, r_from_cm(2), v_from_cm(2), dist_from_cm 
        !     real(kind=8), intent(inout) :: ab(2)
        !     real(kind=8) :: drag_factor, aux_real
        !     real(kind=8) :: Gmcm, v2, acc_radial, vel_radial, mean_movement 

        !     drag_factor = uno2 * (uno + tanh(1.d1 * (uno - t / drag_charac_time)))
        !     Gmcm = G * mcm
        !     v2 = dot_product(v_from_cm, v_from_cm)

        !     ! Debemos chequear que la partícula no esté "desligada"
        !     aux_real = dos * Gmcm / dist_from_cm - v2
        !     if (aux_real < cero) return ! No se puede calcular

        !     mean_movement = aux_real**(1.5d0) / Gmcm ! n
        !     vel_radial = dot_product(v_from_cm, r_from_cm) / dist_from_cm 
        !     acc_radial = - drag_coefficient * mean_movement * vel_radial
            
        !     ab = ab + acc_radial * r_from_cm / dist_from_cm * drag_factor
        ! end subroutine naive_stokes_acceleration

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!! GEO-POTENTIAL !!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! subroutine J2_acceleration(m0, r_from_primary, distance_from_primary, ab)
        !     implicit none
        !     real(kind=8), intent(in) :: m0, r_from_primary(2), distance_from_primary
        !     real(kind=8), intent(inout) :: ab(2)
            
        !     ab = ab - G * m0 * J2_effective * r_from_primary / distance_from_primary**5
        ! end subroutine J2_acceleration     

        end function dydt  
    
end module derivates
