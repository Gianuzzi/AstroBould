module stokes
    use constants, only: cero, uno, uno2, dos, G
    use parameters, only: stokes_C, stokes_alpha, drag_coefficient, drag_charac_time, stokes_charac_time
    implicit none
    procedure (stokes_template), pointer :: apply_stokes => null ()
    procedure (stokes_template), pointer :: apply_naive_stokes => null ()
    
    abstract interface
        subroutine stokes_template(mcm, t, r_from_cm, v_from_cm, dist_from_cm, ab)
            implicit none
            real(kind=8), intent(in) :: mcm, t, r_from_cm(2), v_from_cm(2), dist_from_cm
            real(kind=8), intent(inout) :: ab(2)
        end subroutine stokes_template
    end interface
    
    contains
    
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
        end subroutine naive_stokes_acceleration
        
end module stokes
