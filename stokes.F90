module stokes
    use const
    use parameters
    implicit none
    real(kind=8) :: tau_a=inf, tau_e=inf, C_stk, a_stk
    real(kind=8) :: t_stokes=inf, f_stk=uno

    contains

        subroutine Candalpha(tau_a,tau_e,C,alpha)
            implicit none
            real(kind=8), intent(in)  :: tau_a, tau_e
            real(kind=8), intent(out) :: C, alpha
            C    =  uno / (dos * tau_a) + uno / tau_e
            alpha = (dos * tau_a) / ((dos * tau_a) + tau_e)
        end subroutine Candalpha

        subroutine getfact_stok(t,f_stk)
            implicit none
            real(kind=8), intent(in)    :: t
            real(kind=8), intent(out) :: f_stk
            f_stk = uno2 * (uno + tanh(1.d1 * (uno - t / t_stokes)))
        end subroutine getfact_stok

        subroutine accsto(t,rb,vb,rdd)
            implicit none
            real(kind=8), intent(in) :: t, rb(2), vb(2)
            real(kind=8), intent(inout) :: rdd(2)
            real(kind=8) :: d, vc(2)

            if (t_stokes < tini) return
            call getfact_stok(t,f_stk)
            d = sqrt(rb(1)*rb(1) + rb(2)*rb(2))
            vc  = sqrt(GM / (d*d*d)) * (/-rb(2),rb(1)/)
            rdd = rdd - C_stk * (vb - a_stk * vc) * f_stk ! -C(v * - alpha*vc)
        end subroutine accsto
end module stokes