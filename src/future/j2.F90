module j2
    use constants, only: G
    use parameters, only: J2_effective
    implicit none
    procedure (j2_template), pointer :: apply_J2 => null ()
    
    abstract interface
        subroutine j2_template(m0, r, dist, acc)
            implicit none
            real(kind=8), intent(in) :: m0, r(2), dist
            real(kind=8), intent(inout) :: acc(2)
        end subroutine j2_template
    end interface
    contains
        
        subroutine J2_acceleration(m0, r_from_0, dist_from_0, acc)
            implicit none
            real(kind=8), intent(in) :: m0, r_from_0(2), dist_from_0
            real(kind=8), intent(inout) :: acc(2)
            
            acc = acc - G * m0 * J2_effective * r_from_0 / dist_from_0**5
        end subroutine J2_acceleration
        
        subroutine ignore_J2(m0, r_from_0, dist_from_0, acc)
            implicit none
            real(kind=8), intent(in) :: m0, r_from_0(2), dist_from_0
            real(kind=8), intent(inout) :: acc(2)
        end subroutine ignore_J2
end module j2
