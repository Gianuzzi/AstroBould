module constants
    implicit none
    !!! Constantes.
    real(kind=8), parameter :: twopi      = 8.0 * atan(1.0d0)
    real(kind=8), parameter :: cero       = 0.0d0
    real(kind=8), parameter :: uno        = 1.0d0
    real(kind=8), parameter :: dos        = 2.0d0
    real(kind=8), parameter :: uno2       = 0.5d0
    real(kind=8), parameter :: uno3       = uno / 3.0d0
    real(kind=8), parameter :: pi         = uno2 * twopi
    real(kind=8), parameter :: radian     = twopi / 360.0d0
    real(kind=8), parameter :: infinity   = 1.0d50  ! Infinito
    real(kind=8), parameter :: epsilon    = 1.0d-14 ! Precisión / Error [Usado en elem y coord]
    real(kind=8), parameter :: segundo    = uno / 86400.d0
    real(kind=8), parameter :: tini       = 1.d-30 ! Lowest value
    real(kind=8), parameter :: unit_mass  = 1.d0   ! [kg]
    real(kind=8), parameter :: unit_dist  = 1.d0   ! [km]
    real(kind=8), parameter :: unit_time  = 1.d0   ! [day]
    real(kind=8), parameter :: unit_vel   = unit_dist / unit_time ! [km day⁻¹]
    real(kind=8), parameter :: unit_acc   = unit_vel / unit_time ! [km day⁻²]
    real(kind=8), parameter :: unit_angm  = unit_mass * unit_dist / unit_vel ! [kg km² day⁻¹]
    real(kind=8), parameter :: unit_ener  = unit_angm / unit_time ! [kg km² day⁻²]
    real(kind=8), parameter :: G_aux      = 4.9823394d-10 ! [km³ kg⁻¹ day⁻²]
    real(kind=8), parameter :: G          = G_aux * (unit_dist**3) / unit_mass / unit_time ! [unit_r³ unit_m⁻¹ unit_t⁻²]
end module constants
