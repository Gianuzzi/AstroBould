module const
    implicit none
    !!! Constantes.
    real(kind=8), parameter :: twopi  = 8.0*atan(1.0d0)
    real(kind=8), parameter :: cero   = 0.0d0
    real(kind=8), parameter :: uno    = 1.0d0
    real(kind=8), parameter :: dos    = 2.0d0
    real(kind=8), parameter :: uno2   = 0.5d0
    real(kind=8), parameter :: uno3   = uno/3.0d0
    real(kind=8), parameter :: pi     = uno2*twopi
    real(kind=8), parameter :: rad    = twopi/360.0d0
    real(kind=8), parameter :: inf    = 1.0d30        ! Infinito
    real(kind=8), parameter :: eps    = 1.0d-13       ! Precisión / Error [Usado en elem y coord]
    real(kind=8), parameter :: tini   = 2.2204d-16
    real(kind=8), parameter :: unit_m = 1.d-15 ! [Kg]
    real(kind=8), parameter :: unit_r = 1.d0   ! [Km]
    real(kind=8), parameter :: unit_t = 1.d0   ! [Dia]
    real(kind=8), parameter :: unit_v = unit_r / unit_t ! [Km/dia]
    real(kind=8), parameter :: unit_a = unit_v / unit_t ! [Km/dia^2]
    real(kind=8), parameter :: G_aux  = 4.9823394d-10 ! [km^3 kg^(-1) dia^(-2)]
    real(kind=8), parameter :: G      = G_aux * (unit_r**3) / unit_m / unit_t ! [unit_d^3 unit_m^(-1) unit_t^(-2)]
end module const