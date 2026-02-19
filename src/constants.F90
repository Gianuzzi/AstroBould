!> Module with constant values.
module constants
    use iso_fortran_env, only: real64, real32  ! or int32, etc.
    implicit none

    ! Default precision
#ifdef WP
    integer, parameter :: wp = WP
#else
    integer, parameter :: wp = real64
#endif

    ! Constantes
    real(wp), parameter :: twopi = 8.0e0_wp*atan(1.0e0_wp)
    real(wp), parameter :: cero = 0.0e0_wp
    real(wp), parameter :: uno = 1.0e0_wp
    real(wp), parameter :: dos = 2.0e0_wp
    real(wp), parameter :: tres = 3.0e0_wp
    real(wp), parameter :: uno2 = 0.5e0_wp
    real(wp), parameter :: uno3 = uno/tres
    real(wp), parameter :: pi = uno2*twopi
    real(wp), parameter :: radian = twopi/360.0e0_wp
    real(wp), parameter :: degree = uno/radian
    real(wp), parameter :: infinito = 1.0e50_wp  ! Infinito
    real(wp), parameter :: myepsilon = epsilon(1.0_wp) ! Precisión / Error [Usado en elem y coord]
    real(wp), parameter :: sqepsilon = sqrt(myepsilon) ! Error usado para checkeos
    real(wp), parameter :: tini = tiny(1.0_wp) ! Lowest value to compare
    real(wp), parameter :: unit_mass = 1.e0_wp   ! in [kg]
    real(wp), parameter :: unit_dist = 1.e0_wp   ! in [km]
    real(wp), parameter :: unit_time = 1.e0_wp   ! in [day]
    real(wp), parameter :: unit_vel = unit_dist/unit_time ! [km day⁻¹]
    real(wp), parameter :: unit_acc = unit_vel/unit_time ! [km day⁻²]
    real(wp), parameter :: unit_angm = unit_mass*unit_dist*unit_vel ! [kg km² day⁻¹]
    real(wp), parameter :: unit_ener = unit_angm/unit_time ! [kg km² day⁻²]
    real(wp), parameter :: hora = unit_time/24.e0_wp ! [day]
    real(wp), parameter :: minuto = hora/60.e0_wp ! [day]
    real(wp), parameter :: segundo = minuto/60.e0_wp ! [day]
    real(wp), parameter :: metro = unit_dist/1.e3_wp ! [km]
    real(wp), parameter :: G_aux = 4.9823394e-10_wp ! [km³ kg⁻¹ day⁻²]
    real(wp), parameter :: G = G_aux*(unit_dist**3)/unit_mass/unit_time ! [unit_r³ unit_m⁻¹ unit_t⁻²]
    real(wp), parameter :: megno_factor = uno !
end module constants
