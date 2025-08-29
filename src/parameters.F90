module parameters
    use constants
    use auxiliar
    use bins
    use omp_lib

    implicit none

    
    logical :: existe_configfile = .False.  ! Wether if config file exists
    
    ! This contains only the input parameters
    type :: params
        ! Times for the integration - 
        real(kind=8) :: initial_time = cero
        real(kind=8) :: final_time = cero
        real(kind=8) :: output_timestep = cero
        real(kind=8) :: output_number = cero
        integer(kind=4) :: case_output_type = 2
        ! Parameters for the Integration - 
        logical :: use_parallel = .False.
        integer(kind=4) :: requested_threads = 1
        ! Adaptive step integrations - 
        integer(kind=4) :: error_digits = 12
        real(kind=8) :: learning_rate = uno
        ! The primary -
        real(kind=8) :: mass_primary = cero
        real(kind=8) :: radius_primary = cero
        ! Rotation -
        real(kind=8) :: lambda_kep = cero
        real(kind=8) :: asteroid_rotational_period = cero
        ! Moons - 
        logical :: use_moonsfile = .False.
        character(30) :: moonsfile = ""
        ! Particles - 
        logical :: use_particlesfile = .False.
        character(30) :: particlesfile = ""
        ! Extra forces/effects - 
        real(kind=8) :: use_stokes = .False.
        real(kind=8) :: stokes_a_damping_time = infinity
        real(kind=8) :: stokes_e_damping_time = infinity
        real(kind=8) :: stokes_charac_time = cero
        logical :: use_naive_stokes = .False.
        real(kind=8) :: drag_coefficient = cero
        real(kind=8) :: drag_charac_time = cero
        logical :: use_omega_damping = .False.
        real(kind=8) :: omega_linear_damping_time = infinity
        real(kind=8) :: omega_exp_damping_time = infinity
        real(kind=8) :: omega_exp_poly_A = cero
        real(kind=8) :: omega_exp_poly_B = cero
        real(kind=8) :: omega_charac_time = cero
        real(kind=8) :: mass_exp_damping_time = infinity
        real(kind=8) :: J2_coefficient = cero
        ! Binned forces/effects - [NOT AVAILABLE YET. STILL UNDER DEVELOPMENT]
        logical :: use_self_gravity = .False.
        integer(kind=4) :: Norder_self_gravity = 3
        integer(kind=8) :: Nbins = 0
        integer(kind=4) :: binning_method = 1
        real(kind=8) :: rmin_bins = - uno ! -1 means first particle (initial)
        real(kind=8) :: rmax_bins = - uno ! -1 means last particle (initial)
        ! Conditions for Collision/Escape - 
        real(kind=8) :: min_distance = cero
        real(kind=8) :: max_distance = infinity
        logical :: use_merge = .False.
        ! Manual |(t)imes omega(t) mass_add(t)| file -
        logical :: use_tomfile = .False.
        character(30) :: tomfile = ""
        ! Output -
        logical :: use_screen = .False.
        logical :: use_datafile = .False.
        character(30) :: datafile = ""
        logical :: use_multiple_outputs = .False.
        character(30) :: multfile = ""
        logical :: use_chaosfile = .False.
        character(30) :: chaosfile = ""
        logical :: use_datascreen = .True.
        logical :: use_percentage = .False.
        logical :: use_elements_output = .True.
        logical :: use_potential_map = .False.
        character(30) :: mapfile = ""
        integer(kind=4) :: map_grid_size_x = 500
        integer(kind=4) :: map_grid_size_y = 500
        integer(kind=4) :: map_min_x = -300
        integer(kind=4) :: map_max_x = 300
        integer(kind=4) :: map_min_y = -300
        integer(kind=4) :: map_max_y = 300
        ! < Number of bodies >
        integer(kind=4) :: Nboulders = 0
        integer(kind=4) :: Nmoons = 0
        integer(kind=4) :: Nparticles = 0
        !! < Arrays with of bodies >
        real(kind=8), dimension(:,:), allocatable :: boulders_in !! mass, radius, theta
        real(kind=8), dimension(:,:), allocatable :: moons_in !! mass, a, e, M, w, Radius, MRR
        real(kind=8), dimension(:,:), allocatable :: particles_in !! a, e, M, w, MRR
    end type params
    
    ! Simulation
    integer(kind=4) :: simulation_number

    ! Integration dynamical variables
    real(kind=8) :: time  ! Actual time of the integration
    real(kind=8) :: timestep  ! This timestep 
    real(kind=8) :: adaptive_timestep  ! This adaptive timestep
    integer(kind=4) :: checkpoint_number  ! Number of actual timestep
    real(kind=8) :: min_timestep  ! Minimum timestep

    integer(kind=4) :: available_threads  ! Avaliable threads to use
    real(kind=8) :: my_threads  ! Amount of threads actually used
    logical :: compiled_with_openmp = .False.  ! Flag of compiple with OpenMP

    ! Bodies
    integer(kind=4) :: Ntotal  ! Amount of all bodies (asteroid is just 1)
    integer(kind=4) :: Nactive ! Active bodies
    
    integer(kind=4) :: first_particle
    integer(kind=4) :: last_particle
    integer(kind=4), dimension(:), allocatable :: particles_outcome
    integer(kind=4), dimension(:,:), allocatable :: ij_to_swap


    ! I/O
    !! Pantalla
    logical :: use_screen  ! Imprimir mensajes en pantalla
    logical :: use_datascreen  ! Imprimir salida de datos en pantalla
    logical :: use_percentage  ! Imprimir porcentaje en pantalla
    !! Archivos de salida
    !!! Logicals
    logical :: use_datafile  ! Escribir archivo de salida
    logical :: use_chaosfile  ! Escribir archivo de caos
    logical :: use_multiple_outputs  ! Múltiples salidas. 1 por partículas (body)
    logical :: use_potential_map  ! Escribir archio de potencial
    !!! Nombres
    character(30) :: datafile  ! Nombre de archivo de salida
    character(30) :: chaosfile  ! Nombre de archivo de caos
    character(30) :: multfile  ! Nombre base de archivo de múltiples
    character(30) :: mapfile  ! Nombre de archivo mapa
    !! Configuración de salida
    logical :: use_elements_output  ! Salida de datos en elementos orbitales
    logical :: use_flush_output  ! Hacer flush en los archivos manualmente
    logical :: use_update_chaos  ! Actualziar valores de chaos en cada timestep
    logical :: only_potential_map = .False.  ! Calcular solamente el potencial
    !! Archivos de entrada
    !!! TOM
    logical :: use_tomfile  ! Usar archivo TOM
    character(30) :: tomfile  ! Nombre de archivo TOM
    !!!! Auxiliares
    integer(kind=4) :: tom_index_number  ! Índice para contar que TOM row está activa
    integer(kind=4) :: tom_total_number  ! Total de líneas en TOM
    real(kind=8), dimension(:), allocatable :: tom_times  ! Tiempos en TOM
    real(kind=8), dimension(:), allocatable :: tom_deltaomega  ! Delta Omega en TOM
    real(kind=8), dimension(:), allocatable :: tom_deltamass  ! Delta Mass en TOM
    real(kind=8) :: tom_mass_growth_param
    logical, dimension(:), allocatable :: checkpoint_is_tom  ! Si este checkpoint es TOM
    logical, dimension(:), allocatable :: checkpoint_is_output ! Si este checkpoint es output
    !!! Lunas
    logical :: use_moonsfile  ! Usar archivo con partículas
    character(30) :: moonsfile  ! Nombre de archivo de partículas
    !!! Partículas
    logical :: use_particlesfile  ! Usar archivo con partículas
    character(30) :: particlesfile  ! Nombre de archivo de partículas
    !! Mapa
    integer(kind=4) :: map_grid_size_x
    integer(kind=4) :: map_grid_size_y
    real(kind=8) :: map_min_x
    real(kind=8) :: map_max_x
    real(kind=8) :: map_min_y
    real(kind=8) :: map_max_y

    ! Forces
    !! Internal forces
    logical :: use_torque
    ! External forces
    !! Stokes
    logical :: use_stokes
    real(kind=8) :: stokes_a_damping_time
    real(kind=8) :: stokes_e_damping_time
    real(kind=8) :: stokes_charac_time
    real(kind=8) :: stokes_C
    real(kind=8) :: stokes_alpha
    !! Naive-Stokes (Drag)
    logical :: use_naive_stokes
    real(kind=8) :: drag_coefficient
    real(kind=8) :: drag_charac_time
    !! Mass and omega damping
    logical :: use_omega_damping
    real(kind=8) :: omega_charac_time
    real(kind=8) :: omega_exp_damping_time
    real(kind=8) :: omega_exp_poly_A
    real(kind=8) :: omega_exp_poly_B
    real(kind=8) :: omega_exp_poly_AB
    real(kind=8) :: omega_linear_damping_time
    real(kind=8) :: omega_linear_damping_slope
    real(kind=8) :: mass_exp_damping_time
    !! Geo-potential
    logical :: use_J2
    real(kind=8) :: J2_coefficient
    real(kind=8) :: J2_effective
    !! Self-Gravity [Use BINS] [Not available yet]
    logical :: use_self_gravity
    integer(kind=4) :: Norder_self_gravity
    !! Viscosity [Use BINS] [Not available yet]
    logical :: use_viscosity
    real(kind=8) :: viscosity
    
    !! BINS [Not available yet]
    logical :: use_bins
    integer(kind=4) :: Nbins
    integer(kind=4) :: binning_method
    real(kind=8) :: rmin_bins
    real(kind=8) :: rmax_bins
    type(my_bins) :: disk_bins
    logical :: update_rmin_bins
    logical :: update_rmax_bins
    logical :: update_bins
    integer(kind=4), dimension(:), allocatable :: particles_bins
    

    ! Chaos
    logical :: use_chaos  ! Calcular valores para chaos
    real(kind=8), dimension(:), allocatable :: particles_max_a
    real(kind=8), dimension(:), allocatable :: particles_max_e
    real(kind=8), dimension(:), allocatable :: particles_min_a
    real(kind=8), dimension(:), allocatable :: particles_min_e
    real(kind=8), dimension(:), allocatable :: particles_times
    integer(kind=4) :: discarded_particles
    integer(kind=4) :: staying_particles
    

    ! Extra/Auxiliar variables
    integer(kind=4) :: i
    integer(kind=4) :: j
    integer(kind=4) :: k
    !! Arguments
    integer(kind=4) :: arguments_number
    !! Auxiliar
    logical :: aux_logical
    logical :: is_number
    real(kind=8) :: aux_real
    real(kind=8) :: dummy_real
    real(kind=8) :: dummy_real2
    integer(kind=4) :: aux_integer
    character(1) :: aux_character1
    character(2) :: aux_character2
    character(20) :: aux_character20
    character(30) :: aux_character30
    real(kind=8), dimension(2) :: distance_to_primary
    real(kind=8), dimension(2) :: aux_real_arr2
    real(kind=8), dimension(6) :: aux_real_arr6
    logical :: hard_center
    logical :: use_elements


    ! PARAMETERS ARRAY
    integer, parameter :: equation_size = 4
    real(kind=8), dimension(:), allocatable :: parameters_arr
    real(kind=8), dimension(:), allocatable :: parameters_arr_new
    real(kind=8), dimension(:), allocatable :: parameters_der
    
    !!! Hard Exit
    logical :: is_premature_exit
    integer(kind=4), dimension(:), allocatable, target :: particles_hexit !  Hard Exit integer
    integer(kind=4), pointer :: particles_hexitptr ! pointer to Hard Exit
    procedure (check_continue_template), pointer :: check_continue_ptr => null ()
    
    !!! INITIAL CONDITIONS (Here we store the initial conditions of everything)
    real(kind=8), dimension(:), allocatable :: asteroid_initial_conditions !! mass, radius, pos, vel, theta, omega, inertia, angmom, Prot, acorot
    real(kind=8), dimension(:,:), allocatable :: m0_and_boulders_initial_conditions !! mass, radius, theta
    real(kind=8), dimension(:,:), allocatable :: moons_initial_conditions !! mass, a, e, M, w, Radius, MRR
    real(kind=8), dimension(:,:), allocatable :: particles_initial_conditions !! a, e, M, w, MRR
    

    !!! GENERAL element PTRs
    procedure (int_i_template), pointer :: get_elements_i => null ()
    procedure (int_i_template), pointer :: get_chaos_i => null ()

    !!! MERGE
    logical :: use_merge
    integer :: merged_particles
    real(kind=8) :: mass_to_merge, angular_momentum_to_merge
    procedure (do_merge_template), pointer :: resolve_merge => null ()

    !!! WRITE
    procedure (write_i_to_unit_template), pointer :: write_i_to_general => null ()
    procedure (write_i_to_unit_template), pointer :: write_i_to_screen => null ()
    procedure (write_i_to_unit_template), pointer :: write_i_to_individual => null ()
    procedure (write_i_to_unit_template), pointer :: write_b_to_general => null ()
    procedure (write_i_to_unit_template), pointer :: write_b_to_screen => null ()
    procedure (write_i_to_unit_template), pointer :: write_b_to_individual => null ()
    procedure (int_i_template), pointer :: flush_chaos => null ()
    procedure (int_i_template), pointer :: flush_output => null ()

    !!!!!!!!    FORMATS    !!!!!
    character(19), parameter :: f12    = "(22(A, 1X, I7, 1X))"
    character(45), parameter :: f125131 = "(2(A, 1X, I7, 1X), 2(A, 1X, 1PE22.15, 1X), A)"
    character(29), parameter :: f1233  = "(A, 1X, I7, 22(1X, 1PE22.15))"
    character(25), parameter :: f13    = "(22(A, 1X, 1PE22.15, 1X))"
    character(21), parameter :: f133   = "(A, 22(1X, 1PE22.15))"
    character(27), parameter :: f1331  = "(A, 2(1X, 1PE22.15), 1X, A)"
    character(26), parameter :: f2233  = "(I7, I7, 22(1X, 1PE22.15))"
    character(26), parameter :: f23    = "(22(I7, 1X, 1PE22.15, 1X))"
    character(22), parameter :: f233   = "(I7, 22(1X, 1PE22.15))"


    abstract interface
        subroutine simple_ptr_template ()
            implicit none
        end subroutine simple_ptr_template

        subroutine int_i_template (i)
            implicit none
            integer(kind=4), intent(in) :: i
        end subroutine int_i_template

        subroutine do_merge_template (i, mass, angmom)
            implicit none
            integer(kind=4), intent(in) :: i
            real(kind=8), intent(inout) :: mass, angmom
        end subroutine do_merge_template

        subroutine write_i_to_unit_template (i, unit)
            implicit none
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(in) :: unit
        end subroutine write_i_to_unit_template
        
        function check_continue_template (y) result(keep_going)
            implicit none
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8) :: r0b(2)
            real(kind=8) :: rb(2), dist_to_m0, r_from_m0(2)
            integer(kind=4) :: i, particle_i
            logical :: keep_going
        end function

    end interface

    contains
    
        ! 0. Init default config
        subroutine init_default_config()
            implicit none
            ! Times for the integration - 
            initial_time = cero
            final_time = cero
            output_timestep = cero
            output_number = cero
            case_output_type = 2
            ! Parameters for the Integration - 
            use_parallel = .False.
            requested_threads = 1
            ! Adaptive step integrations - 
            error_digits = 12
            learning_rate = uno
            ! The primary -
            mass_primary = cero
            radius_primary = cero
            ! Rotation - 
            lambda_kep = cero
            asteroid_rotational_period = cero
            ! Moons - 
            use_moonsfile = .False.
            moonsfile = ""
            ! Particles - 
            use_particlesfile = .False.
            particlesfile = ""
            ! Extra forces/effects - 
            use_stokes = .False.
            stokes_a_damping_time = infinity
            stokes_e_damping_time = infinity
            stokes_charac_time = cero
            use_naive_stokes = .False.
            drag_coefficient = cero
            drag_charac_time = cero
            use_omega_damping = .False.
            omega_linear_damping_time = infinity
            omega_exp_damping_time = infinity
            omega_exp_poly_A = cero
            omega_exp_poly_B = cero
            omega_charac_time = cero
            mass_exp_damping_time = infinity
            J2_coefficient = cero
            ! Binned forces/effects - [NOT AVAILABLE YET. STILL UNDER DEVELOPMENT]
            use_self_gravity = .False.
            Norder_self_gravity = 3
            Nbins = 0
            binning_method = 1
            rmin_bins = - uno ! -1 means first particle (initial)
            rmax_bins = - uno ! -1 means last particle (initial)
            ! Conditions for Collision/Escape - 
            min_distance = cero
            max_distance = infinity
            use_merge = .False.
            ! Manual |(t)imes omega(t) mass_add(t)| file -
            use_tomfile = .False.
            tomfile = ""
            ! Output -
            use_screen = .False.
            use_datafile = .False.
            datafile = ""
            use_multiple_outputs = .False.
            multfile = ""
            use_chaosfile = .False.
            chaosfile = ""
            use_datascreen = .True.
            use_percentage = .False.
            use_elements_output = .True.
            use_potential_map = .False.
            mapfile = ""
            map_grid_size_x = 500
            map_grid_size_y = 500
            map_min_x = -300
            map_max_x = 300
            map_min_y = -300
            map_max_y = 300

            ! < Number of bodies >
            Nboulders = 0
            Nmoons = 0
            Nparticles = 0
        end subroutine init_default_config
        
        ! 1. Read config file
        subroutine read_config_file(file_name, existe)
            implicit none
            character(LEN=*), intent(in) :: file_name
            logical, intent(inout) :: existe
            character(256) :: line, param_str, value_str
            character(1) :: auxch1
            character(2) :: auxch2
            character(15) :: auxch15
            integer :: colonPos
            integer :: commentPos
            integer(kind=4) :: nlines, io, len_val, j

            if (.not. existe) return

            inquire(file=trim(file_name), exist=existe)
            if (.not. existe) then
                if (use_screen) write (*,*) "File ", trim(file_name), " does not exists."
                return
            end if

            nlines = 0
            open (unit=10, file=trim(file_name), status='old', action='read')
            do
                ! Read the line from the file
                read (10, '(A)', END=98) line
                nlines = nlines + 1

                if ((len_trim(line) == 0) .or. (auxch15(:2) == "c ") .or. (auxch15(:2) == "! ")) cycle

                colonPos = index(line, ':') !this should be 50

                commentPos = index(line, '!')
                if (commentPos == 0) commentPos = len_trim(line)+1 !this should be 80 at least
                len_val = len_trim(adjustl(line(colonPos + 1 : commentPos - 1)))
                if (colonPos > 0) then
                    param_str = trim(adjustl(line(: colonPos)))
                    value_str = trim(adjustl(line(colonPos + 1 : commentPos - 1)))
                    auxch1 = to_lower(value_str(:1))
                    auxch2 = to_lower(value_str(:2))
                    auxch15 = param_str(:15)
                    select case (auxch15)
                        case("initial time fo")
                            read (value_str, *) initial_time
                        case("total integrati")
                            read (value_str, *) final_time
                        case("output time int")
                            read (value_str, *) output_timestep
                        case("number of outpu")
                            read (value_str, *) output_number
                        case("output distribu")
                            read (value_str, *) case_output_type
                        case("use parallel th")
                            read (value_str,'(i10)',iostat=aux_integer) requested_threads
                            if (aux_integer /= 0) then
                                if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                    use_parallel = .True.
                                    requested_threads = -1
                                else
                                    use_parallel = .False.
                                    requested_threads = 1
                              endif
                            else
                                use_parallel = .True.
                            end if
                        case("precision (digi")
                            read (value_str, *) error_digits
                        case("learning rate")
                            read (value_str, *) learning_rate
                        case("mass of primary")
                            read (value_str, *) mass_primary
                        case("radius of prima")
                            read (value_str, *) radius_primary
                        case("ratio of spin t")
                            read (value_str, *) lambda_kep
                        case("rotational peri")
                            read (value_str, *) asteroid_rotational_period
                        case("moons input fil")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                use_moonsfile = .False.
                                moonsfile = ""
                            else if (len(trim(value_str)) > 0) then
                                use_moonsfile = .True.
                                moonsfile = trim(value_str)
                            else
                                use_moonsfile = .False.
                                moonsfile = ""
                            end if
                        case("particles input")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                use_particlesfile = .False.
                                particlesfile = ""
                            else if (len(trim(value_str)) > 0) then
                                use_particlesfile = .True.
                                particlesfile = trim(value_str)
                            else
                                use_particlesfile = .False.
                                particlesfile = ""
                            end if
                        case("include stokes-")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                use_stokes = .True.
                            else
                                use_stokes = .False.
                            end if
                        case("a damping chara")
                            read (value_str, *) stokes_a_damping_time
                        case("e damping chara")
                            read (value_str, *) stokes_e_damping_time
                        case("F damping chara")
                            read (value_str, *) stokes_charac_time
                        case("include naive-s")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                use_naive_stokes = .True.
                            else
                                use_naive_stokes = .False.
                            end if
                        case("drag force coef")
                            read (value_str, *) drag_coefficient
                        case("drag force char")
                            read (value_str, *) drag_charac_time
                        case("include rotatio")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                use_omega_damping = .True.
                            else
                                use_omega_damping = .False.
                            end if
                        case("rotation linear")
                            read (value_str, *) omega_linear_damping_time
                        case("rotation expone")
                            read (value_str, *) omega_exp_damping_time
                        case("rotation A poly")
                            read (value_str, *) omega_exp_poly_A
                        case("rotation B poly")
                            read (value_str, *) omega_exp_poly_B
                        case("omega damping c")
                            read (value_str, *) omega_charac_time
                        case("mass exponentia")
                            read (value_str, *) mass_exp_damping_time
                        case("geo-potential J")
                            read (value_str, *) J2_coefficient
                        case("include self-gr")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                use_self_gravity = .True.
                            else
                                use_self_gravity = .False.
                            end if
                        case("max Legendre ex")
                            read (value_str, *) Norder_self_gravity
                        case("total bins used")
                            read (value_str, *) Nbins
                        case("binning method ")
                            read (value_str, *) binning_method
                        case("inner bin edge")
                            read (value_str, *) rmin_bins
                        case("outer bin edge")
                            read (value_str, *) rmax_bins
                        case("min distance fr")
                            read (value_str, *) min_distance
                        case("max distance fr")
                            read (value_str, *) max_distance
                        case("merge collision")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                use_merge = .True.
                            else
                                use_merge = .False.
                            end if
                        case("input time-omeg")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                use_tomfile = .False.
                                tomfile = ""
                            else if (len(trim(value_str)) > 0) then
                                use_tomfile = .True.
                                tomfile = trim(value_str)
                            else
                                use_tomfile = .False.
                                tomfile = ""
                            end if
                        case("information on")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                use_screen = .True.
                            else
                                use_screen = .False.
                            end if
                        case("general output")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                use_datafile = .False.
                                datafile = ""
                            else
                                use_datafile = .True.
                                datafile = trim(value_str)
                            end if
                        case("individual file")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                use_multiple_outputs = .False.
                                multfile = ""
                            else
                                use_multiple_outputs = .True.
                                multfile = trim(value_str)
                            end if
                        case("chaos indicator")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                use_chaosfile = .False.
                                chaosfile = ""
                            else
                                use_chaosfile = .True.
                                chaosfile = trim(value_str)
                                use_chaos = .True.
                            end if
                        case("output data on")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                use_datascreen = .True.
                                use_percentage = .False.
                            else if (auxch1 == "%") then
                                use_datascreen = .False.
                                use_percentage = .True.
                            else
                                use_datascreen = .False.
                                use_percentage = .False.
                            end if
                        case("output variable")
                            if (auxch1 == "c") then
                                use_elements_output = .False.
                            else
                                use_elements_output = .True.
                            end if
                        case("create map file")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                use_potential_map = .False.
                                mapfile = ""
                            else
                                use_potential_map = .True.
                                mapfile = trim(value_str)
                            end if
                        case("number x cells")
                            read (value_str, *) map_grid_size_x
                        case("number y cells")
                            read (value_str, *) map_grid_size_y
                        case("lower x bound (")
                            read (value_str, *) map_min_x
                        case("upper x bound (")
                            read (value_str, *) map_max_x
                        case("lower y bound (")
                            read (value_str, *) map_min_y
                        case("upper y bound (")
                            read (value_str, *) map_max_y
                        case default
                            write (*,*) "WARNING: Parámetro no reconocido: ", trim(param_str)
                    end select
                else  ! Read boulders, moons or particles
                    param_str = trim(adjustl(line)) 
                    if (param_str(1:7) .eq. "mass_m0") then  ! Leeemos boulders
                        aux_integer = nlines  ! Save the nlines where to read from (later)
                        Nboulders = 0 ! This will be the boulder count
                        io = 0  ! This is be the reading status
                        do while (io == 0)
                            nlines = nlines + 1  ! Don't forget to count the lines....
                            read (10, '(A)') line
                            if ((line(:2) == "c ") .or. (line(:2) == "! ")) cycle  ! Skip commented
                            if (io == 0) Nboulders = Nboulders + 1
                        end do
                        ! Alocate
                        call allocate_asteriod_arrays(Nboulders)
                        ! Backspace to read again
                        do j = aux_integer, nlines
                            backspace (10)
                        end do
                        ! Read again and store
                        j = 1  ! This is be the amount of boulders read
                        do while (j <= Nboulders)
                            read (10, '(A)') line
                            if ((line(:2) == "c ") .or. (line(:2) == "! ")) cycle  ! Skip commented
                            read (10, *, iostat=io) mu_from_primary(j), radius_ast_arr(j), theta_from_primary(j)
                            if (io /= 0) then
                                write (*,*) "ERROR: Al leer boulder:", j
                                stop 1
                            end if
                            j = j + 1
                        end do
                    end if
                    if (param_str(1:9) .eq. "mass_mAst") then  ! Leeemos las lunas
                        aux_integer = nlines  ! Save the nlines where to read from (later)
                        Nmoons = 0 ! This will be the moon count
                        io = 0  ! This is be the reading status
                        do while (io == 0)
                            nlines = nlines + 1  ! Don't forget to count the lines....
                            read (10, '(A)') line
                            if ((line(:2) == "c ") .or. (line(:2) == "! ")) cycle  ! Skip commented
                            if (io == 0) Nmoons = Nmoons + 1
                        end do
                        ! Alocate
                        call allocate_moons_arrays(Nmoons)
                        ! Backspace to read again
                        do j = aux_integer, nlines
                            backspace (10)
                        end do
                        ! Read again and store
                        j = 1  ! This is be the amount of moons read
                        do while (j <= Nmoons)
                            read (10, '(A)') line
                            if ((line(:2) == "c ") .or. (line(:2) == "! ")) cycle  ! Skip commented
                            read (10, *, iostat=io) &
                                & moons_mass_ratios(j), &
                                & moons_elem(j,1), moons_elem(j,2), &
                                & moons_elem(j,3), moons_elem(j,4), &
                                & moons_radii(j), &
                                & moons_MMR(j)
                            if (io /= 0) then
                                write (*,*) "ERROR: Al leer luna:", j
                                stop 1
                            end if
                            j = j + 1
                        end do
                    end if
                    if (param_str(1:5) .eq. "a(km)") then  ! Leeemos las partículas
                        aux_integer = nlines  ! Save the nlines where to read from (later)
                        Nparticles = 0 ! This will be the particles count
                        io = 0  ! This is be the reading status
                        do while (io == 0)
                            nlines = nlines + 1  ! Don't forget to count the lines....
                            read (10, '(A)') line
                            if ((line(:2) == "c ") .or. (line(:2) == "! ")) cycle  ! Skip commented
                            if (io == 0) Nparticles = Nparticles + 1
                        end do
                        ! Alocate
                        call allocate_particles_arrays(Nparticles)
                        ! Backspace to read again
                        do j = aux_integer, nlines
                            backspace (10)
                        end do
                        ! Read again and store
                        j = 1  ! This is be the amount of particles read
                        do while (j <= Nparticles)
                            read (10, '(A)') line
                            if ((line(:2) == "c ") .or. (line(:2) == "! ")) cycle  ! Skip commented
                            read (10, *, iostat=io) &
                                & particles_elem(j,1), particles_elem(j,2), &
                                & particles_elem(j,3), particles_elem(j,4), &
                                & particles_MMR(j)
                            if (io /= 0) then
                                write (*,*) "ERROR: Al leer partícula:", j
                                stop 1
                            end if
                            j = j + 1
                        end do
                    end if
                end if
            end do
            98 close (10)
        end subroutine read_config_file

        ! 2. Init default hidden (and command line) default parameters
        subroutine init_default_parameters()
            implicit none

            ! Simulation
            simulation_number = 1

            ! Times
            min_timestep = cero
            
            ! OmenMP
            !$ compiled_with_openmp = .True. !! Parece comentado, pero es asi

            ! Forces
            !! Viscosity [unsed yet]
            use_viscosity = .False.
            viscosity = - uno


            ! Extras, fuera de Config.ini. Se editan solo acá
            !! Update chaos at every output timestep, instead of waiting to the end of the run
            use_update_chaos = .True.
            !! Flush output at every output timestep (even though the buffer is not full) 
            use_flush_output = .True.
        end subroutine init_default_parameters

        ! 3 Leer entradas por linea de comandos
        subroutine load_command_line_arguments(existe_configfile)
            implicit none
            logical :: existe_configfile

            ! Leemos de la línea de comandos
            arguments_number = command_argument_count()
            !!!!! PARTÍCULA (en caso de entrada por terminal)
            aux_logical = .False. ! Leí los parámetros numéricos?
            aux_integer = 0
            if (arguments_number > 0) then
                do i = 1, arguments_number
                    if (aux_integer /= 0) then
                        aux_integer = aux_integer - 1
                        cycle
                    end if
                    call get_command_argument(i, aux_character30)
                    select case (trim(aux_character30))
                    case ("-nsim")
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) simulation_number
                        aux_integer = 1
                    case ("-datafile")
                        use_datafile = .True.
                        call get_command_argument(i+1, datafile)
                        aux_integer = 1
                    case ("--nodataf")
                        use_datafile = .False.
                        datafile = ""
                    case ("-chaosfile")
                        use_chaosfile = .True.
                        use_chaos = .True.
                        call get_command_argument(i+1, chaosfile)
                        aux_integer = 1
                    case ("--nochaosf")
                        use_chaosfile = .False.
                        chaosfile = ""
                    case ("--screen")
                        use_screen = .True.
                    case ("--noscreen")
                        use_screen = .False.
                    case ("--perc")
                        use_percentage = .True.
                    case ( "--noperc")
                        use_percentage = .False.
                    case ("--datascr")
                        use_datascreen = .True.
                    case ( "--nodatascr")
                        use_datascreen = .False.
                    case ("-multifile")
                        use_multiple_outputs = .True.
                        call get_command_argument(i+1, multfile)
                        aux_integer = 1
                    case ("--nomultif")
                        use_multiple_outputs = .False.
                        multfile = ""
                    case ("-mapfile")
                        use_potential_map = .True.
                        call get_command_argument(i+1, mapfile)
                        aux_integer = 1
                    case ("--nomapf")
                        use_potential_map = .False.
                        mapfile = ""
                    case ("--elem")
                        use_elements_output = .True.
                    case ("--noelem")
                        use_elements_output = .False.
                    case ("-tomfile")
                        use_tomfile = .True.
                        call get_command_argument(i+1, tomfile)
                        aux_integer = 1
                    case ("--notomfile")
                        use_tomfile = .False.
                        tomfile = ""
                    case ("-moonfile")
                        use_moonsfile = .True.
                        call get_command_argument(i+1, moonsfile)
                        aux_integer = 1
                    case ("--moonfile")
                        use_moonsfile = .False.
                        moonsfile = ""
                    case ("-partfile")
                        use_particlesfile = .True.
                        call get_command_argument(i+1, particlesfile)
                        aux_integer = 1
                    case ("--nopartfile")
                        use_particlesfile = .False.
                        particlesfile = ""
                    case ("--noconfig")
                        existe_configfile = .False.
                    case ("--merge")
                        use_merge = .True.
                    case ("--nomerge")
                        use_merge = .False.
                    case ("-parallel")
                        use_parallel = .True.
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) requested_threads
                        aux_integer = 1
                    case ("--parallel")
                        use_parallel = .True.
                        requested_threads = -1
                    case ("--noparallel")
                        use_parallel = .False.
                        requested_threads = 1
                    case ("--help")
                        call get_command_argument(0, aux_character30)
                        write (*,*) "Uso: " // trim(aux_character30) // " <ea> <ee> <eM> <ew> <eR> [args]"
                        write (*,*) "    ea  : Elemento a de la partícula (km)"
                        write (*,*) "    ee  : Elemento e de la partícula"
                        write (*,*) "    eM  : Elemento M de la partícula (deg)"
                        write (*,*) "    ew  : Elemento w de la partícula (deg)"
                        write (*,*) "    eR  : Elemento R de la partícula [Optional]"
                        write (*,*) "    -mpart       : Masa de la partícula individual"
                        write (*,*) "    -nsim        : Asignar como número de simulación al 'int' que sigue"
                        write (*,*) "    -datafile    : Guardar datos en el archivo que sigue"
                        write (*,*) "    --nodataf    : No guardar datos"
                        write (*,*) "    -chaosfile   : Guardar caos en el archivo que sigue"
                        write (*,*) "    --nochaosf   : No guardar caos"
                        write (*,*) "    --screen     : Imprimir información en pantalla"
                        write (*,*) "    --noscreen   : No imprimir en pantalla"
                        write (*,*) "    --perc       : Imprimir porcentaje de integración"
                        write (*,*) "    --noperc     : No imprimir porcentaje de integración"
                        write (*,*) "    --datascr    : Imprimir datos en pantalla"
                        write (*,*) "    --nodatascr  : No imprimir datos en pantalla"
                        write (*,*) "    -multifile   : Guardar datos en archivos individuales"
                        write (*,*) "    --nomultif   : No guardar datos en archivos individuales"
                        write (*,*) "    -mapfile     : Guardar mapas de potencial en el archivo que sigue"
                        write (*,*) "    --nomapf     : No guardar mapas de potencial"
                        write (*,*) "    --implicit   : Usar método implícito (integra) [default]"
                        write (*,*) "    --explicit   : Usar método explícito (cos, sen)"
                        write (*,*) "    --elem       : Imprimir elementos orbitales (lunas/partículas) [default]"
                        write (*,*) "    --noelem     : Imprimir coordenadas baricéntricas"
                        write (*,*) "    -tomfile     : Utilizar archivo de (t)iempos|omega|masa que sigue"
                        write (*,*) "    --notomfile  : No utilizar archivo de (t)iempos|omega|masa"
                        write (*,*) "    -moonfile    : Utilizar archivo de lunas que sigue"
                        write (*,*) "    --nomoonfile : No utilizar archivo de lunas"
                        write (*,*) "    -partfile    : Utilizar archivo de partículas que sigue"
                        write (*,*) "    --nopartfile : No utilizar archivo de partículas"
                        write (*,*) "    --noconfig   : No leer archivo de configuración"
                        write (*,*) "    --merge      : Incluir asociar colisiones de lunas/partículas a asteroide"
                        write (*,*) "    --nomerge    : No asociar colisiones de lunas/partículas"
                        write (*,*) "    -parallel    : Paralelizar usando la cantida de thread que sique"                        
                        write (*,*) "    --parallel   : Paralelizar usando todos los threads disponibles"
                        write (*,*) "    --noparallel : No usar paralelización para lunas/partículas"
                        write (*,*) "    --help       : Mostrar esta ayuda"
                        stop 0
                    case default  ! Si no es un argumento reconocido...
                        is_number = .False.
                        call get_command_argument(i, aux_character30)
                        do j = 0, 9
                            if (aux_character30(1:1) == char(48 + j)) then ! check if it's a number
                                is_number = .True.
                                exit
                            end if
                        end do
                        if (.not. is_number) then ! No es un número
                            write (*, '(A, I0, A, A)') "ERROR: No se reconoce el argumento ", i, ": ", trim(aux_character30)
                            call get_command_argument(0, aux_character30)
                            write (*,*) "Para ayuda, ejecute: ", trim(aux_character30), " --help"
                            write (*,*) "Saliendo."
                            stop 1
                        end if
                        if (.not. aux_logical)  then ! No leí los parámetros aún
                            if (arguments_number < i+3) then
                                write (*,*) "ERROR: No se pudo leer los elementos orbitales."
                                write (*,'(A,I0,A)') "        Faltan ", 4-arguments_number, " elementos."
                                write (*,*) "Saliendo."
                                stop 1
                            else ! Leo los argumentos numéricos. Considero que es 1 sola partícula
                                Nparticles = 1
                                call allocate_particles_arrays(Nparticles)
                                call get_command_argument(i, aux_character20)
                                read (aux_character20,*) particles_elem(1,1)
                                call get_command_argument(i+1, aux_character20)
                                read (aux_character20,*) particles_elem(1,2)
                                call get_command_argument(i+2, aux_character20)
                                read (aux_character20,*) particles_elem(1,3)
                                call get_command_argument(i+3, aux_character20)
                                read (aux_character20,*) particles_elem(1,4)
                                aux_integer = 3
                                aux_logical = .True.
                                particles_MMR(1) = cero
                            end if
                        else ! Ya leí los numéricos. Falta leer eR
                            call get_command_argument(i, aux_character20)
                            read (aux_character20,*) particles_MMR(1)
                        end if
                    end select
                end do
            else if (existe_configfile) then
                print*, "WARNING: No hay argumentos de entrada."
                print*, "¿Está seguro que desea continuar? [y/N]"
                read (*,*) aux_character1
                if (index("YySs", aux_character1) == 0) then
                    print*, "Exiting."
                    stop
                endif
            end if  
        end subroutine load_command_line_arguments

        ! 4. Setear los parámetros derivados
        subroutine set_derived_parameters()
            implicit none


            !! Times
            if ((case_output_type < 0) .or. (case_output_type > 2)) then
                write (*,*) "ERROR: Tipo de distribución para tiempos de salida no reconocido."
                stop 1
            end if


            !! Boulders
            if (Nboulders < 0) then
                write (*,*) "ERROR: Nboulders < 0"
                stop 1
            end if
            use_boulders = (Nboulders > 0)

            ! -----------------------------------------------------------------
            ! Acá hay que configurar moons y particles en caso de moons sin masa
            ! -----------------------------------------------------------------

            !! Moons
            if (Nmoons < 0) then
                write (*,*) "ERROR: Nmoons < 0"
                stop 1
            end if
            use_moons = (Nmoons > 0)

            !! Moons
            if (Nparticles < 0) then
                write (*,*) "ERROR: Nparticles < 0"
                stop 1
            end if
            use_particles = (Nparticles > 0)
                       
            
            !! Forces
            !!! stokes_a_damping_time y stokes_e_damping_time
            if ((abs(stokes_a_damping_time) < tini) .or. .not. use_stokes) stokes_a_damping_time = cero
            if ((abs(stokes_e_damping_time) < tini) .or. .not. use_stokes) stokes_e_damping_time = cero
            if ((abs(stokes_charac_time) < tini) .or. .not. use_stokes) stokes_charac_time = cero
            if (stokes_charac_time .le. cero) stokes_charac_time = infinity
            if ((abs(stokes_a_damping_time) < tini) .and. (abs(stokes_e_damping_time) < tini)) use_stokes = .False.

            !!! Naive-stokes
            if ((abs(drag_coefficient) < tini) .or. .not. use_naive_stokes) drag_coefficient = cero
            if ((abs(drag_charac_time) < tini) .or. .not. use_naive_stokes) drag_charac_time = cero
            if (drag_charac_time .le. cero) drag_charac_time = infinity
            if (abs(drag_coefficient) < tini) use_naive_stokes = .False.

            !!! tau_m y tau_o
            if ((abs(omega_charac_time) < tini) .or. .not. use_omega_damping) omega_charac_time = cero
            if (omega_charac_time .le. cero) omega_charac_time = infinity
            if ((abs(omega_linear_damping_time) < tini) .or. .not. use_omega_damping) omega_linear_damping_time = infinity
            if ((abs(omega_exp_damping_time) < tini) .or. .not. use_omega_damping) omega_exp_damping_time = infinity
            if ((abs(omega_exp_poly_A) < tini) .or. .not. use_omega_damping) omega_exp_poly_A = cero
            if ((abs(omega_exp_poly_B) < tini) .or. .not. use_omega_damping) omega_exp_poly_B = cero
            if (.not. ((abs(omega_linear_damping_time) > tini) .and. (omega_linear_damping_time < infinity)) .and. &
              & .not. ((abs(omega_exp_damping_time) > tini) .and. (omega_exp_damping_time < infinity)) .and. &
              & .not. ((abs(omega_exp_poly_A) > tini) .and. (abs(omega_exp_poly_B) > tini))) &
              & use_omega_damping = .False.
            if (abs(mass_exp_damping_time) < tini) mass_exp_damping_time = infinity

            !!! Geo-Potential
            if (abs(J2_coefficient) > tini) use_J2 = .True.
            
            ! Self-Gravity or Viscosity
            if (use_self_gravity .or. use_viscosity) then
                use_bins = .True.
                !! Check
                if (binning_method < 1 .or. binning_method > 3) then
                    write (*,*) "ERROR: binning_method must be 1 (equal dr), 2 (equal dA) or 3 (equal Npart)."
                    stop 1
                end if
                if (Nbins < 2) then
                    write (*,*) "ERROR: No se puede crear menos de 2 bines."
                    stop 1
                end if
                !! Set default update values
                update_rmin_bins = rmin_bins < - tini
                update_rmax_bins = rmax_bins < - tini
                update_bins = update_rmin_bins .or. update_rmax_bins .or. binning_method == 3
            end if         

            !! Error
            if (error_digits < 1) then
                write (*,*) "ERROR: El número de dígitos de presición debe ser mayor que 0."
                stop 1
            end if
            error_tolerance = 10.d0**(-error_digits)
            
            !! Cuerpo central
            !!!! Masa
            ! if (mass_primary <= 0) then
            !     write (*,*) "ERROR: La masa del cuerpo central debe ser positiva."
            !     stop 1
            ! end if
            !!!! Radio
            if (radius_primary <= 0) then
                write (*,*) "ERROR: El radio del cuerpo central debe ser positivo."
                stop 1
            end if

            ! Output
            if (use_datascreen .and. use_percentage) then
                write (*,*) "ERROR: No se puede imprimir porcentaje y datos en pantalla."
                stop 1
            end if
            
            ! Auxiliares lógicos
            use_chaos = use_chaosfile .or. use_chaos
            use_elements = use_elements_output .or. use_chaos
            if (.not. use_chaosfile) use_update_chaos = .False.
            if (.not. use_datafile) use_flush_output = .False.

            ! Parallel
            !$ compiled_with_openmp = .True. !! Parece comentado, pero es asi
            if (use_parallel) then
                if ((requested_threads == 0) .or. (requested_threads == 1)) then
                    requested_threads = 0
                    my_threads = 1
                    use_parallel = .False.
                else if (.not. compiled_with_openmp) then
                    write (*,*) "ERROR: No se puede usar paralelización sin OpenMP."
                    stop 1
                end if
                !$ available_threads = OMP_GET_MAX_THREADS()
                !$ if (requested_threads .eq. -1) requested_threads = available_threads
                !$ my_threads = min(available_threads, max(requested_threads,1))
                !$ call OMP_SET_NUM_THREADS(my_threads)
            else
                requested_threads = 0
                my_threads = 1
            end if

        end subroutine set_derived_parameters
        
        ! 5.x Subroutines to allocate arrays

        ! 5.1 Alocatar arrays asteroide
        subroutine allocate_asteriod_arrays(Nboulders)
            implicit none
            integer(kind=4), intent(in) :: Nboulders

            if (allocated(mass_ast_arr)) then
                write (*,*) "ERROR: Number of boulder already defined."
                STOP 1
            end if

            allocate(mass_ast_arr(0:Nboulders))
            allocate(Gmass_ast_arr(0:Nboulders))
            allocate(mu_ast_arr(0:Nboulders))
            allocate(radius_ast_arr(0:Nboulders))
            allocate(pos_ast_arr(0:Nboulders,2))
            allocate(vel_ast_arr(0:Nboulders,2))
            allocate(acc_ast_arr(0:Nboulders,2))
            allocate(theta_ast_arr(0:Nboulders))
            allocate(dist_ast_arr(0:Nboulders))
            allocate(inertia_ast_arr(0:Nboulders))
            
            if (Nboulders > 0) then
                allocate(theta_from_primary(Nboulders))
                allocate(mu_from_primary(Nboulders))
                allocate(pos_from_primary(Nboulders,2))
                allocate(vel_from_primary(Nboulders,2))
                allocate(acc_from_primary(Nboulders,2))
            end if
        end subroutine allocate_asteriod_arrays

        ! 5.2 Alocatar arrays lunas
        subroutine allocate_moons_arrays(Nmoons)
            implicit none
            integer(kind=4), intent(in) :: Nmoons
            
            if (Nmoons > 0) then
                if (allocated(moons_index)) then
                    write (*,*) "ERROR: Number of moons already defined."
                    write (*,*) "        Use moons file or config file; not both."
                    STOP 1
                end if

                allocate(moons_index(1:Nmoons))
                allocate(sorted_moons_index(1:Nmoons))
                allocate(moons_mass(1:Nmoons))
                allocate(moons_mass_ratios(1:Nmoons))
                allocate(moons_radii(1:Nmoons))
                allocate(moons_coord(1:Nmoons,4))
                allocate(moons_elem(1:Nmoons,4))
                allocate(moons_acc(1:Nmoons,2))
                allocate(moons_MMR(1:Nmoons))
                allocate(moons_dist(1:Nmoons))
            end if
        end subroutine allocate_moons_arrays

        ! 5.3 Alocatar particulas
        subroutine allocate_particles_arrays(Nparticles)
            implicit none
            integer(kind=4), intent(in) :: Nparticles

            if (Nparticles > 0) then
                if (allocated(particles_index)) then
                    write (*,*) "ERROR: Number of particles already defined."
                    write (*,*) "        Use particles file or config file; not both."
                    STOP 1
                end if
                !! Parámetros y auxiliares
                allocate(particles_index(1:Nparticles))
                allocate(sorted_particles_index(1:Nparticles))
                allocate(particles_coord(1:Nparticles,4))
                allocate(particles_elem(1:Nparticles,4))
                allocate(particles_acc(1:Nparticles,2))
                allocate(particles_MMR(1:Nparticles))
                allocate(particles_dist(1:Nparticles))
            end if
                ! allocate(particles_outcome(1:Nparticles))
                ! allocate(ij_to_swap(1:Nparticles,2))
                ! if (use_chaos) then
                !     allocate(particles_times(1:Nparticles))
                !     allocate(particles_min_a(1:Nparticles), particles_max_a(1:Nparticles))
                !     allocate(particles_min_e(1:Nparticles), particles_max_e(1:Nparticles))
                ! end if
                ! allocate(particles_hexit(0:Nparticles))
                ! particles_hexitptr => particles_hexit(0) ! Puntero a Hard Exit (NO TOCAR) !!! Ya está en default
                ! if (use_bins) allocate(particles_bins(1:Nparticles))  ! Unused if use_bins is False
        end subroutine allocate_particles_arrays

        ! 5.4 Liberar arrays asteroide
        subroutine free_asteroid_arrays()
            implicit none

            deallocate(mass_ast_arr)
            deallocate(Gmass_ast_arr)
            deallocate(mu_ast_arr)
            deallocate(radius_ast_arr)
            deallocate(pos_ast_arr)
            deallocate(vel_ast_arr)
            deallocate(acc_ast_arr)
            deallocate(theta_ast_arr)
            deallocate(dist_ast_arr)
            deallocate(inertia_ast_arr)
            
            if (use_boulders) then
                deallocate(theta_from_primary)
                deallocate(mu_from_primary)
                deallocate(pos_from_primary)
                deallocate(vel_from_primary)
                deallocate(acc_from_primary)
            end if
        end subroutine free_asteroid_arrays

        ! 5.5 Liberar arrays particles
        subroutine free_moons_arrays()
            implicit none

            if (allocated(moons_index)) then
                deallocate(moons_index)
                deallocate(sorted_moons_index)
                deallocate(moons_mass)
                deallocate(moons_mass_ratios)
                deallocate(moons_radii)
                deallocate(moons_coord)
                deallocate(moons_elem)
                deallocate(moons_acc)
                deallocate(moons_MMR)
                deallocate(moons_dist)
            end if
        end subroutine free_moons_arrays

        ! 5.6 Liberar arrays particles
        subroutine free_particles_arrays()
            implicit none

            if (allocated(particles_index)) then
                deallocate(particles_index)
                deallocate(particles_coord)
                deallocate(particles_elem)
                deallocate(particles_acc)
                deallocate(particles_MMR)
                deallocate(particles_dist)
                deallocate(sorted_particles_index)
                ! if (use_chaos) then
                !     deallocate(particles_times)
                !     deallocate(particles_min_a)
                !     deallocate(particles_max_a)
                !     deallocate(particles_min_e)
                !     deallocate(particles_max_e)
                ! end if
                ! deallocate(particles_outcome)
                ! deallocate(ij_to_swap)
                ! deallocate(particles_hexit)
                ! if (use_bins) deallocate(particles_bins)  ! Unused if use_bins is False
            end if
        end subroutine free_particles_arrays

        ! 7. Leer un archivo con datos en columnas
        subroutine read_columns_file(file_name, values_arr, method)
            implicit none
            character(LEN=*), intent(in) :: file_name
            real(kind=8), dimension(:,:), allocatable, intent(out) :: values_arr
            integer(kind=4), optional :: method
            integer(kind=4), parameter :: MAX_COLS = 8, MAX_ROWS = 10000
            real(kind=8), dimension(:,:), allocatable :: aux_real_arr
            integer(kind=4) :: ncols, nrows, i, j, io, my_method
            real(kind=8) :: aux_real
            character(260) :: auxstr
            logical :: existe
           
            if (present(method)) then
                my_method = method
            else
                my_method = 0
            end if

            inquire (file=trim(file_name), exist=existe)
            if (.not. existe) then
                write (*,*) "ERROR: No se encontró el archivo: ", trim(file_name)
                stop 1
            end if

            open (unit=20, file=trim(file_name), status="old", action="read")

            !! Count number of columns
            ncols = 0
            read (20, '(A)') auxstr
            do i = 1,MAX_COLS
                io = 0
                read (auxstr, *, iostat=io) (aux_real, j=1,i)
                if (io .ne. 0) exit
            end do
            if (io .eq. 0) then 
                ncols = i
            else
                ncols = i - 1
            end if
            
            rewind (20) ! Go to the beginning of the file
            
            if (my_method .eq. 0) then

                !! Count number of rows
                nrows = 0
                do ! Count number of (valid) lines 
                    read (20, *, iostat=io) aux_real
                    if (is_iostat_end(io)) exit
                    nrows = nrows + 1
                end do

                ! Allocate arrays
                allocate (values_arr(nrows, ncols))
                rewind (20) ! Go to the beginning of the file
                do i = 1, nrows
                    read (20, *) (values_arr(i,j), j=1,ncols)
                end do

            else
                ! Allocate auxiliar array
                allocate (aux_real_arr(MAX_ROWS, ncols))
                nrows = 0
                do i = 1, MAX_ROWS
                    read (20, *, iostat=io) (aux_real_arr(i,j), j=1,ncols)
                    if (is_iostat_end(io)) exit
                    nrows = nrows + 1
                end do
                ! Allocate arrays
                allocate (values_arr(nrows, ncols))
                values_arr = aux_real_arr(1:nrows, 1:ncols)
                deallocate (aux_real_arr)
            end if
            close (20)
        end subroutine read_columns_file

        ! 8. Leer el archivo tomfile
        subroutine read_tomfile(t0, tf, t_TOM, domega_TOM, dmass_TOM, file_tout)
            implicit none
            real(kind=8), intent(in) :: t0, tf
            real(kind=8), dimension(:), allocatable, intent(out) :: t_TOM, domega_TOM, dmass_TOM
            character(LEN=*), intent(in) :: file_tout
            integer(kind=4) :: n_TOM
            integer(kind=4) :: i, j, io
            integer(kind=4) :: ncols
            real(kind=8) :: t_aux
            character(80) :: auxstr
            logical :: existe

            n_TOM = 2
            
            inquire (file=trim(file_tout), exist=existe)
            if (.not. existe) then
                write (*,*) "ERROR: No se encontró el archivo TOM: ", trim(file_tout)
                stop 1
            end if

            open (unit=30, file=file_tout, status="old", action="read")

            !! Count number of columns
            ncols = 0
            read (30, '(A)') auxstr
            do i = 1,3   ! The very maximum that the string can contain: 3
                io = 0
                read (auxstr, *, iostat=io) (t_aux, j=1,i)
                if (io .ne. 0) exit
            end do
            if (io .eq. 0) then 
                ncols = i
            else
                ncols = i - 1
            end if
            
            rewind (30) ! Go to the beginning of the file
            do ! Count number of (valid) lines 
                read (30, *, iostat=io) t_aux
                if ((io /= 0) .or. (t_aux > tf)) exit
                if (t_aux < t0) cycle
                n_TOM = n_TOM + 1
            end do

            if (n_TOM == 2) then
                write (*,*) "ERROR: No se encontraron datos válidos en el archivo."
                stop 1
            end if
            
            ! Allocate arrays
            allocate (t_TOM(n_TOM))
            t_TOM = -1.d0
            if (ncols > 1) then
                allocate (domega_TOM(n_TOM))
                domega_TOM = -1.d0
            end if
            if (ncols > 2) then
                allocate (dmass_TOM(n_TOM))
                dmass_TOM = -1.d0
            end if
            rewind (30) ! Go to the beginning of the file

            ! Read data
            i = 2
            if (ncols == 1) then
                do
                    read (30, *, iostat=io) t_TOM(i)
                    if ((io /= 0) .or. (t_TOM(i) > tf)) exit
                    if (t_TOM(i) < t0) cycle
                    i = i + 1
                end do
            else if (ncols == 2) then
                do
                    read (30, *, iostat=io) t_TOM(i), domega_TOM(i)
                    if ((io /= 0) .or. (t_TOM(i) > tf)) exit
                    if (t_TOM(i) < t0) cycle
                    i = i + 1
                end do
            else 
                do
                    read (30, *, iostat=io) t_TOM(i), domega_TOM(i), dmass_TOM(i)
                    if ((io /= 0) .or. (t_TOM(i) > tf)) exit
                    if (t_TOM(i) < t0) cycle
                    i = i + 1
                end do
            end if

            close (30)
        end subroutine read_tomfile

        ! 11. Obtener a de corotación (respecto al asteroide)
        function get_a_corot(masstot, omega1) result(acorot)
            implicit none
            real(kind=8), intent(in) :: masstot, omega1
            real(kind=8) :: acorot

            acorot = (G * masstot / (omega1 * omega1))**(1/3.)
        end function get_a_corot

        ! 12. Obtener centro de masas de asteroide
        subroutine get_center_of_mass(m, rib, vib, mcm, rcm, vcm)
            implicit none
            real(kind=8), intent(in) :: m(0:), rib(0:,:), vib(0:,:)
            real(kind=8), intent(out) :: mcm, rcm(2), vcm(2)
            integer(kind=4) :: i

            rcm = cero
            vcm = cero
            do i = 0, Nboulders
                rcm = rcm + m(i) * rib(i,:)
                vcm = vcm + m(i) * vib(i,:)
            end do
            mcm = sum(m(0:Nboulders))
            rcm = rcm / mcm ! rcm = sum_i m_i * r_i / M
            vcm = vcm / mcm ! vcm = sum_i m_i * v_i / M
        end subroutine get_center_of_mass

        ! 13. Calcular omega a partir de posiciones y velocidades
        subroutine get_asteroid_from_boulders(m, rib, vib, rcm, vcm, omega, inertia)
            implicit none
            real(kind=8), intent(in) :: m(0:), rib(0:,:), vib(0:,:)
            real(kind=8), intent(out) :: rcm(2), vcm(2), omega, inertia
            real(kind=8) :: mcm, angmom
            real(kind=8) :: raux(2), vaux(2)

            call get_center_of_mass(m, rib, vib, mcm, rcm, vcm)
            ! Calculate angmom and inertia
            angmom = cero
            inertia = cero
            do i = 0, Nboulders
                raux = rib(i,:) - rcm
                vaux = vib(i,:) - vcm
                angmom = angmom + m(i) * cross2D(raux, vaux) ! Traslacional
                inertia = inertia + m(i) * sum(raux * raux) ! Rotacional
            end do
            !! Get Omega from L/I
            omega = angmom / inertia
        end subroutine get_asteroid_from_boulders

        ! 15. Setear min y max para caos
        subroutine free_initial_conditions()
            implicit none
            
            deallocate(particles_initial_conditions)
            if (use_boulders) deallocate(m0_and_boulders_initial_conditions)  
            deallocate(asteroid_initial_conditions)
        end subroutine free_initial_conditions

        ! 16. Intercambiar partículas j y i
        subroutine swap_particles(i, j, chaos)
            ! Initial conditions are not swapped
            implicit none
            integer(kind=4), intent(in) :: i, j
            logical, intent(in) :: chaos
            integer(kind=4) :: tmp_integer, tmp_integer2
            integer(kind=4) :: aux_integeri, aux_integerj
            real(kind=8), dimension(13) :: tmp_data
            real(kind=8), dimension(5) :: tmp_chaos_data
            real(kind=8), dimension(4) :: tmp_parameters_arr

            if (i == j) return

            tmp_integer = particles_index(i)
            tmp_data = (/particles_mass(i), particles_coord(i,:), particles_elem(i,:), &
                         particles_MMR(i), particles_dist(i)/)
            tmp_integer2 = particles_outcome(i)

            particles_index(i) = particles_index(j)
            particles_mass(i) = particles_mass(j)
            particles_coord(i,:) = particles_coord(j,:)
            particles_elem(i,:) = particles_elem(j,:)
            particles_MMR(i) = particles_MMR(j)
            particles_dist(i) = particles_dist(j)
            particles_acc(i,:) = particles_acc(j,:)
            particles_outcome(i) = particles_outcome(j)

            particles_index(j) = tmp_integer
            particles_mass(j) = tmp_data(1)
            particles_coord(j,:) = tmp_data(2:5)
            particles_elem(j,:) = tmp_data(6:9)
            particles_MMR(j) = tmp_data(10)
            particles_dist(j) = tmp_data(11)
            particles_acc(j,:) = tmp_data(12:13)
            particles_outcome(j) = tmp_integer2

            aux_integeri = first_particle + 4 * (i - 1)
            tmp_parameters_arr = (/parameters_arr(aux_integeri : aux_integeri+3)/)
            aux_integerj = first_particle + 4 * (j - 1)
            parameters_arr(aux_integeri : aux_integeri+3) = parameters_arr(aux_integerj : aux_integerj+3)
            parameters_arr(aux_integerj : aux_integerj+3) = tmp_parameters_arr

            if (chaos) then
                tmp_chaos_data = (/particles_times(i), particles_min_a(i), particles_max_a(i), &
                                  & particles_min_e(i), particles_max_e(i)/)
                particles_times(i) = particles_times(j)
                particles_min_a(i) = particles_min_a(j)
                particles_max_a(i) = particles_max_a(j)
                particles_min_e(i) = particles_min_e(j)
                particles_max_e(i) = particles_max_e(j)
                particles_times(j) = tmp_chaos_data(1)
                particles_min_a(j) = tmp_chaos_data(2)
                particles_max_a(j) = tmp_chaos_data(3)
                particles_min_e(j) = tmp_chaos_data(4)
                particles_max_e(j) = tmp_chaos_data(5)
            end if
            
            tmp_integer = particles_hexit(i)
            particles_hexit(i) = particles_hexit(j)
            particles_hexit(j) = tmp_integer
        end subroutine swap_particles
        
        
        ! 18 Hacer merger de masa (y ang mom) con asteroide
        subroutine merge_into_asteroid(mass, angmom)
            implicit none
            real(kind=8), intent(in) :: mass, angmom
            real(kind=8) :: growth

            ! Update asteroid mass
            growth = uno + (mass / asteroid_mass)
            mass_primary = mass_primary * growth
            mass_ast_arr = mass_ast_arr * growth
            asteroid_mass = asteroid_mass * growth
            asteroid_inertia = asteroid_inertia * growth
            Gasteroid_mass = G * asteroid_mass ! Update
            Gmass_ast_arr = G * mass_ast_arr ! Update
            
            if (use_boulders) then
                ! Update asteroid angular momentum
                asteroid_angmom = asteroid_angmom + angmom
                !! Update asteroid angular velocity
                asteroid_omega = asteroid_angmom / asteroid_inertia
                asteroid_omega2 = asteroid_omega * asteroid_omega
                !! Update asteroid angle correction
                asteroid_theta_correction = asteroid_theta - asteroid_omega * (time - initial_time)
                !! Update parameters array
                if (use_version_1) then
                    do i = 0, Nboulders
                        aux_integer = i * 4
                        parameters_arr(aux_integer+3) = - pos_ast_arr(i,2) * asteroid_omega
                        parameters_arr(aux_integer+4) = pos_ast_arr(i,1) * asteroid_omega
                    end do
                else
                    parameters_arr(2) = asteroid_omega
                end if
            end if
        end subroutine merge_into_asteroid
        
        ! 19.x Subroutines to accumulate (or not) mass and angmom

        ! 19.1 Acumular masa y momento angular de particulas
        subroutine accumulate_mass_and_angmom(i, acum_mass, acum_angmom)
            implicit none
            integer(kind=4), intent(in) :: i
            real(kind=8), intent(inout) :: acum_mass, acum_angmom

            acum_mass = acum_mass + particles_mass(i)
            acum_angmom = acum_angmom + particles_mass(i) * cross2D(particles_coord(i,1:2), particles_coord(i,3:4))
            particles_mass(i) = cero ! Set to 0,as it is now inside asteroid
        end subroutine accumulate_mass_and_angmom

        ! 19.2 NO acumular masa y momento de particulas
        subroutine do_not_accumulate_mass_and_angmom(i, acum_mass, acum_angmom)
            implicit none
            integer(kind=4), intent(in) :: i
            real(kind=8), intent(inout) :: acum_mass, acum_angmom
        end subroutine do_not_accumulate_mass_and_angmom

        ! 20. Guardar las condiciones iniciales de todo
        subroutine save_initial(part_ic, m0boul_ic, ast_ic)
            implicit none
            real(kind=8), dimension(:,:), allocatable, intent(inout) :: part_ic, m0boul_ic
            real(kind=8), dimension(:), allocatable, intent(inout) :: ast_ic

            ! Particles
            !! mass, a, e, M, w, MRR
            allocate(part_ic(1:Nparticles, 6))
            part_ic(:,1) = particles_mass
            part_ic(:,2) = particles_elem(:,1)
            part_ic(:,3) = particles_elem(:,2)
            part_ic(:,4) = particles_elem(:,3)
            part_ic(:,5) = particles_elem(:,4)
            part_ic(:,6) = particles_MMR

            ! Primary
            !! mass, radius

            ! m0 and Boulders (from m0)
            !! mass, radius, theta
            if (use_boulders) then
                allocate(m0boul_ic(0:Nboulders, 3))
                m0boul_ic(:Nboulders,1) = mass_ast_arr(0:Nboulders)
                m0boul_ic(:Nboulders,2) = radius_ast_arr(0:Nboulders)
                m0boul_ic(0,3) = cero
                m0boul_ic(1:Nboulders,3) = theta_from_primary(1:Nboulders)
            end if

            ! Asteroid
            !! mass, radius, pos, vel, theta, omega, inertia, angmom, Prot, acorot
            allocate(ast_ic(12))
            ast_ic(1) = asteroid_mass
            ast_ic(2) = asteroid_radius
            ast_ic(3:4) = asteroid_pos
            ast_ic(5:6) = asteroid_vel
            ast_ic(7) = asteroid_theta
            ast_ic(8) = asteroid_omega
            ast_ic(9) = asteroid_inertia
            ast_ic(10) = asteroid_angmom
            ast_ic(11) = asteroid_rotational_period
            ast_ic(12) = asteroid_a_corot            
        end subroutine save_initial

        ! 21.x Subroutines to write to a file (or screen)

        ! 21.1 Write elements to unit_file
        subroutine write_elements(i, unit_file)
            implicit none
            integer(kind=4), intent(in) :: i, unit_file
            
            aux_integer = sorted_particles_index(i)
            write (unit_file,f233) &
                & particles_index(aux_integer), &
                & time / unit_time, &
                & particles_elem(aux_integer,1) / unit_dist, &
                & particles_elem(aux_integer,2), &
                & particles_elem(aux_integer,3) / radian, &
                & particles_elem(aux_integer,4) / radian, &
                & particles_MMR(aux_integer), &
                & particles_mass(aux_integer) / unit_mass, &
                & particles_dist(aux_integer) / unit_dist, &
                & asteroid_omega * unit_time, &
                & asteroid_mass / unit_mass
        end subroutine write_elements

        ! 21.2 Write coordinates to unit_file (particles)
        subroutine write_coordinates_particle(i, unit_file)
            implicit none
            integer(kind=4), intent(in) :: i, unit_file

            aux_integer = sorted_particles_index(i)
            write (unit_file,f233) &
                & particles_index(aux_integer) + Nboulders, &
                & time / unit_time, &
                & particles_coord(aux_integer,1:2) / unit_dist, &
                & particles_coord(aux_integer,3:4) / unit_vel, &
                & particles_acc(aux_integer,1:2) / unit_acc, &
                & particles_mass(aux_integer) / unit_mass, &
                & particles_dist(aux_integer) / unit_dist
        end subroutine write_coordinates_particle

        ! 21.3 Write coordinates to unit_file (boulders)
        subroutine write_coordinates_boulders(i, unit_file)
            implicit none
            integer(kind=4), intent(in) :: i, unit_file
            
            write (unit_file,f233) &
                & i, &
                & time / unit_time, &
                & pos_ast_arr(i,:) / unit_dist, &
                & vel_ast_arr(i,:) / unit_vel, &
                & acc_ast_arr(i,:) / unit_acc, &
                & mass_ast_arr(i) / unit_mass, &
                & radius_ast_arr(i) / unit_dist
        end subroutine write_coordinates_boulders

        ! 21.4 DO NOT Write
        subroutine do_not_write(i, unit_file)
            implicit none
            integer(kind=4), intent(in) :: i, unit_file
            
        end subroutine do_not_write

        ! 21.5 Write chaos file
        subroutine write_chaos(unit_file)
            implicit none
            integer(kind=4), intent(in) :: unit_file
            integer(kind=4) :: i, aux_integer
            integer(kind=4), dimension(:), allocatable :: my_sorted_particles_index

            ! Sort particles by index
            allocate (my_sorted_particles_index(Nparticles))
            do i = 1, Nparticles
                my_sorted_particles_index = i
            end do
            call quickargsort_int(particles_index, my_sorted_particles_index, 1, Nparticles)
            call fseek(unit_file, 0, 0)       ! move to beginning
            ! Remember that Initial conditions were not swapped
            do i = 1, Nparticles
                aux_integer = my_sorted_particles_index(i)
                write (unit_file,f2233) &
                & particles_index(aux_integer), & ! i
                & particles_outcome(aux_integer), & ! bad
                & final_time / unit_time, & ! total time to integrate
                & asteroid_initial_conditions(10) / (unit_mass * unit_dist * unit_vel), & ! initial (Asteroid): angular momentum
                & particles_initial_conditions(i,1) / unit_mass, & ! initial: mass
                & particles_initial_conditions(i,2) / unit_dist, & ! initial: a
                & particles_initial_conditions(i,3), & ! initial: e
                & particles_initial_conditions(i,4) / radian, & ! initial: M
                & particles_initial_conditions(i,5) / radian, & ! initial: omega
                & particles_initial_conditions(i,6), & ! initial: MMR
                & sqrt(particles_initial_conditions(i,2) * &
                  & (uno - particles_initial_conditions(i,3)**2) * &
                  & G * (particles_initial_conditions(i,1) + asteroid_initial_conditions(1))) / &
                  & (unit_dist * unit_vel), & ! initial: angular momentum per unit mass
                & particles_times(aux_integer) / unit_time, & ! surviving time
                & particles_elem(aux_integer,1) / unit_dist, particles_elem(aux_integer,2), & ! final: a, e
                & particles_elem(aux_integer,3) / radian, particles_elem(aux_integer,4) / radian, & ! final: M, omega
                & particles_MMR(aux_integer), & ! final: MMR
                & sqrt(particles_elem(aux_integer,1) * &
                  & (uno - particles_elem(aux_integer,2)**2) * &
                  & G * (particles_mass(aux_integer) + asteroid_mass)) / &
                  & (unit_dist * unit_vel), & ! final: angular momentum per unit mass
                & particles_min_a(aux_integer) / unit_dist, & ! a_min 
                & particles_max_a(aux_integer) / unit_dist, & ! a_max
                & particles_min_e(aux_integer), & ! e_min 
                & particles_max_e(aux_integer), & ! e_max
                & (particles_max_a(aux_integer) - particles_min_a(aux_integer)) / unit_dist, & ! Delta a
                & (particles_max_e(aux_integer) - particles_min_e(aux_integer)) ! Delta e
            end do
            deallocate (my_sorted_particles_index)
            flush(unit_file)
        end subroutine write_chaos

        ! 22.x Subroutines to get coordinates and elements from a body

        ! 23. Obtener elementos de partícula i
        subroutine calculate_elements_i(i)
            implicit none
            integer(kind=4), intent(in) :: i
            
            call elem(asteroid_mass, &  !+particles_mass(i) podría ser agregado
                        & (/particles_coord(i,1:2) - asteroid_pos, cero, &
                        & particles_coord(i,3:4) - asteroid_vel, cero/), &
                        & particles_elem(i,1), particles_elem(i,2), dummy_real, &
                        & particles_elem(i,3), particles_elem(i,4), dummy_real)
            
            particles_MMR(i) = (particles_elem(i,1) / asteroid_a_corot)**(1.5) ! MMR
        end subroutine calculate_elements_i

        ! 24. Calcular caos de la partícula i
        subroutine calculate_chaos_i(i)
            implicit none
            integer(kind=4), intent(in) :: i

            particles_min_a(i) = min(particles_min_a(i), particles_elem(i,1))
            particles_max_a(i) = max(particles_max_a(i), particles_elem(i,1))
            particles_min_e(i) = min(particles_min_e(i), particles_elem(i,2))
            particles_max_e(i) = max(particles_max_e(i), particles_elem(i,2))
        end subroutine calculate_chaos_i

        ! 25. Nullify pointers
        subroutine nullify_pointers()
            implicit none

            nullify(particles_hexitptr)
            nullify(get_elements_i)
            nullify(get_chaos_i)
            nullify(resolve_merge)
            nullify(write_i_to_general)
            nullify(write_i_to_screen)
            nullify(write_i_to_individual)
            nullify(write_b_to_general)
            nullify(write_b_to_screen)
            nullify(write_b_to_individual)
            nullify(flush_chaos)
        end subroutine nullify_pointers
        
        ! 26.x Subroutines to switch between lower and upper case

        ! 27. Flush to a file
        subroutine flush_to_file(unit_file)
            implicit none
            integer(kind=4), intent(in) :: unit_file

            flush(unit_file)
        end subroutine flush_to_file
        
        ! 28.x Functions to check collision/escape
        
        ! 28.1 Check distances (Version 1)
        function check_continue_v1 (y) result(keep_going)
            implicit none
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8) :: rb(2)
            real(kind=8) :: rcm(2), dist_to_cm, r_from_cm(2)
            real(kind=8) :: rib(0:Nboulders,2), dist_to_i, r_from_i(2)
            integer(kind=4) :: j, particle_j, i
            integer(kind=4), parameter :: neqs = 4
            logical :: keep_going

            ! Re-set particles_hexit
            particles_hexit = 0  ! 0 is the default or non-action value
            
            ! Calculate the center of mass of the asteroid and boulders
            rcm = cero
            do i = 0, Nboulders 
                rib(i,1) = y((i * neqs) + 1)
                rib(i,2) = y((i * neqs) + 2)
                rcm = rcm + mass_ast_arr(i) * rib(i,:)
            end do
            rcm = rcm / asteroid_mass ! rcm = sum_i m_i * r_i / M

            ! Calculate vector and distance to CM
            do j = 1, Nactive
                particle_j = (j + Nboulders) * neqs
                rb(1) = y(particle_j+1)
                rb(2) = y(particle_j+2)
                r_from_cm = rb - rcm
                dist_to_cm = sqrt(r_from_cm(1)*r_from_cm(1) + r_from_cm(2)*r_from_cm(2))
                if (dist_to_cm > max_distance) then
                    particles_hexit(j) = 2
                else if (dist_to_cm < min_distance) then
                    particles_hexit(j) = 1
                else
                    ! Calculate vector and distance to boulders
                    boulders_loop: do i = 0, Nboulders
                        r_from_i = rb - rib(i,:)
                        dist_to_i = sqrt(r_from_i(1)*r_from_i(1) + r_from_i(2)*r_from_i(2))
                        if (dist_to_i < radius_ast_arr(i)) then
                            particles_hexit(j) = 1
                            exit boulders_loop
                        end if
                    end do boulders_loop
                end if
            end do
            
            keep_going = all(particles_hexit(1:Nactive) .eq. 0)
            if (.not. keep_going) particles_hexit(0) = 1
            
        end function check_continue_v1
        
        ! 28.2 Check distances (Version 2) [from cm]
        function check_continue_v2 (y) result(keep_going)
            implicit none
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8) :: rb(2)
            real(kind=8) :: rcm(2), dist_to_cm, r_from_cm(2)
            real(kind=8) :: rib(0:Nboulders,2), dist_to_i, r_from_i(2)
            integer(kind=4) :: j, particle_j, i
            integer(kind=4), parameter :: neqs = 4
            logical :: keep_going

            ! Re-set particles_hexit
            particles_hexit = 0  ! 0 is the default or non-action value
            
            ! Calculate the center of mass of the asteroid and boulders
            rcm(1) = y(3)
            rcm(2) = y(4)
            do i = 0, Nboulders
                rib(i,1) = cos(y(1) + theta_ast_arr(i)) * dist_ast_arr(i)
                rib(i,2) = sin(y(1) + theta_ast_arr(i)) * dist_ast_arr(i)
            end do

            ! Calculate distance and vector to boulders
            do j = 1, Nactive
                particle_j = j * neqs + 2
                rb(1) = y(particle_j+1)
                rb(2) = y(particle_j+2)
                r_from_cm = rb - rcm
                dist_to_cm = sqrt(r_from_cm(1)*r_from_cm(1) + r_from_cm(2)*r_from_cm(2))
                if (dist_to_cm > max_distance) then
                    particles_hexit(j) = 2
                else if (dist_to_cm < min_distance) then
                    particles_hexit(j) = 1
                else 
                    ! Calculate vector and distance to boulders
                    boulders_loop: do i = 0, Nboulders
                        r_from_i = rb - rib(i,:)
                        dist_to_i = sqrt(r_from_i(1)*r_from_i(1) + r_from_i(2)*r_from_i(2))
                        if (dist_to_i < radius_ast_arr(i)) then
                            particles_hexit(j) = 1
                            exit boulders_loop
                        end if
                    end do boulders_loop
                end if
            end do
            
            keep_going = all(particles_hexit(1:Nactive) .eq. 0)
            if (.not. keep_going) particles_hexit(0) = 1
            
        end function check_continue_v2
        
        ! 99 Subrutina para pointer vacío (no hace nada) con input i
        subroutine do_nothing_i(i)
            implicit none
            integer(kind=4), intent(in) :: i
        end subroutine do_nothing_i


        ! 101 Escribir porcentaje
        subroutine percentage(tout, tstop)
            implicit none
            real(kind=8), intent(in) :: tout, tstop
            integer(kind=4) :: iper
            character(len=100) :: guiones
            character(len=4) :: cestado
        
            iper = int(100.0 * tout / tstop)
            guiones = repeat('.', iper)
        
            if (iper < 100) then
                write(cestado, '(i3)') iper
                write(*, '(a)', advance='no') char(13) // trim(guiones) // trim(adjustl(cestado)) // '%'
            else
                write(*, '(a)', advance='no') char(13) // trim(guiones) // '. FIN'
                write(*, *)
            end if
        end subroutine percentage
        

end module parameters
