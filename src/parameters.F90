module parameters
    use constants
    use bins
    use omp_lib

    implicit none

    ! Simulation
    integer(kind=4) :: simulation_number
    integer(kind=4) :: Nactive

    ! Variables globales. Varias dentro de la subrutina read_conf_file
    logical :: existe_configfile
    
    ! Configuración de integración
    !! Tiempos
    real(kind=8) :: initial_time, final_time, output_timestep
    real(kind=8) :: min_timestep
    integer(kind=4) :: output_number
    integer(kind=4) :: case_output_type ! 0 => linear, 1 => log, 2 => combination
    real(kind=8), dimension(:), allocatable :: output_times
    real(kind=8) :: time, timestep
    real(kind=8) :: adaptive_timestep
    real(kind=8), dimension(:), allocatable :: checkpoint_times
    integer(kind=4) :: checkpoint_number
    !! Parámetros
    integer(kind=4) :: error_digits
    real(kind=8) :: error_tolerance
    real(kind=8) :: learning_rate
    real(kind=8) :: min_distance, max_distance
    !! Integración
    logical :: use_explicit_method
    logical :: use_version_1
    logical :: use_parallel
    !!! Parallel
    integer(kind=4) :: available_threads, requested_threads, my_threads
    logical :: compiled_with_openmp = .False.

    ! I/O
    logical :: use_screen, use_datascreen, use_percentage
    logical :: use_elements_output
    logical :: use_datafile, use_chaosfile, use_potential_map
    logical :: use_tomfile, use_particlesfile
    logical :: use_multiple_outputs
    character(30) :: datafile, chaosfile, mapfile, multfile
    character(30) :: tomfile, particlesfile
    logical :: use_flush_output
    logical :: only_potential_map = .False.
    !! Tomfile
    integer(kind=4) :: tom_index_number, tom_total_number
    real(kind=8), dimension(:), allocatable :: tom_times, tom_deltaomega, tom_deltamass
    logical, dimension(:), allocatable :: checkpoint_is_tom, checkpoint_is_output
    real(kind=8) :: tom_mass_growth_param
    !! Mapa
    integer(kind=4) :: map_grid_size_x, map_grid_size_y
    real(kind=8) :: map_min_x, map_max_x, map_min_y, map_max_y

    ! Forces
    !!Internal forces
    logical :: use_torque
    ! External forces
    !! Stokes
    logical :: use_stokes
    real(kind=8) :: stokes_a_damping_time, stokes_e_damping_time, stokes_charac_time
    real(kind=8) :: stokes_C, stokes_alpha
    !! Naive-Stokes (Drag)
    logical :: use_naive_stokes
    real(kind=8) :: drag_coefficient, drag_charac_time
    !! Mass and omega damping
    logical :: use_omega_damping
    real(kind=8) :: mass_exp_damping_time, omega_exp_damping_time
    real(kind=8) :: omega_exp_poly_A, omega_exp_poly_B, omega_exp_poly_AB
    real(kind=8) :: omega_linear_damping_time, omega_linear_damping_slope
    !! Geo-potential
    logical :: use_J2
    real(kind=8) :: J2_coefficient, J2_effective
    !! Self-Gravity [Use BINS]
    logical :: use_self_gravity
    integer(kind=4) :: Norder_self_gravity
    !! Viscosity [Use BINS]
    logical :: use_viscosity
    real(kind=8) :: viscosity
    
    !!! BINS
    logical :: use_bins
    integer(kind=4) :: Nbins, binning_method
    real(kind=8) :: rmin_bins, rmax_bins
    type(my_bins) :: disk_bins
    logical :: update_rmin_bins, update_rmax_bins, update_bins
    integer(kind=4), dimension(:), allocatable :: particles_bins
    

    ! Bodies
    integer(kind=4) :: Ntotal
    !! Asteroid
    real(kind=8) :: asteroid_mass, asteroid_radius, asteroid_inertia, asteroid_angmom
    real(kind=8) :: asteroid_omega, asteroid_rotational_period, lambda_kep, omega_kep
    real(kind=8) :: asteroid_a_corot 
    real(kind=8), dimension(2) :: asteroid_pos, asteroid_vel, asteroid_acc
    real(kind=8) :: asteroid_theta, asteroid_theta_correction, asteroid_torque
    !!! Primary + boulders
    real(kind=8), dimension(:), allocatable :: mass_ast_arr, radius_ast_arr, mu_ast_arr, theta_ast_arr
    real(kind=8), dimension(:), allocatable :: dist_ast_arr, inertia_ast_arr
    real(kind=8), dimension(:,:), allocatable :: pos_ast_arr, vel_ast_arr, acc_ast_arr
    !!! Primary
    real(kind=8) :: mass_primary, radius_primary
    !!! Boulders
    logical :: use_boulders
    integer(kind=4) :: Nboulders
    real(kind=8), dimension(:), allocatable :: theta_from_primary, mu_from_primary
    real(kind=8), dimension(:,:), allocatable :: pos_from_primary
    real(kind=8), dimension(:,:), allocatable :: vel_from_primary
    real(kind=8), dimension(:,:), allocatable :: acc_from_primary
    !! Single particle
    logical :: use_single_particle
    !!! Masa
    real(kind=8) :: single_part_mass
    !!! Elementos orbitales
    real(kind=8) :: single_part_elem_a, single_part_elem_e
    real(kind=8) :: single_part_elem_M, single_part_elem_w
    real(kind=8) :: single_part_MMR
    !!!! Varias partículas
    integer(kind=4) :: Nparticles, Nparticles_type1, Nparticles_type2, Narticles_in_partfile
    !!!!! Indices
    integer(kind=4), dimension(:), allocatable :: particles_index
    integer(kind=4), dimension(:), allocatable :: sorted_particles_index
    !!!!! Masas
    real(kind=8), dimension(:), allocatable :: particles_mass
    !!!!! Coordenadas y elementos orbitales
    real(kind=8), dimension(:,:), allocatable :: particles_coord, particles_elem
    real(kind=8), dimension(:,:), allocatable :: particles_acc
    real(kind=8), dimension(:), allocatable :: particles_MMR, particles_dist
    ! Bodies (AUX)
    !! Asteroid
    real(kind=8) :: Gasteroid_mass
    real(kind=8) :: asteroid_omega2
    real(kind=8) :: pos_ast_from_primary(2), vel_ast_from_primary(2), acc_ast_from_primary(2)
    real(kind=8) :: theta_ast_from_primary
    !!! Primary + boulders
    real(kind=8), dimension(:), allocatable :: Gmass_ast_arr
    !!! Particles
    real(kind=8), dimension(:,:), allocatable :: aux_particles_arr
    integer(kind=4) :: first_particle, last_particle
    integer(kind=4), dimension(:), allocatable :: particles_outcome
    integer(kind=4), dimension(:,:), allocatable :: ij_to_swap

    ! Chaos
    logical :: use_chaos, use_update_chaos
    real(kind=8), dimension(:), allocatable :: particles_max_a, particles_max_e
    real(kind=8), dimension(:), allocatable :: particles_min_a, particles_min_e
    real(kind=8), dimension(:), allocatable :: particles_times
    integer(kind=4) :: discarded_particles, staying_particles
    

    ! Extra/Auxiliar variables
    integer(kind=4) :: i, j, k
    !! Arguments
    integer(kind=4) :: arguments_number
    !! Auxiliar
    logical :: aux_logical, is_number
    real(kind=8) :: aux_real, dummy_real, dummy_real2
    integer(kind=4) :: aux_integer
    character(1)    :: aux_character1
    character(2)    :: aux_character2
    character(20)   :: aux_character20
    character(30)   :: aux_character30
    real(kind=8), dimension(2) :: distance_to_primary
    real(kind=8), dimension(2) :: aux_real_arr2
    real(kind=8), dimension(6) :: aux_real_arr6
    logical :: hard_center
    logical :: use_elements

    !!! User variables (JUST FOR DEBUGGING)
    real(kind=8), dimension(2) :: user_real_var2

    ! PARAMETERS ARRAY
    integer, parameter :: equation_size = 4
    real(kind=8), dimension(:), allocatable :: parameters_arr, parameters_arr_new, parameters_der
    
    !!! Hard Exit
    logical :: is_premature_exit
    integer(kind=4), dimension(:), allocatable, target :: particles_hexit !  Hard Exit integer
    integer(kind=4), pointer :: particles_hexitptr ! pointer to Hard Exit
    procedure (check_continue_template), pointer :: check_continue_ptr => null ()
    
    !!! INITIAL CONDITIONS (Here we store the initial conditions of everything)
    real(kind=8), dimension(:,:), allocatable :: particles_initial_conditions !! mass, a, e, M, w, MRR
    real(kind=8), dimension(:,:), allocatable :: m0_and_boulders_initial_conditions !! mass, radius, theta
    real(kind=8), dimension(:), allocatable :: asteroid_initial_conditions !! mass, radius, pos, vel, theta, omega, inertia, angmom, Prot, acorot
    

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

        ! 1. Inicializar los parámetros por defecto
        subroutine init_default_parameters()
            implicit none

            ! Simulation
            simulation_number = 1
            ! Times
            initial_time = cero
            final_time = cero
            output_timestep = cero
            min_timestep = cero
            output_number = 0
            case_output_type = 0 ! 0 => linear, 1 => log, 2 => combination
            ! Parameters
            use_explicit_method = .False.
            use_version_1 = .False.
            use_parallel = .False.
            !! Parallel
            requested_threads = 1
            my_threads = 1
            !$ compiled_with_openmp = .True. !! Parece comentado, pero es asi
            ! Adaptive configuration
            error_digits = 12
            learning_rate = uno
            ! Bodies
            !! Primary
            mass_primary = cero 
            radius_primary = cero 
            !! Rotation
            lambda_kep = cero ! Ratio of omega to keplerian 
            asteroid_rotational_period = cero ! Rotational period
            !! Boulders
            use_boulders = .False.
            Nboulders = 0
            !! Particles
            Nparticles = 0
            Nparticles_type1 = 0
            Nparticles_type2 = 0
            !!! Single
            use_single_particle = .False.
            single_part_mass = cero ! Masa
            single_part_elem_a = cero ! Semieje mayor
            single_part_elem_e = cero ! Excentricidad
            single_part_elem_M = cero ! Anomalía media
            single_part_elem_w = cero ! Argumento del perihelio
            single_part_MMR = cero ! Mean Motion Resonance!! Particles
            !!! Particulas en archivo
            Narticles_in_partfile = 0
            particlesfile = "" ! Archivo de partículas
            !! ASTEROID
            asteroid_mass = cero
            asteroid_radius = cero
            asteroid_inertia = cero
            asteroid_pos = (/cero, cero/)
            asteroid_vel = (/cero, cero/)
            asteroid_acc = (/cero, cero/)
            asteroid_theta = cero
            asteroid_theta_correction = cero
            ! Forces
            !! InternalForces
            use_torque = .False.
            ! ExternalForces
            !! Stokes
            use_stokes = .False.
            stokes_charac_time = cero
            stokes_a_damping_time = infinity
            stokes_e_damping_time = infinity
            !! Naive-Stokes
            use_naive_stokes = .False.
            drag_coefficient = cero
            drag_charac_time = cero
            !! Dampings
            use_omega_damping = .False.
            omega_linear_damping_time = infinity
            omega_exp_damping_time = infinity
            omega_exp_poly_A = cero
            omega_exp_poly_B = cero
            mass_exp_damping_time = infinity
            !! Geo-Potential
            use_J2 = .False.
            J2_coefficient = cero
            !! Self Gravity
            use_self_gravity = .False.
            Norder_self_gravity = 3
            !! Viscosity
            use_viscosity = .False.
            viscosity = - uno
            !! BINS (only used if Self Gravity or Viscosity)
            Nbins = 3
            binning_method = 1
            rmin_bins = - uno ! -1 means first particle (initial)
            rmax_bins = - uno ! -1 means last particle (initial)
            ! Escape/Collision (and merge)
            min_distance = cero
            max_distance = infinity
            use_merge = .False.
            ! Tomfile
            tomfile = "" ! Archivo de t_TOM, domega(t_TOM), dm(t_TOM)
            ! I/O
            use_screen  = .False. ; use_percentage = .False. ! Pantalla, porcentaje
            use_datascreen = .False. ! Salida de datos en pantalla
            use_elements_output  = .True. ! Salida en elementos orbitales
            use_multiple_outputs = .False. ! Salida en múltiples archivos
            !! Output (Poner "" o "no" para no crear archivo)
            datafile = ""! Archivo de datos
            chaosfile = "" ! Archivo de datos
            multfile = "" ! Prefijo de archivos de salida individuales
            !!! Mapa de potencial
            mapfile = "" ! Archivo de mapas
            map_grid_size_x = 500 ; map_grid_size_y = 500
            map_min_x = -300. ; map_max_x = 300.
            map_min_y = -300. ; map_max_y = 300.


            ! Extras, fuera de Config.ini. Se editan solo acá
            !! Centrar CM?
            hard_center = .True.
            !! Update chaos at every output timestep, instead of waiting to the end of the run
            use_update_chaos = .True.
            !! Flush output at every output timestep (even though the buffer is not full) 
            use_flush_output = .True.
        end subroutine init_default_parameters

        ! 2. Leer el archivo de configuración
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
                if (use_screen) write (*,*) "El archivo ", trim(file_name), " no existe."
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
                        case("use implicit me")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                use_explicit_method = .False.
                            else
                                use_explicit_method = .True.
                            end if
                        case("use version 2 o")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                use_version_1 = .False.
                            else
                                use_version_1 = .True.
                            end if
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
                        case("number of bould")
                            read (value_str, *) Nboulders
                            if (Nboulders < 0) then
                                write (*,*) "El número de boulders debe ser mayor que 0."
                                stop 1
                            end if
                            !! Allocamos
                            if (.not. allocated(mass_ast_arr)) then
                                call allocate_asteriod_arrays(Nboulders)
                            else 
                                write (*,*) "ERROR: Redefinido."
                            end if
                        case("use single part")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                use_single_particle = .True.
                            else
                                use_single_particle = .False.
                            end if
                        case("mass of particl")
                            read (value_str, *) single_part_mass
                        case("semi-major axis")
                            read (value_str, *) single_part_elem_a
                        case("eccentricity of")
                            read (value_str, *) single_part_elem_e
                        case("mean anomaly of")
                            read (value_str, *) single_part_elem_M
                        case("pericenter long")
                            read (value_str, *) single_part_elem_w
                        case("mean motion rat")
                            read (value_str, *) single_part_MMR
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
                        case("include torque")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                use_torque = .True.
                            else
                                use_torque = .False.
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
                        ! case("include tri-axi")
                        !     if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                        !         loTriAx = .True.
                        !     else
                        !         loTriAx = .False.
                        !     end if
                        ! case("first ellipsoid")
                        !     read (value_str, *) tria
                        ! case("second ellipsoi")
                        !     read (value_str, *) trib
                        ! case("third ellipsoid")
                        !     read (value_str, *) tric
                        case("min distance fr")
                            read (value_str, *) min_distance
                        case("max distance fr")
                            read (value_str, *) max_distance
                        case("merge collsions")
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
                        case("idividual file")
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
                else
                    param_str = trim(adjustl(line)) 
                    if (param_str(1:7) .eq. "mass_m0") then
                        if (Nboulders < 0) then
                            write (*,*) "ERROR: El número de boulders debe especificarse antes que los parametros."
                            stop 1
                        else if (Nboulders == 0) then
                            write (*,*) "WARNING: No se leen boulders."
                        end if
                        do j = 1, Nboulders
                            io = 0
                            check: do while (io == 0)
                                nlines = nlines + 1
                                read (10, '(A)') line
                                if ((len_trim(line) == 0) .or. (line(:2) == "c ") .or. (line(:2) == "! ")) cycle
                                backspace (10)
                                read (10, *, iostat=io) mu_from_primary(j), radius_ast_arr(j), theta_from_primary(j)
                                if (io /= 0) then
                                    write (*,*) "ERROR: Al leer boulder. Línea: ", nlines
                                    stop 1
                                end if
                                io = 1
                            end do check
                        end do
                    end if
                end if
            end do
            98 close (10)
        end subroutine read_config_file

        ! 3 Leer entradas por linea de comandos
        subroutine load_command_line_arguments()
            implicit none

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
                    case ("-mpart")
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) single_part_mass
                        aux_integer = 1
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
                    case ("--implicit")
                        use_explicit_method = .False.
                    case ("--explicit")
                        use_explicit_method = .True.
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
                    case ("-partfile")
                        use_particlesfile = .True.
                        call get_command_argument(i+1, particlesfile)
                        aux_integer = 1
                    case ("--nopartfile")
                        use_particlesfile = .False.
                        particlesfile = ""
                    case ("--noconfig")
                        existe_configfile = .False.
                    case ("--version1")
                        use_version_1 = .True.
                    case ("--version2")
                        use_version_1 = .False.
                    case ("--merge")
                        use_merge = .True.
                    case ("--nomerge")
                        use_merge = .False.
                    case ("--torque")
                        use_torque = .True.
                    case ("--notorque")
                        use_torque = .False.
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
                        write (*,*) "    --elem       : Imprimir elementos orbitales (solo partículas) [default]"
                        write (*,*) "    --noelem     : Imprimir coordenadas baricéntricas"
                        write (*,*) "    -tomfile     : Utilizar archivo de (t)iempos|omega|masa que sigue"
                        write (*,*) "    --notomfile  : No utilizar archivo de (t)iempos|omega|masa"
                        write (*,*) "    -partfile    : Utilizar archivo de partículas que sigue"
                        write (*,*) "    --nopartfile : No utilizar archivo de partículas"
                        write (*,*) "    --noconfig   : No leer archivo de configuración"
                        write (*,*) "    --version1   : Usar versión 1: y=(B0, B1, ..., P0, ...)"
                        write (*,*) "    --version2   : Usar versión 2: y=(θ, ω, P0, ...) [default]"
                        write (*,*) "    --merge      : Incluir asociar colisiones de partículas a asteroide"
                        write (*,*) "    --nomerge    : No asociar colisiones de partículas"
                        write (*,*) "    --torque     : Incluir torque de partículas hacia el asteroide"
                        write (*,*) "    --notorque   : No incluir torque de partículas hacia el asteroide"
                        write (*,*) "    -parallel    : Paralelizar usando la cantida de thread que sique"                        
                        write (*,*) "    --parallel   : Paralelizar usando todos los threads disponibles"
                        write (*,*) "    --noparallel : No usar paralelización para partículas"
                        write (*,*) "    --help       : Mostrar esta ayuda"
                        stop 0
                    case default
                        is_number = .False.
                        call get_command_argument(i, aux_character30)
                        do j = 0, 9
                            if (aux_character30(1:1) == char(48 + j)) then !check if it's a number
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
                            else ! Leo los argumentos numéricos
                                call get_command_argument(i, aux_character20)
                                read (aux_character20,*) single_part_elem_a
                                call get_command_argument(i+1, aux_character20)
                                read (aux_character20,*) single_part_elem_e
                                call get_command_argument(i+2, aux_character20)
                                read (aux_character20,*) single_part_elem_M
                                call get_command_argument(i+3, aux_character20)
                                read (aux_character20,*) single_part_elem_w
                                aux_integer = 3
                                aux_logical = .True.
                                use_single_particle = .True.
                                single_part_MMR = cero
                            end if
                        else ! Ya leí los numéricos. Falta leer eR
                            call get_command_argument(i, aux_character20)
                            read (aux_character20,*) single_part_MMR
                        end if
                    end select
                end do
            else
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
            if (abs(omega_linear_damping_time) < tini) omega_linear_damping_time = infinity
            if (abs(omega_exp_damping_time) < tini) omega_exp_damping_time = infinity
            if (abs(omega_exp_poly_A) < tini) omega_exp_poly_A = cero
            if (abs(omega_exp_poly_B) < tini) omega_exp_poly_B = cero
            if (.not. ((abs(omega_linear_damping_time) > tini) .and. (omega_linear_damping_time < infinity)) .and. &
              & .not. ((abs(omega_exp_damping_time) > tini) .and. (omega_exp_damping_time < infinity)) .and. &
              & .not. ((abs(omega_exp_poly_A) > tini) .and. (abs(omega_exp_poly_B) > tini))) &
              & use_omega_damping = .False.
            if (abs(mass_exp_damping_time) < tini) mass_exp_damping_time = infinity
            !!! Geo-Potential
            if (abs(J2_coefficient) > tini) use_J2 = .True.
            ! !!! Triaxial
            ! if ((tria < tini) .or. .not. loTriAx) tria = cero
            ! if ((trib < tini) .or. .not. loTriAx) trib = cero
            ! if ((tric < tini) .or. .not. loTriAx) tric = cero
            ! if ((abs(tria) < tini) .or. (abs(trib) < tini) .or. (abs(tric) < tini)) loTriAx = .False.

            ! if (loTriAx .and. (Nboulders > 0)) then
            !     write (*,*) "ERROR: No se puede usar triaxialidad con boulders."
            !     stop 1
            ! end if
            
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
            
            !! Torque
            if (use_torque .and. use_explicit_method) then
                write (*,*) "ERROR: No se puede usar torque con el método explícito."
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

            ! Single particle
            if (.not. use_single_particle) then
                single_part_mass = cero
                single_part_elem_a = cero
                single_part_elem_e = cero
                single_part_elem_M = cero
                single_part_elem_w = cero
                single_part_MMR = cero
            end if
        end subroutine set_derived_parameters
        
        ! 5.x Subroutines to allocate arrays

        ! 5.1 Alocatar arrays asteroide
        subroutine allocate_asteriod_arrays(Nboulders)
            implicit none
            integer(kind=4), intent(in) :: Nboulders

            !! Parámetros y auxiliares
            allocate(mass_ast_arr(0:Nboulders), radius_ast_arr(0:Nboulders), mu_ast_arr(0:Nboulders), Gmass_ast_arr(0:Nboulders))
            allocate(pos_ast_arr(0:Nboulders,2), vel_ast_arr(0:Nboulders,2), acc_ast_arr(0:Nboulders,2))
            allocate(theta_ast_arr(0:Nboulders), dist_ast_arr(0:Nboulders), inertia_ast_arr(0:Nboulders))
            if (Nboulders > 0) then
                allocate(theta_from_primary(Nboulders), mu_from_primary(Nboulders))
                allocate(pos_from_primary(Nboulders,2), vel_from_primary(Nboulders,2), acc_from_primary(Nboulders,2))
            end if
        end subroutine allocate_asteriod_arrays

        ! 5.2 Alocatar particulas
        subroutine allocate_particles(Nparticles)
            implicit none
            integer(kind=4), intent(in) :: Nparticles

            !! Parámetros y auxiliares
            allocate(particles_index(1:Nparticles))
            allocate(particles_mass(1:Nparticles))
            allocate(particles_coord(1:Nparticles,4), particles_elem(1:Nparticles,4))
            allocate(particles_acc(1:Nparticles,2))
            allocate(particles_MMR(1:Nparticles))
            allocate(particles_dist(1:Nparticles))
            allocate(sorted_particles_index(1:Nparticles))
            allocate(particles_outcome(1:Nparticles))
            allocate(ij_to_swap(1:Nparticles,2))
            if (use_chaos) then
                allocate(particles_times(1:Nparticles))
                allocate(particles_min_a(1:Nparticles), particles_max_a(1:Nparticles))
                allocate(particles_min_e(1:Nparticles), particles_max_e(1:Nparticles))
            end if
            allocate(particles_hexit(0:Nparticles))
            particles_hexitptr => particles_hexit(0) ! Puntero a Hard Exit (NO TOCAR) !!! Ya está en default
            allocate(particles_bins(1:Nparticles))  ! Unused if use_bins is False
        end subroutine allocate_particles

        ! 5.3 Liberar arrays asteroide
        subroutine free_asteroid_arrays()
            implicit none

            deallocate(mass_ast_arr, radius_ast_arr, mu_ast_arr, Gmass_ast_arr)
            if (use_boulders) then
                deallocate(theta_from_primary, mu_from_primary)
                deallocate(pos_from_primary, vel_from_primary, acc_from_primary)
            end if
            deallocate(pos_ast_arr, vel_ast_arr, acc_ast_arr)
            deallocate(theta_ast_arr, dist_ast_arr, inertia_ast_arr)
        end subroutine free_asteroid_arrays

        ! 5.4 Liberar arrays particles
        subroutine free_particles()
            implicit none

            deallocate(particles_index, particles_mass)
            deallocate(particles_coord, particles_elem)
            deallocate(particles_acc)
            deallocate(particles_MMR, particles_dist)
            deallocate(sorted_particles_index)
            deallocate(particles_outcome)
            deallocate(ij_to_swap)
            if (use_chaos) then
                deallocate(particles_times)
                deallocate(particles_min_a, particles_max_a)
                deallocate(particles_min_e, particles_max_e)
            end if
            deallocate(particles_hexit)
            deallocate(particles_bins)  ! Unused if use_bins is False
        end subroutine free_particles

        ! 6. Setear los tiempos de salida
        subroutine set_output_times(t0, tf, n_out, dt_out, case_dist, t_out)
            implicit none
            real(kind=8), intent(in) :: t0, tf
            integer(kind=4), intent(inout) :: n_out
            real(kind=8), intent(inout) :: dt_out
            integer(kind=4), intent(in) :: case_dist
            real(kind=8), dimension(:), allocatable, intent(out) :: t_out
            real(kind=8) :: t_aux, t_add
            integer(kind=4) :: i
            real(kind=8) :: npointsr

            if (dt_out > 0.d0) then
                if (dt_out > tf - t0) then
                    write (*,*) "ERROR: dt_out >= (tf - t0)"
                    stop 1
                end if
                n_out = int((tf - t0) / dt_out) + 1
                npointsr = n_out * 1.d0
            else
                npointsr = n_out * 1.d0
                dt_out = (tf - t0) / (npointsr - 1.d0)
            end if
            if (n_out < 2) then
                write (*,*) "ERROR: n_out < 2"
                stop 1
            end if
            allocate (t_out(n_out))
            select case (case_dist)
            case(0) ! Lineal
                t_aux = (tf - t0) / (npointsr - uno)
                do i = 2, n_out - 1
                    t_out(i) = t0 + t_aux * (i-1)
                end do
            case(1) ! Logarítmico
                t_aux = exp (log (tf - t0 + 1.) / npointsr)
                t_add = t_aux
                do i = 2, n_out - 1
                    t_add    = t_add * t_aux
                    t_out(i) = t0 + t_add - 1.
                end do
            case(2) ! Combinado
                !! Teniendo en cuenta que 0 es t0 y n_out es tf, si se pide
                !!  una cantidad impar n_out de puntos, se van a hacer 
                !!  lin = (n_out-1)/2 puntos lineales y log = lin-1 puntos logarítmicos
                !!  Si se pide una cantidad par, lin = log = n_out/2
                ! Primero lineal
                npointsr = uno * (ceiling(npointsr * uno2 - tini) + 1) ! - tini para evitar errores de redondeo
                t_aux = (tf - t0) / (npointsr - uno)
                do i = 2, ceiling(npointsr - tini) - 1 ! - tini para evitar errores de redondeo
                    t_out(i) = t0 + t_aux * (i-1)
                end do
                ! Luego logarítmico
                if (mod(n_out, 2) .ne. 0) npointsr = npointsr - uno
                t_aux = exp (log (tf - t0 + 1.) / npointsr)
                t_add = t_aux
                do i = n_out - floor(npointsr + tini) + 2, n_out - 1 ! + tini para evitar errores de redondeo
                    t_add    = t_add * t_aux
                    t_out(i) = t0 + t_add - 1.
                end do
                ! Y al final ordenamos
                call quicksort(t_out, 2, n_out-1)
            end select
            t_out(1) = t0
            t_out(n_out) = tf
        end subroutine set_output_times

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

        ! 9. Combinar, ordenar y eliminar duplicados de dos arreglos
        subroutine merge_sort_and_unique(a, b, ina, inb, c, kfin)
            implicit none
            real(kind=8), dimension(:), intent(in) :: a, b ! Arreglos de reales
            logical, dimension(:), allocatable, intent(out) :: ina, inb ! Booleanos de elementos en a y b
            real(kind=8), dimension(:), allocatable, intent(out) :: c ! Arreglo combinado
            integer(kind=4), intent(out) :: kfin ! Número de elementos en c
            logical, dimension(:), allocatable :: ina0, inb0
            real(kind=8), dimension(:), allocatable :: c0
            integer(kind=4) :: i, j, k

            allocate (c0(size(a)+size(b)))
            allocate (ina0(size(c0)))
            allocate (inb0(size(c0)))
            c0 = 0.
            ina0 = .False.
            inb0 = .False.
            i = 1
            j = 1
            do k = 1, size(c0)
                if (i > size(a)) then
                    if (j > size(b)) then
                        kfin = k - 1
                        exit
                    else
                        c0(k) = b(j)
                        j = j + 1
                        inb0(k) = .True.
                    end if
                else if (j > size(b)) then
                    c0(k) = a(i)
                    i = i + 1
                    ina0(k) = .True.
                else if (a(i) < b(j)) then
                    c0(k) = a(i)
                    i = i + 1
                    ina0(k) = .True.
                else if (a(i) > b(j)) then
                    c0(k) = b(j)
                    j = j + 1
                    inb0(k) = .True.
                else
                    c0(k) = a(i)
                    i = i + 1
                    j = j + 1
                    ina0(k) = .True.
                    inb0(k) = .True.
                end if
                kfin = k
            end do

            allocate (c(kfin))
            allocate (ina(kfin))
            allocate (inb(kfin))
            c = c0(1:kfin)
            ina = ina0(1:kfin)
            inb = inb0(1:kfin)
            deallocate (c0)
            deallocate (ina0)
            deallocate (inb0)
        end subroutine merge_sort_and_unique

        ! 10. Calcular productro vectorial 2D
        function cross2D(a, b) result(res)
            implicit none
            real(kind=8), dimension(2), intent(in) :: a, b
            real(kind=8) :: res

            res = a(1) * b(2) - a(2) * b(1) ! Solo la componente z
        end function cross2D    

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

        ! 14. Liberar arrays de tiempos
        subroutine free_times_arrays()
            implicit none

            deallocate(output_times)
            deallocate(checkpoint_times)
            deallocate(checkpoint_is_output)
            deallocate(checkpoint_is_tom)
        end subroutine free_times_arrays

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
        
        ! 17.x Subroutines to sort arrays

        ! 17.1 Devuelve índices de ordenado
        recursive subroutine quickargsort(a, b, first, last)
            implicit none
            real(kind=8), intent(in) :: a(:)          ! Input array (not modified)
            integer(kind=4), intent(inout) :: b(size(a))    ! Indices array to be sorted
            integer(kind=4), intent(in) :: first, last
            real(kind=8) :: pivot
            integer(kind=4) :: i, j, temp

            ! Pivot selection
            pivot = a(b(floor((first + last) * uno2)))
            i = first
            j = last
            do
                do while (a(b(i)) < pivot)
                    i = i + 1
                end do
                do while (a(b(j)) > pivot)
                    j = j - 1
                end do
                if (i >= j) exit

                ! Swap indices in 'b' array
                temp = b(i)
                b(i) = b(j)
                b(j) = temp

                i = i + 1
                j = j - 1
            end do

            ! Recursive calls
            if (first < i - 1) call quickargsort(a, b, first, i - 1)
            if (j + 1 < last)  call quickargsort(a, b, j + 1, last)
        end subroutine quickargsort


        ! 17.2 Devuelve índices de ordenado (entero)
        recursive subroutine quickargsort_int(a, b, first, last)
            implicit none
            integer(kind=4), intent(in) :: a(:)    ! Input array (not modified)
            integer(kind=4), intent(inout) :: b(size(a)) ! Indices array to be sorted
            integer(kind=4), intent(in) :: first, last
            integer(kind=4) :: pivot
            integer(kind=4) :: i, j, temp

            ! Pivot selection
            pivot = a(b(floor((first + last) * uno2)))
            i = first
            j = last
            do
                do while (a(b(i)) < pivot)
                    i = i + 1
                end do
                do while (a(b(j)) > pivot)
                    j = j - 1
                end do
                if (i >= j) exit

                ! Swap indices in 'b' array
                temp = b(i)
                b(i) = b(j)
                b(j) = temp

                i = i + 1
                j = j - 1
            end do

            ! Recursive calls
            if (first < i - 1) call quickargsort_int(a, b, first, i - 1)
            if (j + 1 < last)  call quickargsort_int(a, b, j + 1, last)
        end subroutine quickargsort_int
        
        ! 17.3 Sort an array
        subroutine quicksort(a, first, last)
            implicit none
            real(kind=8), intent(inout) :: a(:)          ! Input array
            integer(kind=4), intent(in) :: first, last   ! First and last indices
            integer(kind=4), allocatable :: b(:)         ! Indices of array to be sorted
            real(kind=8), allocatable :: a_copy(:)       ! Temporary copy of the array
            integer(kind=4) :: n_sort
            
            ! Total amount of values to sort
            n_sort = last - first + 1
            
            ! Allocate and initialize the indices array
            allocate(b(1 : n_sort))
            allocate(a_copy(1 : n_sort))
            
            ! Initialize indices to [first, first+1, ..., last]
            do i = 1, n_sort
                b(i) = i
                a_copy(i) = a(first + i - 1)
            end do

            ! Sort indices using the quickargsort subroutine
            call quickargsort(a_copy, b, 1, n_sort)

            ! Copy the sorted values back to the original array
            a(first:last) = a_copy(b)

            ! Deallocate temporary arrays
            deallocate(b, a_copy)
        end subroutine quicksort
        
        ! 17.4 Sort and integer array
        subroutine quicksort_int(a, first, last)
            implicit none
            integer(kind=4), intent(inout) :: a(:)       ! Input array
            integer(kind=4), intent(in) :: first, last   ! First and last indices
            integer(kind=4), allocatable :: b(:)         ! Indices of array to be sorted
            integer(kind=4), allocatable :: a_copy(:)    ! Temporary copy of the array
            integer(kind=4) :: n_sort

            ! Total amount of values to sort
            n_sort = last - first + 1
            
            ! Allocate and initialize the indices array
            allocate(b(1 : n_sort))
            allocate(a_copy(1 : n_sort))
            
            ! Initialize indices to [first, first+1, ..., last]
            do i = 1, n_sort
                b(i) = i
                a_copy(i) = a(first + i - 1)
            end do

            ! Sort indices using the quickargsort subroutine
            call quickargsort_int(a_copy, b, 1, n_sort)

            ! Copy the sorted values back to the original array
            a(first:last) = a_copy(b)

            ! Deallocate temporary arrays
            deallocate(b, a_copy)
        end subroutine quicksort_int
        
        ! 17.5 Reorder from given order
        subroutine reorder(array,  order, n)
            real(kind=8), dimension(:), intent(inout) :: array
            integer(kind=4), dimension(:), intent(in) :: order
            integer(kind=4), intent(in) :: n
            real(kind=8), dimension(n) :: tmp
            integer(kind=4) :: i

            ! Loop and reorder
            do i = 1, n
                tmp(i) = array(order(i))
            end do
            array(1:n) = tmp
        end subroutine reorder

        ! 17.5 Reorder from given order
        subroutine reorder_int(array,  order, n)
            integer(kind=4), dimension(:), intent(inout) :: array
            integer(kind=4), dimension(:), intent(in) :: order
            integer(kind=4), intent(in) :: n
            integer(kind=4), dimension(n) :: tmp
            integer(kind=4) :: i

            ! Loop and reorder
            do i = 1, n
                tmp(i) = array(order(i))
            end do
            array(1:n) = tmp
        end subroutine reorder_int

        ! 17.6 Reorder 2D from given order
        subroutine reorder2D(array, order, n)
            real(kind=8), dimension(:,:), intent(inout) :: array
            integer(kind=4), dimension(:), intent(in) :: order
            integer(kind=4), intent(in) :: n
            real(kind=8), dimension(n, size(array, 2)) :: tmp
            integer(kind=4) :: i
        
            ! Loop and reorder the first axis based on the order array
            do i = 1, n
                tmp(i, :) = array(order(i), :)
            end do
        
            ! Copy the reordered array back
            array(1:n, :) = tmp
        end subroutine reorder2D
        
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

        ! 22.1 Get coordinates from a body
        subroutine coord(msum, a, e, inc, capm, omega, capom, xc)
            implicit none
            real(kind=8), intent(in)  :: msum, a, e, inc, capm, omega, capom
            real(kind=8), intent(out) :: xc(6)
            real(kind=8) :: sp, cp, so, co, si, ci
            real(kind=8) :: d11, d12, d13, d21, d22, d23
            real(kind=8) :: cape, dummy, scap, ccap, sqe, sqgma
            real(kind=8) :: ri, xfac1, xfac2, vfac1, vfac2
            
            ! Generate rotation matrices (on p. 42 of Fitzpatrick)
            sp = sin(omega)
            cp = cos(omega)
            so = sin(capom)
            co = cos(capom)
            si = sin(inc)
            ci = cos(inc)
            d11 = cp * co - sp * so * ci
            d12 = cp * so + sp * co * ci
            d13 = sp * si
            d21 = -sp * co - cp * so * ci
            d22 = -sp * so + cp * co * ci
            d23 = cp * si
            
            ! Get the other quantities depending on orbit type (i.e. ialpha)
            call aver(capm, e, cape, dummy)
            scap = sin(cape)
            ccap = cos(cape)
            sqe = sqrt(uno - e * e)
            sqgma = sqrt(G * msum * a)
            xfac1 = a * (ccap - e)
            xfac2 = a * sqe * scap
            ri = uno / (a * (uno - e * ccap))
            vfac1 = -ri * sqgma * scap
            vfac2 = ri * sqgma * sqe * ccap
            
            xc(1) = d11 * xfac1 + d21 * xfac2
            xc(2) = d12 * xfac1 + d22 * xfac2
            xc(3) = d13 * xfac1 + d23 * xfac2
            xc(4) = d11 * vfac1 + d21 * vfac2
            xc(5) = d12 * vfac1 + d22 * vfac2
            xc(6) = d13 * vfac1 + d23 * vfac2
        end subroutine coord

        ! 22.2 Get true f and g from a body
        subroutine aver(dm, e, u, f)
            implicit none
            real(kind=8), intent(in)  :: dm, e
            real(kind=8), intent(out) :: u, f
            real(kind=8) :: u0, dif, seno, cose
            
            u0 = dm
            dif = uno
            do while (dif > epsilon)
                u = dm + e * sin(u0)
                dif = abs(u - u0)
                u0 = u
            end do
            seno = sqrt(uno + e) * sin(uno2 * u)
            cose = sqrt(uno - e) * cos(uno2 * u)
            f = dos * atan2(seno, cose)
        end subroutine aver

        ! 22.3 Get elements from a body
        subroutine elem(msum, xc, a, e, inc, capm, omega, capom)
            implicit none
            real(kind=8), intent(in)  :: msum, xc(6)
            real(kind=8), intent(out) :: a, e, inc, capm, omega, capom
            real(kind=8) :: gmsum
            real(kind=8) :: x, y, z, vx, vy, vz
            real(kind=8) :: hx, hy, hz, h2, h, fac, u
            real(kind=8) :: r, v, v2, vdotr, energy
            real(kind=8) :: cape, cw, sw, w, face, capf, tmpf
            integer(kind=4) :: ialpha = 0
                
            gmsum = G * msum
            x = xc(1)
            y = xc(2)
            z = xc(3)
            vx = xc(4)
            vy = xc(5)
            vz = xc(6)
            
            hx = y * vz - z * vy
            hy = z * vx - x * vz
            hz = x * vy - y * vx
            h2 = hx * hx + hy * hy + hz * hz
            h = sqrt(h2)
            inc = acos(hz / h)
            
            fac = sqrt(hx*hx + hy*hy) / h
            if (fac < epsilon) then
                capom = cero
                u = atan2(y, x)
                if (abs(inc - pi) < 10.d0 * epsilon) then
                    u = -u
                end if
            else
                capom = atan2(hx, -hy)
                u = atan2(z / sin(inc), x * cos(capom) + y * sin(capom))
            end if
            
            if (capom < cero) then
                capom = capom + twopi
            end if
            if (u < cero) then
                u = u + twopi
            end if
            
            r = sqrt(x * x + y * y + z * z)
            v2 = vx * vx + vy * vy + vz * vz
            v = sqrt(v2)
            vdotr = x * vx + y * vy + z * vz
            energy = uno2 * v2 - gmsum / r
            
            if (abs(energy * r / gmsum) < sqrt(epsilon)) then
                ialpha = 0
            else
                if (energy < cero) then
                    ialpha = -1
                else if (energy > cero) then
                    ialpha = 1
                end if
            end if
            
            !! Ellipse
            if (ialpha == -1) then
                a = -uno2 * gmsum / energy
                fac = uno - h2 / (gmsum * a)
                if (fac > epsilon) then
                    e = sqrt(fac)
                    face = (a - r) / (a * e)
                    if (face > uno) then
                        cape = cero
                    else
                        if (face > -uno) then
                            cape = acos(face)
                        else
                            cape = pi
                        end if
                    end if
                    if (vdotr < cero) then
                        cape = twopi - cape
                    end if
                    cw = (cos(cape) - e) / (uno - e * cos(cape))
                    sw = sqrt(uno - e * e) * sin(cape) / (uno - e * cos(cape))
                    w = atan2(sw, cw)
                    if (w < cero) then
                        w = w + twopi
                    end if
                else
                    e = cero
                    w = u
                    cape = u
                end if
                capm = cape - e * sin(cape)
                omega = u - w
                if (omega < cero) then
                    omega = omega + twopi
                end if
                omega = omega - int(omega / twopi) * twopi
            end if
            
            !! Hypérbola
            if (ialpha == 1) then
                a = uno2 * gmsum / energy
                fac = h2 / (gmsum * a)
                if (fac > epsilon) then
                    e = sqrt(uno + fac)
                    tmpf = (a + r) / (a * e)
                    if (tmpf < uno) then
                        tmpf = uno
                    end if
                    capf = log(tmpf + sqrt(tmpf * tmpf - uno))
                    if (vdotr < cero) then
                        capf = -capf
                    end if
                    cw = (e - cosh(capf)) / (e * cosh(capf) - uno)
                    sw = sqrt(e * e - uno) * sinh(capf) / (e * cosh(capf) - uno)
                    w = atan2(sw, cw)
                    if (w < cero) then
                        w = w + twopi
                    end if
                else
                    e = uno
                    tmpf = uno2 * h2 / gmsum
                    w = acos(dos * tmpf / r - uno)
                    if (vdotr < cero) then
                        w = twopi - w
                    end if
                    tmpf = (a + r) / (a * e)
                    capf = log(tmpf + sqrt(tmpf * tmpf - uno))
                end if
                capm = e * sinh(capf) - capf
                omega = u - w
                if (omega < cero) then
                    omega = omega + twopi
                end if
                omega = omega - int(omega / twopi) * twopi
            end if

            !! Parábola    
            if (ialpha == 0) then
                a = uno2 * h2 / gmsum
                e = uno
                w = acos(dos * a / r - uno)
                if (vdotr < cero) then
                    w = twopi - w
                end if
                tmpf = tan(uno2 * w)
                capm = tmpf * (uno + tmpf * tmpf / 3.d0)
                omega = u - w
                if (omega < cero) then 
                    omega = omega + twopi
                end if
                omega = omega - int(omega / twopi) * twopi
            end if
        end subroutine elem

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

        ! 26.1 Go from lower to upper case
        function to_upper(strIn) result(strOut)
            ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
            ! Original author: Clive Page
            implicit none
            character(len=*), intent(in) :: strIn
            character(len=len(strIn)) :: strOut
            integer :: i,j
    
            do i = 1, len(strIn)
                j = iachar(strIn(i:i))
                if (j>= iachar("a") .and. j<=iachar("z") ) then
                    strOut(i:i) = achar(iachar(strIn(i:i))-32)
                else
                    strOut(i:i) = strIn(i:i)
                end if
            end do
        end function to_upper

        ! 26.2 Go from lower to upper case
        function to_lower(strIn) result(strOut)
            implicit none
            character(len=*), intent(in) :: strIn
            character(len=len(strIn)) :: strOut
            integer :: i,j
    
            do i = 1, len(strIn)
                j = iachar(strIn(i:i))
                if (j>= iachar("A") .and. j<=iachar("Z") ) then
                    strOut(i:i) = achar(iachar(strIn(i:i))+32)
                else
                    strOut(i:i) = strIn(i:i)
                end if
            end do
        end function to_lower

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
