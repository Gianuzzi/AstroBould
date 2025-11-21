!> Module with main parameters and implementation routines.
module parameters
    use constants
    use auxiliary
    use bodies, only: system_st, asteroid_st, moon_st, particle_st
    use accelerations, only: init_ellipsoid, init_manual_J2, init_stokes, init_drag, init_damping 
    use filtering
    use omp_lib
    use tomodule

    implicit none

    !! ----  <<<<<    GLOBAL     >>>>>   -----
    integer(kind=4) :: simulation_number = 1
    logical :: keep_integrating = .False.
    !! I/O
    integer(kind=4) :: arguments_number = 0  ! Amount of argument in command line
    logical :: configfile_exists = .False.  ! Wether if config file exists
    logical :: use_configfile = .True.  ! Wether to use config file (if exists)
    !! Parallel
    integer(kind=4) :: available_threads = 1  ! Available threads to use
    integer(kind=4) :: my_threads = 1  ! Amount of threads actually used
    logical :: compiled_with_openmp = .False.  ! Flag of compile with OpenMP
    logical :: any_extra_effect = .False.  ! If any extra is loaded. Just for message.
    logical :: use_flush_chaos = .False.  ! Flush chaos at every checkpoint
    logical :: use_flush_output = .False.  ! Flush data at every checkpoint


    !! ----  <<<<<    BODIES (bodies)     >>>>>   -----
    type(particle_st), allocatable :: particles_arr(:)
    type(moon_st), allocatable :: moons_arr(:)
    type(asteroid_st) :: asteroid
    type(system_st) :: initial_system
    type(system_st) :: system
    type(system_st) :: system_filtered
    
    
    !! ----  <<<<<    GENERAL PARAMETERS     >>>>>   -----
    type :: input_params_st  !!! This contains only the input parameters
        ! Something to do?
        logical :: only_print = .False.    
        ! Times for the integration - 
        real(kind=8) :: final_time = cero
        real(kind=8) :: output_timestep = cero
        integer(kind=4) :: output_number = 0
        integer(kind=4) :: case_output_type = 2
        integer(kind=4) :: extra_checkpoints = 0
        ! Extra parameters for the Integration - 
        logical :: use_parallel = .False.
        integer(kind=4) :: requested_threads = 1
        ! Integrator configuration - 
        integer(kind=4) :: integrator_ID = 0
        logical :: use_adaptive = .True.
        integer(kind=4) :: error_digits = 12
        real(kind=8) :: learning_rate = uno
        real(kind=8) :: dt_min = cero
        ! Primary mass -
        real(kind=8) :: mass_primary = cero
        ! Primary shape -
        real(kind=8) :: radius_primary = cero
        logical :: use_triaxial = .False.
        real(kind=8) :: triax_a_primary = cero
        real(kind=8) :: triax_b_primary = cero
        real(kind=8) :: triax_c_primary = cero
        ! Primary rotation -
        real(kind=8) :: lambda_kep = cero
        real(kind=8) :: asteroid_rotational_period = cero
        ! Moons - 
        logical :: use_moonsfile = .False.
        character(30) :: moonsfile = ""
        ! Particles - 
        logical :: use_particlesfile = .False.
        character(30) :: particlesfile = ""
        ! Extra forces/effects - 
        logical :: use_moon_gravity = .True.
        logical :: use_manual_J2 = .False.
        real(kind=8) :: manual_J2 = cero
        logical :: use_stokes = .False.
        real(kind=8) :: stokes_a_damping_time = infinity
        real(kind=8) :: stokes_e_damping_time = infinity
        real(kind=8) :: stokes_active_time = cero
        logical :: use_drag = .False.
        real(kind=8) :: drag_coefficient = cero
        real(kind=8) :: drag_active_time = cero
        logical :: use_omega_damping = .False.
        real(kind=8) :: omega_lin_damping_time = infinity
        real(kind=8) :: omega_exp_damping_time = infinity
        real(kind=8) :: omega_exp_damp_poly_A = cero
        real(kind=8) :: omega_exp_damp_poly_B = cero
        real(kind=8) :: omega_damp_active_time = cero
        real(kind=8) :: mass_exp_damping_time = infinity
        ! Filter
        logical :: use_filter = .False.
        real(kind=8) :: filter_dt = cero
        integer(kind=4) :: filter_nsamples = 0
        integer(kind=4) :: filter_nwindows = 0
        logical :: filter_use_KH = .False.
        integer(kind=4) :: filter_model = 0
        character(30) :: filter_prefix = ""
        ! [BINS] forces/effects - [NOT AVAILABLE YET. STILL UNDER DEVELOPMENT]
        logical :: use_self_gravity = .False.
        integer(kind=4) :: Norder_self_gravity = 3
        integer(kind=4) :: Nbins = 0
        integer(kind=4) :: binning_method = 1
        real(kind=8) :: rmin_bins = - uno ! -1 means first particle (initial)
        real(kind=8) :: rmax_bins = - uno ! -1 means last particle (initial)
        ! Conditions for Collision/Escape - 
        real(kind=8) :: min_distance = cero
        real(kind=8) :: max_distance = cero  ! Means no check
        logical :: use_merge_part_mass = .True.
        logical :: use_merge_massive = .True.
        logical :: use_stop_no_part_left = .True.
        logical :: use_stop_no_moon_left = .True.
        real(kind=8) :: eta_col = uno  ! 0: Elastic, 1: Plastic
        real(kind=8) :: f_col = cero  ! Bounded: Etot < -f |Epot|
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
        logical :: use_geometricfile = .False.
        character(30) :: geometricfile = ""
        logical :: use_datascreen = .False.
        logical :: use_diagnostics = .True.
        logical :: use_percentage = .False.
        logical :: use_elements_output = .True.
        logical :: use_baryc_output = .True.
        logical :: use_potential_map = .False.
        character(30) :: mapfile = ""
        integer(kind=4) :: map_grid_size_x = 500
        integer(kind=4) :: map_grid_size_y = 500
        real(kind=8) :: map_min_x = -500.d0
        real(kind=8) :: map_max_x = 500.d0
        real(kind=8) :: map_min_y = -500.d0
        real(kind=8) :: map_max_y = 500.d0
        ! < Number of bodies >
        integer(kind=4) :: Nboulders = 0
        integer(kind=4) :: Nmoons = 0
        integer(kind=4) :: Nparticles = 0
        !!! Command line body
        logical :: use_command_body = .False.  ! Wether a body was given through command line
    end type input_params_st

    type, extends(input_params_st) :: sim_params_st  !! Extra DERIVED parameters
        ! Times
        integer(kind=4) :: checkpoint_number = 0
        real(kind=8) :: min_period = infinity  ! This helps to create minimum timestep
        ! Bodies
        integer(kind=4) :: Ntotal = 0  ! Amount of all bodies (asteroid is just 1)
        integer(kind=4) :: Npart_active = 0 ! Active particles
        integer(kind=4) :: Nmoon_active = 0 ! Active moons
        integer(kind=4) :: Nactive = 0 ! Active bodies (including asteroid)
        logical :: use_boulders = .False.
        logical :: use_particles = .False.
        logical :: use_moons = .False.
        !! Stokes
        !! Mass and omega damping
        logical :: use_lin_omega_damp = .False.
        logical :: use_exp_omega_damp = .False.
        logical :: use_poly_omega_damp = .False.
        !! Filter
        integer(kind=4) :: filter_size = 0
        real(kind=8) :: filter_half_width = cero
        !! [BINS] ( Not available yet )
        logical :: use_bins = .False.
        !!! Viscosity [Not available yet]
        logical :: use_viscosity = .False.
        real(kind=8) :: viscosity = - uno
        !!! Self-Gravity
        logical :: update_rmin_bins = .False.
        logical :: update_rmax_bins = .False.
        logical :: update_bins = .False.
        ! Merges
        logical :: use_any_merge = .False.
        ! Stops
        logical :: use_any_stop = .False.
        ! Chaos
        logical :: use_chaos = .False.  ! Need to output chaos values
        character(30) :: geomchaosfile = ""
        logical :: use_geomchaosfile = .False.  ! Need to output geometric chaos values
        ! Error
        real(kind=8) :: error_tolerance = cero
        !!! Command line body
        logical :: cl_body_is_moon = .False. ! Wether a moon was given through command line
    end type sim_params_st
    
    type(input_params_st) :: input  ! This is the structure with the default -> input parameters
    type(sim_params_st) :: sim  ! This is the structure with the input -> derived parameters


    ! ----  <<<<<    INITIAL BODIES     >>>>>   -----
    real(kind=8), dimension(:,:), allocatable :: boulders_in !! mu, radius, theta  | (Nb, 3)
    real(kind=8), dimension(:,:), allocatable :: moons_in !! mass, a, e, M, w, MMR, radius  | (Nm, 7)
    real(kind=8), dimension(:,:), allocatable :: particles_in !! a, e, M, w, MMR  | (Np, 5)
    real(kind=8) :: cl_body_in(7) = cero  !! mass, a, e, M, w, MMR, radius  | (7)  | COMMAND LINE Input


    ! ----  <<<<<    BOULDERS for DYDT     >>>>>   -----
    real(kind=8), dimension(:,:), allocatable :: boulders_coords !! (Nb, 4) |x,y,vx,vy|
    real(kind=8), dimension(:,:), allocatable :: boulders_data !! (Nb, 4) |mass,radius,theta_Ast0,dist_Ast|


    ! ----  <<<<<    TIME     >>>>>   -----
    real(kind=8) :: time  ! Actual time of the integration
    real(kind=8) :: timestep  ! This timestep 
    real(kind=8) :: adaptive_timestep  ! This adaptive timestep
    real(kind=8) :: fixed_timestep  ! This fixed timestep
    

    ! ----  <<<<<    TOM     >>>>>   -----
    type(tom_st) :: tom  ! This is the structure with TOM data

    ! ----  <<<<<    CHECKPOINTS     >>>>>   -----
    real(kind=8), dimension(:), allocatable :: output_times  ! Vector con tiempos de salida
    real(kind=8), dimension(:), allocatable :: checkpoint_times  ! Vector con checkpoints
    logical, dimension(:), allocatable :: checkpoint_is_tom  ! Array with whether each checkpoint is TOM
    logical, dimension(:), allocatable :: checkpoint_is_output ! Array with whether each checkpoint is output
        

    ! ----  <<<<<    PARAMETERS ARRAYS     >>>>>   -----
    integer, parameter :: equation_size = 4
    integer(kind=4) :: y_nvalues = 0  ! Will change dynamically
    real(kind=8), dimension(:), allocatable :: m_arr         ! Mass array
    real(kind=8), dimension(:), allocatable :: R_arr         ! Radius array
    real(kind=8), dimension(:), allocatable :: y_arr         ! Coordinates array
    real(kind=8), dimension(:), allocatable :: y_arr_new     ! Coordinates array (2.0)
    real(kind=8), dimension(:), allocatable :: y_der         ! Derivate of coordinates array


    ! ----  <<<<<    FILTERING     >>>>>   -----
    type(filter_st) :: filter
    integer(kind=4) :: last_idx_no_filter = -1  ! Last checkpoint with no filtering
    integer(kind=4) :: first_idx_yes_filter = 0   ! First checkpoint with no filtering
    real(kind=8), dimension(:), allocatable :: elem_filtered  ! Filtered_data
    real(kind=8), dimension(:), allocatable :: y_pre_filter   ! Data pre-filtering
    ! ==========    EXTRA FILTER    ==========
    real(kind=8) :: next_time  ! Next time
    real(kind=8) :: next_t_filt  ! Next filter time
    real(kind=8) :: time_filt  ! Actual time of the filtering integration
    real(kind=8) :: adaptive_timestep_filt  ! Adaptive timestep of the filtering integration
    integer(kind=4) :: last_output  ! Last output checkpoint
    integer(kind=4) :: next_output  ! Next output checkpoint
    integer(kind=4) :: next_checkpoint  ! Next pure checkpoint
    integer(kind=4) :: tmp_j  ! Temporal j checkpoint
    type(sim_params_st) :: tmp_sim  ! Temporal simulation state
    real(kind=8) :: tmp_time  ! Temporal time
    real(kind=8) :: tmp_adaptive_timestep  ! Temporal time
    type(system_st) :: tmp_system  ! Temporal system
    integer(kind=4) :: tmp_y_nvalues  ! Temporal y nvalues
    real(kind=8), dimension(:), allocatable :: tmp_y_arr  ! Temporal coordinates array
    
    

    ! ----  <<<<<    HARD EXIT     >>>>>   -----
    logical :: is_premature_exit = .False.
    logical :: hard_exit = .False. !  Hard Exit logical
    ! procedure (check_continue_template), pointer :: check_continue_ptr => null ()
    

    ! ----  <<<<<    I/O     >>>>>   -----
    character(18), parameter :: s1i1 = "(3(A, 1X, I7, 1X))"
    character(24), parameter :: s1r1 = "(3(A, 1X, 1PE22.15, 1X))"
    character(21), parameter :: i1r5 = "(I7, 5(1X, 1PE22.15))"
    character(21), parameter :: i1r7 = "(I7, 7(1X, 1PE22.15))"


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    POINTERS     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    procedure (write_to_unit_template), pointer :: write_to_general => null ()
    procedure (write_to_unit_template), pointer :: write_to_screen => null ()
    procedure (write_ch_to_unit_template), pointer :: write_to_chaos => null ()
    procedure (write_to_unit_template), pointer :: write_a_to_individual => null ()
    procedure (write_i_to_unit_template), pointer :: write_m_to_individual => null ()
    procedure (write_i_to_unit_template), pointer :: write_p_to_individual => null ()
    procedure (write_to_unit_template), pointer :: write_to_geometric => null ()
    procedure (write_ch_to_unit_template), pointer :: write_to_geomchaos => null ()
    procedure (int_i_template), pointer :: flush_output => null ()
    procedure (int_i_template), pointer :: flush_chaos => null ()
    procedure (int_i_template), pointer :: flush_geometric => null ()
    procedure (int_i_template), pointer :: flush_geomchaos => null ()



    abstract interface
    
        subroutine simple_ptr_template ()
            implicit none
        end subroutine simple_ptr_template

        subroutine int_i_template (i)
            implicit none
            integer(kind=4), intent(in) :: i
        end subroutine int_i_template
        subroutine write_i_to_unit_template (self, i, unit)
            import :: system_st     ! bring type into scope
            implicit none
            type(system_st), intent(in) :: self
            integer(kind=4), intent(in) :: i, unit
        end subroutine write_i_to_unit_template

        subroutine write_to_unit_template (self, unit)
            import :: system_st     ! bring type into scope
            implicit none
            type(system_st), intent(in) :: self
            integer(kind=4), intent(in) :: unit
        end subroutine write_to_unit_template

        subroutine write_ch_to_unit_template (self, other, unit)
            import :: system_st     ! bring type into scope
            implicit none
            type(system_st), intent(in) :: self
            type(system_st), intent(in) :: other
            integer(kind=4), intent(in) :: unit
        end subroutine write_ch_to_unit_template

    end interface

    contains
    
        ! Allocate asteroid params arrays
        subroutine allocate_params_asteroid(Nboulders)
            implicit none
            integer(kind=4), intent(in) :: Nboulders

            if (allocated(boulders_in)) then
                write(*,*) "ERROR: Number of boulder already defined."
                STOP 1
            end if

            allocate(boulders_in(Nboulders,3))
            !! Allocate params arrays
        end subroutine allocate_params_asteroid

        ! Allocate moons params arrays
        subroutine allocate_params_moons(Nmoons)
            implicit none
            integer(kind=4), intent(in) :: Nmoons

            if (allocated(moons_in)) then
                write(*,*) "ERROR: Number of moons already defined."
                STOP 1
            end if

            allocate(moons_in(Nmoons,7))  ! FROM 1 to Nmoons
        end subroutine allocate_params_moons

        ! Allocate particles params arrays
        subroutine allocate_params_particles(Nparticles)
            implicit none
            integer(kind=4), intent(in) :: Nparticles

            if (allocated(particles_in)) then
                write(*,*) "ERROR: Number of particles already defined."
                STOP 1
            end if

            allocate(particles_in(Nparticles,5))  ! FROM 1 to Nparticles
        end subroutine allocate_params_particles
        

        ! 1. Read command line input
        subroutine load_command_line_arguments(params, use_config)
            implicit none
            type(input_params_st), intent(inout) :: params
            logical, intent(inout)  :: use_config
            logical :: is_number
            integer(kind=4) :: i, j
            logical :: aux_logical = .False.
            integer(kind=4) :: aux_integer
            integer(kind=4) :: merge_type = 0, stop_type = 0
            character(20) :: aux_character20
            character(30) :: aux_character30

            ! Leemos de la línea de comandos
            arguments_number = command_argument_count()  ! Global variable

            aux_logical = .False. ! Leí los parámetros numéricos?
            aux_integer = 0
            !!!!! PARTÍCULA (en caso de entrada por terminal)
            if (arguments_number > 0) then
                do i = 1, arguments_number
                    if (aux_integer /= 0) then
                        aux_integer = aux_integer - 1
                        cycle
                    end if
                    call get_command_argument(i, aux_character30)
                    select case (trim(aux_character30))
                    case ("--onlyprint")
                        params%only_print = .True.
                    case ("-nsim")
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) simulation_number  ! Global variable
                        aux_integer = 1
                    case ("-mumoon")
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) cl_body_in(1)  ! Global variable
                        aux_integer = 1
                    case ("-rmoon")
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) cl_body_in(7)  ! Global variable
                        aux_integer = 1
                    case ("-datafile")
                        params%use_datafile = .True.
                        call get_command_argument(i+1, params%datafile)
                        aux_integer = 1
                    case ("--nodataf")
                        params%use_datafile = .False.
                        params%datafile = ""
                    case ("-chaosfile")
                        params%use_chaosfile = .True.
                        call get_command_argument(i+1, params%chaosfile)
                        aux_integer = 1
                    case ("--nochaosf")
                        params%use_chaosfile = .False.
                        params%chaosfile = ""
                    case ("-geomfile")
                        params%use_geometricfile = .True.
                        call get_command_argument(i+1, params%geometricfile)
                        aux_integer = 1
                    case ("--nogeomf")
                        params%use_geometricfile = .False.
                        params%geometricfile = ""
                    case ("--screen")
                        params%use_screen = .True.
                    case ("--noscreen")
                        params%use_screen = .False.
                    case ("--perc")
                        params%use_percentage = .True.
                        params%use_diagnostics = .False.
                        params%use_datascreen = .False.
                    case ( "--noperc")
                        params%use_percentage = .False.
                    case ("--datascr")
                        params%use_datascreen = .True.
                        params%use_percentage = .False.
                        params%use_diagnostics = .False.
                    case ( "--nodatascr")
                        params%use_datascreen = .False.
                    case ("--diagnostic")
                        params%use_diagnostics = .True.
                        params%use_datascreen = .False.
                        params%use_percentage = .False.
                    case ( "--nodiagnostic")
                        params%use_diagnostics = .False.
                    case ("-multifile")
                        params%use_multiple_outputs = .True.
                        call get_command_argument(i+1, params%multfile)
                        aux_integer = 1
                    case ("--nomultif")
                        params%use_multiple_outputs = .False.
                        params%multfile = ""
                    case ("-mapfile")
                        params%use_potential_map = .True.
                        call get_command_argument(i+1, params%mapfile)
                        aux_integer = 1
                    case ("--nomapf")
                        params%use_potential_map = .False.
                        params%mapfile = ""
                    case ("--elem")
                        params%use_elements_output = .True.
                    case ("--noelem")
                        params%use_elements_output = .False.
                    case ("-tomfile")
                        params%use_tomfile = .True.
                        call get_command_argument(i+1, params%tomfile)
                        aux_integer = 1
                    case ("--notomfile")
                        params%use_tomfile = .False.
                        params%tomfile = ""
                    case ("-moonfile")
                        params%use_moonsfile = .True.
                        call get_command_argument(i+1, params%moonsfile)
                        aux_integer = 1
                    case ("--moonfile")
                        params%use_moonsfile = .False.
                        params%moonsfile = ""
                    case ("-partfile")
                        params%use_particlesfile = .True.
                        call get_command_argument(i+1, params%particlesfile)
                        aux_integer = 1
                    case ("--nopartfile")
                        params%use_particlesfile = .False.
                        params%particlesfile = ""
                    case ("--noconfig")
                        use_config = .False.  ! Inout variable
                    case ("-merge")
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) merge_type
                        if (merge_type == 0) then ! No merge
                            params%use_merge_part_mass = .False.
                            params%use_merge_massive = .False.
                        else if (merge_type == 1) then ! Only particle-massive
                            params%use_merge_part_mass = .True.
                            params%use_merge_massive = .False.
                        else if (merge_type == 2) then ! Only massive-massive
                            params%use_merge_part_mass = .False.
                            params%use_merge_massive = .True.
                        else if (merge_type == 3) then ! All merges
                            params%use_merge_part_mass = .True.
                            params%use_merge_massive = .True.
                        else
                            write(*,*) "ERROR: Merge type not recognized: ", merge_type
                            stop 1
                        end if
                        aux_integer = 1
                    case ("-stopif")
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) stop_type
                        if (stop_type == 0) then ! No stop
                            params%use_stop_no_moon_left = .False.
                            params%use_stop_no_part_left = .False.
                        else if (stop_type == 1) then ! Only if no more moons
                            params%use_stop_no_moon_left = .True.
                            params%use_stop_no_part_left = .False.
                        else if (stop_type == 2) then ! Only if no more particles
                            params%use_stop_no_moon_left = .False.
                            params%use_stop_no_part_left = .True.
                        else if (stop_type == 3) then ! Stop if any deactivation
                            params%use_stop_no_moon_left = .True.
                            params%use_stop_no_part_left = .True.
                        else
                            write(*,*) "ERROR: Stop criteria not recognized: ", stop_type
                            stop 1
                        end if
                        aux_integer = 1
                    case ("-parallel")
                        params%use_parallel = .True.
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) params%requested_threads
                        aux_integer = 1
                    case ("--parallel")
                        params%use_parallel = .True.
                        params%requested_threads = -1
                    case ("--noparallel")
                        params%use_parallel = .False.
                        params%requested_threads = 1
                    case ("--help")
                        call get_command_argument(0, aux_character30)
                        write(*,*) "Uso: " // trim(aux_character30) // " <ea> <ee> <eM> <ew> <eR> [args]"
                        write(*,*) "    ea  : Elemento a de la partícula/luna (km)"
                        write(*,*) "    ee  : Elemento e de la partícula/luna"
                        write(*,*) "    eM  : Elemento M de la partícula/luna (deg)"
                        write(*,*) "    ew  : Elemento w de la partícula/luna (deg)"
                        write(*,*) "    eR  : Elemento R de la partícula/luna [Opcional]"
                        write(*,*) "    --onlyprint   : No integrar; solo imprimir configuraciones"
                        write(*,*) "    -nsim         : Número de simulación [int]"
                        write(*,*) "    -mumoon       : Cociente de masa entre la luna individual y el asteroide"
                        write(*,*) "    -rmoon        : Radio de la luna individual (km). Solo si mumoon > 0"
                        write(*,*) "    -datafile     : Nombre de archivo de salida de datos"
                        write(*,*) "    --nodataf     : No guardar datos de salida"
                        write(*,*) "    -chaosfile    : Nombre de archivo de salida caos"
                        write(*,*) "    --nochaosf    : No guardar salida de caos"
                        write(*,*) "    -geomfile     : Nombre de archivo de salida de elementos geométricos"
                        write(*,*) "    --nogeomf     : No guardar salida de elementos geométricos"
                        write(*,*) "    --screen      : Imprimir información en pantalla"
                        write(*,*) "    --noscreen    : No imprimir en pantalla"
                        write(*,*) "    --perc        : Imprimir porcentaje de integración"
                        write(*,*) "    --noperc      : No imprimir porcentaje de integración"
                        write(*,*) "    --datascr     : Imprimir datos de salida en pantalla"
                        write(*,*) "    --nodatascr   : No imprimir datos de salida en pantalla"
                        write(*,*) "    --diagnostic  : Imprimir datos de diagnóstico en pantalla"
                        write(*,*) "    --nodiagnostic: No imprimir datos de diagnóstico en pantalla"
                        write(*,*) "    -multifile    : Nombre base de archivo de salida de datos individuales"
                        write(*,*) "    --nomultif    : No guardar datos en archivos individuales"
                        write(*,*) "    -mapfile      : Nombre de archivo de mapa"
                        write(*,*) "    --nomapf      : No guardar mapa de potencial"
                        write(*,*) "    --elem        : Imprimir elementos orbitales (lunas/partículas) [default]"
                        write(*,*) "    --noelem      : Imprimir coordenadas baricéntricas"
                        write(*,*) "    -tomfile      : Nombre de archivo de (t)iempos|omega|masa a utilizar"
                        write(*,*) "    --notomfile   : No utilizar archivo de (t)iempos|omega|masa"
                        write(*,*) "    -moonfile     : Nombre de archivo de lunas a utilizar"
                        write(*,*) "    --nomoonfile  : No utilizar archivo de lunas"
                        write(*,*) "    -partfile     : Nombre de archivo de partículas a utilizar"
                        write(*,*) "    --nopartfile  : No utilizar archivo de partículas"
                        write(*,*) "    --noconfig    : No leer archivo de configuración"
                        write(*,*) "    -merge        : Tipo de colisiones (merges) permitidas [int]: "
                        write(*,*) "                    0: Ninguno, 1: Partícula-Masivo, 2: Masivo-Masivo, 3: Todos"
                        write(*,*) "    -stopif       : Detener la integración si no quedan más objetos del tipo [int]:"
                        write(*,*) "                    0: No detener, 1: Luna, 2: Partícula, 3: Ambos"
                        write(*,*) "    -parallel     : Cantida de thread a utilizar en paralelo [int]"
                        write(*,*) "    --parallel    : Paralelizar usando todos los threads disponibles"
                        write(*,*) "    --noparallel  : No usar paralelización para lunas/partículas"
                        write(*,*) "    --help        : Mostrar esta ayuda"
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
                            write(*, '(A, I0, A, A)') "ERROR: Argument not recognized ", i, ": ", trim(aux_character30)
                            call get_command_argument(0, aux_character30)
                            write(*,*) "   For more help run: ", trim(aux_character30), " --help"
                            write(*,*) "Exiting."
                            stop 1
                        end if
                        if (.not. aux_logical)  then ! No leí los parámetros aún
                            if (arguments_number < i+3) then
                                write(*,*) "ERROR: Could not read all particle orbital elements."
                                write(*,'(A,I0,A)') "        Missing ", 4-arguments_number, " elements."
                                write(*,*) "Exiting."
                                stop 1
                            else ! Leo los argumentos numéricos. Considero que es 1 sola partícula
                                ! Para generar prioridad, reallocatamos de ser necesario
                                params%use_command_body = .True.  ! Switch
                                ! Read
                                call get_command_argument(i, aux_character20)
                                read (aux_character20,*) cl_body_in(2)
                                call get_command_argument(i+1, aux_character20)
                                read (aux_character20,*) cl_body_in(3)
                                call get_command_argument(i+2, aux_character20)
                                read (aux_character20,*) cl_body_in(4)
                                call get_command_argument(i+3, aux_character20)
                                read (aux_character20,*) cl_body_in(5)
                                aux_integer = 3
                                aux_logical = .True.
                                cl_body_in(6) = cero
                            end if
                        else ! Ya leí los numéricos. Falta leer eR
                            call get_command_argument(i, aux_character20)
                            read (aux_character20,*) cl_body_in(6)
                        end if
                    end select
                end do
            end if  
        end subroutine load_command_line_arguments

        ! 2. Read config file
        subroutine read_config_file(params, file_name, existe)
            implicit none
            type(input_params_st), intent(inout) :: params
            character(LEN=*), intent(in) :: file_name
            logical, intent(inout) :: existe
            character(256) :: line, param_str, value_str
            character(1) :: auxch1
            character(2) :: auxch2
            character(15) :: auxch15
            character(15) :: auxch30
            integer(kind=4)  :: colonPos, commentPos
            integer(kind=4) :: nlines
            integer(kind=4) :: io
            integer(kind=4) :: len_val
            integer(kind=4) :: j
            integer(kind=4) :: aux_integer
            integer(kind=4) :: Nboulders, Nmoons, Nparticles

            ! CHECK IF 
            aux_integer = command_argument_count()
            do j = 1, aux_integer
                call get_command_argument(j, auxch30)
                if (trim(auxch30) == "--noconfig") then
                    use_configfile = .False.
                    return
                end if
            end do

            ! Init
            Nboulders = 0
            Nmoons = 0
            Nparticles = 0

            inquire(file=trim(file_name), exist=existe)
            if (.not. existe) then
                if (params%use_screen) write(*,*) "File ", trim(file_name), " does not exists."  ! use_screen is the default...
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
                        case("total integrati")
                            read (value_str, *) params%final_time
                        case("output time int")
                            read (value_str, *) params%output_timestep
                        case("number of outpu")
                            read (value_str, *) params%output_number
                        case("number of extra")
                            read (value_str, *) params%extra_checkpoints
                        case("output distribu")
                            read (value_str, *) params%case_output_type
                        case("use parallel th")
                            read (value_str,'(i10)',iostat=aux_integer) params%requested_threads
                            if (aux_integer /= 0) then
                                if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                    params%use_parallel = .True.
                                    params%requested_threads = -1
                                else
                                    params%use_parallel = .False.
                                    params%requested_threads = 1
                              endif
                            else
                                params%use_parallel = .True.
                            end if
                        case("integrator ID t")
                            read (value_str, *) params%integrator_ID
                        case("use adaptive ti")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_adaptive = .True.
                            else
                                params%use_adaptive = .False.
                            end if
                        case("precision (digi")
                            read (value_str, *) params%error_digits
                        case("learning rate")
                            read (value_str, *) params%learning_rate
                        case("minimum timeste")
                            read (value_str, *) params%dt_min
                        case("mass of primary")
                            read (value_str, *) params%mass_primary
                        case("radius of prima")
                            read (value_str, *) params%radius_primary
                        case("use triaxial mo")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_triaxial = .True.
                            else
                                params%use_triaxial = .False.
                            end if
                        case("semi-axis a (km")
                            read (value_str, *) params%triax_a_primary
                        case("semi-axis b (km")
                            read (value_str, *) params%triax_b_primary
                        case("semi-axis c (km")
                            read (value_str, *) params%triax_c_primary
                        case("use manual J2 (")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_manual_J2 = .True.
                            else
                                params%use_manual_J2 = .False.
                            end if
                        case("manual J2 value")
                            read (value_str, *) params%manual_J2
                        case("ratio of spin t")
                            read (value_str, *) params%lambda_kep
                        case("rotational peri")
                            read (value_str, *) params%asteroid_rotational_period
                        case("moons input fil")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                params%use_moonsfile = .False.
                                params%moonsfile = ""
                            else if (len(trim(value_str)) > 0) then
                                params%use_moonsfile = .True.
                                params%moonsfile = trim(value_str)
                            else
                                params%use_moonsfile = .False.
                                params%moonsfile = ""
                            end if
                        case("particles input")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                params%use_particlesfile = .False.
                                params%particlesfile = ""
                            else if (len(trim(value_str)) > 0) then
                                params%use_particlesfile = .True.
                                params%particlesfile = trim(value_str)
                            else
                                params%use_particlesfile = .False.
                                params%particlesfile = ""
                            end if
                        case("deactivate grav")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_moon_gravity = .False.
                            else
                                params%use_moon_gravity = .True.
                            end if
                        case("include stokes-")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_stokes = .True.
                            else
                                params%use_stokes = .False.
                            end if
                        case("a damping chara")
                            read (value_str, *) params%stokes_a_damping_time
                        case("e damping chara")
                            read (value_str, *) params%stokes_e_damping_time
                        case("stokes force da")
                            read (value_str, *) params%stokes_active_time
                        case("include naive-s")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_drag = .True.
                            else
                                params%use_drag = .False.
                            end if
                        case("drag force coef")
                            read (value_str, *) params%drag_coefficient
                        case("drag force acti")
                            read (value_str, *) params%drag_active_time
                        case("include asteroi")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_omega_damping = .True.
                            else
                                params%use_omega_damping = .False.
                            end if
                        case("rotation linear")
                            read (value_str, *) params%omega_lin_damping_time
                        case("rotation expone")
                            read (value_str, *) params%omega_exp_damping_time
                        case("rotation A poly")
                            read (value_str, *) params%omega_exp_damp_poly_A
                        case("rotation B poly")
                            read (value_str, *) params%omega_exp_damp_poly_B
                        case("omega damping a")
                            read (value_str, *) params%omega_damp_active_time
                        case("mass exponentia")
                            read (value_str, *) params%mass_exp_damping_time
                        case("use filtering (")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_filter = .True.
                            else
                                params%use_filter = .False.
                            end if
                        case("time period to")
                            read (value_str, *) params%filter_dt
                        case("filter oversamp")
                            read (value_str, *) params%filter_nsamples
                        case("filter windows")
                            read (value_str, *) params%filter_nwindows
                        case("apply filter to")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%filter_use_KH = .True.
                            else
                                params%filter_use_KH = .False.
                            end if
                        case("kernel weights")
                            read (value_str, *) params%filter_model
                        case("prefix for filt")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                params%filter_prefix = ""
                            else
                                params%filter_prefix = trim(value_str)
                            end if
                        case("include self-gr")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_self_gravity = .True.
                            else
                                params%use_self_gravity = .False.
                            end if
                        case("max Legendre ex")
                            read (value_str, *) params%Norder_self_gravity
                        case("total bins used")
                            read (value_str, *) params%Nbins
                        case("binning method ")
                            read (value_str, *) params%binning_method
                        case("inner bin edge")
                            read (value_str, *) params%rmin_bins
                        case("outer bin edge")
                            read (value_str, *) params%rmax_bins
                        case("min distance fr")
                            read (value_str, *) params%min_distance
                        case("max distance fr")
                            read (value_str, *) params%max_distance
                        case("particle-massiv")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_merge_part_mass = .True.
                            else
                                params%use_merge_part_mass = .False.
                            end if
                        case("massive merges ")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_merge_massive = .True.
                            else
                                params%use_merge_massive = .False.
                            end if
                        case("collisional eta")
                            read (value_str, *) params%eta_col
                        case("collisional f p")
                            read (value_str, *) params%f_col
                        case("stop if no part")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_stop_no_part_left = .True.
                            else
                                params%use_stop_no_part_left = .False.
                            end if
                        case("stop if no moon")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_stop_no_moon_left = .True.
                            else
                                params%use_stop_no_moon_left = .False.
                            end if
                        case("input time-omeg")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                params%use_tomfile = .False.
                                params%tomfile = ""
                            else if (len(trim(value_str)) > 0) then
                                params%use_tomfile = .True.
                                params%tomfile = trim(value_str)
                            else
                                params%use_tomfile = .False.
                                params%tomfile = ""
                            end if
                        case("information on")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                params%use_screen = .True.
                            else
                                params%use_screen = .False.
                            end if
                        case("general output")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                params%use_datafile = .False.
                                params%datafile = ""
                            else
                                params%use_datafile = .True.
                                params%datafile = trim(value_str)
                            end if
                        case("individual file")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                params%use_multiple_outputs = .False.
                                params%multfile = ""
                            else
                                params%use_multiple_outputs = .True.
                                params%multfile = trim(value_str)
                            end if
                        case("chaos indicator")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                params%use_chaosfile = .False.
                                params%chaosfile = ""
                            else
                                params%use_chaosfile = .True.
                                params%chaosfile = trim(value_str)
                            end if
                        case("geometric eleme")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                params%use_geometricfile = .False.
                                params%geometricfile = ""
                            else
                                params%use_geometricfile = .True.
                                params%geometricfile = trim(value_str)
                            end if
                        case("output data on")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                params%use_datascreen = .True.
                                params%use_percentage = .False.
                                params%use_diagnostics = .False.
                            else if (auxch1 == "%") then
                                params%use_datascreen = .False.
                                params%use_percentage = .True.
                                params%use_diagnostics = .False.
                            else if (auxch1 == "d") then
                                params%use_diagnostics = .True.
                                params%use_datascreen = .False.
                                params%use_percentage = .False.
                            else
                                params%use_diagnostics = .False.
                                params%use_percentage = .False.
                                params%use_datascreen = .False.
                            end if
                        case("output variable")
                            if (auxch1 == "c") then
                                params%use_elements_output = .False.
                            else 
                                params%use_elements_output = .True.
                            end if
                        case("output elements")
                            if (auxch1 == "b") then
                                params%use_baryc_output = .True.
                            else
                                params%use_baryc_output = .False.
                            end if
                        case("create map file")
                            if ((to_lower(trim(value_str)) == "n") .or. &
                              & (to_lower(trim(value_str)) == "no")) then
                                params%use_potential_map = .False.
                                params%mapfile = ""
                            else
                                params%use_potential_map = .True.
                                params%mapfile = trim(value_str)
                            end if
                        case("number x cells")
                            read (value_str, *) params%map_grid_size_x
                        case("number y cells")
                            read (value_str, *) params%map_grid_size_y
                        case("lower x bound (")
                            read (value_str, *) params%map_min_x
                        case("upper x bound (")
                            read (value_str, *) params%map_max_x
                        case("lower y bound (")
                            read (value_str, *) params%map_min_y
                        case("upper y bound (")
                            read (value_str, *) params%map_max_y
                        case default
                            write(*,*) "WARNING: Parameter not recongnized: ", trim(param_str)
                    end select
                else  ! Read boulders, moons or particles
                    param_str = trim(adjustl(line)) 
                    if (param_str(1:7) == "mass_m0") then  ! Leeemos boulders
                        aux_integer = nlines  ! Save the nlines where to read from (later)
                        Nboulders = 0 ! This will be the boulder count
                        io = 0  ! This is be the reading status
                        do while (io == 0)
                            nlines = nlines + 1  ! Don't forget to count the lines....
                            read (10, '(A)', iostat=io) line
                            auxch1 = trim(line)  ! get just first char
                            if (len(trim(line)) == 0) exit ! Stop at blank line
                            if ((auxch1 == "c") .or. (auxch1 == "!")) cycle  ! Skip commented
                            if (io == 0) Nboulders = Nboulders + 1
                        end do
                        ! Alocate and set
                        call allocate_params_asteroid(Nboulders)
                        params%Nboulders = Nboulders
                        ! Backspace to read again
                        do j = aux_integer, nlines - 1  ! -1 por ser both included
                            backspace (10)
                        end do
                        ! Read again and store
                        j = 1  ! This is be the amount of boulders read
                        do while (j <= params%Nboulders)
                            read (10, '(A)') line
                            auxch1 = adjustl(line)  ! get just first char
                            if ((auxch1 == "c") .or. (auxch1 == "!")) cycle  ! Skip commented
                            read (line, *, iostat=io) &
                                & boulders_in(j,1), &
                                & boulders_in(j,2), &
                                & boulders_in(j,3)
                            if (io /= 0) then
                                write(*,*) "ERROR: Al leer boulder:", j
                                stop 1
                            end if
                            j = j + 1
                        end do
                    end if
                    if (param_str(1:9) == "mass_mAst") then  ! Leeemos las lunas
                        aux_integer = nlines  ! Save the nlines where to read from (later)
                        Nmoons = 0 ! This will be the moon count
                        io = 0  ! This is be the reading status
                        do while (io == 0)
                            read (10, '(A)', iostat=io) line
                            nlines = nlines + 1  ! Don't forget to count the lines....
                            auxch1 = adjustl(line)  ! get just first char
                            if (len(trim(line)) == 0) exit ! Stop at blank line
                            if ((auxch1 == "c") .or. (auxch1 == "!")) cycle  ! Skip commented
                            if (io == 0) Nmoons = Nmoons + 1
                        end do
                        params%Nmoons = Nmoons
                        if (Nmoons > 0) then
                            ! Alocate
                            call allocate_params_moons(Nmoons)
                            ! Backspace to read again
                            do j = aux_integer, nlines - 1  ! -1 por ser both included
                                backspace (10)
                            end do
                            ! Read again and store
                            j = 1  ! This is be the amount of moons read
                            do while (j <= params%Nmoons)
                                read (10, '(A)') line
                                auxch1 = adjustl(line)  ! get just first char
                                if ((auxch1 == "c") .or. (auxch1 == "!")) cycle  ! Skip commented
                                ! Set MMR to 0
                                moons_in(j, 6) = cero
                                read (line, *, iostat=io) &
                                    & moons_in(j, 1), &
                                    & moons_in(j, 2), &
                                    & moons_in(j, 3), &
                                    & moons_in(j, 4), &
                                    & moons_in(j, 5), &
                                    & moons_in(j, 6), &
                                    & moons_in(j, 7)
                                if (io /= 0) then
                                    write(*,*) "ERROR: Al leer luna:", j
                                    stop 1
                                end if
                                j = j + 1
                            end do
                        end if
                    end if
                    if (param_str(1:5) == "a(km)") then  ! Leeemos las partículas
                        aux_integer = nlines  ! Save the nlines where to read from (later)
                        Nparticles = 0 ! This will be the particles count
                        io = 0  ! This is be the reading status
                        do while (io == 0)
                            read (10, '(A)', iostat=io) line
                            nlines = nlines + 1  ! Don't forget to count the lines....
                            auxch1 = adjustl(line)  ! get just first char
                            if (len(trim(line)) == 0) exit ! Stop at blank line
                            if ((auxch1 == "c") .or. (auxch1 == "!")) cycle  ! Skip commented
                            if (io == 0) Nparticles = Nparticles + 1
                        end do
                        params%Nparticles = Nparticles
                        if (Nparticles > 0) then
                            ! Alocate
                            call allocate_params_particles(Nparticles)
                            ! Backspace to read again
                            do j = aux_integer, nlines - 1  ! -1 por ser both included
                                backspace (10)
                            end do
                            ! Read again and store
                            j = 1  ! This is be the amount of particles read
                            do while (j <= Nparticles)
                                read (10, '(A)') line
                                auxch1 = adjustl(line)  ! get just first char
                                if ((auxch1 == "c") .or. (auxch1 == "!")) cycle  ! Skip commented
                                ! Set MMR to 0
                                particles_in(j, 5) = cero
                                read (line, *, iostat=io) &
                                    & particles_in(j, 1), &
                                    & particles_in(j, 2), &
                                    & particles_in(j, 3), &
                                    & particles_in(j, 4), &
                                    & particles_in(j, 5) 
                                if (io /= 0) then
                                    write(*,*) "ERROR: Al leer partícula:", j
                                    stop 1
                                end if
                                j = j + 1
                            end do
                        end if
                    end if
                end if
            end do
            98 close (10)
        end subroutine read_config_file

        ! 3. Set derived parameters
        subroutine set_derived_parameters(derived)
            implicit none
            type(sim_params_st), intent(inout) :: derived
            integer(kind=4) :: i
            character(:), allocatable :: aux_ch1, aux_ch2

            ! ----------------------- GLOBALS ---------------------------------
            if (derived%use_command_body) then
                if (abs(cl_body_in(1)) > cero)  then  ! Check the mass
                    derived%cl_body_is_moon = .True.
                    derived%Nmoons = 1
                else 
                    derived%cl_body_is_moon = .False.
                    derived%Nparticles = 1
                end if
            end if

            !! Times
            if ((derived%case_output_type < 0) .or. (derived%case_output_type > 2)) then
                write(*,*) "ERROR: Output times distribution method not recognized:", derived%case_output_type
                stop 1
            end if


            !! Expand checkpoints
            if (derived%extra_checkpoints < 0) then
                write(*,*) "ERROR: Extra checkpoint number can not be lower than 0."
                stop 1
            end if


            !! Boulders
            if (derived%Nboulders < 0) then
                write(*,*) "ERROR: Nboulders can not be negative."
                stop 1
            end if
            derived%use_boulders = derived%Nboulders > 0


            ! Primary Shape
            if (derived%use_triaxial) then
                ! Check order and values
                if (derived%triax_c_primary < cero) then
                    write(*,*) "ERROR: Semi-axis 'c' can not be negative"
                    stop 1
                else if (derived%triax_b_primary < cero) then
                    write(*,*) "ERROR: Semi-axis 'b' can not be negative"
                    stop 1
                else if  (derived%triax_a_primary .le. cero) then
                    write(*,*) "ERROR: Semi-axis 'a' must be positive."
                    stop 1
                else if (derived%triax_b_primary > derived%triax_a_primary) then
                    write(*,*) "ERROR: Semi-axis 'b' can not be greater than semi-axis 'a'."
                    stop 1
                else if (derived%triax_c_primary > derived%triax_b_primary) then
                    write(*,*) "ERROR: Semi-axis 'c' can not be greater than semi-axis 'b'."
                    stop 1
                end if
                ! Set values of not defined
                if (derived%triax_b_primary == cero) derived%triax_b_primary = derived%triax_a_primary
                if (derived%triax_c_primary == cero) derived%triax_c_primary = derived%triax_b_primary
                ! Check NO BOULDERS
                if (derived%use_boulders) then
                    write(*,*) "ERROR: Can not set tri-axial model with boulders."
                    stop 1
                end if
                if (derived%triax_a_primary - derived%triax_c_primary < tini) then
                    write(*,*) "WARNING: Tri-axial model with equal semi-axis."
                    write(*,*) "         Switching to single central sphere."
                    derived%radius_primary = (derived%triax_a_primary * &
                                            & derived%triax_b_primary * &
                                            & derived%triax_c_primary)**uno3
                    derived%use_triaxial = .False.
                end if
            else 
                ! Set semi to a sphere
                derived%triax_a_primary = derived%radius_primary
                derived%triax_b_primary = derived%radius_primary
                derived%triax_c_primary = derived%radius_primary
            end if


            !! NMoons
            if (derived%Nmoons < 0) then
                write(*,*) "ERROR: Nmoons can not be negative."
                stop 1
            end if
            derived%use_moons = derived%Nmoons > 0


            !! NParticles
            if (derived%Nparticles < 0) then
                write(*,*) "ERROR: Nparticles can not be negative."
                stop 1
            end if
            derived%use_particles = derived%Nparticles > 0
                       
            
            !! Forces
            !!! stokes_a_damping_time y stokes_e_damping_time
            if (abs(derived%stokes_a_damping_time) < tini) derived%stokes_a_damping_time = cero
            if (abs(derived%stokes_e_damping_time) < tini) derived%stokes_e_damping_time = cero
            if (derived%stokes_active_time .le. cero) derived%stokes_active_time = infinity
            if ((abs(derived%stokes_a_damping_time) < tini) .and. &
              & (abs(derived%stokes_e_damping_time) < tini)) derived%use_stokes = .False.
            if (.not. derived%use_stokes) then
                derived%stokes_a_damping_time = cero
                derived%stokes_e_damping_time = cero
                derived%stokes_active_time = cero
            end if
            
            !!! Naive-stokes
            if (abs(derived%drag_coefficient) < tini) then
                derived%use_drag = .False.
                derived%drag_coefficient = cero
            end if
            if (derived%drag_active_time .le. cero) derived%drag_active_time = infinity
            if (.not. derived%use_drag) then
                derived%drag_coefficient = cero
                derived%drag_active_time = cero
            end if
            
            !!! tau_m y tau_o
            if (derived%omega_damp_active_time .le. cero) derived%omega_damp_active_time = infinity
            derived%use_lin_omega_damp = abs(derived%omega_lin_damping_time) > tini
            derived%use_exp_omega_damp = abs(derived%omega_exp_damping_time) > tini
            derived%use_poly_omega_damp = (abs(derived%omega_exp_damp_poly_A) > tini .or. &
                                           abs(derived%omega_exp_damp_poly_B) > tini)
            if ((.not. derived%use_lin_omega_damp) .and. &
              & (.not. derived%use_exp_omega_damp) .and. &
              & (.not. derived%use_poly_omega_damp)) then
                derived%use_omega_damping = .False.
            end if
            if (.not. derived%use_omega_damping) then
                derived%omega_lin_damping_time = cero
                derived%omega_exp_damping_time = cero
                derived%omega_exp_damp_poly_A = cero
                derived%omega_exp_damp_poly_B = cero
                derived%omega_damp_active_time = cero
            else if ((derived%use_lin_omega_damp .and. (derived%use_exp_omega_damp .or. derived%use_poly_omega_damp)) .or. &
                     (derived%use_exp_omega_damp .and. derived%use_poly_omega_damp)) then
                write(*,*) "ERROR: Can not use multiple omega dampings."
                stop 1
            end if

            !!! manual J2
            if (derived%use_manual_J2) then
                if (abs(derived%manual_J2) < tini) then
                    derived%manual_J2 = cero
                    derived%use_manual_J2 = .False.
                end if
            else
                derived%manual_J2 = cero
            end if
            !!! Check no J2 and triaxial
            if (derived%use_manual_J2 .and. derived%use_triaxial) then
                write(*,*) "ERROR: Can not use both manual J2 and triaxial object at the same time."
                stop 1
            end if
            

            ! ! Mass Damping [UNAVAILABLE YET]
            ! if (abs(derived%mass_exp_damping_time) < tini) derived%mass_exp_damping_time = infinity


            ! Filter (fast check)
            if (derived%use_filter) then
                if (derived%filter_dt == 0) then
                    write(*,*) "ERROR: Filter dt must be different than 0."
                    stop 1
                end if
                if (derived%filter_nsamples .le. 0) then
                    write(*,*) "ERROR: Filter n_samples must be greater than 0."
                    stop 1
                end if
                if (derived%filter_nwindows .le. 0) then
                    write(*,*) "ERROR: Filter n_windows must be greater than 0."
                    stop 1
                end if
                if (to_lower(trim(derived%filter_prefix)) == "") derived%filter_prefix = "filt"
                if ((derived%filter_model < 0) .or. (derived%filter_model > 4)) then
                    write(*,*) "ERROR: Filter model must be between 0 and 4 (included)."
                    stop 1
                end if
            end if
            

            ! ! [BINS]
            ! !! Self-Gravity or Viscosity
            ! if (derived%use_self_gravity .or. derived%use_viscosity) then
            !     derived%use_bins = .True.
            !     !! Check
            !     if (derived%binning_method < 1 .or. derived%binning_method > 3) then
            !         write(*,*) "ERROR: binning_method must be 1 (equal dr), 2 (equal dA) or 3 (equal Npart)."
            !         stop 1
            !     end if
            !     if (derived%Nbins < 2) then
            !         write(*,*) "ERROR: Can not create less than 2 bins."
            !         stop 1
            !     end if
            !     !! Set default update values
            !     derived%update_rmin_bins = derived%rmin_bins < - tini
            !     derived%update_rmax_bins = derived%rmax_bins < - tini
            !     derived%update_bins = derived%update_rmin_bins .or. derived%update_rmax_bins .or. derived%binning_method == 3
            ! end if
            

            ! Moons gravity
            if (derived%Nmoons == 0) derived%use_moon_gravity = .False.  ! Deactivate it


            ! Merges
            if (derived%Nparticles == 0) derived%use_merge_part_mass = .False.  ! Deactivate it
            if (derived%Nmoons == 0) derived%use_merge_massive = .False.  ! Deactivate it
            derived%use_any_merge = derived%use_merge_part_mass .or. derived%use_merge_massive


            ! Stops
            if (derived%Nparticles == 0) derived%use_stop_no_part_left = .False.  ! Deactivate it
            if (derived%Nmoons == 0) derived%use_stop_no_moon_left = .False.  ! Deactivate it
            derived%use_any_stop = derived%use_stop_no_part_left .or. derived%use_stop_no_moon_left


            ! Eta and f collitions values
            if ((derived%eta_col < cero) .or. (derived%eta_col > uno)) then
                write(*,*) "ERROR: Collisional eta must be between 0 and 1."
                stop 1
            end if
            if ((derived%f_col < cero) .or. (derived%f_col > uno)) then
                write(*,*) "ERROR: Collisional f must be between 0 and 1."
                stop 1
            end if


            ! Error
            if (derived%error_digits < 1) then
                write(*,*) "ERROR: Number of presition digits must be greater than 0."
                stop 1
            end if
            derived%error_tolerance = 10.d0**(-derived%error_digits)
            

            ! Primary Radius
            if (derived%radius_primary <= 0) then
                write(*,*) "ERROR: Primary radius must be positive."
                stop 1
            end if


            ! Parallel  (Many here are GLOBAL)
            !$ compiled_with_openmp = .True. !! Parece comentado, pero es asi
            if (derived%use_parallel) then
                if ((derived%requested_threads == 0) .or. (derived%requested_threads == 1)) then
                    derived%requested_threads = 1
                    my_threads = 1
                    derived%use_parallel = .False.
                else if (.not. compiled_with_openmp) then
                    write(*,*) "ERROR: Can not use paralelization withut OpenMP."
                    stop 1
                end if
                !$ available_threads = OMP_GET_MAX_THREADS()
                !$ if (derived%requested_threads == -1) derived%requested_threads = available_threads
                !$ my_threads = min(available_threads, max(derived%requested_threads,1))
                !$ call OMP_SET_NUM_THREADS(my_threads)
            else
                derived%requested_threads = 1
                my_threads = 1
            end if

            
            ! Output
            if (derived%use_datascreen .and. derived%use_percentage) then
                write(*,*) "ERROR: Can not print both percentage and data on screen."
                stop 1
            end if
            if (derived%use_datascreen .and. derived%use_diagnostics) then
                write(*,*) "ERROR: Can not print both data and diagnostics data."
                stop 1
            end if
            if (derived%use_percentage .and. derived%use_diagnostics) then
                write(*,*) "ERROR: Can not print both percentage and diagnostics data."
                stop 1
            end if

            !! Geometric
            if (derived%use_geometricfile .and. .not. derived%use_triaxial) then
                write(*,*) "WARNING: No geometric file created as not triaxial body is used."
                derived%use_geometricfile = .False.
            end if
            
            !! Chaos
            derived%use_chaos = derived%use_chaosfile .or. (derived%extra_checkpoints > 0)
            if (derived%use_chaos .and. derived%use_geometricfile) then
                derived%use_geomchaosfile = .True.  ! Need to output geometric chaos values
                aux_ch1 = trim(adjustl(derived%geometricfile))
                i = index(aux_ch1, ".", back=.True.)
                if (i == 0) i = len_trim(aux_ch1) + 1 ! No "." in geometricfile
                aux_ch2 = aux_ch1(1:i-1)
                if (len_trim(derived%chaosfile) > 0) then 
                    derived%geomchaosfile = trim(aux_ch2) // "_" // trim(adjustl(derived%chaosfile))
                else
                    derived%geomchaosfile = trim(aux_ch2) // "_chaos.out"
                end if
                deallocate(aux_ch1)
                deallocate(aux_ch2)
            else
                derived%use_geomchaosfile = .False.  ! Do not output geometric chaos values
            end if

        
            !! Names
            if ((trim(derived%datafile) /= "") .and. (trim(derived%datafile) /= "no")) then
                derived%use_datafile = .True.
            else 
                derived%use_datafile = .False.
            end if
            if ((trim(derived%chaosfile) /= "") .and. (trim(derived%chaosfile) /= "no")) then
                derived%use_chaosfile = .True.
                derived%use_chaos = .True.
            else 
                derived%use_chaosfile = .False.
            end if
            if ((trim(derived%multfile) /= "") .and. (trim(derived%multfile) /= "no")) then
                derived%use_multiple_outputs = .True.
            else 
                derived%use_multiple_outputs = .False.
            end if
            if ((trim(derived%mapfile) /= "") .and. (trim(derived%mapfile) /= "no")) then
                derived%use_potential_map = .True.
            else 
                derived%use_potential_map = .False.
            end if
            if ((trim(derived%tomfile) /= "") .and. (trim(derived%tomfile) /= "no")) then
                derived%use_tomfile = .True.
            else 
                derived%use_tomfile = .False.
            end if
            if ((trim(derived%moonsfile) /= "") .and. (trim(derived%moonsfile) /= "no")) then
                derived%use_moonsfile = .True.
            else 
                derived%use_moonsfile = .False.
            end if
            if ((trim(derived%particlesfile) /= "") .and. (trim(derived%particlesfile) /= "no")) then
                derived%use_particlesfile = .True.
            else 
                derived%use_particlesfile = .False.
            end if


            ! Check no TOM and Filter
            if (derived%use_tomfile .and. derived%use_filter) then
                write(*,*) "ERROR: Can not use both filtering and TOM at the same time."
                stop 1
            end if


            ! Check no Extra checkpoints and Filter
            if (derived%use_tomfile .and. (derived%extra_checkpoints > 0)) then
                write(*,*) "ERROR: Can not use both filtering and extra checkpoints at the same time (yet)."
                write(*,*) "       Set all checkpoints to outputs."
                stop 1
            end if

        end subroutine set_derived_parameters
        

        ! Define IO pointers
        subroutine define_writing_pointers(simu)
            use bodies, only: write_chaos, write_elem, write_coor, &
                            & write_ast_elem, write_moon_i_elem, write_particle_i_elem, &
                            & write_ast_coor, write_moon_i_coor, write_particle_i_coor, &
                            & write_geom, write_geomchaos
            implicit none
            type (sim_params_st), intent(in) :: simu

            ! Screen
            if (simu%use_datascreen) then
                if (simu%use_elements_output) then
                    write_to_screen => write_elem
                else
                    write_to_screen => write_coor
                end if
            else if (simu%use_diagnostics) then
                write_to_screen => do_not_write
            else
                write_to_screen => do_not_write
            end if

            ! Datafile
            if (simu%use_datafile) then
                if (simu%use_elements_output) then
                    write_to_general => write_elem
                else
                    write_to_general => write_coor
                end if
                if (use_flush_output) then
                    flush_output => flush_to_file
                else 
                    flush_output => do_nothing_i
                end if  
            else 
                write_to_general => do_not_write
                flush_output => do_nothing_i
            end if

            ! Multiple Files
            if (simu%use_multiple_outputs) then
                if (simu%use_elements_output) then
                    write_a_to_individual => write_ast_elem
                    write_m_to_individual => write_moon_i_elem
                    write_p_to_individual => write_particle_i_elem
                else
                    write_a_to_individual => write_ast_coor
                    write_m_to_individual => write_moon_i_coor
                    write_p_to_individual => write_particle_i_coor
                end if
            else 
                write_a_to_individual => do_not_write
                write_m_to_individual => do_not_write_i
                write_p_to_individual => do_not_write_i
            end if

            ! Chaos file
            if (simu%use_chaosfile) then
                write_to_chaos => write_chaos
                ! Flush chaos
                if (use_flush_chaos) then
                    flush_chaos => flush_and_rewind
                else 
                    flush_chaos => rewind_a_file
                end if
            else 
                write_to_chaos => do_not_write_ch
                flush_chaos => do_nothing_i
            end if

            ! Geometric file
            if (simu%use_geometricfile) then
                write_to_geometric => write_geom
                if (use_flush_output) then
                    flush_geometric => flush_to_file
                else 
                    flush_geometric => do_nothing_i
                end if
                if (simu%use_geomchaosfile) then
                    write_to_geomchaos => write_chaos
                    if (use_flush_chaos) then
                        flush_geomchaos => flush_and_rewind
                    else 
                        flush_geomchaos => rewind_a_file
                    end if
                else
                    write_to_geomchaos => do_not_write_ch
                    flush_geomchaos => do_nothing_i
                end if
            else
                write_to_geometric => do_not_write
                write_to_geomchaos => do_not_write_ch
                flush_geometric => do_nothing_i
                flush_geomchaos => do_nothing_i
            end if

        end subroutine define_writing_pointers
            

        ! Check if continue function. This will be passed to the integ caller    
        function check_func (y) result(keep_going)
            implicit none
            real(kind=8), dimension(:), intent(in) :: y
            logical :: keep_going
            ! ! Here, we use the global hexit_arr  (Old Code)
            ! keep_going = all(hexit_arr == 0)
            ! if (.not. keep_going) hexit_arr(1) = 1  ! Set the PROXY
            keep_going = .not. hard_exit
        end function check_func


        ! Deallocate initial arrays
        subroutine free_initial_arays()
            implicit none
            if (allocated(moons_arr)) deallocate(moons_arr)
            if (allocated(particles_arr)) deallocate(particles_arr)
            if (allocated(boulders_in)) deallocate(boulders_in)
            if (allocated(moons_in)) deallocate(moons_in)
            if (allocated(particles_in)) deallocate(particles_in)
        end subroutine free_initial_arays


        ! Deallocate all params arrays
        subroutine free_parameters_arays()
            implicit none
            call free_initial_arays() !! Just in case...
            if (allocated(boulders_coords)) deallocate(boulders_coords)
            if (allocated(boulders_data)) deallocate(boulders_data)
            if (allocated(m_arr)) deallocate(m_arr)
            if (allocated(R_arr)) deallocate(R_arr)
            if (allocated(y_arr)) deallocate(y_arr)
            if (allocated(y_arr_new)) deallocate(y_arr_new)
            if (allocated(y_der)) deallocate(y_der)
            if (allocated(output_times)) deallocate(output_times)
            if (allocated(checkpoint_times)) deallocate(checkpoint_times)
            if (allocated(checkpoint_is_tom)) deallocate(checkpoint_is_tom)
            if (allocated(checkpoint_is_output)) deallocate(checkpoint_is_output)
            if (allocated(y_pre_filter)) deallocate(y_pre_filter)
            if (allocated(elem_filtered)) deallocate(elem_filtered)
            if (allocated(tmp_y_arr)) deallocate(tmp_y_arr)
        end subroutine free_parameters_arays


        ! Nullify pointers
        subroutine nullify_pointers()
            implicit none
            nullify(write_to_general)
            nullify(write_to_screen)
            nullify(write_to_chaos)
            nullify(write_to_geometric)
            nullify(write_to_geomchaos)
            nullify(write_a_to_individual)
            nullify(write_m_to_individual)
            nullify(write_p_to_individual)
            nullify(flush_output)
            nullify(flush_chaos)
            nullify(flush_geometric)
            nullify(flush_geomchaos)
        end subroutine nullify_pointers


        ! Update sim Nactive values
        pure subroutine update_sim_Nactive(simu, Npart_active, Nmoon_active)
            implicit none
            type(sim_params_st), intent(inout) :: simu
            integer(kind=4), intent(in) :: Npart_active, Nmoon_active
            simu%Npart_active = Npart_active  ! Update Npart_active
            simu%Nmoon_active = Nmoon_active  ! Update Nmoon_active
            simu%Nactive = 1 + Npart_active + Nmoon_active ! Update Nactive, including asteroid
        end subroutine update_sim_Nactive


        ! Check if keep integrating after collisions
        subroutine check_after_col(sistema, simu, keep)
            implicit none
            type(system_st), intent(in) :: sistema
            type(sim_params_st), intent(in) :: simu
            logical, intent(inout) :: keep

            if ((.not. simu%use_merge_part_mass) .and. sistema%Nparticles_active < sistema%Nparticles_active) then
                if (simu%use_screen) then
                    write(*,*) ACHAR(10)
                    write(*,*) "A particle-massive merge body occured."
                end if
                keep = .False.
            end if
            if ((.not. simu%use_merge_massive) .and. sistema%Nmoons_active < sistema%Nmoons_active) then
                if (simu%use_screen) then
                    write(*,*) ACHAR(10)
                    write(*,*) "A massive-massive merge body occured."
                end if
                keep = .False.
            end if
        end subroutine check_after_col


        ! Check if keep integrating after escapes (and collisions)
        subroutine check_after_esc(sistema, simu, keep)
            implicit none
            type(system_st), intent(in) :: sistema
            type(sim_params_st), intent(in) :: simu
            logical, intent(inout) :: keep
            
            if (simu%use_stop_no_part_left .and. sistema%Nparticles_active == 0 .and. simu%Nparticles > 0) then
                if (simu%use_screen) then
                    write(*,*) ACHAR(10)
                    write(*,*) "No more particles left."
                end if
                keep = .False.
            end if
            if (simu%use_stop_no_moon_left .and. sistema%Nmoons_active == 0 .and. simu%Nmoons > 0) then
                if (simu%use_screen) then
                    write(*,*) ACHAR(10)
                    write(*,*) "No more moons left."
                end if
                keep = .False.
            end if
        end subroutine check_after_esc


        ! IO subroutines (and dummies)

        subroutine flush_to_file(unit_file)
            implicit none
            integer(kind=4), intent(in) :: unit_file

            flush(unit_file)
        end subroutine flush_to_file

        subroutine flush_and_rewind(unit_file)
            implicit none
            integer(kind=4), intent(in) :: unit_file

            flush(unit_file)
            rewind(unit_file)
        end subroutine flush_and_rewind

        subroutine rewind_a_file(unit_file)
            implicit none
            integer(kind=4), intent(in) :: unit_file

            rewind(unit_file)
        end subroutine rewind_a_file

        subroutine do_nothing_i(i)
            implicit none
            integer(kind=4), intent(in) :: i
        end subroutine do_nothing_i

        subroutine do_not_write_i(self, i, unit_file)
            implicit none
            type(system_st), intent(in) :: self
            integer(kind=4), intent(in) :: i, unit_file
        end subroutine do_not_write_i

        subroutine do_not_write(self, unit_file)
            implicit none
            type(system_st), intent(in) :: self
            integer(kind=4), intent(in) :: unit_file
        end subroutine do_not_write

        subroutine do_not_write_ch(self, other, unit_file)
            implicit none
            type(system_st), intent(in) :: self
            type(system_st), intent(in) :: other
            integer(kind=4), intent(in) :: unit_file
        end subroutine do_not_write_ch

        ! FILTERING
        subroutine apply_filter(y_nvalues, syst_pre_filter, syst_filtered, el_filtered)
            use bodies, only: get_Nactive, copy_objects, update_system_from_array, update_elements, moon_st, particle_st
            implicit none
            integer(kind=4), intent(in) :: y_nvalues
            type(system_st), intent(in) :: syst_pre_filter
            type(system_st), intent(inout) :: syst_filtered
            real(kind=8), dimension(y_nvalues), intent(inout) :: el_filtered
            real(kind=8) :: cos_th, sin_th
            real(kind=8) :: weigth, e_times_fil
            real(kind=8), dimension(:,:), allocatable :: cos_an, sin_an
            integer(kind=4) :: nbodies, i, j, aux_i, idx

            call get_Nactive(syst_pre_filter, nbodies)
            allocate(cos_an(nbodies, 2))
            allocate(sin_an(nbodies, 2))

            cos_th = cero
            sin_th = cero
            cos_an = cero
            sin_an = cero
            
            el_filtered = cero
            call copy_objects(syst_pre_filter, syst_filtered)
            do i = 1, filter%size
                weigth = filter%kernel(i)
                call update_system_from_array(syst_filtered, filter%tmp_times(i), filter%tmp_values(:y_nvalues, i))
                call update_elements(syst_filtered, sim%use_baryc_output)

                ! Theta
                cos_th = cos_th + cos(syst_filtered%asteroid%theta) * weigth
                sin_th = sin_th + sin(syst_filtered%asteroid%theta) * weigth
                ! Omega
                el_filtered(2) = el_filtered(2) + syst_filtered%asteroid%omega * weigth
                ! a
                el_filtered(3) = el_filtered(3) + syst_filtered%asteroid%elements(1) * weigth
                ! e
                e_times_fil = syst_filtered%asteroid%elements(2) * weigth
                el_filtered(4) = el_filtered(4) + e_times_fil
                ! M
                cos_an(1,1) = cos_an(1,1) + cos(syst_filtered%asteroid%elements(3)) * weigth
                sin_an(1,1) = sin_an(1,1) + sin(syst_filtered%asteroid%elements(3)) * weigth
                ! w -> K & H
                cos_an(1,2) = cos_an(1,2) + cos(syst_filtered%asteroid%elements(4)) * e_times_fil
                sin_an(1,2) = sin_an(1,2) + sin(syst_filtered%asteroid%elements(4)) * e_times_fil

                do j = 2, syst_filtered%Nmoons_active + 1
                    aux_i = j - 1
                    idx =  4 * j - 1  ! From derivates

                    ! a
                    el_filtered(idx) = el_filtered(idx) + syst_filtered%moons(aux_i)%elements(1) * weigth
                    ! e
                    e_times_fil = syst_filtered%moons(aux_i)%elements(2) * weigth
                    el_filtered(idx+1) = el_filtered(idx+1) + e_times_fil
                    ! M
                    cos_an(j,1) = cos_an(j,1) + cos(syst_filtered%moons(aux_i)%elements(3)) * weigth
                    sin_an(j,1) = sin_an(j,1) + sin(syst_filtered%moons(aux_i)%elements(3)) * weigth
                    ! w -> K & H
                    cos_an(j,2) = cos_an(j,2) + cos(syst_filtered%moons(aux_i)%elements(4)) * e_times_fil
                    sin_an(j,2) = sin_an(j,2) + sin(syst_filtered%moons(aux_i)%elements(4)) * e_times_fil

                end do

                do j = syst_filtered%Nmoons_active + 2, nbodies
                    aux_i = j - syst_filtered%Nmoons_active - 1
                    idx =  4 * j - 1  ! From derivates

                    ! a
                    el_filtered(idx) = el_filtered(idx) + syst_filtered%particles(aux_i)%elements(1) * weigth
                    ! e
                    e_times_fil = syst_filtered%particles(aux_i)%elements(2) * weigth
                    el_filtered(idx+1) = el_filtered(idx+1) + e_times_fil
                    ! M
                    cos_an(j,1) = cos_an(j,1) + cos(syst_filtered%particles(aux_i)%elements(3)) * weigth
                    sin_an(j,1) = sin_an(j,1) + sin(syst_filtered%particles(aux_i)%elements(3)) * weigth
                    ! w -> K & H
                    cos_an(j,2) = cos_an(j,2) + cos(syst_filtered%particles(aux_i)%elements(4)) * e_times_fil
                    sin_an(j,2) = sin_an(j,2) + sin(syst_filtered%particles(aux_i)%elements(4)) * e_times_fil
                end do

            end do

            ! theta
            el_filtered(1) = modulo(atan2(sin_th, cos_th), twopi)
            do j = 1, nbodies
                idx =  4 * j - 1  ! From derivates
                ! e ?
                if (sim%filter_use_KH) then
                    el_filtered(idx+1) = sqrt(cos_an(j,2) * cos_an(j,2) + sin_an(j,2) * sin_an(j,2))
                end if
                ! M
                el_filtered(idx+2) = modulo(atan2(sin_an(j,1), cos_an(j,1)), twopi)
                ! w
                el_filtered(idx+3) = modulo(atan2(sin_an(j,2), cos_an(j,2)), twopi)
            end do

            deallocate(cos_an)
            deallocate(sin_an)
            
        end subroutine apply_filter


        ! SPECIFIC subroutines to avoid writing the same lines of code

        subroutine generate_output(syst, filtered)
            use bodies, only: write_diagnostics
            implicit none
            type(system_st), intent(in) :: syst
            logical, intent(in), optional :: filtered
            integer(kind=4) :: i, pind, pout

            if (present(filtered)) then
                pind = 3000
                pout = 10
            else
                pind = 0
                pout = 0
            end if

            call write_a_to_individual(syst, 100 + pind)
            !$OMP PARALLEL DEFAULT(SHARED) &
            !$OMP PRIVATE(i)
            !$OMP DO
            do i = 1, syst%Nmoons_active
                call write_m_to_individual(syst, i, 100+syst%moons(i)%id + pind)
            end do 
            !$OMP END DO NOWAIT
            !$OMP DO
            do i = 1, syst%Nparticles_active
                call write_p_to_individual(syst, i, 100+syst%particles(i)%id + pind)
            end do 
            !$OMP END DO NOWAIT
            !$OMP SECTIONS
            !$OMP SECTION
            if (sim%use_diagnostics) then
                call write_diagnostics(initial_system, syst, 6)
            else if (pout == 0) then
                call write_to_screen(syst, 6)
                flush(6)
            end if
            !$OMP SECTION
            call write_to_general(syst, 21 + pout)
            call flush_output(21 + pout)
            !$OMP SECTION
            call write_to_chaos(syst, initial_system, 22 + pout)
            call flush_chaos(22 + pout) ! Update chaos
            !$OMP SECTION
            if (pout == 0) then
                call write_to_geometric(syst, 23)
                call flush_geometric(23)
                call write_to_geomchaos(syst, initial_system, 24)
                call flush_geomchaos(24)
            end if
            !$OMP END SECTIONS
            !$OMP END PARALLEL        
            
        end subroutine generate_output

        subroutine check_esc_and_col(syst, unit_out)
            use bodies, only: resolve_collisions, resolve_escapes
            implicit none
            type(system_st), intent(inout) :: syst
            integer(kind=4), intent(in) :: unit_out

            ! Apply colissions and check
            call resolve_collisions(syst, sim%min_distance, unit_out)
            call check_after_col(syst, sim, keep_integrating)

            ! Apply escapes and check
            call resolve_escapes(syst, sim%max_distance, unit_out)
            call check_after_esc(syst, sim, keep_integrating)
            
        end subroutine check_esc_and_col

end module parameters
