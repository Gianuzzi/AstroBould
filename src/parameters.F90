module parameters
    use constants
    use auxiliar
    use bodies
    use omp_lib

    implicit none

    ! --------------------------- GLOBAL CONFIGURATION ---------------------------
    integer(kind=4) :: simulation_number = 1
    !! I/O
    integer(kind=4) :: arguments_number = 0  ! Amount of argument in command line
    logical :: configfile_exists = .False.  ! Wether if config file exists
    logical :: use_configfile = .True.  ! Wether to use config file (if exists)
    !!! Command line particle
    logical :: use_command_particle = .False.  ! Wether a particle was given through command line
    !! Parallel
    integer(kind=4) :: available_threads = 1  ! Available threads to use
    integer(kind=4) :: my_threads = 1  ! Amount of threads actually used
    logical :: compiled_with_openmp = .False.  ! Flag of compile with OpenMP

    ! Bodies (from bodies)
    type(particle_st), allocatable :: particles_arr(:)
    type(moon_st), allocatable :: moons_arr(:)
    type(asteroid_st) :: asteroid
    type(system_st) :: initial_system
    type(system_st) :: system
    
    ! This contains only the input parameters
    type :: input_params_st
        ! Times for the integration - 
        real(kind=8) :: initial_time = cero
        real(kind=8) :: final_time = cero
        real(kind=8) :: output_timestep = cero
        integer(kind=4) :: output_number = 0
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
        logical :: use_stokes = .False.
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
        ! [BINS] forces/effects - [NOT AVAILABLE YET. STILL UNDER DEVELOPMENT]
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
        real(kind=8) :: map_min_x = -500.d0
        real(kind=8) :: map_max_x = 500.d0
        real(kind=8) :: map_min_y = -500.d0
        real(kind=8) :: map_max_y = 500.d0
        ! < Number of bodies >
        integer(kind=4) :: Nboulders = 0
        integer(kind=4) :: Nmoons = 0
        integer(kind=4) :: Nparticles = 0

        !!!!! Extra DERIVED parameters
        ! Bodies
        integer(kind=4) :: Ntotal = 0  ! Amount of all bodies (asteroid is just 1)
        integer(kind=4) :: Nactive = 0 ! Active bodies
        logical :: use_boulders = .False.
        logical :: use_particles = .False.
        logical :: use_moons = .False.
        ! Forces
        logical :: use_torque = .False.  ! May disapear...
        !! Stokes
        real(kind=8) :: stokes_C = cero  ! Parameter C in Stokes
        real(kind=8) :: stokes_alpha = uno  ! Parameter alpha in Stokes
        !! Mass and omega damping
        real(kind=8) :: omega_exp_poly_AB = cero
        real(kind=8) :: omega_linear_damping_slope = cero
        !! Geo-potential
        logical :: use_J2 = .False.
        real(kind=8) :: J2_effective = cero
        !! [BINS] ( Not available yet )
        logical :: use_bins = .False.
        !!! Viscosity [Not available yet]
        logical :: use_viscosity = .False.
        real(kind=8) :: viscosity = - uno
        !!! Self-Gravity
        logical :: update_rmin_bins = .False.
        logical :: update_rmax_bins = .False.
        logical :: update_bins = .False.
        ! Chaos
        logical :: use_chaos = .False.  ! Need to output chaos values
        ! logical :: use_update_chaos = .True.  ! GLOBAL
        ! logical :: use_flush_output = .True.  ! GLOBAL
        ! Error
        real(kind=8) :: error_tolerance = cero
        ! I/O
        logical :: use_elements = .True.
    end type input_params_st

    !! < Arrays with input data from bodies >
    real(kind=8), dimension(:,:), allocatable :: boulders_in !! mass, radius, theta  | (Nb, 3)
    real(kind=8), dimension(:,:), allocatable :: moons_in !! mass, a, e, M, w, MMR, radius  | (Nm, 7)
    real(kind=8), dimension(:,:), allocatable :: particles_in !! a, e, M, w, MMR  | (Np, 5)

    !!! < Initial boulders angles array > (THis will serve to calculate their positions)
    real(kind=8), dimension(:), allocatable :: boulders_theta_ini !! (Nb)

    ! Time Integration dynamical variables
    real(kind=8) :: time  ! Actual time of the integration
    real(kind=8) :: timestep  ! This timestep 
    real(kind=8) :: adaptive_timestep  ! This adaptive timestep
    integer(kind=4) :: checkpoint_number  ! Number of actual timestep
    real(kind=8) :: min_timestep = cero  ! Minimum timestep

    ! I/O
    logical :: use_flush_output  ! Flush output at every checkpoint
    logical :: use_update_chaos  ! Update values (flush) at every checkpoint
    
    ! TOM
    integer(kind=4) :: tom_index_number  ! Index to count which TOM row is active
    integer(kind=4) :: tom_total_number  ! Total amount of lines in TOM
    real(kind=8), dimension(:), allocatable :: tom_times  ! Times in TOM
    real(kind=8), dimension(:), allocatable :: tom_deltaomega  ! Delta Omega in TOM
    real(kind=8), dimension(:), allocatable :: tom_deltamass  ! Delta Mass in TOM
    real(kind=8) :: tom_mass_growth_param
    logical, dimension(:), allocatable :: checkpoint_is_tom  ! Array with whether each checkpoint is TOM
    logical, dimension(:), allocatable :: checkpoint_is_output ! Array with whether each checkpoint is output
        
    
    ! !! Auxiliar
    ! logical :: aux_logical
    ! real(kind=8) :: aux_real
    ! real(kind=8) :: dummy_real
    ! real(kind=8) :: dummy_real2
    ! integer(kind=4) :: aux_integer
    ! character(1) :: aux_character1
    ! character(2) :: aux_character2
    ! character(20) :: aux_character20
    ! character(30) :: aux_character30
    ! real(kind=8), dimension(2) :: aux_real_arr2
    ! real(kind=8), dimension(6) :: aux_real_arr6
    ! logical :: hard_center
    ! logical :: use_elements


    ! ----  <<<<<    PARAMETERS ARRAYS     >>>>>   -----
    integer, parameter :: equation_size = 4
    real(kind=8), dimension(:), allocatable :: m_arr
    real(kind=8), dimension(:), allocatable :: y_arr
    real(kind=8), dimension(:), allocatable :: y_arr_new
    real(kind=8), dimension(:), allocatable :: y_der
    
    !!! Hard Exit
    logical :: is_premature_exit
    integer(kind=4), dimension(:), allocatable, target :: particles_hexit !  Hard Exit integer
    integer(kind=4), pointer :: particles_hexitptr ! pointer to Hard Exit
    procedure (check_continue_template), pointer :: check_continue_ptr => null ()
    

    !!! GENERAL element PTRs
    procedure (int_i_template), pointer :: get_elements_i => null ()
    procedure (int_i_template), pointer :: get_chaos_i => null ()

    !!! MERGE
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

        ! 0. Init GLOBAL hidden (and command line) default parameters
        subroutine init_default_parameters()
            implicit none

            ! Simulation
            simulation_number = 1

            ! Times
            min_timestep = cero
            
            ! OmenMP
            !$ compiled_with_openmp = .True. !! Parece comentado, pero es asi

            ! Extras, fuera de Config.ini. Se editan solo acá
            !! Update chaos at every output timestep, instead of waiting to the end of the run
            use_update_chaos = .True.
            !! Flush output at every output timestep (even though the buffer is not full) 
            use_flush_output = .True.
        end subroutine init_default_parameters
    
        ! Allocate params arrays
        subroutine allocate_params_asteroid(Nboulders)
            implicit none
            integer(kind=4), intent(in) :: Nboulders

            if (allocated(boulders_in)) then
                write (*,*) "ERROR: Number of boulder already defined."
                STOP 1
            end if

            allocate(boulders_in(Nboulders,3))  ! FROM 1 to Nboulders
        end subroutine allocate_params_asteroid

        subroutine allocate_params_moons(Nmoons)
            implicit none
            integer(kind=4), intent(in) :: Nmoons

            if (allocated(moons_in)) then
                write (*,*) "ERROR: Number of moons already defined."
                STOP 1
            end if

            allocate(moons_in(Nmoons,7))  ! FROM 1 to Nmoons
        end subroutine allocate_params_moons

        subroutine allocate_params_particles(Nparticles)
            implicit none
            integer(kind=4), intent(in) :: Nparticles

            if (allocated(particles_in)) then
                write (*,*) "ERROR: Number of particles already defined."
                STOP 1
            end if

            allocate(particles_in(Nparticles,5))  ! FROM 1 to Nparticles
        end subroutine allocate_params_particles
        
        ! 1. Leer entradas por linea de comandos
        subroutine load_command_line_arguments(params, use_config)
            implicit none
            type(input_params_st), intent(inout) :: params
            logical, intent(inout)  :: use_config
            logical :: is_number
            integer(kind=4) :: i, j
            logical :: aux_logical = .False.
            integer(kind=4) :: aux_integer = 0
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
                    case ("-nsim")
                        call get_command_argument(i+1, aux_character20)
                        read (aux_character20,*) simulation_number  ! Global variable
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
                    case ("--screen")
                        params%use_screen = .True.
                    case ("--noscreen")
                        params%use_screen = .False.
                    case ("--perc")
                        params%use_percentage = .True.
                    case ( "--noperc")
                        params%use_percentage = .False.
                    case ("--datascr")
                        params%use_datascreen = .True.
                    case ( "--nodatascr")
                        params%use_datascreen = .False.
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
                    case ("--merge")
                        params%use_merge = .True.
                    case ("--nomerge")
                        params%use_merge = .False.
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
                                write (*,*) "ERROR: Could not read all particle orbital elements."
                                write (*,'(A,I0,A)') "        Missing ", 4-arguments_number, " elements."
                                write (*,*) "Exiting."
                                stop 1
                            else ! Leo los argumentos numéricos. Considero que es 1 sola partícula
                                ! Para generar prioridad, reallocatamos de ser necesario
                                use_command_particle = .True.  ! GLOBAL parameter
                                if (allocated(particles_in)) deallocate(particles_in)
                                call allocate_params_particles(1)
                                params%Nparticles = 1
                                ! Read
                                call get_command_argument(i, aux_character20)
                                read (aux_character20,*) particles_in(1,1)
                                call get_command_argument(i+1, aux_character20)
                                read (aux_character20,*) particles_in(1,2)
                                call get_command_argument(i+2, aux_character20)
                                read (aux_character20,*) particles_in(1,3)
                                call get_command_argument(i+3, aux_character20)
                                read (aux_character20,*) particles_in(1,4)
                                aux_integer = 3
                                aux_logical = .True.
                                particles_in(1,5) = cero
                                
                            end if
                        else ! Ya leí los numéricos. Falta leer eR
                            call get_command_argument(i, aux_character20)
                            read (aux_character20,*) particles_in(1,5)
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
            integer(kind=4)  :: colonPos, commentPos
            integer(kind=4) :: nlines
            integer(kind=4) :: io
            integer(kind=4) :: len_val
            integer(kind=4) :: j
            integer(kind=4) :: aux_integer
            integer(kind=4) :: Nboulders = 0, Nmoons = 0, Nparticles = 0

            inquire(file=trim(file_name), exist=existe)
            if (.not. existe) then
                if (params%use_screen) write (*,*) "File ", trim(file_name), " does not exists."  ! use_screen is the default...
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
                            read (value_str, *) params%initial_time
                        case("total integrati")
                            read (value_str, *) params%final_time
                        case("output time int")
                            read (value_str, *) params%output_timestep
                        case("number of outpu")
                            read (value_str, *) params%output_number
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
                        case("precision (digi")
                            read (value_str, *) params%error_digits
                        case("learning rate")
                            read (value_str, *) params%learning_rate
                        case("mass of primary")
                            read (value_str, *) params%mass_primary
                        case("radius of prima")
                            read (value_str, *) params%radius_primary
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
                        case("F damping chara")
                            read (value_str, *) params%stokes_charac_time
                        case("include naive-s")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_naive_stokes = .True.
                            else
                                params%use_naive_stokes = .False.
                            end if
                        case("drag force coef")
                            read (value_str, *) params%drag_coefficient
                        case("drag force char")
                            read (value_str, *) params%drag_charac_time
                        case("include asteroi")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_omega_damping = .True.
                            else
                                params%use_omega_damping = .False.
                            end if
                        case("rotation linear")
                            read (value_str, *) params%omega_linear_damping_time
                        case("rotation expone")
                            read (value_str, *) params%omega_exp_damping_time
                        case("rotation A poly")
                            read (value_str, *) params%omega_exp_poly_A
                        case("rotation B poly")
                            read (value_str, *) params%omega_exp_poly_B
                        case("omega damping c")
                            read (value_str, *) params%omega_charac_time
                        case("mass exponentia")
                            read (value_str, *) params%mass_exp_damping_time
                        case("geo-potential J")
                            read (value_str, *) params%J2_coefficient
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
                        case("merge collision")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                params%use_merge = .True.
                            else
                                params%use_merge = .False.
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
                        case("output data on")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                params%use_datascreen = .True.
                                params%use_percentage = .False.
                            else if (auxch1 == "%") then
                                params%use_datascreen = .False.
                                params%use_percentage = .True.
                            else
                                params%use_datascreen = .False.
                                params%use_percentage = .False.
                            end if
                        case("output variable")
                            if (auxch1 == "c") then
                                params%use_elements_output = .False.
                            else
                                params%use_elements_output = .True.
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
                                    write (*,*) "ERROR: Al leer luna:", j
                                    stop 1
                                end if
                                j = j + 1
                            end do
                        end if
                    end if
                    if (param_str(1:5) .eq. "a(km)") then  ! Leeemos las partículas
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
                                    write (*,*) "ERROR: Al leer partícula:", j
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

        ! 3. Setear los parámetros derivados
        subroutine set_derived_parameters(derived)
            implicit none
            type(input_params_st), intent(inout) :: derived

            !! Times
            if ((derived%case_output_type < 0) .or. (derived%case_output_type > 2)) then
                write (*,*) "ERROR: Checkpoint times distribution method not recognized:", derived%case_output_type
                stop 1
            end if


            !! Boulders
            if (derived%Nboulders < 0) then
                write (*,*) "ERROR: Nboulders can not be negative 0"
                stop 1
            end if
            derived%use_boulders = (derived%Nboulders > 0)

            ! -----------------------------------------------------------------
            ! Acá hay que configurar moons y particles en caso de moons sin masa
            ! -----------------------------------------------------------------

            !! Moons
            if (derived%Nmoons < 0) then
                write (*,*) "ERROR: Nmoons can not be negative 0"
                stop 1
            end if
            derived%use_moons = (derived%Nmoons > 0)

            !! Moons
            if (derived%Nparticles < 0) then
                write (*,*) "ERROR: Nparticles can not be negative 0"
                stop 1
            end if
            derived%use_particles = (derived%Nparticles > 0)
                       
            
            !! Forces
            !!! stokes_a_damping_time y stokes_e_damping_time
            if ((abs(derived%stokes_a_damping_time) < tini) .or. &
              & .not. derived%use_stokes) derived%stokes_a_damping_time = cero
            if ((abs(derived%stokes_e_damping_time) < tini) .or. &
              & .not. derived%use_stokes) derived%stokes_e_damping_time = cero
            if ((abs(derived%stokes_charac_time) < tini) .or. &
              & .not. derived%use_stokes) derived%stokes_charac_time = cero
            if (derived%stokes_charac_time .le. cero) derived%stokes_charac_time = infinity
            if ((abs(derived%stokes_a_damping_time) < tini) .and. &
              & (abs(derived%stokes_e_damping_time) < tini)) derived%use_stokes = .False.

            !!! Naive-stokes
            if ((abs(derived%drag_coefficient) < tini) .or. &
              & .not. derived%use_naive_stokes) derived%drag_coefficient = cero
            if ((abs(derived%drag_charac_time) < tini) .or. &
              & .not. derived%use_naive_stokes) derived%drag_charac_time = cero
            if (derived%drag_charac_time .le. cero) derived%drag_charac_time = infinity
            if (abs(derived%drag_coefficient) < tini) derived%use_naive_stokes = .False.

            !!! tau_m y tau_o
            if ((abs(derived%omega_charac_time) < tini) .or. &
              & .not. derived%use_omega_damping) derived%omega_charac_time = cero
            if (derived%omega_charac_time .le. cero) derived%omega_charac_time = infinity
            if ((abs(derived%omega_linear_damping_time) < tini) .or. &
              & .not. derived%use_omega_damping) derived%omega_linear_damping_time = infinity
            if ((abs(derived%omega_exp_damping_time) < tini) .or. &
              & .not. derived%use_omega_damping) derived%omega_exp_damping_time = infinity
            if ((abs(derived%omega_exp_poly_A) < tini) .or. &
              & .not. derived%use_omega_damping) derived%omega_exp_poly_A = cero
            if ((abs(derived%omega_exp_poly_B) < tini) .or. &
              & .not. derived%use_omega_damping) derived%omega_exp_poly_B = cero
            if (.not. ((abs(derived%omega_linear_damping_time) > tini) .and. &
              & (derived%omega_linear_damping_time < infinity)) .and. &
              & .not. ((abs(derived%omega_exp_damping_time) > tini) .and. &
              & (derived%omega_exp_damping_time < infinity)) .and. &
              & .not. ((abs(derived%omega_exp_poly_A) > tini) .and. &
              & (abs(derived%omega_exp_poly_B) > tini))) &
              & derived%use_omega_damping = .False.
            if (abs(derived%mass_exp_damping_time) < tini) derived%mass_exp_damping_time = infinity

            !!! Geo-Potential
            if (abs(derived%J2_coefficient) > tini) derived%use_J2 = .True.
            
            ! ! [BINS]
            ! !! Self-Gravity or Viscosity
            ! if (derived%use_self_gravity .or. derived%use_viscosity) then
            !     derived%use_bins = .True.
            !     !! Check
            !     if (derived%binning_method < 1 .or. derived%binning_method > 3) then
            !         write (*,*) "ERROR: binning_method must be 1 (equal dr), 2 (equal dA) or 3 (equal Npart)."
            !         stop 1
            !     end if
            !     if (derived%Nbins < 2) then
            !         write (*,*) "ERROR: Can not create less than 2 bins."
            !         stop 1
            !     end if
            !     !! Set default update values
            !     derived%update_rmin_bins = derived%rmin_bins < - tini
            !     derived%update_rmax_bins = derived%rmax_bins < - tini
            !     derived%update_bins = derived%update_rmin_bins .or. derived%update_rmax_bins .or. derived%binning_method == 3
            ! end if         

            !! Error
            if (derived%error_digits < 1) then
                write (*,*) "ERROR: Number of presition digist must be greater than 0."
                stop 1
            end if
            derived%error_tolerance = 10.d0**(-derived%error_digits)
            
            ! Primary Radius
            if (derived%radius_primary <= 0) then
                write (*,*) "ERROR: Primary radius must be positive."
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
                    write (*,*) "ERROR: Can not use paralelization withut OpenMP."
                    stop 1
                end if
                !$ available_threads = OMP_GET_MAX_THREADS()
                !$ if (derived%requested_threads .eq. -1) derived%requested_threads = available_threads
                !$ my_threads = min(available_threads, max(derived%requested_threads,1))
                !$ call OMP_SET_NUM_THREADS(my_threads)
            else
                derived%requested_threads = 1
                my_threads = 1
            end if

            
            ! Output
            if (derived%use_datascreen .and. derived%use_percentage) then
                write (*,*) "ERROR: Can not print both percentage and data on screen."
                stop 1
            end if
            
            derived%use_chaos = derived%use_chaosfile .or. derived%use_chaos
            derived%use_elements = derived%use_elements_output .or. derived%use_chaos

            if (.not. derived%use_chaosfile) use_update_chaos = .False.
            if (.not. derived%use_datafile) use_flush_output = .False.

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
        end subroutine set_derived_parameters
        
        ! 4. Leer el archivo tomfile
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


        ! Deallocate
        subroutine free_parameters_arays()
            implicit none
            if (allocated(moons_arr)) deallocate(moons_arr)
            if (allocated(particles_arr)) deallocate(particles_arr)
            if (allocated(boulders_in)) deallocate(boulders_in)
            if (allocated(moons_in)) deallocate(moons_in)
            if (allocated(particles_in)) deallocate(particles_in)
            if (allocated(boulders_theta_ini)) deallocate(boulders_theta_ini)
            if (allocated(m_arr)) deallocate(m_arr)
            if (allocated(y_arr)) deallocate(y_arr)
            if (allocated(y_arr_new)) deallocate(y_arr_new)
            if (allocated(y_der)) deallocate(y_der)
            if (allocated(tom_times)) deallocate(tom_times)
            if (allocated(tom_deltaomega)) deallocate(tom_deltaomega)
            if (allocated(tom_deltamass)) deallocate(tom_deltamass)
            if (allocated(checkpoint_is_tom)) deallocate(checkpoint_is_tom)
            if (allocated(checkpoint_is_output)) deallocate(checkpoint_is_output)
        end subroutine free_parameters_arays

        ! 21.x Subroutines to write to a file (or screen)

        ! ! 21.1 Write elements to unit_file
        ! subroutine write_elements(i, unit_file)
        !     implicit none
        !     integer(kind=4), intent(in) :: i, unit_file
            
        !     aux_integer = sorted_particles_index(i)
        !     write (unit_file,f233) &
        !         & particles_index(aux_integer), &
        !         & time / unit_time, &
        !         & particles_elem(aux_integer,1) / unit_dist, &
        !         & particles_elem(aux_integer,2), &
        !         & particles_elem(aux_integer,3) / radian, &
        !         & particles_elem(aux_integer,4) / radian, &
        !         & particles_MMR(aux_integer), &
        !         & particles_mass(aux_integer) / unit_mass, &
        !         & particles_dist(aux_integer) / unit_dist, &
        !         & asteroid_omega * unit_time, &
        !         & asteroid_mass / unit_mass
        ! end subroutine write_elements

        ! ! 21.2 Write coordinates to unit_file (particles)
        ! subroutine write_coordinates_particle(i, unit_file)
        !     implicit none
        !     integer(kind=4), intent(in) :: i, unit_file

        !     aux_integer = sorted_particles_index(i)
        !     write (unit_file,f233) &
        !         & particles_index(aux_integer) + Nboulders, &
        !         & time / unit_time, &
        !         & particles_coord(aux_integer,1:2) / unit_dist, &
        !         & particles_coord(aux_integer,3:4) / unit_vel, &
        !         & particles_acc(aux_integer,1:2) / unit_acc, &
        !         & particles_mass(aux_integer) / unit_mass, &
        !         & particles_dist(aux_integer) / unit_dist
        ! end subroutine write_coordinates_particle

        ! ! 21.3 Write coordinates to unit_file (boulders)
        ! subroutine write_coordinates_boulders(i, unit_file)
        !     implicit none
        !     integer(kind=4), intent(in) :: i, unit_file
            
        !     write (unit_file,f233) &
        !         & i, &
        !         & time / unit_time, &
        !         & pos_ast_arr(i,:) / unit_dist, &
        !         & vel_ast_arr(i,:) / unit_vel, &
        !         & acc_ast_arr(i,:) / unit_acc, &
        !         & mass_ast_arr(i) / unit_mass, &
        !         & radius_ast_arr(i) / unit_dist
        ! end subroutine write_coordinates_boulders

        ! ! 21.4 DO NOT Write
        ! subroutine do_not_write(i, unit_file)
        !     implicit none
        !     integer(kind=4), intent(in) :: i, unit_file
            
        ! end subroutine do_not_write

        ! ! 21.5 Write chaos file
        ! subroutine write_chaos(unit_file)
        !     implicit none
        !     integer(kind=4), intent(in) :: unit_file
        !     integer(kind=4) :: i, aux_integer
        !     integer(kind=4), dimension(:), allocatable :: my_sorted_particles_index

        !     ! Sort particles by index
        !     allocate (my_sorted_particles_index(Nparticles))
        !     do i = 1, Nparticles
        !         my_sorted_particles_index = i
        !     end do
        !     call quickargsort_int(particles_index, my_sorted_particles_index, 1, Nparticles)
        !     call fseek(unit_file, 0, 0)       ! move to beginning
        !     ! Remember that Initial conditions were not swapped
        !     do i = 1, Nparticles
        !         aux_integer = my_sorted_particles_index(i)
        !         write (unit_file,f2233) &
        !         & particles_index(aux_integer), & ! i
        !         & particles_outcome(aux_integer), & ! bad
        !         & final_time / unit_time, & ! total time to integrate
        !         & asteroid_initial_conditions(10) / (unit_mass * unit_dist * unit_vel), & ! initial (Asteroid): angular momentum
        !         & particles_initial_conditions(i,1) / unit_mass, & ! initial: mass
        !         & particles_initial_conditions(i,2) / unit_dist, & ! initial: a
        !         & particles_initial_conditions(i,3), & ! initial: e
        !         & particles_initial_conditions(i,4) / radian, & ! initial: M
        !         & particles_initial_conditions(i,5) / radian, & ! initial: omega
        !         & particles_initial_conditions(i,6), & ! initial: MMR
        !         & sqrt(particles_initial_conditions(i,2) * &
        !           & (uno - particles_initial_conditions(i,3)**2) * &
        !           & G * (particles_initial_conditions(i,1) + asteroid_initial_conditions(1))) / &
        !           & (unit_dist * unit_vel), & ! initial: angular momentum per unit mass
        !         & particles_times(aux_integer) / unit_time, & ! surviving time
        !         & particles_elem(aux_integer,1) / unit_dist, particles_elem(aux_integer,2), & ! final: a, e
        !         & particles_elem(aux_integer,3) / radian, particles_elem(aux_integer,4) / radian, & ! final: M, omega
        !         & particles_MMR(aux_integer), & ! final: MMR
        !         & sqrt(particles_elem(aux_integer,1) * &
        !           & (uno - particles_elem(aux_integer,2)**2) * &
        !           & G * (particles_mass(aux_integer) + asteroid_mass)) / &
        !           & (unit_dist * unit_vel), & ! final: angular momentum per unit mass
        !         & particles_min_a(aux_integer) / unit_dist, & ! a_min 
        !         & particles_max_a(aux_integer) / unit_dist, & ! a_max
        !         & particles_min_e(aux_integer), & ! e_min 
        !         & particles_max_e(aux_integer), & ! e_max
        !         & (particles_max_a(aux_integer) - particles_min_a(aux_integer)) / unit_dist, & ! Delta a
        !         & (particles_max_e(aux_integer) - particles_min_e(aux_integer)) ! Delta e
        !     end do
        !     deallocate (my_sorted_particles_index)
        !     flush(unit_file)
        ! end subroutine write_chaos

        ! ! 22.x Subroutines to get coordinates and elements from a body

        
        ! ! 25. Nullify pointers
        ! subroutine nullify_pointers()
        !     implicit none

        !     nullify(particles_hexitptr)
        !     nullify(get_elements_i)
        !     nullify(get_chaos_i)
        !     nullify(resolve_merge)
        !     nullify(write_i_to_general)
        !     nullify(write_i_to_screen)
        !     nullify(write_i_to_individual)
        !     nullify(write_b_to_general)
        !     nullify(write_b_to_screen)
        !     nullify(write_b_to_individual)
        !     nullify(flush_chaos)
        ! end subroutine nullify_pointers
        
        ! ! 26.x Subroutines to switch between lower and upper case

        ! ! 27. Flush to a file
        ! subroutine flush_to_file(unit_file)
        !     implicit none
        !     integer(kind=4), intent(in) :: unit_file

        !     flush(unit_file)
        ! end subroutine flush_to_file
        
        ! ! 28.x Functions to check collision/escape
        
        ! ! 28.1 Check distances (Version 1)
        ! function check_continue_v1 (y) result(keep_going)
        !     implicit none
        !     real(kind=8), dimension(:), intent(in) :: y
        !     real(kind=8) :: rb(2)
        !     real(kind=8) :: rcm(2), dist_to_cm, r_from_cm(2)
        !     real(kind=8) :: rib(0:Nboulders,2), dist_to_i, r_from_i(2)
        !     integer(kind=4) :: j, particle_j, i
        !     integer(kind=4), parameter :: neqs = 4
        !     logical :: keep_going

        !     ! Re-set particles_hexit
        !     particles_hexit = 0  ! 0 is the default or non-action value
            
        !     ! Calculate the center of mass of the asteroid and boulders
        !     rcm = cero
        !     do i = 0, Nboulders 
        !         rib(i,1) = y((i * neqs) + 1)
        !         rib(i,2) = y((i * neqs) + 2)
        !         rcm = rcm + mass_ast_arr(i) * rib(i,:)
        !     end do
        !     rcm = rcm / asteroid_mass ! rcm = sum_i m_i * r_i / M

        !     ! Calculate vector and distance to CM
        !     do j = 1, Nactive
        !         particle_j = (j + Nboulders) * neqs
        !         rb(1) = y(particle_j+1)
        !         rb(2) = y(particle_j+2)
        !         r_from_cm = rb - rcm
        !         dist_to_cm = sqrt(r_from_cm(1)*r_from_cm(1) + r_from_cm(2)*r_from_cm(2))
        !         if (dist_to_cm > max_distance) then
        !             particles_hexit(j) = 2
        !         else if (dist_to_cm < min_distance) then
        !             particles_hexit(j) = 1
        !         else
        !             ! Calculate vector and distance to boulders
        !             boulders_loop: do i = 0, Nboulders
        !                 r_from_i = rb - rib(i,:)
        !                 dist_to_i = sqrt(r_from_i(1)*r_from_i(1) + r_from_i(2)*r_from_i(2))
        !                 if (dist_to_i < radius_ast_arr(i)) then
        !                     particles_hexit(j) = 1
        !                     exit boulders_loop
        !                 end if
        !             end do boulders_loop
        !         end if
        !     end do
            
        !     keep_going = all(particles_hexit(1:Nactive) .eq. 0)
        !     if (.not. keep_going) particles_hexit(0) = 1
            
        ! end function check_continue_v1
        
        ! ! 28.2 Check distances (Version 2) [from cm]
        ! function check_continue_v2 (y) result(keep_going)
        !     implicit none
        !     real(kind=8), dimension(:), intent(in) :: y
        !     real(kind=8) :: rb(2)
        !     real(kind=8) :: rcm(2), dist_to_cm, r_from_cm(2)
        !     real(kind=8) :: rib(0:Nboulders,2), dist_to_i, r_from_i(2)
        !     integer(kind=4) :: j, particle_j, i
        !     integer(kind=4), parameter :: neqs = 4
        !     logical :: keep_going

        !     ! Re-set particles_hexit
        !     particles_hexit = 0  ! 0 is the default or non-action value
            
        !     ! Calculate the center of mass of the asteroid and boulders
        !     rcm(1) = y(3)
        !     rcm(2) = y(4)
        !     do i = 0, Nboulders
        !         rib(i,1) = cos(y(1) + theta_ast_arr(i)) * dist_ast_arr(i)
        !         rib(i,2) = sin(y(1) + theta_ast_arr(i)) * dist_ast_arr(i)
        !     end do

        !     ! Calculate distance and vector to boulders
        !     do j = 1, Nactive
        !         particle_j = j * neqs + 2
        !         rb(1) = y(particle_j+1)
        !         rb(2) = y(particle_j+2)
        !         r_from_cm = rb - rcm
        !         dist_to_cm = sqrt(r_from_cm(1)*r_from_cm(1) + r_from_cm(2)*r_from_cm(2))
        !         if (dist_to_cm > max_distance) then
        !             particles_hexit(j) = 2
        !         else if (dist_to_cm < min_distance) then
        !             particles_hexit(j) = 1
        !         else 
        !             ! Calculate vector and distance to boulders
        !             boulders_loop: do i = 0, Nboulders
        !                 r_from_i = rb - rib(i,:)
        !                 dist_to_i = sqrt(r_from_i(1)*r_from_i(1) + r_from_i(2)*r_from_i(2))
        !                 if (dist_to_i < radius_ast_arr(i)) then
        !                     particles_hexit(j) = 1
        !                     exit boulders_loop
        !                 end if
        !             end do boulders_loop
        !         end if
        !     end do
            
        !     keep_going = all(particles_hexit(1:Nactive) .eq. 0)
        !     if (.not. keep_going) particles_hexit(0) = 1
            
        ! end function check_continue_v2
        
        ! 99 Subrutina para pointer vacío (no hace nada) con input i
        subroutine do_nothing_i(i)
            implicit none
            integer(kind=4), intent(in) :: i
        end subroutine do_nothing_i


        

end module parameters
