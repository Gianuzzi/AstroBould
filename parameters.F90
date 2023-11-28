module parameters
    use const
    use run
    use integrators
    implicit none

    !!! Parámetros y auxiliares
    integer(kind=4), parameter :: Npart = 1
    integer(kind=4)  :: Nboul = 0, Ntot !!! Boulders, Particles, Total
    integer(kind=4), parameter :: neqs = 5 ! Número de eqs. a integrar (m, x, y, vx, vy) !! Omega va aparte
    integer(kind=4)  :: FP = 0, NP = 0 !! First particle, y posicion de primer partícula
    
    !!! Variables
    real(kind=8), allocatable :: m(:), mu(:), Gmi(:), mucm(:)
    real(kind=8) :: GM = cero, mcm = cero
    real(kind=8) :: m0 = cero, R0 = cero, GM0 = cero
    real(kind=8), allocatable :: radius(:)
    real(kind=8), allocatable :: rib(:,:), ria(:,:)
    real(kind=8), allocatable :: vib(:,:), via(:,:)
    real(kind=8), allocatable :: aib(:,:), aia(:,:)
    real(kind=8) :: rcm(2) = cero, vcm(2) = cero, acm(2) = cero
    real(kind=8) :: rb(2) = cero, ra(2) = cero
    real(kind=8) :: vb(2) = cero, va(2) = cero
    real(kind=8) :: ab(2) = cero, aa(2) = cero
    real(kind=8) :: db0 = cero, da0 = cero
    real(kind=8) :: da = cero, de = cero, amax = cero, amin = inf, emax = cero, emin = inf
    real(kind=8) :: xc(6) = cero, xe(6) = cero, xc0(6) = cero, xe0(6) = cero, Res0 = cero
    real(kind=8) :: ea = cero, ee = cero, ei = cero, eM = cero, ew = cero, eO = cero, eR = cero
    real(kind=8), allocatable :: theta_a(:)
    real(kind=8) :: theta_acm = cero
    real(kind=8), allocatable :: r_b(:), theta_b(:)
    real(kind=8) :: lambda = cero, Prot = cero
    real(kind=8) :: omega = cero, wk = cero, a_corot = cero, omega2 = cero, lambda2 = cero
    real(kind=8) :: ang_mom = cero
    real(kind=8) :: rmax = cero, rmin = inf

    !! Parámetros de main
    integer(kind=4) :: i = 0, j = 0, ineqs = 0
    integer(kind=4) :: nin = 0, nsim = 0, bad = 0
    integer(kind=4) :: ngx = 0, ngy = 0
    real(kind=8)    :: xmin = inf, xmax = -inf, ymin = inf, ymax = -inf
    character(30)   :: chn, cha, che, chM, chw, chR
    character(30)   :: datafile, chaosfile, infofile, map_file
    integer(kind=4) :: auxin
    character(30)   :: auxch
    character(1)    :: auxch1
    character(2)    :: auxch2
    logical         :: auxlo, is_number, loini
    logical         :: screen = .True., perc = .False., datas = .False., eleout = .False.
    logical         :: map_pot = .False., infoo = .False., datao = .False., chaos = .False.
    logical         :: exact = .True. 
    integer(kind=4) :: dig_err = 13


    !!!! FORCES
    real(kind=8) :: raux(2) = cero ! rb - rb0
    !!! Stokes
    logical :: lostokes = .false.
    real(kind=8) :: t_stokes = cero, f_stk = uno
    real(kind=8) :: tau_a = inf, tau_e = inf, C_stk = cero, a_stk = uno
    !! Dampings (not stokes)
    real(kind=8) :: tau_o = inf ! Omega
    real(kind=8) :: tau_m = inf ! M0
    !! Geo-Potential
    logical :: loJ2 = .False.
    real(kind=8) :: J2 = cero
    ! !! TriAxial
    ! logical :: loTriAx = .False.
    ! real(kind=8) :: tria = cero, trib = cero, tric = cero
    ! real(kind=8) :: C20 = cero, C22 = cero, Re = cero


    !!! Punteros
    integer, target :: hexit = 0 !  Hard Exit integer
    integer, pointer :: hexitptr ! pointer to Hard Exit

    !! Vectores
    !y = /omega, mi, xi, yi, vxi, vyi, .../
    real(kind=8), allocatable :: yb(:), ybnew(:), ya(:), yanew(:)


    contains
        
        !! Modificar algún valor si se quiere otro Default
        subroutine init_params_default() 
            implicit none
            !!!! Particles
            nsim = 0
            ea = cero ; ee = cero ; eM = cero ; ew = cero ; eR = cero
            !!!! Rotation
            Prot = cero
            lambda = cero
            !!!! Asteroid
            m0 = cero
            R0 = cero
            !!!! Integración
            dig_err = 12 ! Dígitos de presición (error)
            exact   = .True.  ! Método exacto (sin | cos)
            !!!! Forces
            !! Stokes
            lostokes = .False.
            t_stokes = cero
            tau_a = inf
            tau_e = inf
            !! Geo-Potential
            loJ2 = .False.
            J2 = cero
            !!!! Times
            t0       = cero    ! Initial time [days]
            tf       = cero    ! Final time [days]
            dt_min   = cero    ! Min timestep [days] ! Almost unused
            logsp    = .FALSE. ! LogSpaced outputs
            n_points = 0       ! Number of outputs (if logsp=.TRUE. or dt_out=0)
            dt_out   = cero    ! Output timestep [days] (if logsp = .False.)
            !!!! Salida
            screen  = .False. ; perc = .False. ! Pantalla, porcentaje
            datas   = .False. ! Salida de datos en pantalla
            eleout  = .False. ! Salida en elementos orbitales
            !!!!!! Información (Poner "" o "no" para no crear archivo)
            infofile = "" ! Archivo de info
            datafile = ""! Archivo de datos
            chaosfile = "chaos.dat" ! Archivo de datos
            !!!! Mapa de potencial
            map_file = "" ! Archivo de mapas
            ngx = 500   ; ngy = 500
            xmin = -300 ; xmax = 300 ; ymin = -300 ; ymax = 300
            !!!! Particle elements
            ea = cero                  ! Element a of the particle (km)
            ee = cero !0.1d0           ! ecc
            eM = cero                  ! Mean anomaly (deg)
            ew = cero                  ! Pericenter argument (deg)
            eR = cero !1.001d0 !3.3d0  ! resonancia nominal correspondiente
        end subroutine init_params_default

        subroutine read_conf_file(file_name, existe)
            character(LEN=*), intent(in) :: file_name
            logical, intent(inout) :: existe
            character(256) :: line, param_str, value_str
            character(15) :: auxch15
            integer :: colonPos
            integer :: commentPos
            integer(kind=4) :: nlines, io, len_val

            inquire(file=trim(file_name), exist=existe)
            if (.not. existe) then
                if (screen) write(*,*) "El archivo ", trim(file_name), " no existe."
                return
            end if
            nlines = 0
            open(1, file=trim(file_name), status='old', action='read')
            do
                ! Read the line from the file
                read(1, '(A)', END=98) line
                nlines = nlines + 1

                if ((len_trim(line) == 0) .or. (auxch15(:2) == "c ") .or. (auxch15(:2) == "! ")) cycle

                colonPos = index(line, ':') !this should be 50

                commentPos = index(line, '!')
                if (commentPos == 0) commentPos = len_trim(line)+1 !this should be 80 at least
                len_val = len_trim(adjustl(line(colonPos + 1 : commentPos - 1)))
                if (colonPos > 0) then
                    param_str = trim(adjustl(line(: colonPos)))
                    value_str = trim(adjustl(line(colonPos + 1 : commentPos - 1)))
                    auxch1 = value_str(:1)
                    auxch2 = value_str(:2)
                    auxch15 = param_str(:15)
                    select case (auxch15)
                        case("initial time fo")
                            read(value_str, *) t0
                        case("total integrati")
                            read(value_str, *) tf
                        case("output time int")
                            read(value_str, *) dt_out
                        case("number of outpu")
                            read(value_str, *) n_points
                        case("logspaced outpu")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                logsp = .true.
                            else
                                logsp = .false.
                            end if
                        case("precision (digi")
                            read(value_str, *) dig_err
                        case("learning rate")
                            read(value_str, *) beta
                        case("mass of primary")
                            read(value_str, *) m0
                        case("radius of prima")
                            read(value_str, *) R0
                        case("ratio of spin t")
                            read(value_str, *) lambda
                        case("rotational peri")
                            read(value_str, *) Prot
                        case("number of bould")
                            read(value_str, *) Nboul
                            if (Nboul < 1) then
                                write(*,*) "El número de boulders debe ser mayor que 0."
                                stop 1
                            end if
                            !! Parámetros y auxiliares
                            Ntot = Nboul + Npart
                            !! Allocamos
                            if (.not. allocated(m)) then
                                call allocate_all(Nboul)
                            else 
                                write(*,*) "ERROR: Redefinido."
                            end if
                        case("include stokes-")
                            if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                                lostokes = .True.
                            else
                                lostokes = .False.
                            end if
                        case("a damping chara")
                            read(value_str, *) tau_a
                        case("e damping chara")
                            read(value_str, *) tau_e
                        case("F damping chara")
                            read(value_str, *) t_stokes
                        case("mass damping ch")
                            read(value_str, *) tau_m
                        case("rotation dampin")
                            read(value_str, *) tau_O
                        case("geo-potential J")
                            read(value_str, *) J2
                        ! case("include tri-axi")
                        !     if (((auxch1 == "y") .or. (auxch1 == "s"))) then
                        !         loTriAx = .True.
                        !     else
                        !         loTriAx = .False.
                        !     end if
                        ! case("first ellipsoid")
                        !     read(value_str, *) tria
                        ! case("second ellipsoi")
                        !     read(value_str, *) trib
                        ! case("third ellipsoid")
                        !     read(value_str, *) tric
                        case("min distance fr")
                            read(value_str, *) rmin
                        case("max distance fr")
                            read(value_str, *) rmax
                        case("information on ")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                screen = .true.
                            else
                                screen = .false.
                            end if
                        case("information sum")
                            if (auxch2 == "no") then
                                infoo    = .False.
                                infofile = ""
                            else if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                infoo = .true.
                                infofile = "info.dat"
                            else
                                infoo = .true.
                                infofile = trim(value_str)                            
                            end if
                        case("general output")
                            if (auxch2 == "no") then
                                datao    = .False.
                                datafile = ""
                            else if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                datao = .true.
                                datafile = "salida.dat"
                            else
                                datao = .true.
                                datafile = trim(value_str)                            
                            end if
                        case("chaos indicator")
                            if (auxch2 == "no") then
                                chaos    = .False.
                                chaosfile = ""
                            else if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                chaos = .true.
                                chaosfile = "chaos.dat"
                            else
                                chaos = .true.
                                chaosfile = trim(value_str)                            
                            end if
                        case("output on scree")
                            if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                datas = .true.
                                perc = .false.
                            else if (auxch1 == "%") then
                                datas = .true.
                                perc = .true.
                            else
                                datas = .false.
                                perc = .false.
                            end if
                        case("output variable")
                            if (auxch1 == "c") then
                                eleout = .false.
                            else
                                eleout = .true.
                            end if
                        case("create map file")
                            if (auxch2 == "no") then
                                map_pot    = .False.
                                map_file = ""
                            else if ((auxch1 == "y") .or. (auxch1 == "s")) then
                                map_pot = .true.
                                map_file = "mapas.dat"
                            else
                                map_pot = .true.
                                map_file = trim(value_str)                            
                            end if
                        case("number x cells")
                            read(value_str, *) ngx
                        case("number y cells")
                            read(value_str, *) ngy
                        case("lower x bound (")
                            read(value_str, *) xmin
                        case("upper x bound (")
                            read(value_str, *) xmax
                        case("lower y bound (")
                            read(value_str, *) ymin
                        case("upper y bound (")
                            read(value_str, *) ymax
                        case default
                            write(*,*) "WARNING: Parámetro no reconocido: ", trim(param_str)
                    end select
                else
                    param_str = trim(adjustl(line)) 
                    if (param_str(1:7) .eq. "mass_m0") then
                        if (Nboul < 1) then
                            write(*,*) "ERROR: El número de boulders debe especificarse antes que los parametros."
                            stop 1
                        end if
                        do j = 1, Nboul
                            io = 0
                            check: do while (io == 0)
                                nlines = nlines + 1
                                read(1, '(A)') line
                                if ((len_trim(line) == 0) .or. (line(:2) == "c ") .or. (line(:2) == "! ")) cycle
                                backspace(1)
                                read(1, *, iostat=io) mu(j), radius(j), theta_a(j)
                                if (io /= 0) then
                                    write(*,*) "ERROR: Al leer boulder. Línea: ", nlines
                                    stop 1
                                end if
                                io = 1
                            end do check
                        end do
                    end if
                end if
            end do
            98 close(1)
        end subroutine read_conf_file

        subroutine set_derived_params()
            if (Nboul < 1) then
                write(*,*) "ERROR: Nboul < 1"
                stop 1
            end if
            !! Indices
            FP = Nboul + 1
            NP = FP * neqs
            
            !! Forces
            !!! tau_m y tau_o
            if (abs(tau_m) < tini) tau_m = cero
            if (abs(tau_o) < tini) tau_o = cero
            !!! tau_a y tau_e
            if ((abs(tau_a) < tini) .or. .not. lostokes) tau_a = cero
            if ((abs(tau_e) < tini) .or. .not. lostokes) tau_e = cero
            if ((abs(t_stokes) < tini) .or. .not. lostokes) t_stokes = cero
            if (((abs(tau_a) < tini) .and. (abs(tau_e) < tini)) .or. (abs(t_stokes) < tini)) lostokes = .False.
            !!! Geo-Potential
            if (abs(J2) > tini) loJ2 = .True.
            ! !!! Triaxial
            ! if ((tria < tini) .or. .not. loTriAx) tria = cero
            ! if ((trib < tini) .or. .not. loTriAx) trib = cero
            ! if ((tric < tini) .or. .not. loTriAx) tric = cero
            ! if ((abs(tria) < tini) .or. (abs(trib) < tini) .or. (abs(tric) < tini)) loTriAx = .False.

            ! if (loTriAx .and. (Nboul > 0)) then
            !     write(*,*) "ERROR: No se puede usar triaxialidad con boulders."
            !     stop 1
            ! end if

            !! Error
            if (dig_err < 1) then
                write(*,*) "ERROR: El número de dígitos de presición debe ser mayor que 0."
                stop 1
            end if
            e_tol = 10.d0**(-dig_err)
            
            !! Cuerpo central
            !!!! Masa
            if (m0 <= 0) then
                write(*,*) "ERROR: La masa del cuerpo central debe ser positiva."
                stop 1
            end if
            m(0) = m0
            m0 = m0 * unit_m
            !!!! Radio
            if (R0 <= 0) then
                write(*,*) "ERROR: El radio del cuerpo central debe ser positivo."
                stop 1
            end if
            radius(0) = R0
            R0 = R0 * unit_r
            !!! GM0
            GM0 = G * m0
            !! NO TOCAR
            hexitptr => hexit ! Puntero a Hard Exit (NO TOCAR)
        end subroutine set_derived_params

        subroutine allocate_all(Nboul)
            integer(kind=4), intent(in) :: Nboul
            !! Parámetros y auxiliares
            Ntot = Nboul + Npart
            allocate(m(0:Ntot), mu(Ntot), Gmi(0:Ntot), mucm(0:Nboul))
            allocate(radius(0:Ntot))
            allocate(rib(0:Ntot,2), ria(1:Ntot,2))
            allocate(vib(0:Ntot,2), via(1:Ntot,2))
            allocate(aib(0:Ntot,2), aia(1:Ntot,2))
            allocate(theta_a(1:Nboul))
            allocate(r_b(0:Nboul), theta_b(0:Nboul))
            allocate(yb(1 + neqs * (Ntot+1)), ybnew(1 + neqs * (Ntot+1)))
            allocate(ya(1 + neqs * (Ntot+1)), yanew(1 + neqs * (Ntot+1)))
            m = cero ; mu = cero ; Gmi = cero ; mucm = cero
            radius = cero
            rib = cero ; ria = cero
            vib = cero ; via = cero
            aib = cero ; aia = cero
            theta_a = cero
            r_b = cero ; theta_b = cero
            yb = cero ; ybnew = cero
            ya = cero ; yanew = cero
        end subroutine allocate_all

        subroutine deallocate_all()
            implicit none
            deallocate(m, mu, Gmi, mucm)
            deallocate(radius)
            deallocate(rib, ria)
            deallocate(vib, via)
            deallocate(aib, aia)
            deallocate(theta_a)
            deallocate(r_b, theta_b)
            deallocate(yb, ybnew)
            deallocate(ya, yanew)
        end subroutine deallocate_all

        subroutine rcmfromast(ria,via,aia,m,rcm,vcm,acm)
            implicit none
            real(kind=8), intent(in)  :: ria(Nboul,2), via(Nboul,2), aia(Nboul,2)
            real(kind=8), intent(in)  :: m(0:Nboul)
            real(kind=8), intent(out) :: rcm(2), vcm(2), acm(2)
            real(kind=8) :: mucm(Nboul)
            integer(kind=4) :: i

            rcm = cero ; vcm = cero ; acm = cero
            mucm = m(1:Nboul) / sum(m)
            do i = 1, Nboul!size(ria,1)
                rcm = rcm + mucm(i) * ria(i,:) ! Posición del centro de masas
                vcm = vcm + mucm(i) * via(i,:) ! Velocidad del centro de masas
                acm = acm + mucm(i) * aia(i,:) ! Aceleración del centro de masas
            end do
        end subroutine rcmfromast

        subroutine ast2bar(ra,va,aa,rcm,vcm,acm,rb,vb,ab)
            implicit none
            real(kind=8), intent(in)  :: ra(2), va(2), aa(2)
            real(kind=8), intent(in)  :: rcm(2), vcm(2), acm(2)
            real(kind=8), intent(out) :: rb(2), vb(2), ab(2)
            rb = ra - rcm ! Posición
            vb = va - vcm ! Velocidad
            ab = aa - acm ! Aceleración
        end subroutine ast2bar

        subroutine bar2ast(rb,vb,ab,rcm,vcm,acm,ra,va,aa)
            implicit none
            real(kind=8), intent(in)  :: rb(2), vb(2), ab(2)
            real(kind=8), intent(in)  :: rcm(2), vcm(2), acm(2)
            real(kind=8), intent(out) :: ra(2), va(2), aa(2)
            ra = rb + rcm ! Posición
            va = vb + vcm ! Velocidad
            aa = ab + acm ! Aceleración
        end subroutine bar2ast

        subroutine ang_mom_bar(rib, vib, mi, radius, omega, ang_mom)
            implicit none
            real(kind=8), intent(in)  :: rib(0:Ntot,2), vib(0:Ntot,2), mi(0:Ntot), radius(0:Ntot)
            real(kind=8), intent(in)  :: omega
            real(kind=8), intent(out) :: ang_mom
            integer(kind=4) :: i
            ang_mom = cero
            do i = 0, Nboul
                ang_mom = ang_mom + mi(i) * (rib(i,1) * vib(i,2) - rib(i,2) * vib(i,1)) ! Traslacional
                ang_mom = ang_mom + 0.4d0 * mi(i) * radius(i) * radius(i) * omega ! Rotacional
            end do
            do i = Nboul+1, Ntot
                ang_mom = ang_mom + mi(i) * (rib(i,1) * vib(i,2) - rib(i,2) * vib(i,1)) ! Traslacional
            end do
        end subroutine ang_mom_bar

end module parameters
