!> Module with TimesOmegaMass routines.
module tomodule
    use constants, only: cero, uno
    implicit none

    type :: tom_st  !!! This contains only the input parameters
        integer(kind=4) :: index_number = -1  ! Index to count which TOM row is active
        integer(kind=4) :: total_number = 0  ! Total amount of lines in TOM
        real(kind=8), dimension(:), allocatable :: times  ! Times in TOM
        real(kind=8), dimension(:), allocatable :: deltaomega  ! Delta Omega in TOM
        real(kind=8), dimension(:), allocatable :: deltamass  ! Delta Mass in TOM
        logical :: use_dmass = .False.
        logical :: use_domega = .False.
        real(kind=8) :: mass_growth_param = cero
    end type tom_st

    contains

        ! Read TOMfile
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
                write(*,*) "ERROR: No se encontró el archivo TOM: ", trim(file_tout)
                stop 1
            end if

            open (unit=12, file=file_tout, status="old", action="read")

            !! Count number of columns
            ncols = 0
            read (12, '(A)') auxstr
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
            
            rewind (12) ! Go to the beginning of the file
            do ! Count number of (valid) lines 
                read (12, *, iostat=io) t_aux
                if ((io /= 0) .or. (t_aux > tf)) exit
                if (t_aux < t0) cycle
                n_TOM = n_TOM + 1
            end do

            if (n_TOM == 2) then
                write(*,*) "ERROR: No se encontraron datos válidos en el archivo."
                stop 1
            end if
            
            ! Allocate arrays
            allocate (t_TOM(n_TOM))
            t_TOM = -uno
            if (ncols > 1) then
                allocate (domega_TOM(n_TOM))
                domega_TOM = -uno
            end if
            if (ncols > 2) then
                allocate (dmass_TOM(n_TOM))
                dmass_TOM = -uno
            end if
            rewind (12) ! Go to the beginning of the file

            ! Read data
            i = 2
            if (ncols == 1) then
                do
                    read (12, *, iostat=io) t_TOM(i)
                    if ((io /= 0) .or. (t_TOM(i) > tf)) exit
                    if (t_TOM(i) < t0) cycle
                    i = i + 1
                end do
            else if (ncols == 2) then
                do
                    read (12, *, iostat=io) t_TOM(i), domega_TOM(i)
                    if ((io /= 0) .or. (t_TOM(i) > tf)) exit
                    if (t_TOM(i) < t0) cycle
                    i = i + 1
                end do
            else 
                do
                    read (12, *, iostat=io) t_TOM(i), domega_TOM(i), dmass_TOM(i)
                    if ((io /= 0) .or. (t_TOM(i) > tf)) exit
                    if (t_TOM(i) < t0) cycle
                    i = i + 1
                end do
            end if

            close (12)
        end subroutine read_tomfile

        ! Set TOM times and arrays
        subroutine setup_TOM(self, tomfile, final_time, use_screen)
            implicit none
            type(tom_st), intent(inout) :: self
            character(LEN=*), intent(in) :: tomfile
            real(kind=8), intent(in) :: final_time
            logical, intent(in) :: use_screen

            !! En este caso, leeremos los tiempos desde un archivo
            if (use_screen) write(*,*) "Reading times from TOM file: ", trim(tomfile)
            call read_tomfile(cero, final_time, & 
                            & self%times, self%deltaomega, self%deltamass, tomfile) ! Read LOOP checkpoints

            self%total_number = size(self%times, 1)

            !!! Unidades
            self%times = self%times
            if (allocated(self%deltaomega)) self%deltaomega = self%deltaomega
            if (allocated(self%deltamass)) self%deltamass = self%deltamass

            !!! Condicion inicial (y final)
            self%times(1) = cero
            self%times(self%total_number) = final_time
            if (allocated(self%deltaomega)) then
                self%use_domega = .True.
                self%deltaomega(1) = cero
                self%deltaomega(self%total_number) = cero
            end if
            if (allocated(self%deltamass)) then
                self%use_dmass = .True.
                self%deltamass(1) = cero
                self%deltamass(self%total_number) = cero
            end if

            !!! Mensaje !
            if (use_screen) then
                if (self%use_dmass) then
                    write(*,*) "  - 3 columns read: t, Delta_omega, Delta_m"
                else if (self%use_domega) then
                    write(*,*) "  - 2 columns read: t, Delta_omega"
                else
                    write(*,*) "  - 1 column read: t"
                end if
            end if

        end subroutine setup_TOM


        ! Deallocate all params arrays
        pure subroutine free_tom(self)
            implicit none
            type(tom_st), intent(inout) :: self
            if (allocated(self%times)) deallocate(self%times)
            if (allocated(self%deltaomega)) deallocate(self%deltaomega)
            if (allocated(self%deltamass)) deallocate(self%deltamass)
        end subroutine free_tom
    
end module tomodule