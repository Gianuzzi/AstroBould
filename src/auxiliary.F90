!> Module with array, string, ordering, ... routines.
module auxiliary
    use constants, only: wp

    implicit none

contains

    ! From lower to upper case
    pure function to_upper(strIn) result(strOut)
        ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
        ! Original author: Clive Page
        implicit none
        character(len=*), intent(in) :: strIn
        character(len=len(strIn)) :: strOut
        integer (kind=4) :: i, j

        do i = 1, len(strIn)
            j = iachar(strIn(i:i))
            if (j >= iachar("a") .and. j <= iachar("z")) then
                strOut(i:i) = achar(iachar(strIn(i:i)) - 32)
            else
                strOut(i:i) = strIn(i:i)
            end if
        end do
    end function to_upper

    ! From lower to upper case
    pure function to_lower(strIn) result(strOut)
        implicit none
        character(len=*), intent(in) :: strIn
        character(len=len(strIn)) :: strOut
        integer (kind=4) :: i, j

        do i = 1, len(strIn)
            j = iachar(strIn(i:i))
            if (j >= iachar("A") .and. j <= iachar("Z")) then
                strOut(i:i) = achar(iachar(strIn(i:i)) + 32)
            else
                strOut(i:i) = strIn(i:i)
            end if
        end do
    end function to_lower

    ! Get indices to order a real array
    pure recursive subroutine quickargsort(a, b, first, last)
        implicit none
        real(wp), intent(in) :: a(:)          ! Input array (not modified)
        integer(kind=4), intent(inout) :: b(size(a))    ! Indices array to be sorted
        integer(kind=4), intent(in) :: first, last
        real(wp) :: pivot
        integer(kind=4) :: i, j, temp

        ! Pivot selection
        pivot = a(b(floor((first + last)*0.5e0_wp)))
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
        if (j + 1 < last) call quickargsort(a, b, j + 1, last)
    end subroutine quickargsort

    ! Get indices to order an integer array
    pure recursive subroutine quickargsort_int(a, b, first, last)
        implicit none
        integer(kind=4), intent(in) :: a(:)    ! Input array (not modified)
        integer(kind=4), intent(inout) :: b(size(a)) ! Indices array to be sorted
        integer(kind=4), intent(in) :: first, last
        integer(kind=4) :: pivot
        integer(kind=4) :: i, j, temp

        ! Pivot selection
        pivot = a(b(floor((first + last)*0.5e0_wp)))
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
        if (j + 1 < last) call quickargsort_int(a, b, j + 1, last)
    end subroutine quickargsort_int

    ! Order a real array
    pure subroutine quicksort(a, first, last)
        implicit none
        real(wp), intent(inout) :: a(:)          ! Input array
        integer(kind=4), intent(in) :: first, last   ! First and last indices
        integer(kind=4), allocatable :: b(:)         ! Indices of array to be sorted
        real(wp), allocatable :: a_copy(:)       ! Temporary copy of the array
        integer(kind=4) :: n_sort, i

        ! Total amount of values to sort
        n_sort = last - first + 1

        ! Allocate and initialize the indices array
        allocate (b(1:n_sort))
        allocate (a_copy(1:n_sort))

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
        deallocate (b, a_copy)
    end subroutine quicksort

    ! Order an integer array
    pure subroutine quicksort_int(a, first, last)
        implicit none
        integer(kind=4), intent(inout) :: a(:)       ! Input array
        integer(kind=4), intent(in) :: first, last   ! First and last indices
        integer(kind=4), allocatable :: b(:)         ! Indices of array to be sorted
        integer(kind=4), allocatable :: a_copy(:)    ! Temporary copy of the array
        integer(kind=4) :: n_sort, i

        ! Total amount of values to sort
        n_sort = last - first + 1

        ! Allocate and initialize the indices array
        allocate (b(1:n_sort))
        allocate (a_copy(1:n_sort))

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
        deallocate (b, a_copy)
    end subroutine quicksort_int

    ! Reorder a real array from given order
    pure subroutine reorder(array, order, n)
        real(wp), dimension(:), intent(inout) :: array
        integer(kind=4), dimension(:), intent(in) :: order
        integer(kind=4), intent(in) :: n
        real(wp), dimension(n) :: tmp
        integer(kind=4) :: i

        ! Loop and reorder
        do i = 1, n
            tmp(i) = array(order(i))
        end do
        array(1:n) = tmp
    end subroutine reorder

    ! Reorder an integer array from given order
    pure subroutine reorder_int(array, order, n)
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

    ! Reorder a 2D real array from given order, along 1st axis
    pure subroutine reorder2D(array, order, n)
        real(wp), dimension(:, :), intent(inout) :: array
        integer(kind=4), dimension(:), intent(in) :: order
        integer(kind=4), intent(in) :: n
        real(wp), dimension(n, size(array, 2)) :: tmp
        integer(kind=4) :: i

        ! Loop and reorder the first axis based on the order array
        do i = 1, n
            tmp(i, :) = array(order(i), :)
        end do

        ! Copy the reordered array back
        array(1:n, :) = tmp
    end subroutine reorder2D

    ! Combine, order, and remove duplicated values between 2 real arrays
    pure subroutine merge_sort_and_unique(a, b, ina, inb, c, kfin)
        implicit none
        real(wp), dimension(:), intent(in) :: a, b ! Arreglos de reales
        logical, dimension(:), allocatable, intent(out) :: ina, inb ! Booleanos de elementos en a y b
        real(wp), dimension(:), allocatable, intent(out) :: c ! Arreglo combinado
        integer(kind=4), intent(out) :: kfin ! Número de elementos en c
        logical, dimension(:), allocatable :: ina0, inb0
        real(wp), dimension(:), allocatable :: c0
        integer(kind=4) :: i, j, k

        allocate (c0(size(a) + size(b)))
        allocate (ina0(size(c0)))
        allocate (inb0(size(c0)))
        c0 = 0.e0_wp
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

    ! Calculate z component of 2D cross product
    pure function cross2D_z(a, b) result(res)
        implicit none
        real(wp), dimension(2), intent(in) :: a, b
        real(wp) :: res

        res = a(1)*b(2) - a(2)*b(1) ! Solo la componente z
    end function cross2D_z

    ! Calculate z component of 3D cross product
    pure function cross3D_z(a, b) result(res)
        implicit none
        real(wp), dimension(3), intent(in) :: a, b
        real(wp) :: res

        res = a(1)*b(2) - a(2)*b(1) ! Solo la componente z
    end function cross3D_z

    ! Rotate a 2D vector an angle theta counter-clockwise
    pure function rotate2D(a, theta) result(res)
        implicit none
        real(wp), dimension(2), intent(in) :: a
        real(wp), intent(in) :: theta
        real(wp), dimension(2) :: res

        res = (/a(1)*cos(theta) - a(2)*sin(theta), & ! x cos(th) - y sin(th),
                a(1)*sin(theta) + a(2)*cos(theta)/)  ! x sin(th) + y cos(th)
    end function rotate2D

    ! Read a file with data structured in columns
    subroutine read_columns_file(file_name, values_arr, method)
        implicit none
        character(LEN=*), intent(in) :: file_name
        real(wp), dimension(:, :), allocatable, intent(out) :: values_arr
        integer(kind=4), optional :: method
        integer(kind=4), parameter :: MAX_COLS = 8, MAX_ROWS = 10000
        real(wp), dimension(:, :), allocatable :: aux_real_arr
        integer(kind=4) :: ncols, nrows
        integer(kind=4) :: i, j, io, my_method
        real(wp) :: aux_real
        character(260) :: auxstr
        logical :: existe

        ! Init
        ncols = 0
        nrows = 0

        if (present(method)) then
            my_method = method
        else
            my_method = 0
        end if

        inquire (file=trim(file_name), exist=existe)
        if (.not. existe) then
            write (*, *) "ERROR: File not found: ", trim(file_name)
            stop 1
        end if

        open (unit=11, file=trim(file_name), status="old", action="read")

        !! Count number of columns
        read (11, '(A)') auxstr
        do i = 1, MAX_COLS
            io = 0
            read (auxstr, *, iostat=io) (aux_real, j=1, i)
            if (io /= 0) exit
        end do
        if (io == 0) then
            ncols = i
        else
            ncols = i - 1
        end if

        rewind (11) ! Go to the beginning of the file

        if (my_method == 0) then

            !! Count number of rows
            do ! Count number of (valid) lines
                read (11, *, iostat=io) aux_real
                if (is_iostat_end(io)) exit
                nrows = nrows + 1
            end do

            ! Allocate arrays
            allocate (values_arr(nrows, ncols))
            rewind (11) ! Go to the beginning of the file
            do i = 1, nrows
                read (11, *) (values_arr(i, j), j=1, ncols)
            end do

        else
            ! Allocate auxiliar array
            allocate (aux_real_arr(MAX_ROWS, ncols))
            do i = 1, MAX_ROWS
                read (11, *, iostat=io) (aux_real_arr(i, j), j=1, ncols)
                if (is_iostat_end(io)) exit
                nrows = nrows + 1
            end do
            ! Allocate arrays
            allocate (values_arr(nrows, ncols))
            values_arr = aux_real_arr(1:nrows, 1:ncols)
            deallocate (aux_real_arr)
        end if
        close (11)
    end subroutine read_columns_file

    ! Write percentage to file unit or std out
    subroutine percentage(tout, tstop, file_unit)
        implicit none
        real(wp), intent(in) :: tout, tstop
        integer(kind=4), intent(in), optional :: file_unit
        integer(kind=4) :: funit = 6  ! 6 is STD_OUTPUT
        integer(kind=4) :: iper
        character(len=100) :: guiones
        character(len=4) :: cestado

        if (present(file_unit)) funit = file_unit

        iper = int(100.0e0_wp*tout/tstop)
        guiones = repeat('.', iper)

        if (iper < 100) then
            write (cestado, '(i3)') iper
            write (funit, '(a)', advance='no') char(13)//trim(guiones)//trim(adjustl(cestado))//'%'
        else
            write (funit, '(a)', advance='no') char(13)//trim(guiones)//'. FIN'
            write (funit, *)
        end if
        flush (funit)
    end subroutine percentage

end module auxiliary
