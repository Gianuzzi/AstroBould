!> Module with array, string, ordering, ... routines.
module auxiliary
    use constants, only: wp
    use iso_fortran_env, only: int64

    implicit none
    public

    !> Bodies sorted by square cell (row-major key), for neighbour searches that are linear in the
    !> number of bodies with any cell size: only occupied cells are stored (no dense grid), so the
    !> cells can be as small as the contact distance whatever the extent of the system.
    !> Neighbours of the body in sorted slot q, in the cell row dcy (-1, 0, 1), have keys in
    !> [key(q) + dcy*Wx - 1, key(q) + dcy*Wx + 1]: a contiguous range of slots.
    type :: cell_index_t
        integer(kind=4) :: n = 0
        real(wp) :: L = 0.0_wp                    ! Cell size
        integer(int64) :: Wx = 0_int64            ! Key stride of a row (with a spare column per side)
        integer(kind=4), allocatable :: ord(:)    ! ord(q): input position (1..n) of the body in slot q
        integer(int64), allocatable :: key(:)     ! Sorted keys
        real(wp), allocatable :: px(:), py(:)     ! Positions, in sorted order
        integer(int64), allocatable :: ktmp(:)    ! Work (radix sort)
        integer(kind=4), allocatable :: otmp(:)   ! Work (radix sort)
    end type cell_index_t

contains

    !> Sort the n positions pos(1:2, 1:n) by cells of size L (stable: equal keys keep their input
    !> order, so the result is deterministic)
    subroutine cell_index_build(ci, n, pos, L)
        implicit none
        type(cell_index_t), intent(inout) :: ci
        integer(kind=4), intent(in) :: n
        real(wp), intent(in) :: pos(:, :)
        real(wp), intent(in) :: L
        integer(kind=4), parameter :: BITS = 11, NBUCK = 2**BITS
        integer(kind=4) :: cnt(0:NBUCK-1), i, b, shift
        integer(int64) :: cx, cy, cxmin, cxmax, cymin, cymax, kmax
        real(wp) :: Lc

        ci%n = n
        if (allocated(ci%key)) then
            if (size(ci%key) < n) deallocate (ci%ord, ci%key, ci%px, ci%py, ci%ktmp, ci%otmp)
        end if
        if (.not. allocated(ci%key)) then
            allocate (ci%ord(max(n, 1)), ci%key(max(n, 1)), ci%px(max(n, 1)), ci%py(max(n, 1)), &
                    & ci%ktmp(max(n, 1)), ci%otmp(max(n, 1)))
        end if
        if (n < 1) return

        ! Cell range. Enlarge the cells only if the keys would not fit in 62 bits (not in practice)
        Lc = L
        do
            cxmin = huge(1_int64); cxmax = -huge(1_int64)
            cymin = huge(1_int64); cymax = -huge(1_int64)
            do i = 1, n
                cx = floor(pos(1, i)/Lc, int64)
                cy = floor(pos(2, i)/Lc, int64)
                cxmin = min(cxmin, cx); cxmax = max(cxmax, cx)
                cymin = min(cymin, cy); cymax = max(cymax, cy)
            end do
            if (real(cxmax - cxmin + 3, wp)*real(cymax - cymin + 3, wp) < 4.0e18_wp) exit
            Lc = 2.0_wp*Lc
        end do
        ci%L = Lc
        ci%Wx = cxmax - cxmin + 3_int64  ! Spare column on each side: neighbour keys never wrap rows

        do i = 1, n
            cx = floor(pos(1, i)/Lc, int64)
            cy = floor(pos(2, i)/Lc, int64)
            ci%key(i) = (cy - cymin + 1_int64)*ci%Wx + (cx - cxmin + 1_int64)
            ci%ord(i) = i
        end do

        ! LSD radix sort of the keys (stable), carrying the permutation
        kmax = maxval(ci%key(1:n))
        shift = 0
        do while (shiftr(kmax, shift) > 0_int64)
            cnt = 0
            do i = 1, n
                b = int(iand(shiftr(ci%key(i), shift), int(NBUCK - 1, int64)))
                cnt(b) = cnt(b) + 1
            end do
            do b = 1, NBUCK - 1
                cnt(b) = cnt(b) + cnt(b - 1)
            end do
            do i = n, 1, -1
                b = int(iand(shiftr(ci%key(i), shift), int(NBUCK - 1, int64)))
                ci%ktmp(cnt(b)) = ci%key(i)
                ci%otmp(cnt(b)) = ci%ord(i)
                cnt(b) = cnt(b) - 1
            end do
            ci%key(1:n) = ci%ktmp(1:n)
            ci%ord(1:n) = ci%otmp(1:n)
            shift = shift + BITS
        end do

        ! Positions in sorted order (sequential memory access in the searches)
        do i = 1, n
            ci%px(i) = pos(1, ci%ord(i))
            ci%py(i) = pos(2, ci%ord(i))
        end do
    end subroutine cell_index_build

    !> Move ptr to the first sorted slot with key >= kval. The first call (ok = .False.) uses a
    !> binary search; later calls only move forward, which is valid when kval does not decrease
    !> (as when the slots q are visited in increasing order). Then ptr..n with key <= kval + 2
    !> is the cell row range.
    pure subroutine cell_index_advance(ci, kval, ptr, ok)
        implicit none
        type(cell_index_t), intent(in) :: ci
        integer(int64), intent(in) :: kval
        integer(kind=4), intent(inout) :: ptr
        logical, intent(inout) :: ok
        integer(kind=4) :: lo, hi, mid

        if (.not. ok) then
            lo = 1
            hi = ci%n + 1
            do while (lo < hi)
                mid = (lo + hi)/2
                if (ci%key(mid) < kval) then
                    lo = mid + 1
                else
                    hi = mid
                end if
            end do
            ptr = lo
            ok = .True.
        else
            do while (ptr <= ci%n)
                if (ci%key(ptr) >= kval) exit
                ptr = ptr + 1
            end do
        end if
    end subroutine cell_index_advance

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
        integer(kind=4), intent(inout) :: b(:)    ! Indices array to be sorted (may be shorter than a: only b(first:last) is used)
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
        integer(kind=4), intent(inout) :: b(:) ! Indices array to be sorted (only b(first:last) is used)
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
        real(wp) :: cos_theta, sin_theta
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        res = (/a(1)*cos_theta - a(2)*sin_theta, & ! x cos(th) - y sin(th),
                a(1)*sin_theta + a(2)*cos_theta/)  ! x sin(th) + y cos(th)
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

    subroutine percentage(tout, tstop, file_unit)
        implicit none

        real(wp), intent(in) :: tout, tstop
        integer(kind=4), intent(in), optional :: file_unit

        integer(kind=4) :: funit, iper
        character(len=100) :: bar
        character(len=4)   :: cestado

        logical, save :: first_call = .True.

        funit = 6
        if (present(file_unit)) funit = file_unit

        iper = min(100, nint(100.0_wp*tout/tstop))
        bar  = repeat('.', iper)

        ! Move cursor up and clear previous output
        if (.not. first_call) then
            write(funit,'(a)',advance='no') achar(27)//'[2F' ! up 2 lines
            write(funit,'(a)',advance='no') achar(27)//'[2K' ! clear line
            write(funit,'(a)',advance='no') achar(27)//'[1B' ! down 1 line
            write(funit,'(a)',advance='no') achar(27)//'[2K' ! clear line
            write(funit,'(a)',advance='no') achar(27)//'[1A' ! back up 1 line
        end if

        ! First line
        write(funit,'(a,1x,1PE15.7," / ",1PE15.7)') "Done:", tout, tstop

        ! Second line
        if (iper < 100) then
            write(cestado,'(i3)') iper
            write(funit,'(a,1x,a,a)') &
                trim(bar), trim(adjustl(cestado)), '%'
        else
            write(funit,'(a)') trim(bar)//' FIN'
        end if

        flush(funit)

        first_call = .False.
    end subroutine

end module auxiliary
