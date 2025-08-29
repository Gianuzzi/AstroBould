!> Module with array, string, ordering, ... routines
module auxiliar
    implicit none
    
    contains


        ! From lower to upper case
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

        ! From lower to upper case
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

        ! Get indices to order a real array
        recursive subroutine quickargsort(a, b, first, last)
            implicit none
            real(kind=8), intent(in) :: a(:)          ! Input array (not modified)
            integer(kind=4), intent(inout) :: b(size(a))    ! Indices array to be sorted
            integer(kind=4), intent(in) :: first, last
            real(kind=8) :: pivot
            integer(kind=4) :: i, j, temp

            ! Pivot selection
            pivot = a(b(floor((first + last) * 0.5d0)))
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

        ! Get indices to order an integer array
        recursive subroutine quickargsort_int(a, b, first, last)
            implicit none
            integer(kind=4), intent(in) :: a(:)    ! Input array (not modified)
            integer(kind=4), intent(inout) :: b(size(a)) ! Indices array to be sorted
            integer(kind=4), intent(in) :: first, last
            integer(kind=4) :: pivot
            integer(kind=4) :: i, j, temp

            ! Pivot selection
            pivot = a(b(floor((first + last) * 0.5d0)))
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
        
        ! Order a real array
        subroutine quicksort(a, first, last)
            implicit none
            real(kind=8), intent(inout) :: a(:)          ! Input array
            integer(kind=4), intent(in) :: first, last   ! First and last indices
            integer(kind=4), allocatable :: b(:)         ! Indices of array to be sorted
            real(kind=8), allocatable :: a_copy(:)       ! Temporary copy of the array
            integer(kind=4) :: n_sort, i
            
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
        
        ! Order an integer array
        subroutine quicksort_int(a, first, last)
            implicit none
            integer(kind=4), intent(inout) :: a(:)       ! Input array
            integer(kind=4), intent(in) :: first, last   ! First and last indices
            integer(kind=4), allocatable :: b(:)         ! Indices of array to be sorted
            integer(kind=4), allocatable :: a_copy(:)    ! Temporary copy of the array
            integer(kind=4) :: n_sort, i

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
        
        ! Reorder a real array from given order
        subroutine reorder(array, order, n)
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

        ! Reorder an integer array from given order
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

        ! Reorder a 2D real array from given order, along 1st axis
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

        ! Combine, order, and remove duplicated values between 2 real arrays
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

        ! Calculate 2D cross product
        function cross2D(a, b) result(res)
            implicit none
            real(kind=8), dimension(2), intent(in) :: a, b
            real(kind=8) :: res

            res = a(1) * b(2) - a(2) * b(1) ! Solo la componente z
        end function cross2D
        

end module auxiliar