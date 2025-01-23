module bins
    use constants, only: pi, tini, uno2, G, cero
    implicit none
    type :: my_bins
        real(kind=8), allocatable :: edges(:)   ! Bin edges
        real(kind=8), allocatable :: centers(:) ! Bin centers
        real(kind=8) :: rmin, rmax              ! Minimum and maximum radii
        real(kind=8), allocatable :: areas(:)   ! Bin areas
        real(kind=8), allocatable :: mass(:)    ! Bin mass
        integer(kind=4) :: Nbins                ! Number of bins
        integer(kind=4) :: binning_method       ! Binning method
        real(kind=8), dimension(2) :: center    ! Center of the disk

        contains
            procedure :: allocate_bins => allocate_bins_impl
            procedure :: set_center => set_center_impl
            procedure :: set_bins => set_bins_impl
            procedure :: get_particles_bins => get_particles_bins_impl
            procedure :: calculate_mass => calculate_mass_impl
            procedure :: calculate_force => calculate_force_impl
            procedure :: free => free_impl
            procedure :: duplicate => duplicate_impl
    end type my_bins
    real(kind=8), dimension(0:15), parameter :: P2_2k = (/ &  ! P_{2k}^2(0)
    & 1.d0, 0.25d0, 0.140625d0, 0.0976563d0, 0.0747681d0, 0.0605621d0, 0.050889d0, 0.0438788d0, &
    & 0.0385653d0, 0.0343993d0, 0.0310454d0, 0.0282872d0, 0.0259791d0, 0.0240191d0, 0.0223341d0, 0.02087d0 &
    &/)

contains

    ! Initialize bins
    subroutine allocate_bins_impl(this, Nbins, rmin, rmax, binning_method)
        implicit none
        class(my_bins), intent(inout) :: this
        integer(kind=4), intent(in) :: Nbins
        real(kind=8), intent(in) :: rmin, rmax
        integer(kind=4), intent(in) :: binning_method

        this%Nbins = Nbins
        this%rmin = rmin
        this%rmax = rmax
        this%binning_method = binning_method
        
        allocate(this%edges(Nbins+1))
        allocate(this%centers(Nbins))
        allocate(this%areas(Nbins))
        allocate(this%mass(Nbins))
    end subroutine allocate_bins_impl
    
    ! Set disk center
    subroutine set_center_impl(this, center)
        implicit none
        class(my_bins), intent(inout) :: this
        real(kind=8), intent(in) :: center(2)
        
        this%center = center
    end subroutine set_center_impl

    ! Set bin edges and centers
    subroutine set_bins_impl(this, particles_dist_sorted)
        implicit none
        class(my_bins), intent(inout) :: this
        real(kind=8), dimension(:), optional, intent(in) :: particles_dist_sorted  ! Out: Sorted
        integer(kind=4) :: i
        real(kind=8) :: dr, da

        this%edges(1) = this%rmin
        this%edges(this%Nbins+1) = this%rmax

        select case (this%binning_method)
        case (1)  ! Equal width bins
            dr = (this%rmax - this%rmin) / this%Nbins
            do i = 2, this%Nbins+1
                this%edges(i) = this%edges(i-1) + dr
                this%areas(i-1) = pi * (this%edges(i)**2 - this%edges(i-1)**2)
                this%centers(i-1) = (this%edges(i-1) + this%edges(i)) * uno2
            end do

        case (2)  ! Equal area bins
            da = (this%rmax**2 - this%rmin**2) / this%Nbins
            do i = 2, this%Nbins+1
                this%edges(i) = sqrt(this%edges(i-1)**2 + da / pi)
                this%areas(i-1) = da
                this%centers(i-1) = (this%edges(i-1) + this%edges(i)) * uno2
            end do

        case (3)  ! Equal particle number bins
            if (.not. present(particles_dist_sorted)) then
                write(*,*) "Error: particles_dist_sorted is required for binning method = 3"
                stop
            end if
            call set_bins_equal_particles(this, particles_dist_sorted)

        case default
            write(*,*) "Error: Invalid binning method:", this%binning_method
            stop
        end select
    end subroutine set_bins_impl

    ! Helper for equal-particle-number binning
    subroutine set_bins_equal_particles(this, particles_dist_sorted)
        implicit none
        class(my_bins), intent(inout) :: this
        real(kind=8), dimension(:), intent(in) :: particles_dist_sorted
        integer(kind=4) :: first_part, last_part
        integer(kind=4) :: npart_per_bin, nbins_with_extra
        integer(kind=4) :: i
        
        ! SORT

        ! Error checks for input range
        if (particles_dist_sorted(1) > this%rmax) then
            write(*, *) "Error: rmax", this%rmax, "below innermost particle", particles_dist_sorted
            stop 3
        end if

        if (particles_dist_sorted(size(particles_dist_sorted)) < this%rmin) then
            write(*, *) "Error: rmin", this%rmin, "further than furthest particle", &
                    & particles_dist_sorted(size(particles_dist_sorted))
            stop 3
        end if

        ! Find the range of particles within [rmin, rmax]
        do first_part = 1, size(particles_dist_sorted)
            ! Move to the next part as long as the particle is outside the current range
            if (particles_dist_sorted(first_part) >= this%rmin) exit
        end do
        do last_part = size(particles_dist_sorted), first_part, -1
            ! Move to the previuos part as long as the particle is outside the current range
            if (particles_dist_sorted(last_part) <= this%rmax) exit
        end do

        npart_per_bin = (last_part - first_part) / this%Nbins
        if (npart_per_bin == 0) then
            write(*, *) "Error: Too many bins. Try rescaling Nbins, or modifying rmax, rmin."
            stop 3
        end if

        nbins_with_extra = mod(last_part - first_part, this%Nbins)

        ! Assign edges and areas
        this%edges(1) = particles_dist_sorted(first_part)
        this%edges(this%Nbins + 1) = particles_dist_sorted(last_part)

        do i = 2, nbins_with_extra + 1
            this%edges(i) = particles_dist_sorted(first_part + (npart_per_bin + 1) * (i - 1))
            this%areas(i - 1) = pi * (this%edges(i)**2 - this%edges(i - 1)**2)
            this%centers(i - 1) = uno2 * (this%edges(i - 1) + this%edges(i))
        end do

        do i = nbins_with_extra + 2, this%Nbins
            this%edges(i) = particles_dist_sorted(first_part + nbins_with_extra + npart_per_bin * (i - 1))
            this%areas(i - 1) = pi * (this%edges(i)**2 - this%edges(i - 1)**2)
            this%centers(i - 1) = uno2 * (this%edges(i - 1) + this%edges(i))
        end do
        
        this%areas(this%Nbins) = pi * (this%edges(this%Nbins + 1)**2 - this%edges(this%Nbins)**2)
        this%centers(this%Nbins) = uno2 * (this%edges(this%Nbins) + this%edges(this%Nbins + 1))
    end subroutine set_bins_equal_particles
    
    ! Get each particle bin
    subroutine get_particles_bins_impl(this, particles_dist, particles_bins, dist_is_sorted)
        implicit none
        class(my_bins), intent(inout) :: this
        real(kind=8), dimension(:), intent(in) :: particles_dist
        integer(kind=4), dimension(:), intent(inout) :: particles_bins
        logical, intent(in) :: dist_is_sorted
        integer(kind=4) :: i, bin_i, aux
        
        particles_bins = -1
        if (dist_is_sorted) then
            bin_i = -1
            aux = 1
            do while (bin_i < 1)  ! Define first no null particle bin
                bin_i = find_bin(this, particles_dist(aux))
                particles_bins(aux) = bin_i
                aux = aux + 1
            end do
            do i = aux, size(particles_dist)
                do while (bin_i <= this%Nbins)
                    if (particles_dist(i) > this%edges(bin_i + 1)) then
                        bin_i = bin_i + 1  ! Separate ifs because impossible edges(Nbins+2)
                    else
                        exit
                    end if
                end do
                particles_bins(i) = bin_i
            end do
        else
            do i = 1, size(particles_dist)
                particles_bins(i) = find_bin(this, particles_dist(i))
            end do
        end if
    end subroutine get_particles_bins_impl

    ! Calculate mass in each bin
    subroutine calculate_mass_impl(this, particles_dist, particles_mass, particles_bins)
        implicit none
        class(my_bins), intent(inout) :: this
        real(kind=8), dimension(:), intent(in) :: particles_dist, particles_mass
        integer(kind=4), dimension(:), intent(in) :: particles_bins
        integer(kind=4) :: i, bin_i

        this%mass = cero
        do i = 1, size(particles_dist)
            bin_i = particles_bins(i)
            if (0 < bin_i .and. bin_i < this%Nbins+1) this%mass(bin_i) = this%mass(bin_i) + particles_mass(i)
        end do
    end subroutine calculate_mass_impl
    
    ! Calculate self-gravity from bins
    subroutine calculate_force_impl(this, particle_bin, particle_pos, particle_dist, order_legendre, particle_acc)
        implicit none
        class(my_bins), intent(in) :: this
        integer(kind=4), intent(in) :: particle_bin   ! Particle position (x, y)
        real(kind=8), intent(in) :: particle_pos(2), particle_dist   ! Particle position (x, y)
        integer(kind=4), intent(in) :: order_legendre
        real(kind=8), intent(inout) :: particle_acc(2)  ! Acceleration components (ax, ay)
        integer(kind=4) :: i, k
        real(kind=8) :: Rring, mring, r_ratio, legendre_sum
        
        ! If particle_bin == 0, only the second loop is used
        ! If particle_bin == Nbins + 1, only the first loop is used
        
        ! Loop over all previous bins
        do i = 1, particle_bin - 1
            Rring = this%centers(i)
            mring = this%mass(i)
            
            r_ratio = Rring / particle_dist

            ! Calculate Legendre polynomial sum up to ORDER
            legendre_sum = cero
            do k = 0, order_legendre
                legendre_sum = legendre_sum + (2 * k + 1) * P2_2k(k) * (r_ratio**(2 * k))
            end do

            ! Acceleration components
            particle_acc(1) = particle_acc(1) + G * mring / particle_dist**3 * legendre_sum * particle_pos(1)
            particle_acc(2) = particle_acc(2) + G * mring / particle_dist**3 * legendre_sum * particle_pos(2)
        end do
        ! Loop over all next bins
        do i = particle_bin + 1, this%Nbins
            Rring = this%centers(i)
            mring = this%mass(i)
            
            r_ratio = particle_dist / Rring

            ! Calculate Legendre polynomial sum up to ORDER
            legendre_sum = cero
            do k = 0, order_legendre
                legendre_sum = legendre_sum + (2 * k + 1) * P2_2k(k) * (r_ratio**(2 * k))
            end do

            ! Acceleration components
            particle_acc(1) = particle_acc(1) + G * mring / particle_dist**3 * legendre_sum * particle_pos(1)
            particle_acc(2) = particle_acc(2) + G * mring / particle_dist**3 * legendre_sum * particle_pos(2)
        end do
    end subroutine calculate_force_impl

    ! Helper function: Find bin index
    function find_bin(this, dist) result(bin)
        class(my_bins), intent(in) :: this
        real(kind=8), intent(in) :: dist
        integer(kind=4) :: bin, low, high, mid

        bin = -1
        if (dist < this%edges(1)) then ! Out of bounds
            bin = 0
            return
        else if (dist > this%edges(this%Nbins+1)) then ! Out of bounds
            bin = this%Nbins+1
            return
        else if (this%edges(this%Nbins+1) - dist < tini) then ! If exactly rightmost edge, bin = Nbins
            bin = this%Nbins
            return
        end if

        low = 1
        high = this%Nbins
        do while (low <= high)
            mid = (low + high) / 2
            if (dist >= this%edges(mid) .and. dist < this%edges(mid + 1)) then
                bin = mid
                exit
            else if (dist < this%edges(mid)) then
                high = mid - 1
            else
                low = mid + 1
            end if
        end do
    end function find_bin
    
    ! Duplicate the struct
    subroutine duplicate_impl(source, dest)
        implicit none
        class(my_bins), intent(in) :: source
        class(my_bins), intent(inout) :: dest

        ! Allocate destination bins if not already allocated
        if (.not. allocated(dest%edges)) allocate(dest%edges(size(source%edges)))
        if (.not. allocated(dest%centers)) allocate(dest%centers(size(source%centers)))
        if (.not. allocated(dest%areas)) allocate(dest%areas(size(source%areas)))
        if (.not. allocated(dest%mass)) allocate(dest%mass(size(source%mass)))

        ! Copy scalar values
        dest%rmin = source%rmin
        dest%rmax = source%rmax
        dest%Nbins = source%Nbins
        dest%binning_method = source%binning_method
        dest%center = source%center

        ! Copy array values
        dest%edges = source%edges
        dest%centers = source%centers
        dest%areas = source%areas
        dest%mass = source%mass
    end subroutine duplicate_impl
    
    ! Free memory for bins
    subroutine free_impl(this)
        implicit none
        class(my_bins), intent(inout) :: this

        if (allocated(this%edges)) deallocate(this%edges)
        if (allocated(this%centers)) deallocate(this%centers)
        if (allocated(this%areas)) deallocate(this%areas)
        if (allocated(this%mass)) deallocate(this%mass)
    end subroutine free_impl
end module bins
