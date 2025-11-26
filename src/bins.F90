!> Module with bins routines. Unused yet.
module bins
    use constants, only: wp, cero, uno2, G, pi, tini

    implicit none

    type :: my_bins
        real(wp), allocatable :: edges(:)   ! Bin edges
        real(wp), allocatable :: centers(:) ! Bin centers
        real(wp) :: rmin, rmax              ! Minimum and maximum radii
        real(wp), allocatable :: areas(:)   ! Bin areas
        real(wp), allocatable :: mass(:)    ! Bin mass
        integer(kind=4) :: Nbins                ! Number of bins
        integer(kind=4) :: binning_method       ! Binning method
        real(wp), dimension(2) :: center    ! Center of the disk
    end type my_bins
    
    real(wp), dimension(0:15), parameter :: P2_2k = (/ &  ! P_{2k}^2(0)
    & 1.e0_wp, 0.25e0_wp, 0.140625e0_wp, 0.0976563e0_wp, 0.0747681e0_wp, 0.0605621e0_wp, 0.050889e0_wp, 0.0438788e0_wp, &
    & 0.0385653e0_wp, 0.0343993e0_wp, 0.0310454e0_wp, 0.0282872e0_wp, 0.0259791e0_wp, 0.0240191e0_wp, 0.0223341e0_wp, 0.02087e0_wp &
    &/)

    contains

        ! Initialize bins
        pure subroutine allocate_bins_impl(self, Nbins, rmin, rmax, binning_method)
            implicit none
            type(my_bins), intent(inout) :: self
            integer(kind=4), intent(in) :: Nbins
            real(wp), intent(in) :: rmin, rmax
            integer(kind=4), intent(in) :: binning_method
            
            self%Nbins = Nbins
            self%rmin = rmin
            self%rmax = rmax
            self%binning_method = binning_method
            
            allocate(self%edges(Nbins+1))
            allocate(self%centers(Nbins))
            allocate(self%areas(Nbins))
            allocate(self%mass(Nbins))
        end subroutine allocate_bins_impl
        
        ! Set disk center
        pure subroutine set_center_impl(self, center)
            implicit none
            type(my_bins), intent(inout) :: self
            real(wp), intent(in) :: center(2)
            
            self%center = center
        end subroutine set_center_impl

        ! Set bin edges and centers
        subroutine set_bins_impl(self, particles_dist_sorted)
            implicit none
            type(my_bins), intent(inout) :: self
            real(wp), dimension(:), optional, intent(in) :: particles_dist_sorted  ! Out: Sorted
            integer(kind=4) :: i
            real(wp) :: dr, da

            self%edges(1) = self%rmin
            self%edges(self%Nbins+1) = self%rmax

            select case (self%binning_method)
            case (1)  ! Equal width bins
                dr = (self%rmax - self%rmin) / self%Nbins
                do i = 2, self%Nbins+1
                    self%edges(i) = self%edges(i-1) + dr
                    self%areas(i-1) = pi * (self%edges(i)**2 - self%edges(i-1)**2)
                    self%centers(i-1) = (self%edges(i-1) + self%edges(i)) * uno2
                end do

            case (2)  ! Equal area bins
                da = (self%rmax**2 - self%rmin**2) / self%Nbins
                do i = 2, self%Nbins+1
                    self%edges(i) = sqrt(self%edges(i-1)**2 + da / pi)
                    self%areas(i-1) = da
                    self%centers(i-1) = (self%edges(i-1) + self%edges(i)) * uno2
                end do

            case (3)  ! Equal particle number bins
                if (.not. present(particles_dist_sorted)) then
                    write(*,*) "Error: particles_dist_sorted is required for binning method = 3"
                    stop
                end if
                call set_bins_equal_particles(self, particles_dist_sorted)

            case default
                write(*,*) "Error: Invalid binning method:", self%binning_method
                stop
            end select
        end subroutine set_bins_impl

        ! Helper for equal-particle-number binning
        subroutine set_bins_equal_particles(self, particles_dist_sorted)
            implicit none
            type(my_bins), intent(inout) :: self
            real(wp), dimension(:), intent(in) :: particles_dist_sorted
            integer(kind=4) :: first_part, last_part
            integer(kind=4) :: npart_per_bin, nbins_with_extra
            integer(kind=4) :: i
            
            ! SORT

            ! Error checks for input range
            if (particles_dist_sorted(1) > self%rmax) then
                write(*, *) "Error: rmax", self%rmax, "below innermost particle", particles_dist_sorted
                stop 3
            end if

            if (particles_dist_sorted(size(particles_dist_sorted)) < self%rmin) then
                write(*, *) "Error: rmin", self%rmin, "further than furthest particle", &
                        & particles_dist_sorted(size(particles_dist_sorted))
                stop 3
            end if

            ! Find the range of particles within [rmin, rmax]
            do first_part = 1, size(particles_dist_sorted)
                ! Move to the next part as long as the particle is outside the current range
                if (particles_dist_sorted(first_part) >= self%rmin) exit
            end do
            do last_part = size(particles_dist_sorted), first_part, -1
                ! Move to the previuos part as long as the particle is outside the current range
                if (particles_dist_sorted(last_part) <= self%rmax) exit
            end do

            npart_per_bin = (last_part - first_part) / self%Nbins
            if (npart_per_bin == 0) then
                write(*, *) "Error: Too many bins. Try rescaling Nbins, or modifying rmax, rmin."
                stop 3
            end if

            nbins_with_extra = mod(last_part - first_part, self%Nbins)

            ! Assign edges and areas
            self%edges(1) = particles_dist_sorted(first_part)
            self%edges(self%Nbins + 1) = particles_dist_sorted(last_part)

            do i = 2, nbins_with_extra + 1
                self%edges(i) = particles_dist_sorted(first_part + (npart_per_bin + 1) * (i - 1))
                self%areas(i - 1) = pi * (self%edges(i)**2 - self%edges(i - 1)**2)
                self%centers(i - 1) = uno2 * (self%edges(i - 1) + self%edges(i))
            end do

            do i = nbins_with_extra + 2, self%Nbins
                self%edges(i) = particles_dist_sorted(first_part + nbins_with_extra + npart_per_bin * (i - 1))
                self%areas(i - 1) = pi * (self%edges(i)**2 - self%edges(i - 1)**2)
                self%centers(i - 1) = uno2 * (self%edges(i - 1) + self%edges(i))
            end do
            
            self%areas(self%Nbins) = pi * (self%edges(self%Nbins + 1)**2 - self%edges(self%Nbins)**2)
            self%centers(self%Nbins) = uno2 * (self%edges(self%Nbins) + self%edges(self%Nbins + 1))
        end subroutine set_bins_equal_particles
        
        ! Get each particle bin
        pure subroutine get_particles_bins_impl(self, particles_dist, particles_bins, dist_is_sorted)
            implicit none
            type(my_bins), intent(inout) :: self
            real(wp), dimension(:), intent(in) :: particles_dist
            integer(kind=4), dimension(:), intent(inout) :: particles_bins
            logical, intent(in) :: dist_is_sorted
            integer(kind=4) :: i, bin_i, aux
            
            particles_bins = -1
            if (dist_is_sorted) then
                bin_i = -1
                aux = 1
                do while (bin_i < 1)  ! Define first no null particle bin
                    bin_i = find_bin(self, particles_dist(aux))
                    particles_bins(aux) = bin_i
                    aux = aux + 1
                end do
                do i = aux, size(particles_dist)
                    do while (bin_i <= self%Nbins)
                        if (particles_dist(i) > self%edges(bin_i + 1)) then
                            bin_i = bin_i + 1  ! Separate ifs because impossible edges(Nbins+2)
                        else
                            exit
                        end if
                    end do
                    particles_bins(i) = bin_i
                end do
            else
                do i = 1, size(particles_dist)
                    particles_bins(i) = find_bin(self, particles_dist(i))
                end do
            end if
        end subroutine get_particles_bins_impl

        ! Calculate mass in each bin
        pure subroutine calculate_mass_impl(self, particles_dist, particles_mass, particles_bins)
            implicit none
            type(my_bins), intent(inout) :: self
            real(wp), dimension(:), intent(in) :: particles_dist, particles_mass
            integer(kind=4), dimension(:), intent(in) :: particles_bins
            integer(kind=4) :: i, bin_i

            self%mass = cero
            do i = 1, size(particles_dist)
                bin_i = particles_bins(i)
                if (0 < bin_i .and. bin_i < self%Nbins+1) self%mass(bin_i) = self%mass(bin_i) + particles_mass(i)
            end do
        end subroutine calculate_mass_impl
        
        ! Calculate self-gravity from bins
        pure subroutine calculate_force_impl(self, particle_bin, particle_pos, particle_dist, order_legendre, particle_acc)
            implicit none
            type(my_bins), intent(in) :: self
            integer(kind=4), intent(in) :: particle_bin   ! Particle position (x, y)
            real(wp), intent(in) :: particle_pos(2), particle_dist   ! Particle position (x, y)
            integer(kind=4), intent(in) :: order_legendre
            real(wp), intent(inout) :: particle_acc(2)  ! Acceleration components (ax, ay)
            integer(kind=4) :: i, k
            real(wp) :: Rring, mring, r_ratio, legendre_sum
            
            ! If particle_bin == 0, only the second loop is used
            ! If particle_bin == Nbins + 1, only the first loop is used
            
            ! Loop over all previous bins
            do i = 1, particle_bin - 1
                Rring = self%centers(i)
                mring = self%mass(i)
                
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
            do i = particle_bin + 1, self%Nbins
                Rring = self%centers(i)
                mring = self%mass(i)
                
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
        pure function find_bin(self, dist) result(bin)
            type(my_bins), intent(in) :: self
            real(wp), intent(in) :: dist
            integer(kind=4) :: bin, low, high, mid

            bin = -1
            if (dist < self%edges(1)) then ! Out of bounds
                bin = 0
                return
            else if (dist > self%edges(self%Nbins+1)) then ! Out of bounds
                bin = self%Nbins+1
                return
            else if (self%edges(self%Nbins+1) - dist < tini) then ! If exactly rightmost edge, bin = Nbins
                bin = self%Nbins
                return
            end if

            low = 1
            high = self%Nbins
            do while (low <= high)
                mid = (low + high) / 2
                if (dist >= self%edges(mid) .and. dist < self%edges(mid + 1)) then
                    bin = mid
                    exit
                else if (dist < self%edges(mid)) then
                    high = mid - 1
                else
                    low = mid + 1
                end if
            end do
        end function find_bin
        
        ! Free memory for bins
        pure subroutine free_bins_impl(self)
            implicit none
            type(my_bins), intent(inout) :: self
            if (allocated(self%edges)) deallocate(self%edges)
            if (allocated(self%centers)) deallocate(self%centers)
            if (allocated(self%areas)) deallocate(self%areas)
            if (allocated(self%mass)) deallocate(self%mass)
        end subroutine free_bins_impl
end module bins
