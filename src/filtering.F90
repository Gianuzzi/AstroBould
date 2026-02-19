!> Module with filtering methods.
module filtering
    use constants, only: wp, pi, twopi, cero, uno, uno2, dos, myepsilon

    implicit none

    type :: filter_st
        ! filter parameters
        real(wp) :: dt = cero
        real(wp) :: half_width = cero
        real(wp) :: total_dt = cero
        integer(kind=4) :: size = 0
        integer(kind=4) :: half_size = 0
        integer(kind=4) :: n_windows = 0
        integer(kind=4) :: n_samples = 0
        integer(kind=4) :: weight_model = -1
        real(wp) :: omega_pass = cero
        logical :: low_pass = .True.
        real(wp), dimension(:), allocatable :: kernel  ! convolution kernel
        real(wp), dimension(:), allocatable :: tmp_times  ! times to use for filtering
        real(wp), dimension(:, :), allocatable :: tmp_values  ! values to use for filtering
    end type filter_st

contains

    subroutine create_filter(self, cutoff_period, samples_per_cutoff, support_factor, low_pass, model)
        implicit none
        type(filter_st), intent(inout) :: self
        real(wp), intent(in) :: cutoff_period
        integer(kind=4), intent(in) :: samples_per_cutoff, support_factor
        logical, intent(in) :: low_pass
        integer(kind=4), intent(in) :: model
        real(wp) :: t_j
        real(wp) :: val, weight, N
        integer(kind=4) :: j

        !------------------------
        ! Assign parameters
        !------------------------
        if (cutoff_period < myepsilon) then
            write (*, *) "ERROR: Filter dt can not be too low."
            stop 1
        end if

        self%total_dt = cutoff_period*support_factor

        self%n_samples = samples_per_cutoff
        self%n_windows = support_factor
        if (mod(support_factor, 2) == 1) self%n_windows = support_factor + 1  ! even

        self%size = self%n_samples*self%n_windows + 1  ! +1 for the first condition  ! odd
        self%half_size = int((self%size - 1)/2, 4)
        self%dt = real(self%total_dt/(self%size - 1), kind=wp)
        self%half_width = self%total_dt*uno2

        self%omega_pass = twopi/cutoff_period  ! Frequency allowed

        ! Allocate kernel
        allocate (self%kernel(self%size))

        !------------------------
        ! Build sinc low-pass kernel
        ! kernel_j = sin(omega_pass * t_j) / (pi * t_j)
        ! For j = 0, use limit value 1
        !------------------------
        N = real(self%size - 1, kind=wp)
        do j = -self%half_size, self%half_size
            t_j = real(j, kind=wp)*self%dt
            if (abs(t_j) < myepsilon) then
                val = self%omega_pass/pi
            else
                val = sin(self%omega_pass*t_j)/(pi*t_j)
            end if
            modelo:select case(model)
            case (0)  ! boxcar
            weight = uno
            case (1)  ! Hann
            weight = uno2 - &
                   & uno2*cos(twopi*(j + self%half_size)/N)
            case (2)  ! Hamming
            weight = 0.54e0_wp - &
                   & 0.46e0_wp*cos(twopi*(j + self%half_size)/N)
            case (3)  ! Blackman
            weight = 0.42e0_wp - &
                   & uno2*cos(twopi*(j + self%half_size)/N) + &
                   & 0.08e0_wp*cos(dos*twopi*(j + self%half_size)/N)
            case (4)  ! Flat top
            weight = uno - &
                   & 1.93e0_wp*cos(twopi*(j + self%half_size)/N) + &
                   & 1.29e0_wp*cos(dos*twopi*(j + self%half_size)/N) - &
                   & 0.388e0_wp*cos(3.e0_wp*twopi*(j + self%half_size)/N) + &
                   & 0.032e0_wp*cos(4.e0_wp*twopi*(j + self%half_size)/N)
            case default  ! boxcar
            weight = uno
            end select modelo
            self%kernel(j + self%half_size + 1) = val*weight
        end do

        !------------------------
        ! Normalize kernel to unity sum
        !------------------------
        self%kernel = self%kernel/sum(self%kernel)

        if (.not. low_pass) then
            !------------------------
            ! Reverse to obtain high-pass kernel
            !------------------------
            self%kernel = -self%kernel
            self%kernel(self%half_size + 1) = uno + self%kernel(self%half_size + 1)
            self%low_pass = .False.
        else
            self%low_pass = .True.
        end if

    end subroutine create_filter

    pure subroutine allocate_filter(self, data_size)
        implicit none
        type(filter_st), intent(inout) :: self
        integer(kind=4), intent(in) :: data_size

        allocate (self%tmp_values(data_size, self%size))
        allocate (self%tmp_times(self%size))

    end subroutine allocate_filter

    subroutine setup_filter(self, cutoff_period, samples_per_cutoff, support_factor, low_pass, model, data_size)
        implicit none
        type(filter_st), intent(inout) :: self
        real(wp), intent(in) :: cutoff_period
        integer(kind=4), intent(in) :: samples_per_cutoff, support_factor, model, data_size
        logical, intent(in) :: low_pass

        call create_filter(self, cutoff_period, samples_per_cutoff, support_factor, low_pass, model)

        call allocate_filter(self, data_size)

    end subroutine setup_filter

    pure subroutine store_to_filter(self, time, state_vector, data_size, filter_index)
        implicit none
        type(filter_st), intent(inout) :: self
        real(wp), intent(in) :: time
        real(wp), dimension(:), intent(in) :: state_vector
        integer(kind=4), intent(in) :: data_size, filter_index

        self%tmp_values(1:data_size, filter_index) = state_vector(1:data_size)
        self%tmp_times(filter_index) = time

    end subroutine store_to_filter

    pure subroutine free_filter(self)
        implicit none
        type(filter_st), intent(inout) :: self

        if (allocated(self%kernel)) deallocate (self%kernel)
        if (allocated(self%tmp_values)) deallocate (self%tmp_values)
        if (allocated(self%tmp_times)) deallocate (self%tmp_times)
    end subroutine free_filter

end module filtering
