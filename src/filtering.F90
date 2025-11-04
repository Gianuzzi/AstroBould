!> Module with filtering methods.
module filtering
    use constants, only: pi, twopi, cero, uno, uno2, dos, tini
    
    implicit none

    type :: filter_st
        ! filter parameters
        real(kind=8) :: dt = cero
        real(kind=8) :: half_width = cero
        integer(kind=4) :: size = 0
        integer(kind=4) :: half_size = 0
        integer(kind=4) :: n_windows = 0
        integer(kind=4) :: n_samples = 0
        integer(kind=4) :: weight_model = -1
        real(kind=8) :: omega_pass = cero
        logical :: low_pass = .True.
        real(kind=8), dimension(:), allocatable :: kernel  ! convolution kernel
        real(kind=8), dimension(:), allocatable :: tmp_times  ! times to use for filtering
        real(kind=8), dimension(:,:), allocatable :: tmp_values  ! values to use for filtering
    end type filter_st
    
contains

    subroutine create_filter(self, dt_cutoff, oversample, window_factor, low_pass, model)
        implicit none
        type(filter_st), intent(inout) :: self
        real(kind=8), intent(in) :: dt_cutoff
        integer(kind=4), intent(in) :: oversample, window_factor
        logical, intent(in) :: low_pass
        integer(kind=4), intent(in) :: model
        real(kind=8) :: filter_long_dt
        real(kind=8) :: t_j
        real(kind=8) :: val, weight, N
        integer(kind=4) :: j

        !------------------------
        ! Assign parameters
        !------------------------
        if (dt_cutoff < tini) then
            write(*,*) "ERROR: Filter dt can not be too low."
            stop 1
        end if

        filter_long_dt = dt_cutoff * window_factor

        self%n_samples = oversample
        self%n_windows = window_factor
        if (mod(window_factor, 2) == 1) self%n_windows = window_factor + 1  ! even

        self%size = self%n_samples * self%n_windows + 1  ! +1 for the first condition  ! odd
        self%half_size = int((self%size - 1) / 2, 4)
        self%dt = dble(filter_long_dt / (self%size - 1))
        self%half_width = filter_long_dt * uno2
        
        self%omega_pass = twopi / dt_cutoff  ! Frequency allowed

        ! Allocate kernel
        allocate(self%kernel(self%size))

        !------------------------
        ! Build sinc low-pass kernel
        ! kernel_j = sin(omega_pass * t_j) / (pi * t_j)
        ! For j = 0, use limit value 1
        !------------------------
        N = dble(self%size - 1)
        do j = -self%half_size, self%half_size
            t_j = dble(j) * self%dt
            if (abs(t_j) < tini) then
                val = self%omega_pass / pi
            else
                val = sin(self%omega_pass * t_j) / (pi * t_j)
            end if
            modelo: select case (model)
                case (0)  ! boxcar
                    weight = uno
                case (1)  ! Hann
                    weight = uno2 - &
                           & uno2 * cos(twopi * (j + self%half_size) / N)
                case (2)  ! Hamming
                    weight = 0.54d0 - &
                           & 0.46d0 * cos(twopi * (j + self%half_size) / N)
                case (3)  ! Blackman
                    weight = 0.42d0 - &
                           & uno2 * cos(twopi * (j + self%half_size) / N) + &
                           & 0.08d0 * cos(dos * twopi * (j + self%half_size) / N)
                case (4)  ! Flat top
                    weight = uno - &
                           & 1.93d0 * cos(twopi * (j + self%half_size) / N) + &
                           & 1.29d0 * cos(dos * twopi * (j + self%half_size) / N) - &
                           & 0.388d0 * cos(3.d0 * twopi * (j + self%half_size) / N) + &
                           & 0.032d0 * cos(4.d0 * twopi * (j + self%half_size) / N)
                case default  ! boxcar
                    weight = uno
            end select modelo
            self%kernel(j + self%half_size + 1) = val * weight
        end do

        !------------------------
        ! Normalize kernel to unity sum
        !------------------------
        self%kernel = self%kernel / sum(self%kernel)

        if (.not. low_pass) then
            !------------------------
            ! Reverse to obtain high-pass kernel
            !------------------------
            self%kernel = - self%kernel
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

        allocate(self%tmp_values(data_size, self%size))
        allocate(self%tmp_times(self%size))
        
    end subroutine allocate_filter

    subroutine setup_filter(self, dt_cutoff, oversample, window_factor, low_pass, model, data_size)
        implicit none
        type(filter_st), intent(inout) :: self
        real(kind=8), intent(in) :: dt_cutoff
        integer(kind=4), intent(in) :: oversample, window_factor, model, data_size
        logical, intent(in) :: low_pass

        call create_filter(self, dt_cutoff, oversample, window_factor, low_pass, model)

        call allocate_filter(self, data_size)   
        
    end subroutine setup_filter

    pure subroutine store_to_filter(self, time, state_vector, data_size, filter_index)
        implicit none
        type(filter_st), intent(inout) :: self
        real(kind=8), intent(in) :: time
        real(kind=8), dimension(:), intent(in) :: state_vector
        integer(kind=4), intent(in) :: data_size, filter_index

        self%tmp_values(1:data_size, filter_index) = state_vector(1:data_size)    
        self%tmp_times(filter_index) = time
        
    end subroutine store_to_filter
    
    pure subroutine free_filter(self)
        implicit none
        type(filter_st), intent(inout) :: self

        if (allocated(self%kernel)) deallocate(self%kernel)
        if (allocated(self%tmp_values)) deallocate(self%tmp_values)
        if (allocated(self%tmp_times)) deallocate(self%tmp_times)
    end subroutine free_filter

    
end module filtering