!> Module with filtering methods
module filter
    use constants, only: pi, twopi, cero, uno, uno2, dos, tini
    implicit none
    ! filter parameters
    real(kind=8) :: filt_dt, filt_half_width_dt
    integer(kind=4) :: last_idx_no_filter, first_idx_yes_filter
    integer(kind=4) :: filt_size, filt_half_size
    real(kind=8), dimension(:), allocatable :: filt_kernel  ! convolution kernel
    real(kind=8), dimension(:), allocatable :: elem_filtered    ! Filtered_data
    real(kind=8), dimension(:), allocatable :: y_pre_filter  ! Data pre-filtering
    real(kind=8), dimension(:,:), allocatable :: filt_tmp_values  ! values to use for filtering
    real(kind=8), dimension(:), allocatable :: filt_tmp_times  ! times to use for filtering
    
contains

    subroutine create_filter(dt_cutoff, oversample, window_factor, low_pass)
        implicit none
        real(kind=8), intent(in) :: dt_cutoff
        integer(kind=4), intent(in) :: oversample, window_factor
        logical, intent(in) :: low_pass
        real(kind=8) :: filter_long_dt
        integer(kind=4) :: effective_nwindows
        real(kind=8) :: omega_pass, t_j
        real(kind=8) :: norm, val, weight
        integer(kind=4) :: j

        !------------------------
        ! Assign parameters
        !------------------------
        if (dt_cutoff < tini) then
            write(*,*) "ERROR: Filter dt can not be too low."
            stop 1
        end if

        filter_long_dt = dt_cutoff * window_factor

        effective_nwindows = window_factor
        if (mod(window_factor, 2) .eq. 1) effective_nwindows = window_factor + 1  ! even

        filt_size = oversample * effective_nwindows + 1  ! +1 for the first condition  ! odd
        filt_half_size = int((filt_size - 1) / 2, 4)
        filt_dt = dble(filter_long_dt / (filt_size - 1))
        filt_half_width_dt = filter_long_dt * uno2
        
        omega_pass = twopi / dt_cutoff  ! Frequency allowed

        ! Allocate kernel
        allocate(filt_kernel(filt_size))

        !------------------------
        ! Build sinc low-pass kernel
        ! kernel_j = sin(omega_pass * t_j) / (pi * t_j)
        ! For j = 0, use limit value 1
        !------------------------
        do j = -filt_half_size, filt_half_size
            t_j = dble(j) * filt_dt
            if (abs(t_j) < tini) then
                val = omega_pass / pi
            else
                val = sin(omega_pass * t_j) / (pi * t_j)
            end if
            ! Hann window
            weight = uno2 - uno2 * cos(twopi * (j + filt_half_size) / dble(filt_size - 1))
            !weight = uno
            filt_kernel(j + filt_half_size + 1) = val * weight
        end do

        !------------------------
        ! Normalize kernel to unity sum
        !------------------------
        norm = sum(filt_kernel)
        filt_kernel = filt_kernel / norm

        if (.not. low_pass) then
            !------------------------
            ! Reverse to obtain high-pass kernel
            !------------------------
            filt_kernel = - filt_kernel
            filt_kernel(filt_half_size + 1) = uno + filt_kernel(filt_half_size + 1)
        end if

    end subroutine create_filter

    subroutine allocate_filter(filter_size, data_size)
        implicit none
        integer(kind=4), intent(in) :: filter_size, data_size

        allocate(filt_tmp_values(data_size, filter_size))
        allocate(filt_tmp_times(filter_size))
        allocate(elem_filtered(data_size))
        allocate(y_pre_filter(data_size))
        
    end subroutine allocate_filter

    subroutine setup_filter(dt_cutoff, oversample, window_factor, low_pass, data_size)
        implicit none
        real(kind=8), intent(in) :: dt_cutoff
        integer(kind=4), intent(in) :: oversample, window_factor, data_size
        logical, intent(in) :: low_pass

        call create_filter(dt_cutoff, oversample, window_factor, low_pass)

        call allocate_filter(filt_size, data_size)   
        
    end subroutine setup_filter

    subroutine store_to_filter(time, state_vector, data_size, filter_index)
        implicit none
        real(kind=8), intent(in) :: time
        real(kind=8), dimension(:), intent(in) :: state_vector
        integer(kind=4), intent(in) :: data_size, filter_index

        filt_tmp_values(1:data_size, filter_index) = state_vector(1:data_size)    
        filt_tmp_times(filter_index) = time
        
    end subroutine store_to_filter

    subroutine simple_apply_filter(data_size, filtered)
        implicit none
        integer(kind=4), intent(in) :: data_size
        real(kind=8), dimension(:), intent(out) :: filtered
        real(kind=8) :: cos_th, sin_th
        integer(kind=4) :: j

        ! init to 0
        filtered(1:data_size) = cero
        cos_th = cero
        sin_th = cero

        do j = 1, filt_size
            !! Assumes the first index contains and angle (and the only one)
            cos_th = cos_th + cos(filt_tmp_values(1, j)) * filt_kernel(j)
            sin_th = sin_th + sin(filt_tmp_values(1, j)) * filt_kernel(j)
            filtered(2:data_size) = filtered(2:data_size) + filt_tmp_values(2:data_size, j) * filt_kernel(j)
        end do

        ! Angle filtered
        filtered(1) = modulo(atan2(sin_th, cos_th), twopi)
        
    end subroutine simple_apply_filter
    
    subroutine free_filter_arrays()
        implicit none

        if (allocated(filt_kernel)) deallocate(filt_kernel)
        if (allocated(elem_filtered)) deallocate(elem_filtered)
        if (allocated(y_pre_filter)) deallocate(y_pre_filter)
        if (allocated(filt_tmp_values)) deallocate(filt_tmp_values)
        if (allocated(filt_tmp_times)) deallocate(filt_tmp_times)
        
    end subroutine free_filter_arrays

    
end module filter