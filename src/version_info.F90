module version_info
  implicit none
  private
  public :: print_version

  !-----------------------------
  ! Package name
#ifdef PACKAGE_VERSION
  character(len=*), parameter :: package_name = PACKAGE_NAME
#else
  character(len=*), parameter :: package_name = "UNKNOWN"
#endif

  !-----------------------------
  ! Version string
#ifdef PACKAGE_VERSION
  character(len=*), parameter :: package_version = PACKAGE_VERSION
#else
  character(len=*), parameter :: package_version = "unknown"
#endif

  !-----------------------------
  ! Git branch
#ifdef PACKAGE_BRANCH
  character(len=*), parameter :: package_branch = PACKAGE_BRANCH
#else
  character(len=*), parameter :: package_branch = "n/a"
#endif

  !-----------------------------
  ! Commit date
#ifdef PACKAGE_DATE
  character(len=*), parameter :: package_date = PACKAGE_DATE
#else
  character(len=*), parameter :: package_date = "unknown"
#endif

  !-----------------------------
  ! Debug marker
#ifdef DEBUG
  logical, parameter :: is_debug = .true.
#else
  logical, parameter :: is_debug = .false.
#endif

contains

  subroutine print_version(unit_file, program_name)
    implicit none
    integer(kind=4), intent(in) :: unit_file
    character(len=*), intent(in), optional :: program_name
    character(len=50) :: my_package_name

    if ((trim(package_name) .eq. "UNKNOWN") .and. present(program_name)) then
        my_package_name = program_name
    else
        my_package_name = package_name
    end if

    write(unit_file,'(A)', advance='no') &
      trim(my_package_name)//"  "//trim(package_version)

    if (package_branch /= "n/a") then
      write(unit_file,'(A)', advance='no') &
        " ["//trim(package_branch)//"]"
    end if

    if (package_date /= "unknown") then
      write(unit_file,'(A)', advance='no') &
        " ("//trim(package_date)//")"
    end if

    if (is_debug) then
      write(unit_file,'(A)', advance='no') "  DEBUG"
    end if

    write(unit_file,*)

  end subroutine print_version

end module version_info
