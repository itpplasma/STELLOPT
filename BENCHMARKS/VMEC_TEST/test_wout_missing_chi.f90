program test_wout_missing_chi
    ! A WOUT file without chi must yield the poloidal flux VMEC would write.
    use stel_kinds, only: rprec
    use read_wout_mod, only: read_wout_file, read_wout_deallocate, chi
    implicit none
    character(len=512) :: full, stripped
    real(rprec), allocatable :: expected(:)
    integer :: ierr

    call get_command_argument(1, full)
    call get_command_argument(2, stripped)
    call read_wout_file(trim(full), ierr)
    if (ierr /= 0) error stop 'cannot read WOUT with chi'
    expected = chi
    call read_wout_deallocate
    call read_wout_file(trim(stripped), ierr)
    if (ierr /= 0) error stop 'cannot read WOUT without chi'
    print '(a,es10.3)', 'max |chi - chi_VMEC|/max|chi_VMEC| = ', &
        maxval(abs(chi - expected))/maxval(abs(expected))
    if (maxval(abs(chi - expected)) > 1.0e-13_rprec*maxval(abs(expected))) &
        error stop 'reconstructed chi differs from VMEC'
    print '(a)', 'PASS: missing chi reconstruction'
end program test_wout_missing_chi
