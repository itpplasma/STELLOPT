program test_profile_limits
    ! The ion temperature is clamped at the last point of its own radial
    ! grid, which may differ from the density grid.
    use beams3d_input_mod, only: init_beams3d_input, read_beams3d_input
    use beams3d_globals, only: s_max_ti, s_max_ne
    implicit none
    integer :: iunit, istat
    open(newunit=iunit, file='input.profile_limits', status='replace')
    write(iunit, '(A)') '&BEAMS3D_INPUT'
    write(iunit, '(A)') '  NE_AUX_S = 0.0 0.5 1.0'
    write(iunit, '(A)') '  NE_AUX_F = 1.0E19 0.5E19 0.1E19'
    write(iunit, '(A)') '  TI_AUX_S = 0.0 0.4 0.8 1.2'
    write(iunit, '(A)') '  TI_AUX_F = 1.0E3 0.8E3 0.4E3 0.1E3'
    write(iunit, '(A)') '/'
    close(iunit)
    call init_beams3d_input
    call read_beams3d_input('input.profile_limits', istat)
    if (istat /= 0) error stop 2
    if (abs(s_max_ti - 1.2d0) > 1.0d-12 .or. abs(s_max_ne - 1.0d0) > 1.0d-12) then
        print '(A,2ES16.8)', 'FAIL s_max_ti, s_max_ne:', s_max_ti, s_max_ne
        error stop 1
    end if
    print '(A)', 'PASS test_profile_limits'
end program test_profile_limits
