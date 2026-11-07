program test_bongrid_edge_gradient
    ! bongrid must give the one-sided Te and Ti gradients at the outermost
    ! surface, (T(irup) - T(irup-1))/drho, each from its own profile.
    use parambs
    use read_boozer_mod, only: mboz_b, nboz_b, lasym_b
    implicit none
    integer, parameter :: ians = 21
    integer :: ihere, i
    real(rprec) :: expected_te, expected_ti

    irup = 4
    irdim = irup + 1
    mboz_b = 1
    nboz_b = 1
    lasym_b = .false.
    lasym_bootsj = .false.
    mbuse = 1
    nbuse = 1
    nthetah = 8
    nzetah = 8
    call allocate_radial
    amnfit(:, 0, 0) = 1
    amnfit(:, 1, 0) = 0.1_dp
    d_rho(:irup) = 1.0_dp/irup
    do i = 1, irup
        rhoar(i) = (i - 0.5_dp)/irup
    end do
    tempe1(:irup) = 3 - 2*rhoar
    tempi1(:irup) = 1 - 0.5_dp*rhoar**2
    dense(:irup) = 1
    open (ians, status='scratch')

    ihere = 0
    call bongrid(irup - 1, ians, ihere)
    call bongrid(irup, ians, ihere)
    expected_te = (tempe1(irup) - tempe1(irup - 1))/drho
    expected_ti = (tempi1(irup) - tempi1(irup - 1))/drho
    if (abs(temperho1 - expected_te) > 1e-12_dp*abs(expected_te)) then
        print *, 'FAIL: Te gradient at the outermost surface', temperho1, expected_te
        stop 1
    end if
    if (abs(tempirho1 - expected_ti) > 1e-12_dp*abs(expected_ti)) then
        print *, 'FAIL: Ti gradient at the outermost surface', tempirho1, expected_ti
        stop 1
    end if
    print '(a)', 'PASS: BOOTSJ edge temperature gradients'
end program test_bongrid_edge_gradient
