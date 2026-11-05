program test_coulomb_rates
    ! The BEAMS3D collision coefficients against the NRL Plasma Formulary
    ! (2019, p. 31) energy-loss and transverse rates, written with the
    ! Maxwell integral psi(x) instead of the Chandrasekhar G, and the
    ! zero-flux condition of a Maxwellian at equal temperatures. The NRL
    ! rates assume a constant logarithm; with ln(Lambda(v)) the zero-flux
    ! drift carries the extra term v D dln(Lambda)/dv.
    use beams3d_physics_mod, only: coulomb_rates
    implicit none
    integer, parameter :: dp = kind(1.0d0), ns = 3
    real(dp), parameter :: e = 1.602176634d-19, eps0 = 8.8541878128d-12
    real(dp), parameter :: hbar = 1.054571817d-34, pi = acos(-1.0_dp)
    real(dp), parameter :: amu = 1.66053906660d-27, me = 9.1093837015d-31
    real(dp), parameter :: ma = 2.01410177812_dp*amu, qa = e, temp = 1.0d3
    real(dp) :: nb(ns), tb(ns), mb(ns), zb(ns), v, h, vth
    real(dp) :: dpar, drift, dhh, dpe, dre, dplus, dminus, expected
    real(dp) :: loss, defl, x, nu0, lam, dlam, logflux
    integer :: k, b
    real(dp), parameter :: speeds(4) = [0.3_dp, 1.0_dp, 3.0_dp, 30.0_dp]

    ! Electrons, H+ and C6+ at 1 keV; quasi-neutral.
    mb = [me, 1.00727646688_dp*amu, 12.0_dp*amu]
    zb = [-1.0_dp, 1.0_dp, 6.0_dp]
    nb = [0.0_dp, 4.0d19, 1.0d18]
    nb(1) = nb(2) + 6.0_dp*nb(3)
    tb = temp
    vth = sqrt(2*temp*e/mb(2))
    do k = 1, 4
        v = speeds(k)*vth
        call coulomb_rates(ma, qa, v, ns, nb, tb, mb, zb, dpar, drift, dhh, dpe, dre)
        loss = 0.0_dp
        defl = 0.0_dp
        logflux = 0.0_dp
        do b = 1, ns
            x = mb(b)*v*v/(2*temp*e)
            lam = coulomb_log(b, v)
            nu0 = nb(b)*(qa*zb(b)*e)**2*lam/(4*pi*eps0**2*ma**2*v**3)
            loss = loss + 2*(ma/mb(b)*psi(x) - dpsi(x))*nu0
            defl = defl + 2*((1 - 0.5_dp/x)*psi(x) + dpsi(x))*nu0
            dlam = (coulomb_log(b, v*(1 + 1.0d-6)) - coulomb_log(b, v*(1 - 1.0d-6)))/(2.0d-6*v)
            logflux = logflux + v*(nu0*v*v*psi(x)/(2*x))*dlam/lam
        end do
        ! Ito: d<m v^2/2>/dt = m (v drift + D) = -nu_eps m v^2/2.
        expected = -loss*0.5_dp*v*v + logflux
        if (abs((v*drift + dpar) - expected) > 1.0d-10*abs(expected)) &
            error stop 'Energy-loss rate differs from the NRL rate'
        ! NRL's transverse rate nu_perp is 2 nu_D; dhh is the Lorentz
        ! coefficient nu_D/2.
        if (abs(dhh - 0.25_dp*defl) > 1.0d-10*0.25_dp*defl) &
            error stop 'Pitch scattering differs from the NRL transverse rate'
        h = 1.0d-5*v
        call coulomb_rates(ma, qa, v + h, ns, nb, tb, mb, zb, dplus, drift, dhh, dpe, dre)
        call coulomb_rates(ma, qa, v - h, ns, nb, tb, mb, zb, dminus, drift, dhh, dpe, dre)
        call coulomb_rates(ma, qa, v, ns, nb, tb, mb, zb, dpar, drift, dhh, dpe, dre)
        expected = (dplus - dminus)/(2*h) + dpar*(2/v - ma*v/(temp*e))
        if (abs(drift - expected) > 1.0d-6*abs(drift)) &
            error stop 'Speed operator does not preserve the Maxwellian'
    end do
    print '(a)', 'PASS: BEAMS3D Coulomb rates'

contains

    real(dp) function psi(x)
        real(dp), intent(in) :: x
        psi = erf(sqrt(x)) - 2*sqrt(x/pi)*exp(-x)
    end function psi

    real(dp) function dpsi(x)
        real(dp), intent(in) :: x
        dpsi = 2*sqrt(x/pi)*exp(-x)
    end function dpsi

    real(dp) function coulomb_log(b, v)
        integer, intent(in) :: b
        real(dp), intent(in) :: v
        real(dp) :: debye, v2, mr, bcl, bqm
        debye = sqrt(eps0*temp*e/(e**2*sum(nb*zb**2)))
        v2 = v*v + 2*temp*e/mb(b)
        mr = ma*mb(b)/(ma + mb(b))
        bcl = abs(qa*zb(b)*e)/(4*pi*eps0*mr*v2)
        bqm = hbar/(2*mr*sqrt(v2))
        coulomb_log = log(debye/max(bcl, bqm))
    end function coulomb_log

end program test_coulomb_rates
