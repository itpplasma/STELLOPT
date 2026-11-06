program test_coulomb_log_general
    ! ln(Lambda) = ln(lambda_D / max(b_cl, b_qm)) with the SI classical
    ! distance of closest approach and the quantum (de Broglie) limit,
    ! NRL Plasma Formulary (2019), p. 34. Both limits are checked, and an
    ! electron (Z < 0) must give the same logarithm as a positive charge.
    use collision_operators, only: coulomb_log_nrl_general
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: e = 1.602176634d-19, eps0 = 8.8541878128d-12
    real(dp), parameter :: hbar = 1.054571817d-34, pi = acos(-1.0_dp)
    real(dp), parameter :: mp = 1.67262192369d-27, me = 9.1093837015d-31
    real(dp), parameter :: ld = 1.0d-4
    real(dp) :: u(2), mr, bcl, bqm, expected, got
    integer :: k, nfail

    nfail = 0
    ! Classical limit dominates at 1e5 m/s, quantum limit at 1e7 m/s.
    u = [1.0d5, 1.0d7]
    mr = mp/2
    do k = 1, 2
        bcl = e*e/(4*pi*eps0*mr*u(k)**2)
        bqm = hbar/(2*mr*u(k))
        expected = log(ld/max(bcl, bqm))
        got = coulomb_log_nrl_general(mp, 1.0_dp, mp, 1.0_dp, ld, u(k))
        call check('p-p', got, expected)
    end do
    mr = me*mp/(me + mp)
    bcl = e*e/(4*pi*eps0*mr*u(1)**2)
    bqm = hbar/(2*mr*u(1))
    expected = log(ld/max(bcl, bqm))
    got = coulomb_log_nrl_general(me, -1.0_dp, mp, 1.0_dp, ld, u(1))
    call check('e-p', got, expected)
    if (nfail > 0) error stop 1
    print '(A)', 'PASS test_coulomb_log_general'
contains
    subroutine check(label, got, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected
        if (abs(got - expected) > 1.0d-6*abs(expected)) then
            print '(A,2ES16.8)', 'FAIL '//label//' got, expected:', got, expected
            nfail = nfail + 1
        end if
    end subroutine check
end program test_coulomb_log_general
