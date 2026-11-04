program test_lambda_optional
    use read_wout_mod, only: tosuvspaceLambda
    implicit none
    real(8) :: lam, lu, lv

    ! The invalid-domain return needs no equilibrium or external fixture.
    lam = 1d0
    call tosuvspaceLambda(-1d0, 0d0, 0d0, lam=lam)
    if (lam /= 0d0) error stop 'lambda was not initialized'
    lu = 1d0
    call tosuvspaceLambda(-1d0, 0d0, 0d0, lu=lu)
    if (lu /= 0d0) error stop 'lambda_u was not initialized'
    lv = 1d0
    call tosuvspaceLambda(-1d0, 0d0, 0d0, lv=lv)
    if (lv /= 0d0) error stop 'lambda_v was not initialized'
    call tosuvspaceLambda(-1d0, 0d0, 0d0)
    call tosuvspaceLambda(-1d0, 0d0, 0d0, lam, lu, lv)
    if (any([lam,lu,lv] /= 0d0)) error stop 'output initialization changed'
    print *, 'PASS: optional lambda outputs'
end program
