program test_flx2cyl_c1
    ! The radial interpolation of flx2cyl must be C1 across full-mesh nodes
    ! and cubically accurate for smooth coefficients, so that quantities
    ! built from Rs and Zs (grad u, the vector potential) are continuous.
    use stel_kinds, only: rprec
    use vmec_utils, only: flx2cyl
    implicit none
    integer, parameter :: ns = 17
    real(rprec) :: coeff(ns,0:0,0:1,3), q(3), x(3), rs, zs, s, rsl, rsr, worst
    integer :: j, status

    coeff = 0
    do j = 1, ns
        s = real(j-1,rprec)/(ns-1)
        coeff(j,0,0,1) = exp(s)
        coeff(j,0,1,1) = sqrt(s)*cos(s)
        coeff(j,0,1,2) = sqrt(s)*cos(s)
    end do
    ! Rs from both sides of interior nodes.
    worst = 0
    do j = 3, ns - 2
        s = real(j-1,rprec)/(ns-1)
        q = [s - 1e-9_rprec, 0.7_rprec, 0.0_rprec]
        call flx2cyl(coeff,q,x,ns,0,2,1,.false.,.false.,status,rs=rsl,zs=zs)
        q(1) = s + 1e-9_rprec
        call flx2cyl(coeff,q,x,ns,0,2,1,.false.,.false.,status,rs=rsr,zs=zs)
        worst = max(worst, abs(rsr - rsl))
    end do
    print '(a,es10.3)', 'max Rs jump across nodes: ', worst
    if (worst > 1e-6_rprec) error stop 'Rs is discontinuous at radial nodes'
    ! Interval midpoints away from the axis: interpolation error.
    worst = 0
    do j = 3, ns - 2
        s = (real(j,rprec) - 0.5_rprec)/(ns-1)
        q = [s, 0.7_rprec, 0.0_rprec]
        call flx2cyl(coeff,q,x,ns,0,2,1,.false.,.false.,status,rs=rs,zs=zs)
        worst = max(worst, abs(x(1) - exp(s) - sqrt(s)*cos(s)*cos(q(2))))
    end do
    print '(a,es10.3)', 'max midpoint error: ', worst
    if (worst > 2e-6_rprec) error stop 'radial interpolation is not cubic'
    print *, 'PASS: C1 radial interpolation'
end program test_flx2cyl_c1
