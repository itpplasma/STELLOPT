program test_flx2cyl_fullmesh
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use stel_kinds, only: rprec
    use vmec_utils, only: flx2cyl
    implicit none
    integer, parameter :: ns = 5
    real(rprec) :: coeff(ns,0:0,0:1,3), q(3), x(3), rs, zs, expected, slope
    real(rprec), parameter :: samples(7) = [0.13_rprec, 0.12_rprec, &
        0.249999_rprec, 0.250001_rprec, 0.374999_rprec, 0.375001_rprec, 1.1_rprec]
    real(rprec) :: s, slo, shi, weight
    integer :: j, k, lower, status

    coeff = 0
    do j = 1, ns
        s = real(j-1,rprec)/(ns-1)
        coeff(j,0,0,1) = 4 + 0.2_rprec*s*s
        coeff(j,0,1,1) = sqrt(s)
        coeff(j,0,1,2) = sqrt(s)
    end do
    do k = 1, size(samples)
        s = samples(k); q = [s, 0.7_rprec, 0.0_rprec]
        call flx2cyl(coeff,q,x,ns,0,2,1,.false.,.false.,status,rs=rs,zs=zs)
        if (status /= 0 .or. .not.all(ieee_is_finite(x))) error stop 'forward map'
        lower = min(floor(s*(ns-1)), ns-2)
        slo = real(lower,rprec)/(ns-1); shi = real(lower+1,rprec)/(ns-1)
        weight = (s-slo)/(shi-slo)
        expected = 4 + 0.2_rprec*((1-weight)*slo*slo + weight*shi*shi) &
                   + sqrt(s)*cos(q(2))
        slope = 0.2_rprec*(slo+shi) + cos(q(2))/(2*sqrt(s))
        print '(a,4es25.16)', 'MAP ', s, x(1), expected, rs
        if (abs(x(1)-expected) > 1e-12_rprec) error stop 'wrong full-mesh interval'
        if (abs(rs-slope) > 1e-12_rprec) error stop 'wrong full-mesh derivative'
    end do
    print *, 'PASS: full-mesh radial intervals'
end program
