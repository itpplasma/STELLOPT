program test_vector_potential_flux
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use stel_kinds, only: rprec
    use stel_constants, only: twopi
    use read_wout_mod
    use vmec_utils, only: GetAcyl
    implicit none
    integer :: j, k, idir, ierr
    real(rprec), parameter :: major = 4, torflux = 3, polflux = 0.8_rprec
    real(rprec) :: s, u, radius, r, z, ar, ap, az, lam, cov(3), expected(3)

    ! Synthetic circular surfaces make the flux circulations analytic.
    ns = 9; mpol = 2; ntor = 0; nfp = 1; mnmax = 2; ntmax = 1
    lasym = .false.; lthreed = .false.; lwout_opened = .true.
    allocate(xm(2), xn(2), rmnc(2,ns), zmns(2,ns), lmns(2,ns))
    allocate(phi(ns), chi(ns), phipf(ns))
    xm = [0.0_rprec, 1.0_rprec]; xn = 0
    rmnc = 0; zmns = 0; lmns = 0
    do idir = -1, 1, 2
        ! VMEC writes phi and phipf times signgs, chi without it.
        isigng = idir
        do j = 1, ns
            s = real(j-1,rprec)/(ns-1)
            rmnc(1,j) = major; rmnc(2,j) = sqrt(s); zmns(2,j) = sqrt(s)
            phi(j) = idir*torflux*s; chi(j) = polflux*s; phipf(j) = idir*torflux
            if (j > 1) lmns(2,j) = 0.07_rprec*sqrt((j-1.5_rprec)/(ns-1))
        end do
        do j = 1, 3
            s = 0.15_rprec + 0.2_rprec*j
            radius = sqrt(s)
            do k = 1, 4
                u = 0.3_rprec + k
                r = major + radius*cos(u); z = radius*sin(u)
                call GetAcyl(r,0.0_rprec,z,ar,ap,az,info=ierr)
                if (ierr /= 0) error stop 'synthetic inverse mapping failed'
                cov(1) = (ar*cos(u)+az*sin(u))/(2*radius)
                cov(2) = radius*(-ar*sin(u)+az*cos(u))
                cov(3) = r*ap
                lam = 0.07_rprec*radius*sin(u)
                ! A poloidal/toroidal circuit must enclose the full WOUT flux,
                ! not another factor of 2*pi times that flux.
                expected = [-lam*torflux, torflux*s, -polflux*s]/twopi
                if (.not.all(ieee_is_finite(cov))) error stop 'nonfinite vector potential'
                if (maxval(abs(cov-expected)) > 2e-7_rprec) then
                    print *, 'covariants / expected: ', cov, expected
                    error stop 'vector potential flux normalization'
                end if
            end do
        end do
    end do
    print *, 'PASS: vector potential flux normalization'
end program
