program test_vector_potential_curl
    ! The poloidal field of curl(A) from GetAcyl must reproduce GetBcyl inside
    ! a VMEC equilibrium. In this axisymmetric fixture B_R and B_Z come from
    ! A_phi alone, which carries the poloidal flux chi. The fixture has
    ! signgs = -1, where WOUT phi carries the Jacobian sign but chi does not.
    use stel_kinds, only: rprec
    use read_wout_mod, only: read_wout_file, rmnc, ns, isigng
    use vmec_utils, only: GetAcyl, GetBcyl
    implicit none
    character(len=512) :: path
    integer :: i, k, ierr
    real(rprec) :: r0, minor, r, z, phi, h, dp, b(3), curl(3), worst, bmag
    real(rprec) :: ap(3), am(3), arp(3), arm(3), azp(3), azm(3)

    call get_command_argument(1, path)
    call read_wout_file(trim(path), ierr)
    if (ierr /= 0) error stop 'cannot read WOUT fixture'
    if (isigng >= 0) error stop 'fixture must have signgs = -1'
    r0 = rmnc(1,1); minor = rmnc(2,ns)
    h = 1.0e-4_rprec*minor; dp = 1.0e-4_rprec; phi = 0.3_rprec
    worst = 0
    do i = 1, 3
        do k = 0, 5
            r = r0 + 0.2_rprec*i*minor*cos(1.1_rprec*k + 0.2_rprec)
            z = 0.2_rprec*i*minor*sin(1.1_rprec*k + 0.2_rprec)
            call GetBcyl(r, phi, z, b(1), b(2), b(3), info=ierr)
            if (ierr /= 0) error stop 'GetBcyl failed'
            arp = potential(r + h, phi, z); arm = potential(r - h, phi, z)
            azp = potential(r, phi, z + h); azm = potential(r, phi, z - h)
            ap = potential(r, phi + dp, z); am = potential(r, phi - dp, z)
            curl(1) = (ap(3) - am(3))/(2*dp*r) - (azp(2) - azm(2))/(2*h)
            curl(2) = (azp(1) - azm(1))/(2*h) - (arp(3) - arm(3))/(2*h)
            curl(3) = ((r + h)*arp(2) - (r - h)*arm(2))/(2*h*r) &
                      - (ap(1) - am(1))/(2*dp*r)
            bmag = sqrt(sum(b*b))
            worst = max(worst, max(abs(curl(1) - b(1)), abs(curl(3) - b(3)))/bmag)
        end do
    end do
    print '(a,es10.3)', 'max poloidal |curl A - B|/|B| = ', worst
    if (worst > 1.0e-3_rprec) error stop 'curl of vector potential differs from B'
    print '(a)', 'PASS: vector potential curl'

contains

    function potential(rr, pp, zz) result(a)
        real(rprec), intent(in) :: rr, pp, zz
        real(rprec) :: a(3)
        integer :: info

        call GetAcyl(rr, pp, zz, a(1), a(2), a(3), info=info)
        if (info /= 0) error stop 'GetAcyl failed'
    end function potential

end program test_vector_potential_curl
