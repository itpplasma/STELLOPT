PROGRAM test_part2gc
    ! Particle-to-gyrocentre transformation in a purely toroidal field of
    ! constant magnitude: curvature -Rhat/R, no grad|B| and no torsion, so the
    ! first-order corrections follow from the curvature alone,
    !    vll    = vphi (1 + rho rhohat_R / R),
    !    moment = mu0 - rho m vphi^2 rhohat_R / (B R),
    ! with rho = m vperp/(q B) and rhohat_R = vZ/vperp for a positive charge.
    USE stel_kinds, ONLY: dp => rprec
    USE beams3d_physics_mod, ONLY: beams3d_part2gc
    USE beams3d_grid, ONLY: nr,nphi,nz,raxis,phiaxis,zaxis,phimax, &
        hr,hp,hz,hri,hpi,hzi,BR4D,BPHI4D,BZ4D,MODB4D
    USE beams3d_lines, ONLY: mymass,mycharge,moment
    USE beams3d_runtime, ONLY: pi,pi2
    IMPLICIT NONE
    INTEGER :: i,failures
    REAL(dp) :: q(6),vr,vphi,vz,vperp,rho,mu0,vll_exact,mu_exact
    REAL(dp), PARAMETER :: b0 = 2.5D0

    pi=ACOS(-1D0); pi2=2*pi
    nr=4; nphi=4; nz=4; phimax=pi2
    ALLOCATE(raxis(nr),phiaxis(nphi),zaxis(nz))
    ALLOCATE(hr(nr),hp(nphi),hz(nz),hri(nr),hpi(nphi),hzi(nz))
    ALLOCATE(BR4D(8,nr,nphi,nz),BPHI4D(8,nr,nphi,nz), &
             BZ4D(8,nr,nphi,nz),MODB4D(8,nr,nphi,nz))
    DO i=1,4
        raxis(i)=3D0+i; phiaxis(i)=(i-1)*pi2/3; zaxis(i)=-1D0+2D0*(i-1)/3
    END DO
    hr=1; hp=pi2/3; hz=2D0/3
    hri=1/hr; hpi=1/hp; hzi=1/hz
    BR4D=0; BPHI4D=0; BZ4D=0; MODB4D=0
    BPHI4D(1,:,:,:)=b0; MODB4D(1,:,:,:)=b0
    mymass=1.67262192369D-27; mycharge=1.602176634D-19
    failures=0
    DO i=0,5
        vr=1.0D6-3.0D5*i; vphi=2.0D6-4.5D5*i; vz=0.4D6+2.0D5*i
        q=[4.6D0+0.3D0*i,0.4D0+1.1D0*i,0.2D0,vr,vphi,vz]
        vperp=HYPOT(vr,vz)
        rho=mymass*vperp/(mycharge*b0)
        mu0=0.5D0*mymass*vperp**2/b0
        vll_exact=vphi*(1+rho*(vz/vperp)/q(1))
        mu_exact=mu0-rho*mymass*vphi**2*(vz/vperp)/(b0*q(1))
        CALL beams3d_part2gc(q)
        IF (ABS(q(4)/vll_exact-1) > 1D-12 .or. ABS(moment/mu_exact-1) > 1D-12) THEN
            WRITE(*,*) 'case ',i,' vll/exact-1 ',q(4)/vll_exact-1, &
                       ' moment/exact-1 ',moment/mu_exact-1
            failures=failures+1
        END IF
    END DO
    IF (failures /= 0) THEN
        WRITE(*,*) 'Particle-to-gyrocentre curvature failures: ',failures
        ERROR STOP 1
    END IF
    WRITE(*,*) 'Particle-to-gyrocentre curvature terms: PASS'
END PROGRAM test_part2gc
