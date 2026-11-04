PROGRAM test_gc2fo
    USE stel_kinds, ONLY: dp => rprec
    USE beams3d_physics_mod, ONLY: beams3d_gc2fo
    USE beams3d_grid, ONLY: nr,nphi,nz,raxis,phiaxis,zaxis,phimax, &
        hr,hp,hz,hri,hpi,hzi,BR4D,BPHI4D,BZ4D,MODB4D
    USE beams3d_lines, ONLY: mymass,mycharge
    USE beams3d_runtime, ONLY: pi,pi2
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
    IMPLICIT NONE
    INTEGER :: i,j,k,failures
    REAL(dp) :: q(6),bphi,bz,bnorm,vperp,radius,measured,speed,mu
    REAL(dp), PARAMETER :: tilt(3) = [0D0,.4D0,2.4D0]

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
    mymass=6.69048769476D-27; mycharge=3.204353268D-19
    vperp=1D7
    failures=0
    DO k=1,3
        bphi=2.4D0
        IF (k == 3) bphi=0
        bz=tilt(k)
        bnorm=HYPOT(bphi,bz)
        mu=.5D0*mymass*vperp**2/bnorm
        radius=mymass*vperp/(mycharge*bnorm)
        BR4D=0; BPHI4D=0; BZ4D=0; MODB4D=0
        BPHI4D(1,:,:,:)=bphi; BZ4D(1,:,:,:)=bz; MODB4D(1,:,:,:)=bnorm
        DO j=1,256
            q=[5D0,0D0,0D0,1D6,mu,0D0]
            CALL beams3d_gc2fo(q)
            IF (.NOT.ALL(ieee_is_finite(q))) THEN
                failures=failures+1
                CYCLE
            END IF
            measured=SQRT((q(1)*COS(q(2))-5D0)**2+(q(1)*SIN(q(2)))**2+q(3)**2)
            speed=SQRT(SUM(q(4:6)**2))
            IF (ABS(measured/radius-1) > 5D-12) failures=failures+1
            IF (ABS(speed/HYPOT(1D6,vperp)-1) > 5D-12) failures=failures+1
        END DO
    END DO
    IF (failures /= 0) THEN
        WRITE(*,*) 'GC-to-particle radius or velocity failures: ',failures
        ERROR STOP 1
    END IF
    WRITE(*,*) 'GC-to-particle radius and velocity: PASS'
END PROGRAM test_gc2fo
