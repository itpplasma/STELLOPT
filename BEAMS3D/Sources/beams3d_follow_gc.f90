!-----------------------------------------------------------------------
!     Module:        beams3d_follow_gc
!     Authors:       S. Lazerson (lazerson@pppl.gov), M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   This subroutine follows the particles through
!                    the grid.  The general ODE which must be solved
!                    can be written:
!                        dX        
!                       ----   = fgc(t,q,qdot)
!                        dt        
!                    where X=(R,phi,Z).  This allows particle trajectories
!                    to be parametrized in terms of time.
!                    
!-----------------------------------------------------------------------
SUBROUTINE beams3d_follow_gc
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE beams3d_runtime
    USE beams3d_lines
    USE beams3d_grid, ONLY: tmin, tmax, delta_t, wall_load, wall_shine, &
                            plasma_mass, plasma_Zmean, therm_factor, &
                            rho_fullorbit, rho_help, E_kick, freq_kick, &
                            nr_fida, nphi_fida, nz_fida, nenergy_fida, npitch_fida,raxis
    USE mpi_params ! MPI
    USE beams3d_write_par
    USE beams3d_physics_mod, ONLY: beams3d_calc_dt, coll_t0
    USE collision_operators, ONLY: SET_COULOMB_FACTOR
    USE safe_open_mod, ONLY: safe_open
    USE mpi_inc
    !-----------------------------------------------------------------------
    !     Local Variables
    !          status       MPI stats indicator
    !          mystart      Starting fieldline for each processor
    !          mypace       Number of fieldlines each processor should follow
    !          i            Dummy index
    !          j            Dummy index
    !          sender       Dummy index
    !          ier          Error flag
    !          l            Dummy index
    !          neqs_nag     Number of ODEs to solve (not limited to NAG routines)
    !          l2           Index helper
    !          itol         LSODE tolerance flag
    !          itask        LSODE flag (overshoot and interpolate)
    !          istate       LSODE restart flag
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER :: ier, l, neqs_nag, itol, itask, &
               istate, iopt, lrw, liw, mf, out, iunit
    INTEGER, ALLOCATABLE :: iwork(:), itemp(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: w(:), q(:)
    DOUBLE PRECISION :: tf_nag, eps_temp, t_nag, &
                        tol_nag, rtol, s_fullorbit, &
                        dtmin, dtmax
    DOUBLE PRECISION :: atol(4)
    DOUBLE PRECISION :: rkh_work(4, 2)
    DOUBLE PRECISION :: h_dp, vscale_dp
    CHARACTER*1 :: relab
    DOUBLE PRECISION, PARAMETER :: e_charge      = 1.602176634D-19 !e_c

    !-----------------------------------------------------------------------
    !     External Functions
    !          fgc_eom            RHS of ODE integrator (for NAG)    for BEAMS3D
    !          D02CJF               NAG ODE Solver
    !          D02CJW               NAG Dummy function
    !          jacobian_lsode       Jacobian function (for LSODE, not currently utilized)
    !-----------------------------------------------------------------------
    EXTERNAL D02CJF, D02CJW, fgc_eom, D02CJX, out_beams3d_gc
    EXTERNAL fgc_lsode, jacobian_lsode
    EXTERNAL fgc_rkh68
    !-----------------------------------------------------------------------
    !     Begin Subroutine
    !-----------------------------------------------------------------------
    ! Initializations
    s_fullorbit = SIGN(rho_fullorbit*rho_fullorbit,rho_fullorbit)
    ier = 0
    tol_nag = follow_tol
    neqs_nag = 4
    relab = "M"
    mf = 10
    ALLOCATE(q(neqs_nag))

! Screen output so we know what's happening
    IF (lverb) THEN
       ! Do a calculation of nstep and delta-t
       myline = MAXLOC(ABS(vll_start),1)
       my_end = t_end(myline)
       CALL beams3d_calc_dt(1,q(1),q(2),q(3),dtmin)
       myline = MINLOC(ABS(vll_start),1)
       my_end = t_end(myline)
       CALL beams3d_calc_dt(1,q(1),q(2),q(3),dtmax)
       ! Screen output
       WRITE(6, '(A)') '----- FOLLOWING GYROCENTER TRAJECTORIES -----'
       WRITE(6, '(A,A)')          '       Method: ', TRIM(int_type)
       WRITE(6, '(A,I9)')          '   Particles: ', nparticles
       WRITE(6, '(A,I9,2(A,EN12.3))') '       Steps: ', ndt_max*NPOINC, '   dt_min: ', dtmin,'   dt_max: ', dtmax
       WRITE(6, '(A,I9)')          '      NPOINC: ', NPOINC
       SELECT CASE(TRIM(int_type))
          CASE("NAG")
             WRITE(6, '(A,EN12.3,A,A1)') '         Tol: ', follow_tol, '  Type: ', relab
          CASE("LSODE")
             WRITE(6, '(A,EN12.3,A,I2)') '         Tol: ', follow_tol, '  Type: ', mf
          CASE("DOPRI5")
             WRITE(6, '(A,EN12.3)') '         Tol: ', follow_tol
       END SELECT
       WRITE(6, '(5X,A,I3,A)', ADVANCE = 'no') 'Trajectory Calculation [', 0, ']%'
       CALL FLUSH(6)
    END IF

    ! Follow Trajectories
    IF (mystart_save <= nparticles) THEN
        SELECT CASE (TRIM(int_type))
            CASE ("NAG")
                ALLOCATE(w(neqs_nag * 21 + 28), STAT = ier)
                IF (ier /= 0) CALL handle_err(ALLOC_ERR, 'W', ier)
                DO l = mystart_save, myend_save
                    t_nag = t_last(l)
                    ! Don't do particle if stopped or beyond the full_orbit limit
                    IF (t_nag>t_end(l)) CYCLE
                    ! Particle indicies
                    myline = l
                    mytdex = 1; ndt = 0
                    IF (lbeam) mytdex = 3
                    ! Don't do full_orbit particles
                    IF (sqrt(S_lines(mytdex-1,l))>rho_fullorbit) CYCLE
                    ! Particle Parameters
                    q(1) = R_lines(mytdex-1,l)
                    q(2) = PHI_lines(mytdex-1,l)
                    q(3) = Z_lines(mytdex-1,l)
                    q(4) = vll_lines(mytdex-1,l)
                    xlast = q(1)*cos(q(2))
                    ylast = q(1)*sin(q(2))
                    zlast = q(3)
                    previous_time = t_nag
                    previous_vll = q(4)
                    moment = moment_lines(mytdex-1,l)
                    previous_moment = moment
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    previous_energy = 0.5*mymass*q(4)*q(4) + &
                                      moment*B_lines(mytdex-1,l)
                    E_by_v=mymass*0.5d-3/e_charge
                    mybeam = Beam(l)
                    my_end = t_end(l)
                    ltherm = .false.
                    lneut  = .false.
                    lcollision = lbeam
                    ! Collision parameters
                    fact_kick = 2*E_kick*mycharge/(mymass*pi2*pi2*freq_kick*freq_kick*SQRT(pi*1E-7*plasma_mass))
                    fact_pa   = plasma_mass/(mymass*plasma_Zmean)
                    CALL SET_COULOMB_FACTOR(mymass,myZ,plasma_mass)
                    ! Now calc dt
                    CALL beams3d_calc_dt(1,q(1),q(2),q(3),dt)
                    coll_t0 = t_nag
                    tf_nag = t_nag+dt
                    ndt = 0
                    DO ! Must do it this way becasue lbeam changes q(4) values
#if defined(NAG)
                       CALL D02CJF(t_nag,tf_nag,neqs_nag,q,fgc_eom,tol_nag,relab,out_beams3d_gc,D02CJW,w,ier)
#endif
                       IF (ier < 0) CALL handle_err(D02CJF_ERR, 'beams3d_follow', ier)
                       t_last(l) = tf_nag ! Save the value here in case out_beams3d changes it
                       CALL out_beams3d_gc(tf_nag,q)
                       IF (ABS(tf_nag) > ABS(my_end)+0.5D0*ABS(dt)) EXIT
                    END DO
                END DO
            CASE ("RKH68")
                ier = 0
                DO l = mystart_save, myend_save
                    tf_nag = t_last(l)
                    ! Don't do particle if stopped
                    IF (tf_nag>t_end(l)) CYCLE
                    ! Particle indicies
                    myline = l
                    mytdex = 1; ndt = 0
                    IF (lbeam) mytdex = 3
                    ! Don't do full_orbit particles
                    IF (sqrt(S_lines(mytdex-1,l))>rho_fullorbit) CYCLE
                    ! Particle Parameters
                    q(1) = R_lines(mytdex-1,l)
                    q(2) = PHI_lines(mytdex-1,l)
                    q(3) = Z_lines(mytdex-1,l)
                    q(4) = vll_lines(mytdex-1,l)
                    xlast = q(1)*cos(q(2))
                    ylast = q(1)*sin(q(2))
                    zlast = q(3)
                    moment = moment_lines(mytdex-1,l)
                    t_nag = tf_nag
                    previous_time = t_nag
                    previous_vll = q(4)
                    previous_moment = moment
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    previous_energy = 0.5*mymass*q(4)*q(4) + &
                                      moment*B_lines(mytdex-1,l)
                    E_by_v=mymass*0.5d-3/e_charge
                    mybeam = Beam(l)
                    my_end = t_end(l)
                    ltherm = .false.
                    lneut  = .false.
                    lcollision = lbeam
                    ! Collision parameters
                    fact_kick = 2*E_kick*mycharge/(mymass*pi2*pi2*freq_kick*freq_kick*SQRT(pi*1E-7*plasma_mass))
                    fact_pa   = plasma_mass/(mymass*plasma_Zmean)
                    CALL SET_COULOMB_FACTOR(mymass,myZ,plasma_mass)
                    ! Now calc dt
                    CALL beams3d_calc_dt(1,q(1),q(2),q(3),dt)
                    coll_t0 = t_nag
                    tf_nag = t_nag+dt
                    ndt = 0
                    ! Setup DRKHVG parameters
                    iopt = 0
                    DO
                        CALL drkhvg(t_nag, q, neqs_nag, dt, 2, fgc_rkh68, rkh_work, iopt, ier)
                        IF (ier < 0) CALL handle_err(RKH68_ERR, 'beams3d_follow', ier)
                        q(1:neqs_nag) = rkh_work(1:neqs_nag,2)
                        ! q is now the state at tf_nag; out_beams3d_gc
                        ! advances tf_nag by dt for the next step.
                        t_nag = tf_nag
                        t_last(l) = tf_nag ! Save the value here in case out_beams3d changes it
                        CALL out_beams3d_gc(tf_nag,q)
                        IF (ABS(tf_nag) > ABS(my_end)+0.5D0*ABS(dt)) EXIT
                    END DO
                END DO
            CASE ("LSODE","DLSODE","DOPRI5")
                ! Open a file
                iunit = myworkid+400
                ier = 0
                IF (ldebug) THEN
                   !CALL safe_open(iunit,ier,'lsode_err_'//TRIM(id_string),'replace','formatted')
                   iunit = 6
                   CALL xsetun(iunit)
                   iunit = 1
                   CALL xsetf(iunit)
                ELSE
                   iunit = 0
                   CALL xsetf(iunit)
                END IF
                ! Allocate helpers
                lrw = 20 + 16 * neqs_nag
                liw = 20
                ALLOCATE(w(lrw))
                ALLOCATE(iwork(liw))
                ier = 0
                DO l = mystart_save, myend_save
                    t_nag = t_last(l)
                    ! Don't do particle if stopped
                    IF (t_nag>t_end(l)) CYCLE
                    ! Particle indicies
                    myline = l
                    mytdex = 1
                    IF (lbeam) mytdex = 3
                    ! Don't do full_orbit particles
                    IF (S_lines(mytdex-1,l)>=s_fullorbit) CYCLE
                    ! Particle Parameters
                    q(1) = R_lines(mytdex-1,l)
                    q(2) = PHI_lines(mytdex-1,l)
                    q(3) = Z_lines(mytdex-1,l)
                    q(4) = vll_lines(mytdex-1,l)
                    xlast = q(1)*cos(q(2))
                    ylast = q(1)*sin(q(2))
                    zlast = q(3)
                    previous_time = t_nag
                    previous_vll = q(4)
                    moment = moment_lines(mytdex-1,l)
                    previous_moment = moment
                    mycharge = charge(l)
                    myZ = Zatom(l)
                    mymass = mass(l)
                    previous_energy = 0.5*mymass*q(4)*q(4) + &
                                      moment*B_lines(mytdex-1,l)
                    E_by_v=mymass*0.5d-3/e_charge
                    mybeam = Beam(l)
                    my_end = t_end(l)
                    ltherm = .false.
                    lneut  = .false.
                    lcollision = lbeam .or. lcollision
                    ! Collision parameters
                    fact_kick = 2*E_kick*mycharge/(mymass*pi2*pi2*freq_kick*freq_kick*SQRT(pi*1E-7*plasma_mass))
                    fact_pa   = plasma_mass/(mymass*plasma_Zmean)
                    CALL SET_COULOMB_FACTOR(mymass,myZ,plasma_mass)
                    ! Now calc dt
                    CALL beams3d_calc_dt(1,q(1),q(2),q(3),dt)
                    coll_t0 = t_nag
                    tf_nag = t_nag+dt
                    ndt = 0
                    ! Setup LSODE parameters
                    iopt = 0 ! No optional output
                    w = 0; iwork = 0; itask = 1; istate = 1;
                    itol = 2; rtol = follow_tol; atol(:) = follow_tol
                    h_dp = 0
                    vscale_dp = SQRT(2*previous_energy/mymass)
                    DO
                        IF (TRIM(int_type) == 'DOPRI5') THEN
                           CALL dopri5_advance(t_nag, tf_nag, q, h_dp, follow_tol, vscale_dp)
                           t_last(l) = tf_nag
                           CALL out_beams3d_gc(tf_nag,q)
                           IF ((ABS(tf_nag) > ABS(my_end)+0.5D0*ABS(dt)) &
                               .or. (rho_help > rho_fullorbit)) EXIT
                           CYCLE
                        END IF
                        IF (lcollision) istate = 1
                        CALL FLUSH(6)
                        CALL DLSODE(fgc_lsode, neqs_nag, q, t_nag, tf_nag, itol, rtol, atol, itask, istate, &
                                   iopt, w, lrw, iwork, liw, jacobian_lsode, mf)
                        IF ((istate == -3) .or. (istate == -4)) THEN
                           ! BIG  DEBUG MESSAGE
                           CALL FLUSH(6)
                           WRITE(6,*) '------------------'
                           WRITE(6,*) '     ',myworkid, l, myline,mf
                           WRITE(6,*) '     ',myworkid, neqs_nag, t_nag, tf_nag
                           WRITE(6,*) '     ',myworkid, q, moment
                           WRITE(6,*) '     ',myworkid, itol, rtol, atol
                           WRITE(6,*) '     ',myworkid, itask,istate,iopt
                           WRITE(6,*) '     ',myworkid, w(5),w(6),w(7)
                           WRITE(6,*) '     ',myworkid, iwork(5),iwork(6),iwork(7)
                           WRITE(6,*) '     ',myworkid, w(11:14)
                           WRITE(6,*) '     ',myworkid, iwork(11:18)
                           WRITE(6,*) '------------------'
                           CALL FLUSH(6)
                           CALL handle_err(LSODE_ERR, 'beams3d_follow_gc', istate)
                        END IF
                        iwork(11) = 0; iwork(12) = 0; iwork(13) = 0
                        t_last(l) = tf_nag ! Save the value here in case out_beams3d changes it
                        CALL out_beams3d_gc(tf_nag,q)
                        IF ( (istate == -1) .or. (istate ==-2) &
                                            .or. (ABS(tf_nag) > ABS(my_end)+0.5D0*ABS(dt)) &
                                            .or. (rho_help > rho_fullorbit)) EXIT
                    END DO
                END DO
                IF (ldebug) CLOSE(iunit)
             CASE ('DEBUG')
                DO l = 0, npoinc
                   R_lines(l,mystart_save:myend_save) = REAL(l)
                END DO
                B_lines(0:npoinc,mystart_save:myend_save) = REAL(myworkid)
        END SELECT
    END IF
      
    IF (lverb) THEN
       CALL backspace_out(6,36)
       CALL FLUSH(6)
       WRITE(6,'(36X)',ADVANCE='no')
       CALL FLUSH(6)
       CALL backspace_out(6,36)
       WRITE(6,*)
       CALL FLUSH(6)
    END IF    

    CALL dist5d_flush

    ! Check for crash
    istate = 0
    CALL handle_err(MPI_CHECK,'beams3d_follow',istate)

    !Deallocations
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF (ALLOCATED(w)) DEALLOCATE(w)
    IF (ALLOCATED(iwork)) DEALLOCATE(iwork)
    RETURN

CONTAINS

    SUBROUTINE dopri5_advance(t, tout, y, h, tol, vscale)
    ! Dormand-Prince 5(4) from t to tout with error control. A one-step
    ! method needs no restart after the collision kick at every output time,
    ! and h carries the accepted step size over to the next interval.
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(INOUT) :: t, y(4), h
    DOUBLE PRECISION, INTENT(IN) :: tout, tol, vscale
    DOUBLE PRECISION, PARAMETER :: &
       c2 = 1D0/5, c3 = 3D0/10, c4 = 4D0/5, c5 = 8D0/9, &
       a21 = 1D0/5, &
       a31 = 3D0/40, a32 = 9D0/40, &
       a41 = 44D0/45, a42 = -56D0/15, a43 = 32D0/9, &
       a51 = 19372D0/6561, a52 = -25360D0/2187, a53 = 64448D0/6561, a54 = -212D0/729, &
       a61 = 9017D0/3168, a62 = -355D0/33, a63 = 46732D0/5247, a64 = 49D0/176, a65 = -5103D0/18656, &
       b1 = 35D0/384, b3 = 500D0/1113, b4 = 125D0/192, b5 = -2187D0/6784, b6 = 11D0/84, &
       e1 = 71D0/57600, e3 = -71D0/16695, e4 = 71D0/1920, e5 = -17253D0/339200, &
       e6 = 22D0/525, e7 = -1D0/40
    DOUBLE PRECISION :: k1(4), k2(4), k3(4), k4(4), k5(4), k6(4), k7(4), &
                        ynew(4), yerr(4), sc(4), hs, errn, fac, dir
    LOGICAL :: lclip
    INTEGER :: nrej
    dir = SIGN(1D0, tout - t)
    IF (h <= 0) h = ABS(tout - t)
    CALL fgc_eom(t, y, k1)
    nrej = 0
    DO WHILE (dir*(tout - t) > 0)
       lclip = h >= ABS(tout - t)
       hs = dir*MIN(h, ABS(tout - t))
       CALL fgc_eom(t + c2*hs, y + hs*a21*k1, k2)
       CALL fgc_eom(t + c3*hs, y + hs*(a31*k1 + a32*k2), k3)
       CALL fgc_eom(t + c4*hs, y + hs*(a41*k1 + a42*k2 + a43*k3), k4)
       CALL fgc_eom(t + c5*hs, y + hs*(a51*k1 + a52*k2 + a53*k3 + a54*k4), k5)
       CALL fgc_eom(t + hs, y + hs*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5), k6)
       ynew = y + hs*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
       CALL fgc_eom(t + hs, ynew, k7)
       yerr = hs*(e1*k1 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7)
       ! Positions relative to the major radius, phi in radians, v_par to the speed.
       sc(1) = tol*MAX(ABS(y(1)), ABS(ynew(1)))
       sc(2) = tol
       sc(3) = sc(1)
       sc(4) = tol*vscale
       errn = SQRT(SUM((yerr/sc)**2)/4)
       ! The rejection cap only guards against a NaN error norm.
       IF (errn <= 1 .or. nrej >= 50) THEN
          t = t + hs
          IF (lclip) t = tout
          y = ynew
          k1 = k7
          fac = MIN(5D0, MAX(0.2D0, 0.9D0*MAX(errn, 1D-10)**(-0.2D0)))
          IF (nrej > 0) fac = MIN(fac, 1D0)
          IF (lclip) THEN
             h = MAX(h, ABS(hs)*fac)
          ELSE
             h = ABS(hs)*fac
          END IF
          nrej = 0
       ELSE
          h = ABS(hs)*MAX(0.2D0, 0.9D0*errn**(-0.2D0))
          nrej = nrej + 1
       END IF
    END DO
    END SUBROUTINE dopri5_advance

END SUBROUTINE beams3d_follow_gc
