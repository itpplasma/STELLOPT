!-----------------------------------------------------------------------
!     Module:        beams3d_spline3d_setup
!     Authors:       S. Lazerson
!     Date:          10/22/2025
!     Description:   Automate allocation and deallocation of 3D splines
!-----------------------------------------------------------------------
      SUBROUTINE beams3d_spline3d_setup
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE beams3d_globals, ONLY: nr, nphi, nz, spline_type
      USE beams3d_runtime, ONLY: handle_err, lvac, EZSPLINE_ERR, &
         MPI_BARRIER_ERR
      USE mpi_params, ONLY: master, myid_sharmem, MPI_COMM_SHARMEM
      USE mpi_sharmem, ONLY: mpialloc, mpidealloc
      USE mpi_inc
      USE EZspline_obj, ONLY: EZspline3_r8
      USE EZspline, ONLY: EZspline_init, EZspline_setup, EZspline_free
      USE stel_kinds, ONLY: rprec
      USE beams3d_grid

!-----------------------------------------------------------------------
!     Input parameters
!-----------------------------------------------------------------------
      IMPLICIT NONE


!-----------------------------------------------------------------------
!     Local variables
!-----------------------------------------------------------------------
      INTEGER :: i, ier, nshar
      INTEGER :: bcs1(2), bcs2(2), bcs3(2)

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
!     Each spline is built by one rank of the node, round robin, since
!     the coefficient arrays are shared; the inputs are freed only after
!     all ranks are done with them.
!-----------------------------------------------------------------------
      ier = 0 ; bcs1=(/ 0, 0/); bcs2=(/-1,-1/); bcs3=(/ 0, 0/)
      nshar = 1
      lcubic_grid = (spline_type(1:1) == 'C' .or. spline_type(1:1) == 'c')
#if defined(MPI_OPT)
      CALL MPI_COMM_SIZE(MPI_COMM_SHARMEM, nshar, ier)
#endif

      IF (.not. lvac) THEN
         CALL mpialloc(TE4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TE4D)
         CALL mpialloc(NE4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NE4D)
         CALL mpialloc(TI4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TI4D)
         CALL mpialloc(ZEFF4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_ZEFF4D)
         CALL mpialloc(OMEG4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_OMEG4D)
         CALL mpialloc(NI5D, 8, nr, nphi, nz, NION, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NI5D)
      END IF
      CALL mpialloc(BR4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BR4D)
      CALL mpialloc(BPHI4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BPHI4D)
      CALL mpialloc(BZ4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_BZ4D)
      CALL mpialloc(MODB4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_MODB4D)
      CALL mpialloc(U4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_U4D)
      CALL mpialloc(POT4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_POT4D)
      CALL mpialloc(RHO4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_RHO4D)
      CALL mpialloc(XRHO4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_XRHO4D)
      CALL mpialloc(YRHO4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_YRHO4D)
      CALL mpialloc(S4D, 8, nr, nphi, nz, myid_sharmem, 0, MPI_COMM_SHARMEM, win_S4D)

      IF (.not. lvac) THEN
         CALL build(TE, TE4D, 0)
         CALL build(NE, NE4D, 1)
         CALL build(TI, TI4D, 2)
         CALL build(ZEFF_ARR, ZEFF4D, 3)
         CALL build(OMEG_ARR, OMEG4D, 4)
         DO i = 1, NION
            CALL build(NI(i,:,:,:), NI5D(:,:,:,:,i), 4+i)
         END DO
      END IF
      CALL build(B_R, BR4D, 10)
      CALL build(B_PHI, BPHI4D, 11)
      CALL build(B_Z, BZ4D, 12)
      CALL build(MODB, MODB4D, 13)
      CALL build(U_ARR, U4D, 14)
      CALL build(POT_ARR, POT4D, 15)
      CALL build(RHO_ARR, RHO4D, 16)
      CALL build(XRHO_ARR, XRHO4D, 17)
      CALL build(YRHO_ARR, YRHO4D, 18)
      CALL build_cubic(S_ARR, S4D, 19)
      lpot = ANY(POT_ARR /= 0)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_SHARMEM, ier)
      IF (ier /= 0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_spline3d_setup',ier)
#endif

      IF (.not. lvac) THEN
         CALL mpidealloc(TE,win_TE)
         CALL mpidealloc(NE,win_NE)
         CALL mpidealloc(TI,win_TI)
         CALL mpidealloc(ZEFF_ARR,win_ZEFF_ARR)
         CALL mpidealloc(OMEG_ARR,win_OMEG_ARR)
         CALL mpidealloc(NI,win_NI)
      END IF
      CALL mpidealloc(B_R,win_B_R)
      CALL mpidealloc(B_PHI,win_B_PHI)
      CALL mpidealloc(B_Z,win_B_Z)
      CALL mpidealloc(MODB,win_MODB)
      CALL mpidealloc(U_ARR,win_U_ARR)
      CALL mpidealloc(POT_ARR,win_POT_ARR)
      CALL mpidealloc(RHO_ARR,win_RHO_ARR)
      CALL mpidealloc(XRHO_ARR,win_XRHO_ARR)
      CALL mpidealloc(YRHO_ARR,win_YRHO_ARR)
      CALL mpidealloc(S_ARR,win_S_ARR)

      RETURN

      CONTAINS

      SUBROUTINE build(fin, fout, k)
         DOUBLE PRECISION, INTENT(in)    :: fin(nr,nphi,nz)
         DOUBLE PRECISION, INTENT(inout) :: fout(8,nr,nphi,nz)
         INTEGER, INTENT(in)             :: k
         TYPE(EZspline3_r8) :: F_spl
         INTEGER :: ierr
         IF (lcubic_grid) THEN
            CALL build_cubic(fin, fout, k)
            RETURN
         END IF
         IF (MOD(k,nshar) /= myid_sharmem) RETURN
         ierr = 0
         CALL EZspline_init(F_spl,nr,nphi,nz,bcs1,bcs2,bcs3,ierr)
         IF (ierr /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ierr)
         F_spl%isHermite   = 1
         F_spl%x1   = raxis
         F_spl%x2   = phiaxis
         F_spl%x3   = zaxis
         CALL EZspline_setup(F_spl,fin,ierr,EXACT_DIM=.true.)
         IF (ierr /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ierr)
         fout = F_spl%fspl
         CALL EZspline_free(F_spl,ierr)
      END SUBROUTINE build

      ! The flux label of the collision step: a tricubic spline with natural
      ! ends in R and Z and periodic in phi, the interpolant ASCOT5 (B_STS)
      ! and SIMPLE build from the same grid. Unlike the Hermite fits it is C2.
      SUBROUTINE build_cubic(fin, fout, k)
         DOUBLE PRECISION, INTENT(in)    :: fin(nr,nphi,nz)
         DOUBLE PRECISION, INTENT(inout) :: fout(8,nr,nphi,nz)
         INTEGER, INTENT(in)             :: k
         TYPE(EZspline3_r8) :: F_spl
         INTEGER :: ierr
         IF (MOD(k,nshar) /= myid_sharmem) RETURN
         ierr = 0
         CALL EZspline_init(F_spl,nr,nphi,nz,(/2,2/),bcs2,(/2,2/),ierr)
         IF (ierr /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZsplineinit',ierr)
         F_spl%isHermite = 0
         F_spl%x1 = raxis
         F_spl%x2 = phiaxis
         F_spl%x3 = zaxis
         F_spl%bcval1min = 0; F_spl%bcval1max = 0
         F_spl%bcval3min = 0; F_spl%bcval3max = 0
         CALL EZspline_setup(F_spl,fin,ierr,EXACT_DIM=.true.)
         IF (ierr /=0) CALL handle_err(EZSPLINE_ERR,'beams3d_spline3d_setup: EZspline_setup',ierr)
         fout = F_spl%fspl
         CALL EZspline_free(F_spl,ierr)
      END SUBROUTINE build_cubic

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_spline3d_setup
