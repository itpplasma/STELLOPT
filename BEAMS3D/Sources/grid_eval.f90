!-----------------------------------------------------------------------
!     Subroutine:    grid_eval
!     Description:   R8HERM3FCN with the same arguments, or the C2 tricubic
!                    spline evaluation when the grids hold spline
!                    coefficients (SPLINE_TYPE = 'CUBIC'). ict selects f,
!                    df/dR, df/dphi and df/dZ. External, like the pspline
!                    routines, so that callers may pass grid elements.
!-----------------------------------------------------------------------
      SUBROUTINE grid_eval(ict, ivec, ivecd, fval, ii, jj, kk, xparam, yparam, zparam, &
                           hx, hxi, hy, hyi, hz, hzi, fin, inf2, inf3, nzd)
      USE beams3d_grid, ONLY: lcubic_grid
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ict(8), ivec, ivecd, inf2, inf3, nzd
      INTEGER, INTENT(IN) :: ii(ivec), jj(ivec), kk(ivec)
      REAL*8, INTENT(IN)  :: xparam(ivec), yparam(ivec), zparam(ivec), &
                             hx(ivec), hxi(ivec), hy(ivec), hyi(ivec), hz(ivec), hzi(ivec)
      REAL*8, INTENT(IN)  :: fin(0:7,inf2,inf3,nzd)
      REAL*8, INTENT(OUT) :: fval(ivecd,*)
      INTEGER :: ict10(10)
      IF (lcubic_grid) THEN
         ict10 = 0
         ict10(1:4) = ict(1:4)
         CALL R8FVTRICUB(ict10, ivec, ivecd, fval, ii, jj, kk, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, fin, inf2, inf3, nzd)
      ELSE
         CALL R8HERM3FCN(ict, ivec, ivecd, fval, ii, jj, kk, xparam, yparam, zparam, &
                         hx, hxi, hy, hyi, hz, hzi, fin, inf2, inf3, nzd)
      END IF
      END SUBROUTINE grid_eval
