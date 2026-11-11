PROGRAM test_collision_noise
    USE beams3d_physics_mod, ONLY: collision_noise_normal
    IMPLICIT NONE
    ! Fixture shared with ASCOT5 test_collision_noise.c and SIMPLE
    ! test_collision_noise.f90: the same keys must give the same numbers.
    INTEGER(8), PARAMETER :: key(4,4) = RESHAPE([0_8, 1_8, 0_8, 3_8, &
        20261110_8, 1_8, 0_8, 4_8, 20261110_8, 16384_8, 1000000_8, 3_8, &
        7_8, 42_8, 123456789_8, 0_8], [4,4])
    DOUBLE PRECISION, PARAMETER :: expected(4) = [1.0115071454929696D0, &
        -0.70803017369489951D0, 1.1794644411780411D0, -0.65541438399747554D0]
    INTEGER, PARAMETER :: n = 200000
    DOUBLE PRECISION :: x, a, b, s1, s2, c
    INTEGER :: i, failures
    INTEGER(8) :: k

    failures = 0
    DO i = 1, 4
        x = collision_noise_normal(key(1,i), key(2,i), key(3,i), key(4,i))
        IF (ABS(x - expected(i)) > 4D-15) THEN
            WRITE(*,'(A,I0,2ES25.17)') 'fixture ', i, x, expected(i)
            failures = failures + 1
        END IF
    END DO
    s1 = 0; s2 = 0; c = 0
    DO k = 0, n - 1
        a = collision_noise_normal(99_8, 5_8, k, 3_8)
        b = collision_noise_normal(99_8, 5_8, k, 4_8)
        s1 = s1 + a; s2 = s2 + a*a; c = c + a*b
    END DO
    s1 = s1/n; s2 = s2/n - s1*s1; c = c/n
    IF (ABS(s1) > 1.2D-2 .or. ABS(s2 - 1) > 2D-2 .or. ABS(c) > 1.2D-2) THEN
        WRITE(*,*) 'moments mean, var-1, corr: ', s1, s2 - 1, c
        failures = failures + 1
    END IF
    IF (failures /= 0) ERROR STOP 1
    WRITE(*,*) 'Counter-based collision noise: PASS'
END PROGRAM test_collision_noise
