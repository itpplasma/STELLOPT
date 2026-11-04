program test_wall_grazing
    use mpi
    use wall_mod, only: wall_load_txt, wall_free, collide, get_wall_ik
    implicit none
    integer :: ierr, world, u, shift, direction
    logical :: hit
    double precision :: p0(3), p1(3), point(3), fraction, expected, x

    call MPI_INIT(ierr)
    world = MPI_COMM_WORLD
    p0 = [0.9999999999000001d0, -0.5d0, 0d0]
    p1 = [1.0000000001d0, 0.5d0, 0d0]
    expected = -0.5d0 + (1d0-p0(1))/(p1(1)-p0(1))
    do shift = 0, 1
        x = 1d0-shift
        open(newunit=u, file='grazing_wall.dat', status='replace')
        write(u,'(a)') 'MACHINE: synthetic', 'DATE: synthetic', '6 2'
        write(u,*) x, -1d0, -1d0
        write(u,*) x, 2d0, -1d0
        write(u,*) x, -1d0, 2d0
        write(u,*) x+1d0, -1d0, -1d0
        write(u,*) x+1d0, 2d0, -1d0
        write(u,*) x+1d0, -1d0, 2d0
        write(u,*) 1, 2, 3
        write(u,*) 4, 5, 6
        close(u)
        ierr = 0
        call wall_load_txt('grazing_wall.dat', ierr, .false., world)
        if (ierr /= 0) error stop 'wall reader failed'
        do direction = 1, 2
            if (direction == 1) then
                call collide(p0(1)-shift, p0(2), p0(3), p1(1)-shift, p1(2), p1(3), &
                    point(1), point(2), point(3), hit, fraction)
            else
                call collide(p1(1)-shift, p1(2), p1(3), p0(1)-shift, p0(2), p0(3), &
                    point(1), point(2), point(3), hit, fraction)
            end if
            if (.not. hit .or. get_wall_ik() /= 1) error stop 'grazing wall missed'
            if (abs(point(1)-x) > 1d-12 .or. abs(point(2)-expected) > 1d-12 &
                .or. abs(point(3)) > 1d-12) then
                print *, 'Grazing point:', point, 'expected y:', expected
                error stop 'grazing plane cancellation'
            end if
        end do
        call wall_free(ierr, world)
    end do
    call MPI_FINALIZE(ierr)
end program test_wall_grazing
