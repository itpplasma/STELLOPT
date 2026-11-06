program test_wall_interior
    ! A segment that crosses a large triangle far from all of its edges must
    ! hit it, also when the acceleration blocks are much smaller than the face.
    use mpi
    use wall_mod, only: wall_load_txt, wall_free, collide, get_wall_ik, set_nface
    implicit none
    integer, parameter :: ntiny = 4000
    integer :: ierr, world, u, k
    logical :: hit
    double precision :: point(3), fraction, x, y

    call MPI_INIT(ierr)
    world = MPI_COMM_WORLD
    open(newunit=u, file='interior_wall.dat', status='replace')
    write(u,'(a)') 'MACHINE: synthetic', 'DATE: synthetic'
    write(u,*) 3*(ntiny+1), ntiny+1
    write(u,*) 0d0, -20d0, -20d0
    write(u,*) 0d0, 60d0, -20d0
    write(u,*) 0d0, -20d0, 60d0
    ! Small faces away from the segment refine the block grid to a few metres.
    do k = 1, ntiny
        x = -40d0 + 80d0*dble(mod(k-1, 80))/79d0
        y = 45d0 + 15d0*dble((k-1)/80)/49d0
        write(u,*) x, y, 55d0
        write(u,*) x+1d-2, y, 55d0
        write(u,*) x, y+1d-2, 55d0
    end do
    do k = 0, ntiny
        write(u,*) 3*k+1, 3*k+2, 3*k+3
    end do
    close(u)
    call set_nface(1)
    ierr = 0
    call wall_load_txt('interior_wall.dat', ierr, .false., world)
    if (ierr /= 0) error stop 'interior wall read failed'
    call collide(-1d0, 10d0, -5d0, 1d0, 10d0, -5d0, &
        point(1), point(2), point(3), hit, fraction)
    if (.not. hit .or. get_wall_ik() /= 1) error stop 'interior of large face missed'
    if (maxval(abs(point-[0d0, 10d0, -5d0])) > 1d-12 .or. abs(fraction-0.5d0) > 1d-12) &
        error stop 'interior hit misplaced'
    call wall_free(ierr, world)
    call MPI_FINALIZE(ierr)
end program test_wall_interior
