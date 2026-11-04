program test_wall_planar
    use mpi
    use wall_mod, only: vertex, nface, get_block_size, wall_load_txt, wall_free, collide
    implicit none
    integer :: ierr, world, u
    logical :: hit
    double precision :: size, point(3), fraction

    call MPI_INIT(ierr)
    world = MPI_COMM_WORLD
    allocate(vertex(3,3))
    vertex(1,:) = [1d0, -1d0, -1d0]
    vertex(2,:) = [1d0, 2d0, -1d0]
    vertex(3,:) = [1d0, -1d0, 2d0]
    nface = 1
    call get_block_size(size)
    if (size < 1d0 .or. size > 3.1d0) error stop 'planar grid over-refined'
    deallocate(vertex)
    nullify(vertex)
    open(newunit=u, file='planar_wall.dat', status='replace')
    write(u,'(a)') 'MACHINE: synthetic', 'DATE: synthetic', '3 1'
    write(u,*) 1d0, -1d0, -1d0
    write(u,*) 1d0, 2d0, -1d0
    write(u,*) 1d0, -1d0, 2d0
    write(u,*) 1, 2, 3
    close(u)
    ierr = 0
    call wall_load_txt('planar_wall.dat', ierr, .false., world)
    if (ierr /= 0) error stop 'planar wall read failed'
    call collide(0d0, 0d0, 0d0, 2d0, 0d0, 0d0, &
        point(1), point(2), point(3), hit, fraction)
    if (.not. hit .or. maxval(abs(point-[1d0, 0d0, 0d0])) > 1d-12 &
        .or. abs(fraction-0.5d0) > 1d-12) error stop 'planar wall intersection failed'
    call wall_free(ierr, world)
    call MPI_FINALIZE(ierr)
end program test_wall_planar
