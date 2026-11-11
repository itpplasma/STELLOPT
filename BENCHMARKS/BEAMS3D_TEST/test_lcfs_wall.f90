program test_lcfs_wall
    ! The faceted boundary that a plasma-only run uses as its loss surface must
    ! follow the boundary closely: for a circular torus of W7-X size the facets
    ! built with the default resolution lie within 0.05 mm of the surface. With
    ! the former 120 poloidal and 180 toroidal points they lay 0.9 mm inside.
    use mpi
    use beams3d_globals, only: lcfs_nu, lcfs_nv
    use beams3d_input_mod, only: init_beams3d_input
    use wall_mod, only: vertex, face, nface, wall_load_mn, wall_free
    implicit none
    double precision, parameter :: r0 = 5.5d0, a = 0.53d0
    double precision :: rmn(2), zmn(2), xm(2), xn(2), c(3), dev, coarse
    integer :: ierr, world, i

    call MPI_INIT(ierr)
    world = MPI_COMM_WORLD
    call init_beams3d_input
    rmn = [r0, a]; zmn = [0d0, a]; xm = [0d0, 1d0]; xn = 0d0
    coarse = deviation(120, 180)
    dev = deviation(lcfs_nu, lcfs_nv)
    write(*,'(a,es9.2,a,es9.2,a)') ' Faceted boundary inside the torus by ', dev, &
        ' m (120 x 180 points: ', coarse, ' m)'
    if (coarse < 5d-4) error stop 'the coarse mesh should miss the surface by 0.9 mm'
    if (dev > 5d-5) error stop 'default boundary mesh further than 0.05 mm from the surface'
    write(*,*) 'Faceted boundary of plasma-only runs: PASS'
    call MPI_FINALIZE(ierr)
contains
    double precision function deviation(nu, nv)
        ! Largest distance of a facet centroid from the torus.
        integer, intent(in) :: nu, nv
        call wall_load_mn(rmn, zmn, xm, xn, 2, nu, nv, .false., world)
        deviation = 0
        do i = 1, nface
            c = (vertex(face(i,1),:) + vertex(face(i,2),:) + vertex(face(i,3),:))/3
            deviation = max(deviation, abs(hypot(hypot(c(1), c(2)) - r0, c(3)) - a))
        end do
        call wall_free(ierr, world)
    end function deviation
end program test_lcfs_wall
