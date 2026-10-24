! Existing native reader and exact collision test, no orbit or field.
program native_wall_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mpi
    use wall_mod, only: wall_load_txt, wall_free, collide, get_wall_ik, nface, nvertex
    implicit none
    character(len=512) :: filename, output, argument, cycle_output
    character(len=32) :: suffix
    integer :: status, unit_in, unit_out, n, row, facet, world
    integer :: rank, ranks, repeat, cycle
    logical :: hit
    real(dp) :: p0(3), p1(3), point(3), fraction
    call MPI_INIT(status)
    world=MPI_COMM_WORLD
    call MPI_COMM_RANK(world,rank,status)
    call MPI_COMM_SIZE(world,ranks,status)
    call get_command_argument(1,filename)
    call get_command_argument(2,output)
    repeat=1
    if (command_argument_count() >= 3) then
        call get_command_argument(3,argument)
        read(argument,*) repeat
    end if
    if (repeat < 1 .or. repeat > 3) error stop 'Invalid repeat count'
    if (ranks > 1) then
        write(suffix,'(".rank",i0)') rank
        output=trim(output)//trim(suffix)
    end if
    do cycle=1,repeat
        status=0
        if (storage_size(nface) /= 32) error stop 'Native INTEGER range contract'
        call wall_load_txt(trim(filename),status,.true.,world)
        if (status /= 0) error stop 'Native wall input failed'
        open(newunit=unit_in,file='rays.dat',status='old',action='read')
        read(unit_in,*) n
        cycle_output=output
        if (repeat > 1) then
            write(suffix,'(".cycle",i0)') cycle
            cycle_output=trim(output)//trim(suffix)
        end if
        open(newunit=unit_out,file=trim(cycle_output),status='replace',action='write')
        write(unit_out,'(3i14)') nface,nvertex,n
        do row=1,n
            read(unit_in,*) p0,p1
            call collide(p0(1),p0(2),p0(3),p1(1),p1(2),p1(3), &
                         point(1),point(2),point(3),hit,fraction)
            facet=get_wall_ik()
            write(unit_out,'(3i12,4es26.17)') row,merge(1,0,hit),facet,point,fraction
            flush(unit_out)
        end do
        close(unit_in)
        close(unit_out)
        call wall_free(status,world)
        if(status /= 0) error stop 'Native wall cleanup failed'
        write(*,'(a,i0,a,i0)') 'WALL_CYCLE ',cycle,' RANK ',rank
    end do
    call MPI_FINALIZE(status)
end program native_wall_probe
