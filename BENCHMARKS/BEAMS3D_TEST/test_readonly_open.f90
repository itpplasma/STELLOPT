PROGRAM test_readonly_open
    ! An input file must open without write access to it: a write-protected
    ! file opens with LREADONLY and its contents read back. The default
    ! read-write open of the same file fails, which is also what a second
    ! job meets while a first one holds the file open for writing.
    USE ez_hdf5
    IMPLICIT NONE
    INTEGER :: ier,n
    CHARACTER(LEN=*), PARAMETER :: fname = 'test_readonly_open.h5'

    CALL EXECUTE_COMMAND_LINE('rm -f '//fname)
    CALL open_hdf5(fname,fid,ier,LCREATE=.true.)
    IF (ier /= 0) ERROR STOP 'cannot create test file'
    CALL write_scalar_hdf5(fid,'value',ier,INTVAR=42)
    IF (ier /= 0) ERROR STOP 'cannot write test file'
    CALL close_hdf5(fid,ier)
    CALL EXECUTE_COMMAND_LINE('chmod a-w '//fname)

    n = 0
    CALL open_hdf5(fname,fid,ier,LCREATE=.false.,LREADONLY=.true.)
    IF (ier /= 0) THEN
        CALL EXECUTE_COMMAND_LINE('rm -f '//fname)
        ERROR STOP 'write-protected file did not open read-only'
    END IF
    CALL read_scalar_hdf5(fid,'value',ier,INTVAR=n)
    CALL close_hdf5(fid,ier)
    CALL EXECUTE_COMMAND_LINE('rm -f '//fname)
    IF (n /= 42) ERROR STOP 'wrong value read from read-only file'
    WRITE(*,*) 'Read-only open of input files: PASS'
END PROGRAM test_readonly_open
