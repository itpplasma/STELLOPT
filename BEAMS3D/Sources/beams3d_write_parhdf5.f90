!-----------------------------------------------------------------------
!     Module:        beams3d_write_parhdf5
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/03/2017
!     Description:   Write trajectory info using parallel HDF5. Note,
!                    you must be on a parallel file system for this to
!                    work.
!-----------------------------------------------------------------------
      MODULE beams3d_write_par
      CONTAINS
      SUBROUTINE beams3d_write_parhdf5(n1,n2,m1,m2,mystart,myend,var_name,INTVAR,FLTVAR,DBLVAR,DESCRIPTION)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mpi_params
      USE hdf5
      USE ez_hdf5, ONLY: write_att_hdf5
      USE beams3d_runtime, ONLY: id_string, nprocs_beams, handle_err, &
                                 MPI_BARRIER_ERR
      USE mpi_inc
!-----------------------------------------------------------------------
!     Input Variables
!          var_name     Name of variable
!          n1,n2        Extent of first dimesion of global array
!          m1,m2        Extent of second dimension of global array
!          mystart      Starting point of data (2nd dim)
!          myend        Ending point of data (2nd dim)
!          INTVAR       Array (INTEGER)
!          FLTVAR       Array (REAL)
!          DBLVAR       Array (DOUBLE PRCISION)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER                                :: n1,n2,m1,m2,mystart,myend
      CHARACTER(LEN=*), INTENT(in)           :: var_name
      INTEGER, INTENT(in),  OPTIONAL         :: INTVAR(n1:n2,mystart:myend)
      REAL, INTENT(in),  OPTIONAL            :: FLTVAR(n1:n2,mystart:myend)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: DBLVAR(n1:n2,mystart:myend)
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: DESCRIPTION
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          info         ?
!          plist_id     Property list
!          file_id      HDF5 file ID
!          rank         Data space rank
!          dimsf        Data space dimensions
!          fspace_id    Data space identifier
!          dset_id      Dataset identifier
!          mspace_id    Memory identifier
!          dataspace    Data identifier
!-----------------------------------------------------------------------
      LOGICAL :: livar, lfvar, ldvar
      INTEGER :: ier, info, rank, i
      INTEGER(HID_T) :: file_id, fspace_id, dset_id, mspace_id, &
                        fapl_id, dcpl_id, dxpl_id, driver_id
      INTEGER(HSIZE_T), ALLOCATABLE :: dimsf(:), counts(:), chunk_dims(:),&
                                       offset(:)

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Handle types
      livar = .false.; lfvar=.false.; ldvar=.false.
      IF (PRESENT(INTVAR)) livar=.true.
      IF (PRESENT(FLTVAR)) lfvar=.true.
      IF (PRESENT(DBLVAR)) ldvar=.true.

      ! Setup Helper Arrays
      rank =2;
      ALLOCATE(dimsf(rank),chunk_dims(rank),counts(rank),offset(rank))
      dimsf(1) = n2-n1+1
      dimsf(2) = m2-m1+1
      chunk_dims(1) = n2-n1+1
      chunk_dims(2) = myend-mystart+1
      counts(1) = 1
      counts(2) = 1
      offset(1) = 0
      offset(2) = mystart-1

#if defined(HDF5_PAR)
      ! Do this so we define the chunking correctly
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,chunk_dims(2),1,MPI_INTEGER,MPI_MAX,MPI_COMM_BEAMS,ier)

      ier = 0
      info = MPI_INFO_NULL
      ! Initialize
      CALL h5open_f(ier)
      ! Setup File access
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)
      CALL h5pset_fapl_mpio_f(fapl_id, MPI_COMM_BEAMS, info, ier)
      CALL h5pget_driver_f(fapl_id, driver_id, ier)
      ! Open file
      CALL h5fopen_f('beams3d_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = fapl_id)
!!!!!!! Begin writing

!      WRITE(6,*) myworkid,dimsf,chunk_dims,counts,offset; CALL FLUSH(6)

      ! Create Data Space
      CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)

      ! Enable Chunking
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
      CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)

      ! Create Dataset
      IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier, dcpl_id)
      IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)
      IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)
      IF (PRESENT(DESCRIPTION)) CALL write_att_hdf5(dset_id,'description',DESCRIPTION,ier)

      ! Close the file space
      CALL h5sclose_f(fspace_id, ier)

      ! Create the Memory Space
      chunk_dims(2) = myend-mystart+1 ! So we write correctly
      CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)

      ! Select the Hyperslab in data
      CALL h5dget_space_f(dset_id, fspace_id, ier)
      CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, chunk_dims, ier)

      ! Create Properties
      CALL h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ier)
      CALL h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ier)

      ! Write dataset
      IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)

      ! Close Property list
      CALL h5pclose_f(fapl_id, ier)
      CALL h5pclose_f(dcpl_id, ier)
      CALL h5pclose_f(dxpl_id, ier)
      CALL h5sclose_f(mspace_id, ier)
      CALL h5sclose_f(fspace_id, ier)
      CALL h5dclose_f(dset_id, ier)

!!!!!!!CLOSE FILE
      ! Close the file
      CALL h5fclose_f(file_id, ier)
      ! Close the fortran interface
      CALL h5close_f(ier)

#else

#if defined(MPI_OPT)
      IF (.not. lfvar .and. INT(n2-n1+1,8)*INT(m2-m1+1,8) < HUGE(i)) THEN
         CALL write_gathered(rank, dimsf, n2-n1+1, m1, mystart, myend, &
                             'beams3d_'//TRIM(id_string)//'.h5', var_name, &
                             livar, INTVAR, DBLVAR, DESCRIPTION)
         DEALLOCATE(dimsf,chunk_dims,counts,offset)
         RETURN
      END IF
#endif

      DO i = 0, nprocs_beams-1
         IF (myworkid == i) THEN
            ! Open the fotran interface
            CALL h5open_f(ier)
            !PRINT *,'h5open ',ier

            ! Setup File access
            CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)
            !PRINT *,'h5pcreate_f ',ier

            ! Open file
            CALL h5fopen_f('beams3d_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = H5P_DEFAULT_F)
            !PRINT *,'h5fopen_f ',ier

            ! Open or create the dataset and get/create dataspace identifer
            IF  (myworkid == master) THEN
               CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)
            !PRINT *,'h5screate_simple_f ',ier
               IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier)
               IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
               IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
               IF (PRESENT(DESCRIPTION)) CALL write_att_hdf5(dset_id,'description',DESCRIPTION,ier)
            !PRINT *,'h5dcreate_f ',ier
            ELSE
               CALL h5dopen_f(file_id, TRIM(var_name), dset_id, ier)
            !PRINT *,'h5dopen_f ',ier
               CALL h5dget_space_f(dset_id, fspace_id, ier)
            !PRINT *,'h5dget_space_f ',ier
            END IF

            ! Select Hyperslab
            CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F,offset, chunk_dims, ier) 
            !PRINT *,'h5sselect_hyperslab_f ',ier
 
            ! Get Memory Space
            CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)
            !PRINT *,'h5screate_simple_f ',ier
            
            ! Enable Chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
            !PRINT *,'h5pcreate_f ',ier
            CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)
            !PRINT *,'h5pset_chunk_f ',ier

            ! Write dataset
            IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            !PRINT *,'h5dwrite_f ',ier
            
            ! Close down
            CALL h5sclose_f(fspace_id, ier)
            CALL h5sclose_f(mspace_id, ier)
            CALL h5dclose_f(dset_id, ier)
      
            ! Close the file
            CALL h5fclose_f(file_id, ier)
      
            ! Close the fortran interface
            CALL h5close_f(ier)
         END IF
         CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_write_parhdf5',ierr_mpi)
         !IF (i == 1) STOP
      END DO
#endif

      ! Deallocate Helpers
      DEALLOCATE(dimsf,chunk_dims,counts,offset)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_write_parhdf5

      SUBROUTINE beams3d_write1d_parhdf5(m1,m2,mystart,myend,var_name,INTVAR,FLTVAR,DBLVAR,FILENAME,DESCRIPTION)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE mpi_params
      USE hdf5
      USE ez_hdf5, ONLY: write_att_hdf5
      USE beams3d_runtime, ONLY: id_string, nprocs_beams, handle_err, &
                                 MPI_BARRIER_ERR
      USE mpi_inc
!-----------------------------------------------------------------------
!     Input Variables
!          var_name     Name of variable
!          n1,n2        Extent of first dimesion of global array
!          m1,m2        Extent of second dimension of global array
!          mystart      Starting point of data (2nd dim)
!          myend        Ending point of data (2nd dim)
!          INTVAR       Array (INTEGER)
!          FLTVAR       Array (REAL)
!          DBLVAR       Array (DOUBLE PRCISION)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER                                :: m1,m2,mystart,myend
      CHARACTER(LEN=*), INTENT(in)           :: var_name
      INTEGER, INTENT(in),  OPTIONAL         :: INTVAR(mystart:myend)
      REAL, INTENT(in),  OPTIONAL            :: FLTVAR(mystart:myend)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: DBLVAR(mystart:myend)
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: FILENAME
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: DESCRIPTION
!-----------------------------------------------------------------------
!     Local Variables
!          ier          Error Flag
!          info         ?
!          plist_id     Property list
!          file_id      HDF5 file ID
!          rank         Data space rank
!          dimsf        Data space dimensions
!          fspace_id    Data space identifier
!          dset_id      Dataset identifier
!          mspace_id    Memory identifier
!          dataspace    Data identifier
!-----------------------------------------------------------------------
      LOGICAL :: livar, lfvar, ldvar
      INTEGER :: ier, info, rank, i
      INTEGER(HID_T) :: file_id, fspace_id, dset_id, mspace_id, &
                        fapl_id, dcpl_id, dxpl_id, driver_id
      INTEGER(HSIZE_T), ALLOCATABLE :: dimsf(:), counts(:), chunk_dims(:),&
                                       offset(:)

!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Handle types
      livar = .false.; lfvar=.false.; ldvar=.false.
      IF (PRESENT(INTVAR)) livar=.true.
      IF (PRESENT(FLTVAR)) lfvar=.true.
      IF (PRESENT(DBLVAR)) ldvar=.true.

      ! Setup Helper Arrays
      rank =1
      ALLOCATE(dimsf(rank),chunk_dims(rank),counts(rank),offset(rank))
      dimsf(1) = m2-m1+1
      chunk_dims(1) = myend-mystart+1
      counts(1) = 1
      offset(1) = mystart-1

#if defined(HDF5_PAR)
      ! Do this so we define the chunking correctly
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,chunk_dims(1),1,MPI_INTEGER,MPI_MAX,MPI_COMM_BEAMS,ier)
      ier = 0
      info = MPI_INFO_NULL
      ! Initialize
      CALL h5open_f(ier)
      ! Setup File access
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)
      CALL h5pset_fapl_mpio_f(fapl_id, MPI_COMM_BEAMS, info, ier)
      CALL h5pget_driver_f(fapl_id, driver_id, ier)
      ! Open file
      IF (PRESENT(FILENAME)) THEN
         CALL h5fopen_f(TRIM(FILENAME)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = fapl_id)
      ELSE
         CALL h5fopen_f('ascot5_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = fapl_id)
      END IF
!!!!!!! Begin writing

!      WRITE(6,*) myworkid,dimsf,chunk_dims,counts,offset; CALL FLUSH(6)

      ! Create Data Space
      CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)

      ! Enable Chunking
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
      CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)

      ! Create Dataset
      IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier, dcpl_id)
      IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)
      IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier, dcpl_id)
      IF (PRESENT(DESCRIPTION)) CALL write_att_hdf5(dset_id,'description',DESCRIPTION,ier)

      ! Close the file space
      CALL h5sclose_f(fspace_id, ier)

      ! Create the Memore Space
      chunk_dims(1) = myend-mystart+1
      CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)

      ! Select the Hyperslab in data
      CALL h5dget_space_f(dset_id, fspace_id, ier)
      CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, chunk_dims, ier)

      ! Create Properties
      CALL h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ier)
      CALL h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ier)

      ! Write dataset
      IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)
      IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id, xfer_prp = dxpl_id)

      ! Close Property list
      CALL h5pclose_f(fapl_id, ier)
      CALL h5pclose_f(dcpl_id, ier)
      CALL h5pclose_f(dxpl_id, ier)
      CALL h5sclose_f(mspace_id, ier)
      CALL h5sclose_f(fspace_id, ier)
      CALL h5dclose_f(dset_id, ier)

!!!!!!!CLOSE FILE
      ! Close the file
      CALL h5fclose_f(file_id, ier)
      ! Close the fortran interface
      CALL h5close_f(ier)

#else

#if defined(MPI_OPT)
      IF (.not. lfvar) THEN
         IF (PRESENT(FILENAME)) THEN
            CALL write_gathered(rank, dimsf, 1, m1, mystart, myend, &
                                TRIM(FILENAME)//'.h5', var_name, &
                                livar, INTVAR, DBLVAR, DESCRIPTION)
         ELSE
            CALL write_gathered(rank, dimsf, 1, m1, mystart, myend, &
                                'ascot5_'//TRIM(id_string)//'.h5', var_name, &
                                livar, INTVAR, DBLVAR, DESCRIPTION)
         END IF
         DEALLOCATE(dimsf,chunk_dims,counts,offset)
         RETURN
      END IF
#endif

      DO i = 0, nprocs_beams-1
         IF (myworkid == i) THEN
            ! Open the fotran interface
            CALL h5open_f(ier)
            !PRINT *,'h5open ',ier

            ! Setup File access
            CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ier)
            !PRINT *,'h5pcreate_f ',ier

            ! Open file
            IF (PRESENT(FILENAME)) THEN
               CALL h5fopen_f(TRIM(FILENAME)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = H5P_DEFAULT_F)
            ELSE
               CALL h5fopen_f('ascot5_'//TRIM(id_string)//'.h5', H5F_ACC_RDWR_F, file_id, ier, access_prp = H5P_DEFAULT_F)
            END IF
            !PRINT *,'h5fopen_f ',ier

            ! Open or create the dataset and get/create dataspace identifer
            IF  (myworkid == master) THEN
               CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)
            !PRINT *,'h5screate_simple_f ',ier
               IF (livar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier)
               IF (lfvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
               IF (ldvar) CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
               IF (PRESENT(DESCRIPTION)) CALL write_att_hdf5(dset_id,'description',DESCRIPTION,ier)
            !PRINT *,'h5dcreate_f ',ier
            ELSE
               CALL h5dopen_f(file_id, TRIM(var_name), dset_id, ier)
            !PRINT *,'h5dopen_f ',ier
               CALL h5dget_space_f(dset_id, fspace_id, ier)
            !PRINT *,'h5dget_space_f ',ier
            END IF

            ! Select Hyperslab
            CALL h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F,offset, chunk_dims, ier) 
            !PRINT *,'h5sselect_hyperslab_f ',ier
 
            ! Get Memory Space
            CALL h5screate_simple_f(rank, chunk_dims, mspace_id, ier)
            !PRINT *,'h5screate_simple_f ',ier
            
            ! Enable Chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ier)
            !PRINT *,'h5pcreate_f ',ier
            CALL h5pset_chunk_f(dcpl_id, rank, chunk_dims, ier)
            !PRINT *,'h5pset_chunk_f ',ier

            ! Write dataset
            IF (livar) CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (lfvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, FLTVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            IF (ldvar) CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DBLVAR, dimsf, ier, mem_space_id = mspace_id, file_space_id = fspace_id)
            !PRINT *,'h5dwrite_f ',ier
            
            ! Close down
            CALL h5sclose_f(fspace_id, ier)
            CALL h5sclose_f(mspace_id, ier)
            CALL h5dclose_f(dset_id, ier)
      
            ! Close the file
            CALL h5fclose_f(file_id, ier)
      
            ! Close the fortran interface
            CALL h5close_f(ier)
         END IF
         CALL MPI_BARRIER(MPI_COMM_BEAMS,ierr_mpi)
         IF (ierr_mpi /=0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_write_parhdf5',ierr_mpi)
         !IF (i == 1) STOP
      END DO

#endif

      ! Deallocate Helpers
      DEALLOCATE(dimsf,chunk_dims,counts,offset)
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE beams3d_write1d_parhdf5

#if defined(MPI_OPT)
      SUBROUTINE write_gathered(rank,dimsf,ncol,m1,mystart,myend,fname,var_name, &
                                livar,INTVAR,DBLVAR,DESCRIPTION)
      !-----------------------------------------------------------------
      ! Rank by rank, each rank opening and closing the file in turn,
      ! costs seconds per dataset on network file systems. Every rank
      ! holds the contiguous columns mystart:myend, so gather them on
      ! master and write the dataset with one open.
      !-----------------------------------------------------------------
      USE mpi_params
      USE hdf5
      USE ez_hdf5, ONLY: write_att_hdf5
      USE beams3d_runtime, ONLY: nprocs_beams, handle_err, MPI_BARRIER_ERR
      USE mpi_inc
      IMPLICIT NONE
      INTEGER, INTENT(in)                    :: rank, ncol, m1, mystart, myend
      INTEGER(HSIZE_T), INTENT(in)           :: dimsf(rank)
      CHARACTER(LEN=*), INTENT(in)           :: fname, var_name
      LOGICAL, INTENT(in)                    :: livar
      INTEGER, INTENT(in), OPTIONAL          :: INTVAR(*)
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: DBLVAR(*)
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: DESCRIPTION
      INTEGER :: ier, i, nsend, myrange(2)
      INTEGER, ALLOCATABLE :: ranges(:,:), rcounts(:), displs(:)
      INTEGER, ALLOCATABLE :: ifull(:)
      DOUBLE PRECISION, ALLOCATABLE :: dfull(:)
      INTEGER(HID_T) :: file_id, fspace_id, dset_id

      nsend = MAX(myend-mystart+1,0)*ncol
      myrange = (/ mystart, myend /)
      ALLOCATE(ranges(2,nprocs_beams), rcounts(nprocs_beams), displs(nprocs_beams))
      CALL MPI_GATHER(myrange, 2, MPI_INTEGER, ranges, 2, MPI_INTEGER, master, MPI_COMM_BEAMS, ier)
      rcounts = 0; displs = 0
      IF (myworkid == master) THEN
         DO i = 1, nprocs_beams
            rcounts(i) = MAX(ranges(2,i)-ranges(1,i)+1,0)*ncol
            IF (rcounts(i) > 0) displs(i) = (ranges(1,i)-m1)*ncol
         END DO
      END IF
      IF (livar) THEN
         ALLOCATE(ifull(PRODUCT(dimsf)))
         ifull = 0
         CALL MPI_GATHERV(INTVAR, nsend, MPI_INTEGER, ifull, rcounts, displs, &
                          MPI_INTEGER, master, MPI_COMM_BEAMS, ier)
      ELSE
         ALLOCATE(dfull(PRODUCT(dimsf)))
         dfull = 0
         CALL MPI_GATHERV(DBLVAR, nsend, MPI_DOUBLE_PRECISION, dfull, rcounts, displs, &
                          MPI_DOUBLE_PRECISION, master, MPI_COMM_BEAMS, ier)
      END IF
      IF (myworkid == master) THEN
         CALL h5open_f(ier)
         CALL h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, ier, access_prp = H5P_DEFAULT_F)
         CALL h5screate_simple_f(rank, dimsf, fspace_id, ier)
         IF (livar) THEN
            CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_INTEGER, fspace_id, dset_id, ier)
         ELSE
            CALL h5dcreate_f(file_id, TRIM(var_name), H5T_NATIVE_DOUBLE, fspace_id, dset_id, ier)
         END IF
         IF (PRESENT(DESCRIPTION)) CALL write_att_hdf5(dset_id,'description',DESCRIPTION,ier)
         IF (livar) THEN
            CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ifull, dimsf, ier)
         ELSE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dfull, dimsf, ier)
         END IF
         CALL h5sclose_f(fspace_id, ier)
         CALL h5dclose_f(dset_id, ier)
         CALL h5fclose_f(file_id, ier)
         CALL h5close_f(ier)
      END IF
      CALL MPI_BARRIER(MPI_COMM_BEAMS, ier)
      IF (ier /= 0) CALL handle_err(MPI_BARRIER_ERR,'beams3d_write_parhdf5',ier)
      END SUBROUTINE write_gathered
#endif
END MODULE beams3d_write_par
