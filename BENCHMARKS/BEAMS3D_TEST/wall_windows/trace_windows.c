/* Diagnostic only: forward the native GNU/Open MPI Fortran binding unchanged. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <mpi.h>
#include "ompi/communicator/communicator.h"
#include "ompi/mca/pml/pml.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef void (*allocate_shared_f)(MPI_Aint *, MPI_Fint *, MPI_Fint *, MPI_Fint *,
                                  char *, MPI_Fint *, MPI_Fint *);
typedef int (*comm_dup_c)(ompi_communicator_t *, ompi_communicator_t **);
static uint64_t window_attempt;

int ompi_comm_dup(ompi_communicator_t *comm, ompi_communicator_t **newcomm)
{
    static comm_dup_c original;
    int result;
    if (!original) {
        original = (comm_dup_c)dlsym(RTLD_NEXT, "ompi_comm_dup");
        if (!original) {
            fprintf(stderr, "window profiler could not resolve communicator binding\n");
            exit(126);
        }
    }
    result = original(comm, newcomm);
    if (result != MPI_SUCCESS) {
        opal_pointer_array_t *table = dlsym(RTLD_DEFAULT, "ompi_mpi_communicators");
        mca_pml_base_module_t *pml = dlsym(RTLD_DEFAULT, "mca_pml");
        fprintf(stderr, "COMM_DUP_FAILURE WINDOW_ATTEMPT %llu RETURN %d LOWEST_CID %d CID_LIMIT %u\n",
                (unsigned long long)window_attempt, result,
                table ? table->lowest_free : -1, pml ? pml->pml_max_contextid : 0);
        fflush(stderr);
    }
    return result;
}

void mpi_win_allocate_shared_cptr_(MPI_Aint *size, MPI_Fint *disp, MPI_Fint *info,
                              MPI_Fint *comm, char *base, MPI_Fint *win,
                              MPI_Fint *error)
{
    static allocate_shared_f original;
    static opal_pointer_array_t *communicators;
    static mca_pml_base_module_t *pml;
    if (!original) {
        original = (allocate_shared_f)dlsym(RTLD_NEXT, "mpi_win_allocate_shared_cptr_");
        if (!original) {
            fprintf(stderr, "window profiler could not resolve native binding\n");
            exit(126);
        }
        communicators = dlsym(RTLD_DEFAULT, "ompi_mpi_communicators");
        pml = dlsym(RTLD_DEFAULT, "mca_pml");
    }
    ++window_attempt;
    fprintf(stderr, "WINDOW_ATTEMPT %llu BYTES %lld COMM %d LOWEST_CID %d CID_LIMIT %u\n",
            (unsigned long long)window_attempt, (long long)*size, (int)*comm,
            communicators ? communicators->lowest_free : -1,
            pml ? pml->pml_max_contextid : 0);
    fflush(stderr);
    original(size, disp, info, comm, base, win, error);
    if (*error != MPI_SUCCESS) {
        fprintf(stderr, "WINDOW_RETURN %llu ERROR %d\n",
                (unsigned long long)window_attempt, (int)*error);
        fflush(stderr);
    }
}
