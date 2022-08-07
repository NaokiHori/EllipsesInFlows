#if !defined(PARALLEL_H)
#define PARALLEL_H

#include <stddef.h>
#include <mpi.h>

#include "structure.h"

typedef struct {
  int *sendcounts, *recvcounts;
  int *sdispls, *rdispls;
  MPI_Datatype *sendtypes, *recvtypes, *temptypes;
} parallel_transpose_t;

/* ! definition of a structure parallel_t_ ! 4 ! */
struct parallel_t_ {
  int mpisize, mpirank;
  int ymrank, yprank;
};

extern parallel_t *parallel_init(void);
extern int parallel_finalise(parallel_t *parallel);
extern int parallel_get_size(const int num, const int size, const int rank);
extern int parallel_get_offset(const int num, const int size, const int rank);
extern double parallel_get_wtime(const MPI_Op op);

/* parallel matrix transpose */
extern parallel_transpose_t *parallel_transpose_init(const int g_isize, const int g_jsize, const size_t dtypesize, const MPI_Datatype mpi_dtype);
extern int parallel_transpose_execute(parallel_transpose_t *str, const void *sendbuf, void *recvbuf);
extern int parallel_transpose_finalise(parallel_transpose_t *str);

/* parallel halo communication */
extern int parallel_communicate_halo_with_ymrank(const parallel_t *parallel, const size_t size, const void *sendbuf, void *recvbuf);
extern int parallel_communicate_halo_with_yprank(const parallel_t *parallel, const size_t size, const void *sendbuf, void *recvbuf);

#endif // PARALLEL_H
