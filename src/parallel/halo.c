#include <stddef.h>
#include <mpi.h>
#include "common.h"
#include "parallel.h"


int parallel_communicate_halo_with_ymrank(const parallel_t *parallel, const size_t size, const void *sendbuf, void *recvbuf){
  const int sendtag = 0;
  const int recvtag = 0;
  const int ymrank = parallel->ymrank;
  const int yprank = parallel->yprank;
  MPI_Sendrecv(
      sendbuf, size, MPI_BYTE, yprank, sendtag,
      recvbuf, size, MPI_BYTE, ymrank, recvtag,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE
  );
  return 0;
}

int parallel_communicate_halo_with_yprank(const parallel_t *parallel, const size_t size, const void *sendbuf, void *recvbuf){
  const int sendtag = 0;
  const int recvtag = 0;
  const int ymrank = parallel->ymrank;
  const int yprank = parallel->yprank;
  MPI_Sendrecv(
      sendbuf, size, MPI_BYTE, ymrank, sendtag,
      recvbuf, size, MPI_BYTE, yprank, recvtag,
      MPI_COMM_WORLD, MPI_STATUS_IGNORE
  );
  return 0;
}

