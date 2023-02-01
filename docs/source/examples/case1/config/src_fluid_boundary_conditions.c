#include <mpi.h>
#include "param.h"
#include "parallel.h"
#include "fluid.h"


int fluid_update_boundaries_ux(const param_t *param, const parallel_t *parallel, double *ux){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const int nitems = itot-1;
  const size_t size = sizeof(double)*nitems;
  /* ! update halo values of ux ! 2 ! */
  parallel_communicate_halo_with_ymrank(parallel, size, &UX(2, jsize), &UX(2,       0));
  parallel_communicate_halo_with_yprank(parallel, size, &UX(2,     1), &UX(2, jsize+1));
  /* ! set boundary values of ux ! 4 ! */
  for(int j=0; j<=jsize+1; j++){
    UX(     1, j) = 0.; // impermeable
    UX(itot+1, j) = 0.; // impermeable
  }
  return 0;
}

int fluid_update_boundaries_uy(const param_t *param, const parallel_t *parallel, double *uy){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const int nitems = itot;
  const size_t size = sizeof(double)*nitems;
  /* ! update halo values of uy ! 2 ! */
  parallel_communicate_halo_with_ymrank(parallel, size, &UY(1, jsize), &UY(1,       0));
  parallel_communicate_halo_with_yprank(parallel, size, &UY(1,     1), &UY(1, jsize+1));
  /* ! set boundary values of uy ! 4 ! */
  for(int j=0; j<=jsize+1; j++){
    UY(     0, j) = -0.5; // no-slip
    UY(itot+1, j) = +0.5; // no-slip
  }
  return 0;
}

int fluid_update_boundaries_p(const param_t *param, const parallel_t *parallel, double *p){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const int nitems = itot;
  const size_t size = sizeof(double)*nitems;
  /* ! update halo values of p ! 2 ! */
  parallel_communicate_halo_with_ymrank(parallel, size, &P(1, jsize), &P(1,       0));
  parallel_communicate_halo_with_yprank(parallel, size, &P(1,     1), &P(1, jsize+1));
  /* ! set boundary values of p ! 4 ! */
  for(int j=0; j<=jsize+1; j++){
    P(     0, j) = P(   1, j); // Neumann
    P(itot+1, j) = P(itot, j); // Neumann
  }
  return 0;
}

