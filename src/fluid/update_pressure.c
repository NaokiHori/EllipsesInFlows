#include <math.h>
#include "param.h"
#include "parallel.h"
#include "fluid.h"


int fluid_update_pressure(const param_t *param, const parallel_t *parallel, fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double *psi = fluid->psi;
  double *p = fluid->p;
  /* ! add correction ! 5 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      P(i, j) += PSI(i, j);
    }
  }
  /* ! boundary and halo values are updated ! 1 ! */
  fluid_update_boundaries_p(param, parallel, p);
  return 0;
}

