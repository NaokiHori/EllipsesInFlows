#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "tdm.h"


static buffers_compute_potential_t *init(const param_t *param, const parallel_t *parallel){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  buffers_compute_potential_t *str = common_calloc(1, sizeof(buffers_compute_potential_t));
  // tri-diagonal matrix solver
  str->tdm_solver = tdm_init(jtot, true);
  // buffers
  {
    int sizes[2];
    // x-aligned (y-decomposed)
    sizes[0] = itot;
    sizes[1] = parallel_get_size(jtot, mpisize, mpirank);
    str->qx = common_calloc(sizes[0]*sizes[1], sizeof(double));
    // y-aligned (x-decomposed)
    sizes[0] = parallel_get_size(itot, mpisize, mpirank);
    sizes[1] = jtot;
    str->qy = common_calloc(sizes[0]*sizes[1], sizeof(double));
  }
  /* ! create fftw plans ! 16 ! */
  {
    str->fftw_buf_r = common_calloc(itot, sizeof(double));
    str->fftw_plan_fwrd = fftw_plan_r2r_1d(itot, str->fftw_buf_r, str->fftw_buf_r, FFTW_REDFT10, FFTW_PATIENT);
    str->fftw_plan_bwrd = fftw_plan_r2r_1d(itot, str->fftw_buf_r, str->fftw_buf_r, FFTW_REDFT01, FFTW_PATIENT);
  }
  /* parallel matrix transpose */
  {
    str->transposer_x_to_y = parallel_transpose_init(itot, jtot, sizeof(double), MPI_DOUBLE);
    str->transposer_y_to_x = parallel_transpose_init(jtot, itot, sizeof(double), MPI_DOUBLE);
  }
  return str;
}

#define QX(I, J) (qx[((J)-1)*(itot)+((I)-1)])
#define QY(I, J) (qy[((I)-1)*(jtot)+((J)-1)])

int fluid_compute_potential(const param_t *param, const parallel_t *parallel, const int rkstep, fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int isize = parallel_get_size(itot, mpisize, mpirank);
  const int ioffset = parallel_get_offset(itot, mpisize, mpirank);
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double dx = param->dx;
  const double dy = param->dy;
  double *psi = fluid->psi;
  if(fluid->buffers_compute_potential == NULL){
    fluid->buffers_compute_potential = init(param, parallel);
  }
  buffers_compute_potential_t *buffers = fluid->buffers_compute_potential;
  double *qx = buffers->qx;
  double *qy = buffers->qy;
  /* ! compute right-hand-side ! 17 ! */
  const double gamma = param->rkcoefs[rkstep].gamma;
  const double dt = param->dt;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      double ux_xm = UX(i  , j  );
      double ux_xp = UX(i+1, j  );
      double uy_ym = UY(i  , j  );
      double uy_yp = UY(i  , j+1);
      QX(i, j) =
        1./(gamma*dt)*(
         +(ux_xp-ux_xm)/dx
         +(uy_yp-uy_ym)/dy
        );
    }
  }
  /* ! project to wave space ! 9 ! */
  {
    fftw_plan plan = buffers->fftw_plan_fwrd;
    double *r = buffers->fftw_buf_r;
    for(int j = 1; j <= jsize; j++){
      memcpy(r, &QX(1, j), sizeof(double)*itot);
      fftw_execute(plan);
      memcpy(&QX(1, j), r, sizeof(double)*itot);
    }
  }
  /* ! transpose x-aligned matrix to y-aligned matrix ! 1 ! */
  parallel_transpose_execute(buffers->transposer_x_to_y, qx, qy);
  /* solve linear systems */
  {
    tdm_t *tdm_solver = buffers->tdm_solver;
    double *tdm_l = tdm_solver->l;
    double *tdm_c = tdm_solver->c;
    double *tdm_u = tdm_solver->u;
    for(int i = 1; i <= isize; i++){
      /* ! compute eigenvalue of this i position ! 4 ! */
      double eigenvalue = -4./pow(dx, 2.)*pow(
          sin(M_PI*(i+ioffset-1)/(2.*itot)),
          2.
      );
      /* ! initialise tri-diagonal matrix ! 5 ! */
      for(int j = 0; j < jtot; j++){
        tdm_l[j] = 1./dy/dy;
        tdm_u[j] = 1./dy/dy;
        tdm_c[j] = -tdm_l[j]-tdm_u[j]+eigenvalue;
      }
      /* ! solve linear system ! 7 ! */
      tdm_solve_double(tdm_solver, &QY(i, 1));
    }
  }
  /* ! transpose y-aligned matrix to x-aligned matrix ! 1 ! */
  parallel_transpose_execute(buffers->transposer_y_to_x, qy, qx);
  /* ! project to physical space ! 9 ! */
  {
    fftw_plan plan = buffers->fftw_plan_bwrd;
    double *r = buffers->fftw_buf_r;
    for(int j = 1; j <= jsize; j++){
      memcpy(r, &QX(1, j), sizeof(double)*itot);
      fftw_execute(plan);
      memcpy(&QX(1, j), r, sizeof(double)*itot);
    }
  }
  /* ! normalise and store result ! 5 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      PSI(i, j) = QX(i, j)/(2.*itot);
    }
  }
  // NOTE: psi and p are assumed to have the same array shape
  fluid_update_boundaries_p(param, parallel, psi);
  return 0;
}

#undef QX
#undef QY

