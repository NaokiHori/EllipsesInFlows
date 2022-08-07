#include <math.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "suspensions.h"


static int compute_inertia(const param_t *param, const parallel_t *parallel, const int cnstep, const fluid_t *fluid, suspensions_t *suspensions){
  // \int ux dV
  // \int uy dV
  // \int ( - ry ux + rx uy ) dV
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double ly = param->ly;
  const double *xc = param->xc;
  const double *yf = param->yf;
  const double *yc = param->yc;
  const double dx = param->dx;
  const double dy = param->dy;
  const double grid_size = fmin(dx, dy);
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  const int n_particles = suspensions->n_particles;
  particle_t **particles = suspensions->particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    // constant parameters
    const double pden = p->den;
    const double pa   = p->a;
    const double pb   = p->b;
    const double px   = p->x+p->dx;
    const double py   = p->y+p->dy;
    const double paz  = p->az+p->daz;
    const double pm   = suspensions_compute_mass(pden, pa, pb);
    const double pim  = suspensions_compute_moment_of_inertia(pden, pa, pb);
    // buffers (for simplicity)
    double iux = 0.;
    double iuy = 0.;
    double ivz = 0.;
    // ux contribution
    for(int periodic = -1; periodic <= 1; periodic++){
      double py_ = py+ly*periodic;
      int imin, imax, jmin, jmax;
      suspensions_decide_loop_size(1, itot,  dx, fmax(pa, pb), px,        &imin, &imax);
      suspensions_decide_loop_size(1, jsize, dy, fmax(pa, pb), py_-YF(1), &jmin, &jmax);
      for(int j = jmin; j <= jmax; j++){
        double y = YC(j);
        for(int i = imin; i <= imax; i++){
          double x = XC(i);
          double w = suspensions_v_weight(grid_size, pa, pb, px, py_, paz, x, y);
          double valx = w*0.5*(UX(i  , j  )+UX(i+1, j  ))*(dx*dy);
          double valy = w*0.5*(UY(i  , j  )+UY(i  , j+1))*(dx*dy);
          iux += valx/pm;
          ivz += -(y-py_)*valx/pim;
          iuy += valy/pm;
          ivz += +(x-px)*valy/pim;
        }
      }
    }
    // assign results
    p->iux[cnstep] = iux;
    p->iuy[cnstep] = iuy;
    p->ivz[cnstep] = ivz;
  }
  return 0;
}

static int synchronise_information(const int cnstep, suspensions_t *suspensions){
  const int n_particles = suspensions->n_particles;
  particle_t **particles = suspensions->particles;
  // message buffer
  double *buf = suspensions->buf;
  // pack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    buf[3*n+0] = p->iux[cnstep];
    buf[3*n+1] = p->iuy[cnstep];
    buf[3*n+2] = p->ivz[cnstep];
  }
  // sum up all
  MPI_Allreduce(MPI_IN_PLACE, buf, 3*n_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // unpack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    p->iux[cnstep] = buf[3*n+0];
    p->iuy[cnstep] = buf[3*n+1];
    p->ivz[cnstep] = buf[3*n+2];
  }
  return 0;
}

int suspensions_compute_inertia(const param_t *param, const parallel_t *parallel, const int cnstep, const fluid_t *fluid, suspensions_t *suspensions){
  // update for each particle LOCALLY
  compute_inertia(param, parallel, cnstep, fluid, suspensions);
  // communicate updated information
  synchronise_information(cnstep, suspensions);
  return 0;
}

