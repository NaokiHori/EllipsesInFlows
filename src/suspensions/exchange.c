#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "suspensions.h"


static int kernel_exchange_momentum(const param_t *param, const parallel_t *parallel, const fluid_t *fluid, suspensions_t *suspensions){
  // \int -fx dV
  // \int -fy dV
  // \int ( + ry fx - rx fy ) dV
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double dt = param->dt;
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
  double *dux = suspensions->dux;
  double *duy = suspensions->duy;
  memset(dux, 0, DUX_MEMSIZE);
  memset(duy, 0, DUY_MEMSIZE);
  for(int n = 0; n < n_particles; n++){
    particle_t *p = suspensions->particles[n];
    // constant parameters
    const double pden = p->den;
    const double pa   = p->a;
    const double pb   = p->b;
    const double pm   = suspensions_compute_mass(pden, pa, pb);
    const double pim  = suspensions_compute_moment_of_inertia(pden, pa, pb);
    const double px   = p->x;
    const double py   = p->y;
    const double paz  = p->az;
    const double pux  = p->ux;
    const double puy  = p->uy;
    const double pvz  = p->vz;
    // buffers (for simplicity)
    double fux = 0.;
    double fuy = 0.;
    double tvz = 0.;
    for(int periodic = -1; periodic <= 1; periodic++){
      double py_ = py+ly*periodic;
      int imin, imax, jmin, jmax;
      suspensions_decide_loop_size(1, itot,  dx, fmax(pa, pb), px,        &imin, &imax);
      suspensions_decide_loop_size(1, jsize, dy, fmax(pa, pb), py_-YF(1), &jmin, &jmax);
      for(int j = jmin; j <= jmax; j++){
        double y = YC(j);
        for(int i = imin; i <= imax; i++){
          double x = XC(i);
          double w = suspensions_s_weight(grid_size, pa, pb, px, py_, paz, x, y);
          double ux_p = pux-pvz*(y-py_);
          double uy_p = puy+pvz*(x-px);
          double ux_f = 0.5*(UX(i  , j  )+UX(i+1, j  ));
          double uy_f = 0.5*(UY(i  , j  )+UY(i  , j+1));
          double fx = w*(ux_p-ux_f)/dt;
          double fy = w*(uy_p-uy_f)/dt;
          DUX(i, j) += fx*dt;
          DUY(i, j) += fy*dt;
          fux -= +(fx*dx*dy)/pm;
          tvz -= -(y-py_)*(fx*dx*dy)/pim;
          fuy -= +(fy*dx*dy)/pm;
          tvz -= +(x-px)*(fy*dx*dy)/pim;
        }
      }
    }
    fluid_update_boundaries_p(param, parallel, dux);
    fluid_update_boundaries_p(param, parallel, duy);
    // assign results
    p->fux = fux;
    p->fuy = fuy;
    p->tvz = tvz;
  }
  return 0;
}

static int synchronise_information(suspensions_t *suspensions){
  const int n_particles = suspensions->n_particles;
  particle_t **particles = suspensions->particles;
  // prepare message buffer
  double *buf = suspensions->buf;
  // pack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    buf[3*n+0] = p->fux;
    buf[3*n+1] = p->fuy;
    buf[3*n+2] = p->tvz;
  }
  // sum up all
  MPI_Allreduce(MPI_IN_PLACE, buf, 3*n_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // unpack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    p->fux = buf[3*n+0];
    p->fuy = buf[3*n+1];
    p->tvz = buf[3*n+2];
  }
  return 0;
}

int suspensions_exchange_momentum(const param_t *param, const parallel_t *parallel, const fluid_t *fluid, suspensions_t *suspensions){
  // exchange (translational and angular) momenta with each particle
  kernel_exchange_momentum(param, parallel, fluid, suspensions);
  // communicate updated information
  synchronise_information(suspensions);
  return 0;
}

