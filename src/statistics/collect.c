#include <math.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "suspensions.h"
#include "statistics.h"


static int collect_mean_ux(const param_t *param, const parallel_t *parallel, const fluid_t *fluid, statistics_t *statistics){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double *ux = fluid->ux;
  double *ux1 = statistics->ux1;
  double *ux2 = statistics->ux2;
  /* ! ux and its square are added ! 5 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot+1; i++){
      UX1(i, j) += pow(UX(i, j), 1.);
      UX2(i, j) += pow(UX(i, j), 2.);
    }
  }
  return 0;
}

static int collect_mean_uy(const param_t *param, const parallel_t *parallel, const fluid_t *fluid, statistics_t *statistics){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double *uy = fluid->uy;
  double *uy1 = statistics->uy1;
  double *uy2 = statistics->uy2;
  /* ! uy and its square are added ! 5 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=0; i<=itot+1; i++){
      UY1(i, j) += pow(UY(i, j), 1.);
      UY2(i, j) += pow(UY(i, j), 2.);
    }
  }
  return 0;
}

static int collect_mean_phi(const param_t *param, const parallel_t *parallel, const suspensions_t *suspensions, statistics_t *statistics){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double ly = param->ly;
  const double *xc = param->xc;
  const double *yc = param->yc;
  double *phi = statistics->phi;
  const int n_particles = suspensions->n_particles;
  for(int n = 0; n < n_particles; n++){
    const particle_t *p = suspensions->particles[n];
    const double pa = p->a;
    const double pb = p->b;
    const double px = p->x;
    const double py = p->y;
    const double paz = p->az;
    for(int periodic = -1; periodic <= 1; periodic++){
      double py_ = py+ly*periodic;
      for(int j = 1; j <= jsize; j++){
        const double y = YC(j);
        for(int i = 1; i <= itot; i++){
          const double x = XC(i);
          double dx = x-px;
          double dy = y-py_;
          double c = cos(-paz);
          double s = sin(-paz);
          double x_ = c*dx-s*dy;
          double y_ = s*dx+c*dy;
          // add 1 if this point is inside the ellipse, otherwise do nothing
          PHI(i, j) += 1.-pow(x_/pa, 2.)-pow(y_/pb, 2.) > 0. ? 1. : 0.;
        }
      }
    }
  }
  return 0;
}

int statistics_collect(param_t *param, const parallel_t *parallel, const fluid_t *fluid, const suspensions_t *suspensions, statistics_t *statistics){
  /* ! collect temporally-averaged quantities ! 2 ! */
  collect_mean_ux(param, parallel, fluid, statistics);
  collect_mean_uy(param, parallel, fluid, statistics);
  collect_mean_phi(param, parallel, suspensions, statistics);
  /* ! number of samples is incremented ! 1 ! */
  statistics->num += 1;
  return 0;
}

