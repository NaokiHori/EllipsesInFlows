#include <math.h>
#include <float.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "suspensions.h"


static int increment_particle_velocities(const param_t *param, const int rkstep, suspensions_t *suspensions, double *residual){
  const double Fr = param->Fr;
  const double dt = param->dt;
  const double gamma = param->rkcoefs[rkstep].gamma;
  const int n_particles = suspensions->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = suspensions->particles[n];
    const double den = p->den;
    // store previous increments
    double dux_prev = p->dux;
    double duy_prev = p->duy;
    double dvz_prev = p->dvz;
    // compute new increments
    p->dux =
      +0.5*gamma*dt*(p->cfx[0]+p->cfx[1]) // collision
      +          dt* p->fux               // boundary force
      +(p->iux[1]-p->iux[0])              // internal inertia
      +1.0*gamma*dt*(1.-den)/pow(Fr, 2.)  // buoyancy
      ;
    p->duy =
      +0.5*gamma*dt*(p->cfy[0]+p->cfy[1]) // collision
      +          dt* p->fuy               // boundary force
      +(p->iuy[1]-p->iuy[0])              // internal inertia
      ;
    p->dvz =
      +0.5*gamma*dt*(p->ctz[0]+p->ctz[1]) // collision
      +          dt* p->tvz               // boundary torque
      +(p->ivz[1]-p->ivz[0])              // internal inertia
      ;
    // check convergence
    *residual = fmax(*residual, fabs(p->dux-dux_prev));
    *residual = fmax(*residual, fabs(p->duy-duy_prev));
    *residual = fmax(*residual, fabs(p->dvz-dvz_prev));
  }
  return 0;
}

static int increment_particle_positions(const param_t *param, const int rkstep, suspensions_t *suspensions){
  const double dt = param->dt;
  const double gamma = param->rkcoefs[rkstep].gamma;
  const int n_particles = suspensions->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = suspensions->particles[n];
    double ux0 = p->ux;
    double uy0 = p->uy;
    double vz0 = p->vz;
    double ux1 = p->ux+p->dux;
    double uy1 = p->uy+p->duy;
    double vz1 = p->vz+p->dvz;
    p->dx  = 0.5*gamma*dt*(ux0+ux1);
    p->dy  = 0.5*gamma*dt*(uy0+uy1);
    p->daz = 0.5*gamma*dt*(vz0+vz1);
  }
  return 0;
}

int suspensions_increment_particles(const param_t *param, const int rkstep, suspensions_t *suspensions, double *residual){
  *residual = 0.;
  increment_particle_velocities(param, rkstep, suspensions, residual);
  increment_particle_positions(param, rkstep, suspensions);
  MPI_Allreduce(MPI_IN_PLACE, residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return 0;
}

