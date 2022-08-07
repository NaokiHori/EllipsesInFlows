#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "suspensions.h"


static int update_particle_velocities(suspensions_t *suspensions){
  const int n_particles = suspensions->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = suspensions->particles[n];
    p->ux += p->dux;
    p->uy += p->duy;
    p->vz += p->dvz;
  }
  return 0;
}

static int update_particle_positions(const param_t *param, suspensions_t *suspensions){
  const int n_particles = suspensions->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = suspensions->particles[n];
    p->x  += p->dx;
    p->y  += p->dy;
    p->az += p->daz;
    // correct periodicity
    const double ly = param->ly;
    if(p->y < 0.) p->y += ly;
    if(p->y > ly) p->y -= ly;
    if(p->az <      0.) p->az += 2.*M_PI;
    if(p->az > 2.*M_PI) p->az -= 2.*M_PI;
  }
  return 0;
}

int suspensions_update_momentum_fleid(const param_t *param, const parallel_t *parallel, fluid_t *fluid, const suspensions_t *suspensions){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  double *ux = fluid->ux;
  double *uy = fluid->uy;
  const double *dux = suspensions->dux;
  const double *duy = suspensions->duy;
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= itot; i++){
      UX(i, j) += 0.5*(DUX(i-1, j  )+DUX(i  , j  ));
    }
  }
  fluid_update_boundaries_ux(param, parallel, ux);
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= itot; i++){
      UY(i, j) += 0.5*(DUY(i  , j-1)+DUY(i  , j  ));
    }
  }
  fluid_update_boundaries_uy(param, parallel, uy);
  return 0;
}

int suspensions_update_particles(const param_t *param, suspensions_t *suspensions){
  update_particle_velocities(suspensions);
  update_particle_positions(param, suspensions);
  return 0;
}

