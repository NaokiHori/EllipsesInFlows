#include <string.h>
#include "param.h"
#include "parallel.h"
#include "suspensions.h"


int suspensions_reset_particle_increments(suspensions_t *suspensions){
  const int n_particles = suspensions->n_particles;
  for(int n = 0; n < n_particles; n++){
    suspensions->particles[n]->dux = 0.;
    suspensions->particles[n]->duy = 0.;
    suspensions->particles[n]->dvz = 0.;
    suspensions->particles[n]->dx  = 0.;
    suspensions->particles[n]->dy  = 0.;
    suspensions->particles[n]->daz = 0.;
  }
  return 0;
}

