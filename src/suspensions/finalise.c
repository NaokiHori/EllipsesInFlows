#include "common.h"
#include "suspensions.h"


int suspensions_finalise(suspensions_t *suspensions){
  // particles
  const int n_particles = suspensions->n_particles;
  for(int n = 0; n < n_particles; n++){
    common_free(suspensions->particles[n]);
  }
  common_free(suspensions->particles);
  // Euler variables
  common_free(suspensions->dux);
  common_free(suspensions->duy);
  // buffers to communicate Lagrange info
  common_free(suspensions->buf);
  // main structure
  common_free(suspensions);
  return 0;
}

