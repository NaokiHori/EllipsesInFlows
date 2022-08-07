#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "suspensions.h"
#include "fileio.h"


static int allocate(const param_t *param, const parallel_t *parallel, suspensions_t **suspensions){
  *suspensions = common_calloc(1, sizeof(suspensions_t));
  // response on momentum field
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  (*suspensions)->dux = common_calloc(1, DUX_MEMSIZE);
  (*suspensions)->duy = common_calloc(1, DUY_MEMSIZE);
  return 0;
}

static int init_or_load(const param_t *param, suspensions_t *suspensions){
  char *dirname = NULL;
  int n_particles;
  if(param->load_flow_field){
    // load from given restart directory (param->dirname_restart)
    dirname = common_calloc(strlen(param->dirname_restart)+1, sizeof(char));
    strcpy(dirname, param->dirname_restart);
  }else{
    // load from default directory with initial particle information
    const char dirname_init_particles[] = {"initparticles"};
    dirname = common_calloc(strlen(dirname_init_particles)+1, sizeof(char));
    strcpy(dirname, dirname_init_particles);
  }
  fileio_r_0d_serial(dirname, "n_particles", sizeof(int), &n_particles);
  suspensions->n_particles = n_particles;
  // particles
  suspensions->particles = common_calloc(n_particles, sizeof(particle_t*));
  for(int n = 0; n < n_particles; n++){
    suspensions->particles[n] = common_calloc(1, sizeof(particle_t));
  }
  // prepare buffers
  double *dens = common_calloc(n_particles, sizeof(double));
  double *as   = common_calloc(n_particles, sizeof(double));
  double *bs   = common_calloc(n_particles, sizeof(double));
  double *xs   = common_calloc(n_particles, sizeof(double));
  double *ys   = common_calloc(n_particles, sizeof(double));
  double *azs  = common_calloc(n_particles, sizeof(double));
  double *uxs  = common_calloc(n_particles, sizeof(double));
  double *uys  = common_calloc(n_particles, sizeof(double));
  double *vzs  = common_calloc(n_particles, sizeof(double));
  fileio_r_1d_serial(dirname, "particle_dens", sizeof(double), n_particles, dens);
  fileio_r_1d_serial(dirname, "particle_as",   sizeof(double), n_particles,   as);
  fileio_r_1d_serial(dirname, "particle_bs",   sizeof(double), n_particles,   bs);
  fileio_r_1d_serial(dirname, "particle_xs",   sizeof(double), n_particles,   xs);
  fileio_r_1d_serial(dirname, "particle_ys",   sizeof(double), n_particles,   ys);
  fileio_r_1d_serial(dirname, "particle_azs",  sizeof(double), n_particles,  azs);
  fileio_r_1d_serial(dirname, "particle_uxs",  sizeof(double), n_particles,  uxs);
  fileio_r_1d_serial(dirname, "particle_uys",  sizeof(double), n_particles,  uys);
  fileio_r_1d_serial(dirname, "particle_vzs",  sizeof(double), n_particles,  vzs);
  for(int n = 0; n < n_particles; n++){
    particle_t *p = suspensions->particles[n];
    p->den = dens[n];
    p->a   = as[n];
    p->b   = bs[n];
    p->x   = xs[n];
    p->y   = ys[n];
    p->az  = azs[n];
    p->ux  = uxs[n];
    p->uy  = uys[n];
    p->vz  = vzs[n];
  }
  common_free(dens);
  common_free(as);
  common_free(bs);
  common_free(xs);
  common_free(ys);
  common_free(azs);
  common_free(uxs);
  common_free(uys);
  common_free(vzs);
  common_free(dirname);
  return 0;
}

suspensions_t *suspensions_init(const param_t *param, const parallel_t *parallel){
  suspensions_t *suspensions = NULL;
  allocate(param, parallel, &suspensions);
  init_or_load(param, suspensions);
  // buffers to communicate Lagrange info
  suspensions->buf = common_calloc(3*suspensions->n_particles, sizeof(double));
  return suspensions;
}

