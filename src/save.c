#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "suspensions.h"
#include "save.h"
#include "fileio.h"


static int save_param(const char dirname[], const param_t *param){
  fileio_w_0d_serial(dirname, "step", NPYIO_INT,    sizeof(int),    &(param->step));
  fileio_w_0d_serial(dirname, "itot", NPYIO_INT,    sizeof(int),    &(param->itot));
  fileio_w_0d_serial(dirname, "jtot", NPYIO_INT,    sizeof(int),    &(param->jtot));
  fileio_w_0d_serial(dirname, "time", NPYIO_DOUBLE, sizeof(double), &(param->time));
  fileio_w_0d_serial(dirname, "lx",   NPYIO_DOUBLE, sizeof(double), &(param->lx  ));
  fileio_w_0d_serial(dirname, "ly",   NPYIO_DOUBLE, sizeof(double), &(param->ly  ));
  fileio_w_0d_serial(dirname, "Re",   NPYIO_DOUBLE, sizeof(double), &(param->Re  ));
  fileio_w_1d_serial(dirname, "xf",   NPYIO_DOUBLE, sizeof(double), param->itot+1, param->xf);
  fileio_w_1d_serial(dirname, "xc",   NPYIO_DOUBLE, sizeof(double), param->itot+2, param->xc);
  // GLOBAL y coordinates (NOTE: param->yc and param->yf are LOCAL)
  const int jtot = param->jtot;
  const double dy = param->dy;
  double *yf = common_calloc(jtot+1, sizeof(double));
  double *yc = common_calloc(jtot  , sizeof(double));
  for(int j = 0; j < jtot+1; j++){
    yf[j] = 1.*j*dy;
  }
  for(int j = 0; j < jtot; j++){
    yc[j] = 0.5*(yf[j]+yf[j+1]);
  }
  fileio_w_1d_serial(dirname, "yf", NPYIO_DOUBLE, sizeof(double), jtot+1, yf);
  fileio_w_1d_serial(dirname, "yc", NPYIO_DOUBLE, sizeof(double), jtot  , yc);
  common_free(yf);
  common_free(yc);
  return 0;
}

static int save_fluid(const char dirname[], const param_t *param, const parallel_t *parallel, const fluid_t *fluid){
  fileio_w_ux_like_parallel(dirname, "ux", param, parallel, fluid->ux);
  fileio_w_uy_like_parallel(dirname, "uy", param, parallel, fluid->uy);
  fileio_w_p_like_parallel (dirname, "p",  param, parallel, fluid->p );
  return 0;
}

static int save_particles(const char dirname[], const suspensions_t *suspensions){
  const int n_particles = suspensions->n_particles;
  fileio_w_0d_serial(dirname, "n_particles", NPYIO_INT, sizeof(int), &n_particles);
  // save as SOA (radii, x coordinates, y coordinates etc.) instead of AOS (particle_t **particles)
  // to reduce the number of files
  double *dens = common_calloc(n_particles, sizeof(double));
  double *as   = common_calloc(n_particles, sizeof(double));
  double *bs   = common_calloc(n_particles, sizeof(double));
  double *xs   = common_calloc(n_particles, sizeof(double));
  double *ys   = common_calloc(n_particles, sizeof(double));
  double *azs  = common_calloc(n_particles, sizeof(double));
  double *uxs  = common_calloc(n_particles, sizeof(double));
  double *uys  = common_calloc(n_particles, sizeof(double));
  double *vzs  = common_calloc(n_particles, sizeof(double));
  for(int n = 0; n < n_particles; n++){
    particle_t *p = suspensions->particles[n];
    dens[n] = p->den;
    as[n]   = p->a;
    bs[n]   = p->b;
    xs[n]   = p->x;
    ys[n]   = p->y;
    azs[n]  = p->az;
    uxs[n]  = p->ux;
    uys[n]  = p->uy;
    vzs[n]  = p->vz;
  }
  fileio_w_1d_serial(dirname, "particle_dens", NPYIO_DOUBLE, sizeof(double), n_particles, dens);
  fileio_w_1d_serial(dirname, "particle_as",   NPYIO_DOUBLE, sizeof(double), n_particles,   as);
  fileio_w_1d_serial(dirname, "particle_bs",   NPYIO_DOUBLE, sizeof(double), n_particles,   bs);
  fileio_w_1d_serial(dirname, "particle_xs",   NPYIO_DOUBLE, sizeof(double), n_particles,   xs);
  fileio_w_1d_serial(dirname, "particle_ys",   NPYIO_DOUBLE, sizeof(double), n_particles,   ys);
  fileio_w_1d_serial(dirname, "particle_azs",  NPYIO_DOUBLE, sizeof(double), n_particles,  azs);
  fileio_w_1d_serial(dirname, "particle_uxs",  NPYIO_DOUBLE, sizeof(double), n_particles,  uxs);
  fileio_w_1d_serial(dirname, "particle_uys",  NPYIO_DOUBLE, sizeof(double), n_particles,  uys);
  fileio_w_1d_serial(dirname, "particle_vzs",  NPYIO_DOUBLE, sizeof(double), n_particles,  vzs);
  common_free(dens);
  common_free(as);
  common_free(bs);
  common_free(xs);
  common_free(ys);
  common_free(azs);
  common_free(uxs);
  common_free(uys);
  common_free(vzs);
  return 0;
}

static char *generate_dirname(const int step){
  const char prefix[] = {FILEIO_SAVE "/step"};
  const int nzeros = 10;
  char *dirname = common_calloc(
      strlen(prefix)+nzeros+1, // + NUL
      sizeof(char)
  );
  sprintf(dirname, "%s%0*d", prefix, nzeros, step);
  return dirname;
}

int save(param_t *param, const parallel_t *parallel, const fluid_t *fluid, const suspensions_t *suspensions){
  /* ! create directory from main process ! 2 ! */
  char *dirname = generate_dirname(param->step);
  fileio_mkdir_by_main_process(dirname, parallel);
  /* ! save parameters ! 4 ! */
  const int mpirank = parallel->mpirank;
  if(mpirank == 0){
    save_param(dirname, param);
  }
  /* ! save flow fields ! 1 ! */
  save_fluid(dirname, param, parallel, fluid);
  if(mpirank == 0){
    save_particles(dirname, suspensions);
  }
  common_free(dirname);
  return 0;
}

