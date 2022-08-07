#include <stdio.h>
#include <math.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "suspensions.h"
#include "fileio.h"
#include "logging.h"


static int get_ngidits(const double rate){
  int retval = 1;
  int rateinv = (int)(1./rate);
  while(rateinv /= 10){
    retval++;
  }
  return retval;
}

static int show_progress(const char fname[], const param_t *param, const parallel_t *parallel){
  if(parallel->mpirank == 0){
    FILE *fp = NULL;
    if(param->step == 0){
      fp = fileio_fopen(fname, "w");
    }else{
      fp = fileio_fopen(fname, "a");
    }
    if(fp != NULL){
      int ndigits = get_ngidits(param->log.rate);
      /* ! show progress to standard output and file ! 2 ! */
      fprintf(fp,     "step %8d, time %*.*f, dt %.2e\n", param->step, ndigits+3, ndigits, param->time, param->dt);
      fprintf(stdout, "step %8d, time %*.*f, dt %.2e\n", param->step, ndigits+3, ndigits, param->time, param->dt);
      fileio_fclose(fp);
    }
  }
  return 0;
}

static int check_divergence(const char fname[], const param_t *param, const parallel_t *parallel, const fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double *dxf = param->dxf;
  const double dy = param->dy;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  double divmax = 0.;
  double divsum = 0.;
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      /* ! compute local divergence ! 7 ! */
      double ux_xm = UX(i  , j  );
      double ux_xp = UX(i+1, j  );
      double uy_ym = UY(i  , j  );
      double uy_yp = UY(i  , j+1);
      double div =
        (ux_xp-ux_xm)/DXF(i)
       +(uy_yp-uy_ym)/dy;
      /* ! check local maximum divergence and global divergence ! 2 ! */
      divmax = fmax(divmax, fabs(div));
      divsum += div;
    }
  }
  /* ! collect information among all processes ! 2 ! */
  MPI_Allreduce(MPI_IN_PLACE, &divmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &divsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  /* ! information are output, similar for other functions ! 10 ! */
  if(mpirank == 0){
    FILE *fp = NULL;
    if(param->step == 0){
      fp = fileio_fopen(fname, "w");
    }else{
      fp = fileio_fopen(fname, "a");
    }
    if(fp != NULL){
      int ndigits = get_ngidits(param->log.rate);
      fprintf(fp, "%*.*f % .1e % .1e\n", ndigits+3, ndigits, param->time, divmax, divsum);
      fileio_fclose(fp);
    }
  }
  return 0;
}

static int check_momentum(const char fname[], const param_t *param, const parallel_t *parallel, const fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double *dxf = param->dxf;
  const double *dxc = param->dxc;
  const double dy = param->dy;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  double moms[2] = {0.};
  /* ! compute total x-momentum ! 6 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot+1; i++){
      double cellsize = DXC(i)*dy;
      moms[0] += UX(i, j)*cellsize;
    }
  }
  /* ! compute total y-momentum ! 6 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      double cellsize = DXF(i)*dy;
      moms[1] += UY(i, j)*cellsize;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, moms, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(mpirank == 0){
    FILE *fp = NULL;
    if(param->step == 0){
      fp = fileio_fopen(fname, "w");
    }else{
      fp = fileio_fopen(fname, "a");
    }
    if(fp != NULL){
      int ndigits = get_ngidits(param->log.rate);
      fprintf(fp, "%*.*f % 18.15e % 18.15e\n", ndigits+3, ndigits, param->time, moms[0], moms[1]);
      fileio_fclose(fp);
    }
  }
  return 0;
}

static int check_energy(const char fname[], const param_t *param, const parallel_t *parallel, const fluid_t *fluid){
  /* compute kinetic and thermal energies */
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double *dxf = param->dxf;
  const double *dxc = param->dxc;
  const double dy = param->dy;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  double quantities[2] = {0.};
  /* ! compute quadratic quantity in x direction ! 6 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot+1; i++){
      double cellsize = DXC(i)*dy;
      quantities[0] += 0.5*pow(UX(i, j), 2.)*cellsize;
    }
  }
  /* ! compute quadratic quantity in y direction ! 6 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      double cellsize = DXF(i)*dy;
      quantities[1] += 0.5*pow(UY(i, j), 2.)*cellsize;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, quantities, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(mpirank == 0){
    FILE *fp = NULL;
    if(param->step == 0){
      fp = fileio_fopen(fname, "w");
    }else{
      fp = fileio_fopen(fname, "a");
    }
    if(fp != NULL){
      int ndigits = get_ngidits(param->log.rate);
      fprintf(fp, "%*.*f % 18.15e % 18.15e\n", ndigits+3, ndigits, param->time, quantities[0], quantities[1]);
      fileio_fclose(fp);
    }
  }
  return 0;
}

static int check_particle(const char fname[], const param_t *param, const parallel_t *parallel, const particle_t *p){
  const int mpirank = parallel->mpirank;
  if(mpirank == 0){
    FILE *fp = NULL;
    if(param->step == 0){
      fp = fileio_fopen(fname, "w");
    }else{
      fp = fileio_fopen(fname, "a");
    }
    if(fp != NULL){
      int ndigits = get_ngidits(param->log.rate);
      fprintf(fp, "%*.*f % .7e % .7e % .7e % .7e % .7e % .7e\n",
          ndigits+3, ndigits,
          param->time,
          p->x, p->y, p->az,
          p->ux, p->uy, p->vz
      );
      fileio_fclose(fp);
    }
  }
  return 0;
}

int logging(param_t *param, const parallel_t *parallel, const fluid_t *fluid, const suspensions_t *suspensions){
  show_progress   (FILEIO_LOG "/progress.dat",   param, parallel);
  check_divergence(FILEIO_LOG "/divergence.dat", param, parallel, fluid);
  check_momentum  (FILEIO_LOG "/momentum.dat",   param, parallel, fluid);
  check_energy    (FILEIO_LOG "/energy.dat",     param, parallel, fluid);
  for(int n = 0; n < suspensions->n_particles; n++){
    char fname[128];
    sprintf(fname, "%s/particle%010d.dat", FILEIO_LOG, n);
    check_particle(fname, param, parallel, suspensions->particles[n]);
  }
  return 0;
}

