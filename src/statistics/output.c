#include <string.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "statistics.h"
#include "fileio.h"


static int save_param(const char dirname[], const param_t *param, const int num){
  fileio_w_0d_serial(dirname, "num",  NPYIO_INT,    sizeof(int),    &num);
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

static int save_fluid(const char dirname[], const param_t *param, const parallel_t *parallel, const statistics_t *statistics){
  fileio_w_ux_like_parallel(dirname, "ux1", param, parallel, statistics->ux1);
  fileio_w_uy_like_parallel(dirname, "uy1", param, parallel, statistics->uy1);
  fileio_w_ux_like_parallel(dirname, "ux2", param, parallel, statistics->ux2);
  fileio_w_uy_like_parallel(dirname, "uy2", param, parallel, statistics->uy2);
  return 0;
}

static int save_suspensions(const char dirname[], const param_t *param, const parallel_t *parallel, const statistics_t *statistics){
  fileio_w_p_like_parallel(dirname, "phi", param, parallel, statistics->phi);
  return 0;
}

static char *generate_dirname(const int step){
  const char prefix[] = {FILEIO_STAT "/step"};
  const int nzeros = 10;
  char *dirname = common_calloc(
      strlen(prefix)+nzeros+1, // + NUL
      sizeof(char)
  );
  sprintf(dirname, "%s%0*d", prefix, nzeros, step);
  return dirname;
}

int statistics_output(const param_t *param, const parallel_t *parallel, const statistics_t *statistics){
  /* ! create directory from main process ! 2 ! */
  char *dirname = generate_dirname(param->step);
  fileio_mkdir_by_main_process(dirname, parallel);
  /* ! save parameters ! 4 ! */
  const int mpirank = parallel->mpirank;
  if(mpirank == 0){
    save_param(dirname, param, statistics->num);
  }
  /* ! save 2d statistics ! 1 ! */
  save_fluid(dirname, param, parallel, statistics);
  save_suspensions(dirname, param, parallel, statistics);
  common_free(dirname);
  return 0;
}

