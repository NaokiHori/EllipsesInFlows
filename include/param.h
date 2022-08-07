#if !defined(PARAM_H)
#define PARAM_H

#include <stdbool.h>

#include "structure.h"
#include "arrays/param.h"

// 3-step Runge-Kutta scheme
// set to 1 to use EF1
#define RKSTEPMAX 3

// wall-to-wall distance is fixed to unity
#define LX 1.

typedef struct rkcoef_t_ {
  double alpha;
  double beta;
  double gamma;
} rkcoef_t;

typedef struct schedule_t_ {
  double rate;
  double after;
  double next;
} schedule_t;

/* ! definition of a structure param_t_ ! 24 !*/
struct param_t_ {
  // restart / initialise
  bool load_flow_field;
  char *dirname_restart;
  // domain sizes
  int itot, jtot;
  double lx, ly;
  // grid sizes
  double *xf, *xc;
  double *dxf, *dxc;
  double *yf, *yc;
  double dx, dy;
  // non-dimensional params
  double Re, Fr;
  // external force in y
  double extfrcy;
  // temporal integration
  rkcoef_t rkcoefs[3];
  double time, dt;
  double safefactors[3];
  int step;
  // when to stop, when to write log, etc.
  double timemax, wtimemax;
  schedule_t log, save, stat;
};

extern param_t *param_init(void);
extern int param_finalise(param_t *param);

extern int param_decide_dt(param_t *param, const parallel_t *parallel, const fluid_t *fluid, const suspensions_t *suspensions);

#endif // PARAM_H
