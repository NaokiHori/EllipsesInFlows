#if !defined(FLUID_H)
#define FLUID_H

#include <fftw3.h>
#include "structure.h"
#include "arrays/fluid.h"
#include "parallel.h"
#include "tdm.h"

typedef struct {
  double *qx, *qy;
  fftw_plan fftw_plan_fwrd, fftw_plan_bwrd;
  double *fftw_buf_r;
  tdm_t *tdm_solver;
  parallel_transpose_t *transposer_x_to_y, *transposer_y_to_x;
} buffers_compute_potential_t;

/* ! definition of a structure fluid_t_ ! 6 ! */
struct fluid_t_ {
  double *ux, *uy;
  double *p, *psi;
  double *srcuxa, *srcuxb, *srcuxg;
  double *srcuya, *srcuyb, *srcuyg;
  buffers_compute_potential_t *buffers_compute_potential;
};

extern fluid_t *fluid_init(const param_t *param, const parallel_t *parallel);
extern int fluid_finalise(fluid_t *fluid);

/* boundary condition and halo update */
extern int fluid_update_boundaries_ux(const param_t *param, const parallel_t *parallel, double *ux);
extern int fluid_update_boundaries_uy(const param_t *param, const parallel_t *parallel, double *uy);
extern int fluid_update_boundaries_p(const param_t *param, const parallel_t *parallel, double *p);

extern int fluid_update_velocity(const param_t *param, const parallel_t *parallel, const int rkstep, fluid_t *fluid);
extern int fluid_compute_potential(const param_t *param, const parallel_t *parallel, const int rkstep, fluid_t *fluid);
extern int fluid_correct_velocity(const param_t *param, const parallel_t *parallel, const int rkstep, fluid_t *fluid);
extern int fluid_update_pressure(const param_t *param, const parallel_t *parallel, fluid_t *fluid);

#endif // FLUID_H
