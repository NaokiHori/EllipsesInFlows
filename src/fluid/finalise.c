#include "common.h"
#include "fluid.h"


static int deallocate_buffers_compute_potential(buffers_compute_potential_t *str){
  common_free(str->qx);
  common_free(str->qy);
  common_free(str->fftw_buf_r);
  fftw_destroy_plan(str->fftw_plan_fwrd);
  fftw_destroy_plan(str->fftw_plan_bwrd);
  fftw_cleanup();
  tdm_finalise(str->tdm_solver);
  parallel_transpose_finalise(str->transposer_x_to_y);
  parallel_transpose_finalise(str->transposer_y_to_x);
  common_free(str);
  return 0;
}

int fluid_finalise(fluid_t *fluid){
  common_free(fluid->ux);
  common_free(fluid->uy);
  common_free(fluid->p);
  common_free(fluid->psi);
  common_free(fluid->srcuxa);
  common_free(fluid->srcuxb);
  common_free(fluid->srcuxg);
  common_free(fluid->srcuya);
  common_free(fluid->srcuyb);
  common_free(fluid->srcuyg);
  deallocate_buffers_compute_potential(fluid->buffers_compute_potential);
  common_free(fluid);
  return 0;
}

