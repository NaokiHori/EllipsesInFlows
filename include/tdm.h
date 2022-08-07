#if !defined(TDM_H)
#define TDM_H

#include <stdbool.h>
#include <fftw3.h>

typedef struct {
  // number of elements
  int n;
  // periodic or not
  bool is_periodic;
  // original matrix
  double *l, *c, *u;
  // buffers, both for periodic and non-periodic cases
  double *l0, *c0, *u0;
  // buffers, only for periodic cases
  double *l1, *c1, *u1, *q1;
} tdm_t;

extern tdm_t *tdm_init(const int n, const bool is_periodic);
extern int tdm_solve_double(tdm_t *tdm, double *q);
extern int tdm_solve_fftw_complex(tdm_t *tdm, fftw_complex *q);
extern int tdm_finalise(tdm_t *tdm);

#endif // TDM_H
