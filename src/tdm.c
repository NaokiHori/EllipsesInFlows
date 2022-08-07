#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>
#include "common.h"
#include "tdm.h"


tdm_t *tdm_init(const int n, const bool is_periodic){
  tdm_t *str = common_calloc(1, sizeof(tdm_t));
  str->n = n;
  str->is_periodic = is_periodic;
  str->l = common_calloc(n, sizeof(double));
  str->c = common_calloc(n, sizeof(double));
  str->u = common_calloc(n, sizeof(double));
  if(is_periodic){
    str->l0 = common_calloc(n-1, sizeof(double));
    str->c0 = common_calloc(n-1, sizeof(double));
    str->u0 = common_calloc(n-1, sizeof(double));
    str->l1 = common_calloc(n-1, sizeof(double));
    str->c1 = common_calloc(n-1, sizeof(double));
    str->u1 = common_calloc(n-1, sizeof(double));
    str->q1 = common_calloc(n-1, sizeof(double));
  }else{
    str->l0 = common_calloc(n,   sizeof(double));
    str->c0 = common_calloc(n,   sizeof(double));
    str->u0 = common_calloc(n,   sizeof(double));
    str->l1 = NULL;
    str->c1 = NULL;
    str->u1 = NULL;
    str->q1 = NULL;
  }
  return str;
}

#define MY_GTSV_B(T) \
static int my_gtsv_b_##T(const int n, double *l, double *c, double *u, T *q){ \
  /* ! divide first row by center-diagonal term ! 2 ! */ \
  u[0] = u[0]/c[0]; \
  q[0] = q[0]/c[0]; \
  /* ! forward sweep ! 11 ! */ \
  for(int i=1; i<=n-1; i++){ \
    double val = c[i]-l[i]*u[i-1]; \
    if(fabs(val) < DBL_EPSILON){ \
      u[i] = 0.; \
      q[i] = 0.; \
    }else{ \
      val = 1./val; \
      u[i] = val* u[i]; \
      q[i] = val*(q[i]-l[i]*q[i-1]); \
    } \
  } \
  /* ! backward sweep ! 3 ! */ \
  for(int i=n-2; i>=0; i--){ \
    q[i] -= u[i]*q[i+1]; \
  } \
  return 0; \
}

MY_GTSV_B(double)
MY_GTSV_B(fftw_complex)

#define MY_GTSV_P(T) \
static int my_gtsv_p_##T(tdm_t *str, T *q){ \
  const int n = str->n; \
  double *l = str->l; \
  double *c = str->c; \
  double *u = str->u; \
  /* 1st small linear system */ \
  double *l0 = str->l0; \
  double *c0 = str->c0; \
  double *u0 = str->u0; \
  T      *q0 = q; \
  /* 2nd small linear system */ \
  double *l1 = str->l1; \
  double *c1 = str->c1; \
  double *u1 = str->u1; \
  double *q1 = str->q1; \
  /* ! create 1st small tri-diagonal system ! 6 ! */ \
  for(int i=0; i<n-1; i++){ \
    l0[i] = l[i]; \
    c0[i] = c[i]; \
    u0[i] = u[i]; \
  } \
  /* ! create 2nd small tri-diagonal system ! 9 ! */ \
  for(int i=0; i<n-1; i++){ \
    l1[i] = l[i]; \
    c1[i] = c[i]; \
    u1[i] = u[i]; \
    q1[i] \
      = i ==   0 ? -l[i] \
      : i == n-2 ? -u[i] \
      : 0.; \
  } \
  /* ! solve two small systems ! 2 ! */ \
  my_gtsv_b_##T   (n-1, l0, c0, u0, q0); \
  my_gtsv_b_double(n-1, l1, c1, u1, q1); \
  /* ! compute bottom solution ! 3 ! */ \
  T      num = q[n-1]-u[n-1]*q0[0]-l[n-1]*q0[n-2]; \
  double den = c[n-1]+u[n-1]*q1[0]+l[n-1]*q1[n-2]; \
  q[n-1] = fabs(den) < DBL_EPSILON ? 0. : num / den; \
  /* ! assign answers of the original system ! 3 ! */ \
  for(int i=0; i<=n-2; i++){ \
    q[i] = q0[i]+q[n-1]*q1[i]; \
  } \
  return 0; \
}

MY_GTSV_P(double)
MY_GTSV_P(fftw_complex)

#define tdm_SOLVE(T) \
int tdm_solve_##T(tdm_t *str, T *q){ \
  if(str->is_periodic){ \
    my_gtsv_p_##T(str, q); \
  }else{ \
    const int n = str->n; \
    double *l0 = str->l0; \
    double *c0 = str->c0; \
    double *u0 = str->u0; \
    memcpy(l0, str->l, n*sizeof(double)); \
    memcpy(c0, str->c, n*sizeof(double)); \
    memcpy(u0, str->u, n*sizeof(double)); \
    my_gtsv_b_##T(n, l0, c0, u0, q); \
  } \
  return 0; \
}

tdm_SOLVE(double)
tdm_SOLVE(fftw_complex)

int tdm_finalise(tdm_t *str){
  common_free(str->l);
  common_free(str->c);
  common_free(str->u);
  common_free(str->l0);
  common_free(str->c0);
  common_free(str->u0);
  if(str->is_periodic){
    common_free(str->l1);
    common_free(str->c1);
    common_free(str->u1);
    common_free(str->q1);
  }
  common_free(str);
  return 0;
}

