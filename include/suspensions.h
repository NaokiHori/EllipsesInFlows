#if !defined(SUSPENSIONS_H)
#define SUSPENSIONS_H

#include <stdbool.h>
#include "structure.h"
#include "arrays/suspensions.h"

typedef struct particle_t_ {
  // fixed parameters
  // density (ratio), major and minor axes
  double den, a, b;
  // gravity center,
  // translational (x, y) and angular (az) positions
  double x, y, az;
  double dx, dy, daz;
  // translational and rotational velocities
  double ux, uy, vz;
  double dux, duy, dvz;
  // surface forces and torque at k+1/2 step
  // 1/m * \int_{S} a_i dS or 1/I * \int_{S} \epsilon_{ijk} \omega_j a_k dS
  // a_i^k = \alpha^k ( U_i^k - u_i^* ) / ( \gamma \Delta t )
  double fux, fuy, tvz;
  // internal inertia, e.g.,
  //   0: 1/C * \int_{V^{k  }} u_i^{k  } dV^{k  }
  //   1: 1/C * \int_{V^{k+1}} u_i^{k+1} dV^{k+1}
  // similar to the rotational component
  // C is a pre-factor, mass or moment of inertia
  double iux[2], iuy[2], ivz[2];
  // collision force based on k (0) and k+1 (1) step particle positions and velocities
  double cfx[2], cfy[2], ctz[2];
} particle_t;

struct suspensions_t_ {
  int n_particles;
  particle_t **particles;
  // responses of surface forces and torque on the momentum fields
  double *dux, *duy;
  // buffers to communicate Lagrange information
  // whose size is sizeof(double) * 3*n_particles
  double *buf;
};

/* constructor and destructor */
extern suspensions_t *suspensions_init(const param_t *param, const parallel_t *parallel);
extern int suspensions_finalise(suspensions_t *suspensions);

/* called by main update routine */
extern int suspensions_reset_particle_increments(suspensions_t *suspensions);
extern int suspensions_compute_inertia(const param_t *param, const parallel_t *parallel, const int cnstep, const fluid_t *fluid, suspensions_t *suspensions);
extern int suspensions_compute_collision_force(const param_t *param, const parallel_t *parallel, const int cnstep, suspensions_t *suspensions);
extern int suspensions_exchange_momentum(const param_t *param, const parallel_t *parallel, const fluid_t *fluid, suspensions_t *suspensions);
extern int suspensions_increment_particles(const param_t *param, const int rkstep, suspensions_t *suspensions, double *residual);
extern int suspensions_update_momentum_fleid(const param_t *param, const parallel_t *parallel, fluid_t *fluid, const suspensions_t *suspensions);
extern int suspensions_update_particles(const param_t *param, suspensions_t *suspensions);

/* other supportive functions */
extern int suspensions_decide_loop_size(const int lbound, const int ubound, const double grid_size, const double radius, const double grav_center, int *min, int *max);
extern double suspensions_compute_volume(const double a, const double b);
extern double suspensions_compute_mass(const double den, const double a, const double b);
extern double suspensions_compute_moment_of_inertia(const double den, const double a, const double b);

extern double suspensions_s_weight(const double grid_size, const double pa, const double pb, const double px, const double py, const double paz, const double x, const double y);
extern double suspensions_v_weight(const double grid_size, const double pa, const double pb, const double px, const double py, const double paz, const double x, const double y);

#endif // SUSPENSIONS_H
