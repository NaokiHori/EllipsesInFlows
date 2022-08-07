#include <stdio.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "suspensions.h"
#include "ellipse.h"


typedef struct ellipse_t_ {
  // major and minor axes
  double a, b;
  // gravity center
  double x, y;
  // rotation angle (from major axis)
  double angle;
} ellipse_t;

typedef struct circle_t_ {
  // center
  double x, y;
  // radius
  double r;
} circle_t;

/*
 * Collision pair matrix
 *
 *    | 0  1    2     3 ...  N-2  N-1  <- n1
 * ---+------------------------------
 *  0 |    0    1     2 ...  N-3  N-2
 *  1 |       N-1     N ... 2N-5 2N-4
 *  2 |            2N-3 ... 3N-8 3N-7
 * ...
 *  ^ |
 *  |
 *  n0
 */

static int get_l(const int n_particles, const int n0){
  // left-most index of the collision pair matrix for particle index "n0"
  return (2*n_particles-n0-1)*(n0  )/2  ;
}

static int get_r(const int n_particles, const int n0){
  // right-most index of the collision pair matrix for particle index "n0"
  return (2*n_particles-n0-2)*(n0+1)/2-1;
}

static int get_particle_indices(const int n_particles, const int n, int *n0, int *n1){
  // convert index of collision pair "n" to the corresponding particle indices "n0" and "n1"
  for(*n0 = 0; ; (*n0)++){
    if(get_l(n_particles, *n0) <= n && n <= get_r(n_particles, *n0)){
      break;
    }
  }
  *n1 = n-get_l(n_particles, *n0)+(*n0)+1;
  return 0;
}

static int get_my_range(const parallel_t *parallel, const int n_total, int *n_min, int *n_max){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  *n_min = parallel_get_offset(n_total, mpisize, mpirank);
  *n_max = parallel_get_size(n_total, mpisize, mpirank) + *n_min;
  return 0;
}

static double harmonic_average(const double v0, const double v1){
  double retval = 0.5*(1./v0+1./v1);
  return 1./retval;
}

static int compute_evolute(const ellipse_t *e, const double x, const double y, double *ex, double *ey, double *ix, double *iy, double *r){
  // e: ellipse
  // (x, y): target point
  // (ex, ey): center of evolute
  // (ix, iy): intersection of fitted circle and ellipse
  // forward transformation
  double x_, y_;
  {
    double x__ = x - e->x;
    double y__ = y - e->y;
    x_ = cos(-e->angle) * x__ - sin(-e->angle) * y__;
    y_ = sin(-e->angle) * x__ + cos(-e->angle) * y__;
  }
  double t = find_normal_t(e->a, e->b, x_, y_);
  *r = 1./fmax(compute_curvature(e->a, e->b, t), DBL_EPSILON);
  // find corresponding evolute
  x_ = compute_ex(e->a, e->b, t);
  y_ = compute_ey(e->a, e->b, t);
  // inverse transformation
  {
    *ex = cos(+e->angle) * x_ - sin(+e->angle) * y_;
    *ey = sin(+e->angle) * x_ + cos(+e->angle) * y_;
    *ex = *ex + e->x;
    *ey = *ey + e->y;
  }
  // find intersection
  x_ = e->a*cos(t);
  y_ = e->b*sin(t);
  // inverse transformation
  {
    *ix = cos(+e->angle) * x_ - sin(+e->angle) * y_;
    *iy = sin(+e->angle) * x_ + cos(+e->angle) * y_;
  }
  return 0;
}

static int find_equivalent_circles(const ellipse_t *e0, const ellipse_t *e1, circle_t *c0, circle_t *c1, double *ix0, double *iy0, double *ix1, double *iy1){
  const double small = 1.e-8;
  c0->x = e0->x;
  c0->y = e0->y;
  c1->x = e1->x;
  c1->y = e1->y;
  for(int n = 0; ; n++){
    double c0_x, c0_y;
    double c1_x, c1_y;
    compute_evolute(e0, c1->x, c1->y, &c0_x, &c0_y, ix0, iy0, &c0->r);
    compute_evolute(e1, c0->x, c0->y, &c1_x, &c1_y, ix1, iy1, &c1->r);
    double dc0_x = c0_x-c0->x;
    double dc0_y = c0_y-c0->y;
    double dc1_x = c1_x-c1->x;
    double dc1_y = c1_y-c1->y;
    c0->x += dc0_x;
    c0->y += dc0_y;
    c1->x += dc1_x;
    c1->y += dc1_y;
    if(fmax(fmax(fabs(dc0_x), fabs(dc0_y)), fmax(fabs(dc1_x), fabs(dc1_y))) < small){
      break;
    }
  }
  return 0;
}

static int find_equivalent_circle(const ellipse_t *e, const double wallx, circle_t *c, double *ix, double *iy){
  const double small = 1.e-8;
  c->x = e->x;
  c->y = e->y;
  for(int n = 0; ; n++){
    double c_x, c_y;
    double mirrorx = 2.*wallx-c->x;
    compute_evolute(e, mirrorx, c->y, &c_x, &c_y, ix, iy, &c->r);
    double dc_x = c_x-c->x;
    double dc_y = c_y-c->y;
    c->x += dc_x;
    c->y += dc_y;
    if(fmax(fabs(dc_x), fabs(dc_y)) < small){
      break;
    }
  }
  return 0;
}

static int compute_collision_force_p_p(const param_t *param, const int cnstep, particle_t *p0, particle_t *p1){
  // correct periodicity
  double yoffset = 0.;
  {
    const double ly = param->ly;
    const double p0y = p0->y+p0->dy;
    const double p1y = p1->y+p1->dy;
    double minval = DBL_MAX;
    for(int periodic = -1; periodic <= 1; periodic++){
      double val = fabs((p1y+ly*periodic) - p0y);
      if(val < minval){
        minval = val;
        yoffset = ly*periodic;
      }
    }
  }
  // check collision of larger circles for early return
  {
    double x0 = p0->x+p0->dx;
    double y0 = p0->y+p0->dy;
    double r0 = fmax(p0->a, p0->b);
    double x1 = p1->x+p1->dx;
    double y1 = p1->y+p1->dy+yoffset;
    double r1 = fmax(p1->a, p1->b);
    double nx = x1-x0;
    double ny = y1-y0;
    double norm = fmax(hypot(nx, ny), DBL_EPSILON);
    nx /= norm;
    ny /= norm;
    double overlap_dist = r0+r1-norm;
    if(overlap_dist < 0.){
      return 0;
    }
  }
  // convert ellipse to circles
  ellipse_t e0 = {
    .a = p0->a,
    .b = p0->b,
    .x = p0->x+p0->dx,
    .y = p0->y+p0->dy,
    .angle = p0->az
  };
  ellipse_t e1 = {
    .a = p1->a,
    .b = p1->b,
    .x = p1->x+p1->dx,
    .y = p1->y+p1->dy+yoffset,
    .angle = p1->az
  };
  circle_t c0, c1;
  double ix0, iy0, ix1, iy1;
  find_equivalent_circles(&e0, &e1, &c0, &c1, &ix0, &iy0, &ix1, &iy1);
  // compute force and torque
  const double p0mass = suspensions_compute_mass(p0->den, p0->a, p0->b);
  const double p1mass = suspensions_compute_mass(p1->den, p1->a, p1->b);
  const double p0im   = suspensions_compute_moment_of_inertia(p0->den, p0->a, p0->b);
  const double p1im   = suspensions_compute_moment_of_inertia(p1->den, p1->a, p1->b);
  // k: pre-factor (spring stiffness)
  double k;
  {
    const double mass = harmonic_average(p0mass, p1mass);
    // equivalent radius based on the area
    const double r0 = sqrt(p0->a*p0->b);
    const double r1 = sqrt(p1->a*p1->b);
    const double reft = pow(harmonic_average(r0, r1), 1.5);
    k = mass*pow(M_PI, 2.)/pow(reft, 2.);
  }
  // compute normal vector from particle 0 to 1 (N.B. periodicity in y)
  {
    double nx = c1.x-c0.x;
    double ny = c1.y-c0.y;
    double norm = fmax(hypot(nx, ny), DBL_EPSILON);
    nx /= norm;
    ny /= norm;
    double overlap_dist = c0.r+c1.r-norm;
    // impose force when overlapped
    if(overlap_dist > 0.){
      double cfx = k*overlap_dist*nx;
      double cfy = k*overlap_dist*ny;
      // force in x
      p0->cfx[cnstep] += 1./p0mass*(-cfx);
      p1->cfx[cnstep] += 1./p1mass*(+cfx);
      // force in y
      p0->cfy[cnstep] += 1./p0mass*(-cfy);
      p1->cfy[cnstep] += 1./p1mass*(+cfy);
      // torque in z
      p0->ctz[cnstep] += 1./p0im*(ix0*(-cfy)-iy0*(-cfx));
      p1->ctz[cnstep] += 1./p1im*(ix1*(+cfy)-iy1*(+cfx));
    }
  }
  return 0;
}

static int compute_collision_force_p_w(const double wallx, const int cnstep, particle_t *p){
  // check collision of larger circles for early return
  {
    double x = p->x+p->dx;
    double r = fmax(p->a, p->b);
    double nx = wallx-x;
    double norm = fabs(nx);
    nx /= norm;
    double overlap_dist = r-norm;
    if(overlap_dist < 0.){
      return 0;
    }
  }
  // convert ellipse to circles
  ellipse_t e = {
    .a = p->a,
    .b = p->b,
    .x = p->x+p->dx,
    .y = p->y+p->dy,
    .angle = p->az
  };
  circle_t c;
  double ix, iy;
  find_equivalent_circle(&e, wallx, &c, &ix, &iy);
  // compute force and torque
  const double pmass = suspensions_compute_mass(p->den, p->a, p->b);
  const double pim   = suspensions_compute_moment_of_inertia(p->den, p->a, p->b);
  // k: pre-factor (spring stiffness)
  double k;
  {
    const double mass = pmass;
    // equivalent radius based on the area
    const double r = sqrt(p->a*p->b);
    const double reft = pow(r, 1.5);
    k = mass*pow(M_PI, 2.)/pow(reft, 2.);
  }
  // compute normal vector from particle 0 to the wall
  {
    double nx = wallx-c.x;
    double norm = fabs(nx);
    nx /= norm;
    double overlap_dist = c.r-norm;
    // impose force when overlapped
    if(overlap_dist > 0.){
      double cfx = k*overlap_dist*nx;
      double cfy = 0.;
      // force in x
      p->cfx[cnstep] += 1./pmass*(-cfx);
      // force in y
      p->cfy[cnstep] += 1./pmass*(-cfy);
      // torque in z
      p->ctz[cnstep] += 1./pim*(ix*(-cfy)-iy*(-cfx));
    }
  }
  return 0;
}

int suspensions_compute_collision_force(const param_t *param, const parallel_t *parallel, const int cnstep, suspensions_t *suspensions){
  /*
   * NOTE: only the spring in the normal direction is considered for simplicity
   * Although this is sufficient to avoid over-penetrations between particles,
   *   obviously not collect from a physical perspective
   */
  const double lx = param->lx;
  const int n_particles = suspensions->n_particles;
  particle_t **particles = suspensions->particles;
  // reset collision force at this CN step
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    p->cfx[cnstep] = 0.;
    p->cfy[cnstep] = 0.;
    p->ctz[cnstep] = 0.;
  }
  // particle-particle collisions
  {
    const int n_total = n_particles*(n_particles-1)/2;
    int n_min, n_max;
    get_my_range(parallel, n_total, &n_min, &n_max);
    for(int n = n_min; n < n_max; n++){
      int n0, n1;
      get_particle_indices(n_particles, n, &n0, &n1);
      compute_collision_force_p_p(param, cnstep, particles[n0], particles[n1]);
    }
  }
  // particle-wall collisions
  {
    double wall_locations[2] = {0., lx};
    for(int wall_index = 0; wall_index < 2; wall_index++){
      double wall_location = wall_locations[wall_index];
      int n_min, n_max;
      get_my_range(parallel, n_particles, &n_min, &n_max);
      for(int n = n_min; n < n_max; n++){
        compute_collision_force_p_w(wall_location, cnstep, particles[n]);
      }
    }
  }
  // synchronise computed forcings
  {
    // prepare message buffer
    double *buf = suspensions->buf;
    // pack
    for(int n = 0; n < n_particles; n++){
      particle_t *p = particles[n];
      buf[3*n+0] = p->cfx[cnstep];
      buf[3*n+1] = p->cfy[cnstep];
      buf[3*n+2] = p->ctz[cnstep];
    }
    // sum up all
    MPI_Allreduce(MPI_IN_PLACE, buf, 3*n_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // unpack
    for(int n = 0; n < n_particles; n++){
      particle_t *p = particles[n];
      p->cfx[cnstep] = buf[3*n+0];
      p->cfy[cnstep] = buf[3*n+1];
      p->ctz[cnstep] = buf[3*n+2];
    }
  }
  return 0;
}

