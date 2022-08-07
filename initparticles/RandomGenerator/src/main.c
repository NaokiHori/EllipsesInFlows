#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include "fileio.h"


typedef struct ellipse_t_ {
  double x, y;
  double a, b, angle;
} ellipse_t;

typedef struct circle_t_ {
  double x, y;
  double r;
} circle_t;

static double gen_random(const double min, const double max){
  return (max-min)*rand()/RAND_MAX+min;
}

static int check_stats(double *vfrac, const double lx, const double ly, const int n_particles, const double *as, const double *bs){
  *vfrac = 0.;
  for(int n = 0; n < n_particles; n++){
    *vfrac += M_PI*as[n]*bs[n];
  }
  *vfrac /= lx*ly;
  return 0;
}

static int output(const int n_particles, const double *dens, const double *as, const double *bs, const double *xs, const double *ys, const double *azs, const double *uxs, const double *uys, const double *vzs){
  fileio_w_0d_serial("..", "n_particles",   NPYIO_INT,    sizeof(int),    &n_particles);
  fileio_w_1d_serial("..", "particle_dens", NPYIO_DOUBLE, sizeof(double), n_particles, dens);
  fileio_w_1d_serial("..", "particle_as",   NPYIO_DOUBLE, sizeof(double), n_particles,   as);
  fileio_w_1d_serial("..", "particle_bs",   NPYIO_DOUBLE, sizeof(double), n_particles,   bs);
  fileio_w_1d_serial("..", "particle_xs",   NPYIO_DOUBLE, sizeof(double), n_particles,   xs);
  fileio_w_1d_serial("..", "particle_ys",   NPYIO_DOUBLE, sizeof(double), n_particles,   ys);
  fileio_w_1d_serial("..", "particle_azs",  NPYIO_DOUBLE, sizeof(double), n_particles,  azs);
  fileio_w_1d_serial("..", "particle_uxs",  NPYIO_DOUBLE, sizeof(double), n_particles,  uxs);
  fileio_w_1d_serial("..", "particle_uys",  NPYIO_DOUBLE, sizeof(double), n_particles,  uys);
  fileio_w_1d_serial("..", "particle_vzs",  NPYIO_DOUBLE, sizeof(double), n_particles,  vzs);
  return 0;
}

double compute_ex(const double a, const double b, const double t){
  return a*(1.-pow(b/a, 2.))*pow(cos(t), 3.);
}

double compute_ey(const double a, const double b, const double t){
  return b*(1.-pow(a/b, 2.))*pow(sin(t), 3.);
}

double compute_curvature(const double a, const double b, const double t){
  const double cost = cos(t);
  const double sint = sin(t);
  const double cost2 = pow(cost, 2.);
  const double sint2 = pow(sint, 2.);
  const double a2 = pow(a, 2.);
  const double b2 = pow(b, 2.);
  const double num = a*b*(cost2+sint2);
  const double den = pow(+a2*sint2+b2*cost2, 1.5);
  return num/den;
}

double find_normal_t(const double a, const double b, const double x_, const double y_){
  double px = fabs(x_);
  double py = fabs(y_);
  double t = 0.25*M_PI;
  for(int n = 0; ; n++){
    // point on ellipse
    double x = a*cos(t);
    double y = b*sin(t);
    // evolute of the given t
    double ex = compute_ex(a, b, t);
    double ey = compute_ey(a, b, t);
    // determine dt
    double rx =  x-ex;
    double ry =  y-ey;
    double qx = px-ex;
    double qy = py-ey;
    double r = sqrt(pow(ry, 2.)+pow(rx, 2.));
    double q = sqrt(pow(qy, 2.)+pow(qx, 2.));
    double dc = r*asin((rx*qy-ry*qx)/(r*q));
    double dt = dc/sqrt(a*a+b*b-x*x-y*y);
    t += dt;
    // consider in 1st quadrant
    t = fmin(0.5*M_PI, t);
    t = fmax(0.,       t);
    if(fabs(dt) < 1.e-8){
      break;
    }
    if(n > 100){
      break;
    }
  }
  // for the other quadrants
  if(x_ < 0.) t = M_PI-t;
  if(y_ < 0.) t =     -t;
  return t;
}

static int compute_evolute(const ellipse_t *e, const double x, const double y, double *ex, double *ey, double *r){
  // e: ellipse
  // (x, y): target point
  // (ex, ey): center of evolute
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
  return 0;
}

static int find_equivalent_circles(const ellipse_t *e0, const ellipse_t *e1, circle_t *c0, circle_t *c1){
  const double small = 1.e-8;
  c0->x = e0->x;
  c0->y = e0->y;
  c1->x = e1->x;
  c1->y = e1->y;
  for(int n = 0; ; n++){
    double c0_x, c0_y;
    double c1_x, c1_y;
    compute_evolute(e0, c1->x, c1->y, &c0_x, &c0_y, &c0->r);
    compute_evolute(e1, c0->x, c0->y, &c1_x, &c1_y, &c1->r);
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

int main(void){
  srand(1);
  const double lx = 1.;
  const double ly = 4.;
  const int n_particles = 256;
  double *dens = common_calloc(n_particles, sizeof(double));
  double *as   = common_calloc(n_particles, sizeof(double));
  double *bs   = common_calloc(n_particles, sizeof(double));
  double *xs   = common_calloc(n_particles, sizeof(double));
  double *ys   = common_calloc(n_particles, sizeof(double));
  double *azs  = common_calloc(n_particles, sizeof(double));
  double *uxs  = common_calloc(n_particles, sizeof(double));
  double *uys  = common_calloc(n_particles, sizeof(double));
  double *vzs  = common_calloc(n_particles, sizeof(double));
  // density
  for(int n = 0; n < n_particles; n++){
    dens[n] = 1.;
  }
  // radius
  for(int n = 0; n < n_particles; n++){
    as[n] = gen_random(0.05, 0.05);
    bs[n] = as[n]*gen_random(0.50, 0.50);
    azs[n] = gen_random(0., 2.*M_PI);
  }
  // position, no overlaps
  for(int n0 = 0; n0 < n_particles; n0++){
regen:
    {
      double x0 = gen_random(fmax(as[n0], bs[n0]), lx-fmax(as[n0], bs[n0]));
      double y0 = gen_random(                  0.,                      ly);
      for(int n1 = 0; n1 < n0; n1++){
        // correct periodicity
        double yoffset = 0.;
        {
          double minval = DBL_MAX;
          for(int periodic = -1; periodic <= 1; periodic++){
            double val = fabs((ys[n1]+ly*periodic)-y0);
            if(val < minval){
              minval = val;
              yoffset = ly*periodic;
            }
          }
        }
        // check collision of larger circles for early return
        {
          double r0 = fmax(as[n0], bs[n0]);
          double x1 = xs[n1];
          double y1 = ys[n1]+yoffset;
          double r1 = fmax(as[n1], bs[n1]);
          double nx = x1-x0;
          double ny = y1-y0;
          double norm = fmax(hypot(nx, ny), DBL_EPSILON);
          nx /= norm;
          ny /= norm;
          double dist = norm-r0-r1;
          if(dist < 0.){
            goto regen;
          }
        }
        // convert ellipse to circles and check circles' overlap
        {
          ellipse_t e0 = {
            .a = as[n0],
            .b = bs[n0],
            .x = x0,
            .y = y0,
            .angle = azs[n0]
          };
          ellipse_t e1 = {
            .a = as[n1],
            .b = bs[n1],
            .x = xs[n1],
            .y = ys[n1]+yoffset,
            .angle = azs[n1]
          };
          circle_t c0, c1;
          find_equivalent_circles(&e0, &e1, &c0, &c1);
          double nx = c1.x-c0.x;
          double ny = c1.y-c0.y;
          double norm = fmax(hypot(nx, ny), DBL_EPSILON);
          nx /= norm;
          ny /= norm;
          double dist = norm-c0.r-c1.r;
          if(dist < 0.){
            goto regen;
          }
        }
      }
      xs[n0] = x0;
      ys[n0] = y0;
      printf("particle %*d @ % .3f % .3f\n", 5, n0, x0, y0);
    }
  }
  {
    double vfrac;
    check_stats(&vfrac, lx, ly, n_particles, as, bs);
    printf("volume fraction: % .1e\n", vfrac);
  }
  // velocities, still
  memset(uxs, 0, sizeof(double)*n_particles);
  memset(uys, 0, sizeof(double)*n_particles);
  memset(vzs, 0, sizeof(double)*n_particles);
  // output
  output(n_particles, dens, as, bs, xs, ys, azs, uxs, uys, vzs);
  // clean-up buffers
  common_free(dens);
  common_free(  as);
  common_free(  xs);
  common_free(  ys);
  common_free( azs);
  common_free( uxs);
  common_free( uys);
  common_free( vzs);
  return 0;
}

