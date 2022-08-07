#include <math.h>
#include "common.h"
#include "suspensions.h"
#include "ellipse.h"


// number of grid points outside circles
#define NEXTRA 3
#define POW2(x) ((x)*(x))

int suspensions_decide_loop_size(const int lbound, const int ubound, const double grid_size, const double radius, const double grav_center, int *min, int *max){
  *min = (int)(floor((grav_center-radius)/grid_size)-NEXTRA);
  *max = (int)(ceil ((grav_center+radius)/grid_size)+NEXTRA);
  *min = *min < lbound ? lbound : *min;
  *max = *max > ubound ? ubound : *max;
  return 0;
}

static double compute_signed_dist(const double grid_size, const double pa, const double pb, const double px, const double py, const double paz, const double x, const double y){
  // transform to the origin and cancel rotation
  // an ellise goes to (0, 0) and the rotation leads 0
  // the target point (x, y) moves to (x_, y_)
  double x_, y_, sign;
  {
    double dx = x-px;
    double dy = y-py;
    double c = cos(-paz);
    double s = sin(-paz);
    x_ = c*dx-s*dy;
    y_ = s*dx+c*dy;
  }
  // check whether (x_, y_) is inside an ellipse whose center is at (0, 0) and major/minor axes are pa and pb.
  // inside: sign=+1, outside: sign=-1
  {
    // NOTE: f0 > 0 : inside, f0 < 0 : outside
    double f0 = 1.-POW2(x_/pa)-POW2(y_/pb);
    sign = f0 > 0. ? +1. : -1.;
  }
  double t = find_normal_t(pa, pb, x_, y_);
  // the nearest point on the ellipse: (a cos(t), b sin(t))
  // compute signed distance
  double dist = sign*sqrt(
      +POW2(pa*cos(t)-x_)
      +POW2(pb*sin(t)-y_)
  );
  return dist/grid_size;
}

#define BETA 2.

double suspensions_s_weight(const double grid_size, const double pa, const double pb, const double px, const double py, const double paz, const double x, const double y){
  double dist = compute_signed_dist(grid_size, pa, pb, px, py, paz, x, y);
  return 0.5*BETA*(1.-POW2(tanh(BETA*dist)));
}

double suspensions_v_weight(const double grid_size, const double pa, const double pb, const double px, const double py, const double paz, const double x, const double y){
  double dist = compute_signed_dist(grid_size, pa, pb, px, py, paz, x, y);
  return 0.5*(1.+tanh(BETA*dist));
}

#undef BETA

double suspensions_compute_volume(const double a, const double b){
  return M_PI*a*b;
}

double suspensions_compute_mass(const double den, const double a, const double b){
  double vol = suspensions_compute_volume(a, b);
  return den*vol;
}

double suspensions_compute_moment_of_inertia(const double den, const double a, const double b){
  double mass = suspensions_compute_mass(den, a, b);
  return 0.25*mass*(POW2(a)+POW2(b));
}

#undef NEXTRA
#undef POW2

