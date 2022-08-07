#include <math.h>
#include <float.h>
#include "common.h"
#include "ellipse.h"


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

// find t, with which a vector from (a cos(t), b sin(t)) to (x_, y_)
//   becomes a normal vector to the ellipse
// https://blog.chatfield.io/simple-method-for-distance-to-ellipse/
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

