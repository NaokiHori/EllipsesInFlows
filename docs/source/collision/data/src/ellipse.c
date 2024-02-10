#include <math.h>
#include "others.h"
#include "ellipse.h"


// compute x-coordinate of evolute
double ellipse_compute_xc(const double a, const double b, const double t){
  return a*(1.-pow(b/a, 2.))*pow(cos(t), 3.);
}

// compute y-coordinate of evolute
double ellipse_compute_yc(const double a, const double b, const double t){
  return b*(1.-pow(a/b, 2.))*pow(sin(t), 3.);
}

double ellipse_compute_curvature(const double a, const double b, const double t){
  const double num = a*b;
  const double den = pow(pow(a*sin(t), 2.)+pow(b*cos(t), 2.), 1.5);
  return num/den;
}

// https://blog.chatfield.io/simple-method-for-distance-to-ellipse/
double ellipse_find_normal_t(const double a, const double b, const double xp, const double yp){
  /* ! consider in the 1st quadrant ! 2 ! */
  double abs_xp = fabs(xp);
  double abs_yp = fabs(yp);
  /* ! initialise t ! 1 ! */
  double t = 0.25*M_PI;
  // iterate until converged
  while(1){
    /* ! compute point on the ellipse ! 2 ! */
    double xe = a*cos(t);
    double ye = b*sin(t);
    /* ! compute center of the fitted circle ! 2 ! */
    double xc = ellipse_compute_xc(a, b, t);
    double yc = ellipse_compute_yc(a, b, t);
    /* ! compute residual ! 8 ! */
    double dxe =     xe-xc;
    double dye =     ye-yc;
    double dxp = abs_xp-xc;
    double dyp = abs_yp-yc;
    double norme = hypot(dxe, dye);
    double normp = hypot(dxp, dyp);
    double dc = norme*asin((dxe*dyp-dye*dxp)/(norme*normp));
    double dt = dc/sqrt(a*a+b*b-xe*xe-ye*ye);
    /* ! update t ! 4 ! */
    t += dt;
    // limit to the 1st quadrant
    t = fmin(0.5*M_PI, t);
    t = fmax(      0., t);
    /* ! terminate iteration when the residual is sufficiently small ! 3 ! */
    if(fabs(dt) < 1.e-8){
      break;
    }
  }
  /* ! recover information of the quadrants ! 2 ! */
  if(xp < 0.) t = M_PI-t;
  if(yp < 0.) t =     -t;
  /* ! return final t ! 1 ! */
  return t;
}

