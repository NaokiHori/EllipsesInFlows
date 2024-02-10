#include <stdlib.h>
#include <math.h>
#include "others.h"


// coordinate transformations
int shift_and_rotate(const double dx, const double dy, const double dt, const double x0, const double y0, double *x1, double *y1){
  *x1 = cos(dt)*(x0+dx) - sin(dt)*(y0+dy);
  *y1 = sin(dt)*(x0+dx) + cos(dt)*(y0+dy);
  return 0;
}

int rotate_and_shift(const double dx, const double dy, const double dt, const double x0, const double y0, double *x1, double *y1){
  *x1 = dx + cos(dt)*x0 - sin(dt)*y0;
  *y1 = dy + sin(dt)*x0 + cos(dt)*y0;
  return 0;
}

double get_random_value(const double min, const double max){
  return (max-min)*(1.*rand()/RAND_MAX)+min;
}

