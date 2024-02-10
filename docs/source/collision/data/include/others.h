#if !defined(OTHERS_H)
#define OTHERS_H

#if !defined(M_PI)
#define M_PI 3.1415926535897932
#endif

extern int shift_and_rotate(const double dx, const double dy, const double dt, const double x0, const double y0, double *x1, double *y1);
extern int rotate_and_shift(const double dx, const double dy, const double dt, const double x0, const double y0, double *x1, double *y1);
extern double get_random_value(const double min, const double max);

#endif // OTHERS_H
