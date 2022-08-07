#if !defined(ELLIPSE_H)
#define ELLIPSE_H

extern double compute_ex(const double a, const double b, const double t);
extern double compute_ey(const double a, const double b, const double t);
extern double compute_curvature(const double a, const double b, const double t);
extern double find_normal_t(const double a, const double b, const double x_, const double y_);

#endif // ELLIPSE_H
