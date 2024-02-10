#if !defined(ELLIPSE_H)
#define ELLIPSE_H

extern double ellipse_compute_xc(const double a, const double b, const double t);
extern double ellipse_compute_yc(const double a, const double b, const double t);
extern double ellipse_compute_curvature(const double a, const double b, const double t);
extern double ellipse_find_normal_t(const double a, const double b, const double xp, const double yp);

#endif // ELLIPSE_H
