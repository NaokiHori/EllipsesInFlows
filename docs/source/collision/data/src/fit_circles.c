#include <stdio.h>
#include <stdlib.h>
#include "others.h"
#include "ellipse.h"
#include "example.h"


int fit_circles(void){
  /* ! initialise geometry of ellipses ! 6 ! */
  const double e0_a     = 2.0;
  const double e0_b     = 1.5;
  const double e0_theta = 0.2;
  const double e1_a     = 1.5;
  const double e1_b     = 1.0;
  const double e1_theta = 2.0;
  for(int index = 0; index < 2; index++){
    /* ! initialise centers of ellipses ! 4 ! */
    const double e0_x = index == 0 ? 0.5 : 0.2;
    const double e0_y = index == 0 ? 0.5 : 0.2;
    const double e1_x = index == 0 ? 2.0 : 2.3;
    const double e1_y = index == 0 ? 2.5 : 2.8;
    /* ! initialise evolutes using the centers of ellipses ! 4 ! */
    double e0_xc = e0_x;
    double e0_yc = e0_y;
    double e1_xc = e1_x;
    double e1_yc = e1_y;
    double e0_r, e1_r;
    for(int n = 0; n < 10; n++){
      /* ! transform coordinates, forward ! 4 ! */
      double e0_xc_, e0_yc_;
      double e1_xc_, e1_yc_;
      shift_and_rotate(-e0_x, -e0_y, -e0_theta, e1_xc, e1_yc, &e0_xc_, &e0_yc_);
      shift_and_rotate(-e1_x, -e1_y, -e1_theta, e0_xc, e0_yc, &e1_xc_, &e1_yc_);
      /* ! find desired t ! 2 ! */
      double e0_t = ellipse_find_normal_t(e0_a, e0_b, e0_xc_, e0_yc_);
      double e1_t = ellipse_find_normal_t(e1_a, e1_b, e1_xc_, e1_yc_);
      /* ! update center of fitted circles ! 4 ! */
      e0_xc_ = ellipse_compute_xc(e0_a, e0_b, e0_t);
      e0_yc_ = ellipse_compute_yc(e0_a, e0_b, e0_t);
      e1_xc_ = ellipse_compute_xc(e1_a, e1_b, e1_t);
      e1_yc_ = ellipse_compute_yc(e1_a, e1_b, e1_t);
      /* ! transform coordinates, backward ! 2 ! */
      rotate_and_shift(+e0_x, +e0_y, +e0_theta, e0_xc_, e0_yc_, &e0_xc, &e0_yc);
      rotate_and_shift(+e1_x, +e1_y, +e1_theta, e1_xc_, e1_yc_, &e1_xc, &e1_yc);
      // radii of the fitted circles
      e0_r = 1./ellipse_compute_curvature(e0_a, e0_b, e0_t);
      e1_r = 1./ellipse_compute_curvature(e1_a, e1_b, e1_t);
    }
    char fname[128];
    sprintf(fname, "fit-circles-%d.dat", index);
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "e0_a     = % .15f\n", e0_a    );
    fprintf(fp, "e0_b     = % .15f\n", e0_b    );
    fprintf(fp, "e0_x     = % .15f\n", e0_x    );
    fprintf(fp, "e0_y     = % .15f\n", e0_y    );
    fprintf(fp, "e0_theta = % .15f\n", e0_theta);
    fprintf(fp, "e0_xc    = % .15f\n", e0_xc   );
    fprintf(fp, "e0_yc    = % .15f\n", e0_yc   );
    fprintf(fp, "e0_r     = % .15f\n", e0_r    );
    fprintf(fp, "e1_a     = % .15f\n", e1_a    );
    fprintf(fp, "e1_b     = % .15f\n", e1_b    );
    fprintf(fp, "e1_x     = % .15f\n", e1_x    );
    fprintf(fp, "e1_y     = % .15f\n", e1_y    );
    fprintf(fp, "e1_theta = % .15f\n", e1_theta);
    fprintf(fp, "e1_xc    = % .15f\n", e1_xc   );
    fprintf(fp, "e1_yc    = % .15f\n", e1_yc   );
    fprintf(fp, "e1_r     = % .15f\n", e1_r    );
    fclose(fp);
  }
  return 0;
}



