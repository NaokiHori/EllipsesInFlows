#include <string.h>
#include <math.h>
#include "common.h"
#include "param.h"
#include "parallel.h"
#include "fluid.h"


static int compute_src_ux(const param_t *param, const parallel_t *parallel, fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double *dxf = param->dxf;
  const double *dxc = param->dxc;
  const double dy = param->dy;
  const double Re = param->Re;
  const double * restrict ux = fluid->ux;
  const double * restrict uy = fluid->uy;
  const double * restrict p = fluid->p;
  double * restrict srcuxa = fluid->srcuxa;
  double * restrict srcuxb = fluid->srcuxb;
  double * restrict srcuxg = fluid->srcuxg;
  /* ! previous k-step source term of ux is copied ! 1 ! */
  memcpy(srcuxb, srcuxa, SRCUXA_MEMSIZE);
  // UX(i=1, j) and UX(itot+1, j) are fixed to 0
  for(int j=1; j<=jsize; j++){
    for(int i=2; i<=itot; i++){
      /* ! velocity-gradient tensor L_xx ! 2 ! */
      double duxdx_xm = (-UX(i-1, j  )+UX(i  , j  ))/DXF(i-1);
      double duxdx_xp = (-UX(i  , j  )+UX(i+1, j  ))/DXF(i  );
      /* ! velocity-gradient tensor L_xy ! 2 ! */
      double duxdy_ym = (-UX(i  , j-1)+UX(i  , j  ))/dy;
      double duxdy_yp = (-UX(i  , j  )+UX(i  , j+1))/dy;
      /* advection */
      /* ! x-momentum is transported by ux ! 10 ! */
      double adv1;
      {
        double c_xm = DXF(i-1)/(2.*DXC(i));
        double c_xp = DXF(i  )/(2.*DXC(i));
        double ux_xm = 0.5*UX(i-1, j  )+0.5*UX(i  , j  );
        double ux_xp = 0.5*UX(i  , j  )+0.5*UX(i+1, j  );
        adv1 =
          -c_xm*ux_xm*duxdx_xm
          -c_xp*ux_xp*duxdx_xp;
      }
      /* ! x-momentum is transported by uy ! 10 ! */
      double adv2;
      {
        double c_xm = DXF(i-1)/(2.*DXC(i));
        double c_xp = DXF(i  )/(2.*DXC(i));
        double uy_ym = c_xm*UY(i-1, j  )+c_xp*UY(i  , j  );
        double uy_yp = c_xm*UY(i-1, j+1)+c_xp*UY(i  , j+1);
        adv2 =
          -0.5*uy_ym*duxdy_ym
          -0.5*uy_yp*duxdy_yp;
      }
      /* diffusion */
      /* ! x-momentum is diffused in x ! 5 ! */
      double dif1;
      {
        const double prefactor = 1./Re;
        dif1 = prefactor*(duxdx_xp-duxdx_xm)/DXC(i);
      }
      /* ! x-momentum is diffused in y ! 5 ! */
      double dif2;
      {
        const double prefactor = 1./Re;
        dif2 = prefactor*(duxdy_yp-duxdy_ym)/dy;
      }
      /* pressure gradient */
      /* ! pressure gradient in x ! 6 ! */
      double pre;
      {
        double p_xm = P(i-1, j  );
        double p_xp = P(i  , j  );
        pre = -(p_xp-p_xm)/DXC(i);
      }
      /* summation */
      /* ! summation of ux explicit terms ! 3 ! */
      SRCUXA(i, j) =
        +(adv1+adv2)
        +(dif1+dif2);
      /* ! summation of ux implicit terms ! 1 ! */
      SRCUXG(i, j) = pre;
    }
  }
  return 0;
}

static int compute_src_uy(const param_t *param, const parallel_t *parallel, fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double *dxf = param->dxf;
  const double *dxc = param->dxc;
  const double dy = param->dy;
  const double Re = param->Re;
  const double ext = param->extfrcy;
  const double * restrict ux = fluid->ux;
  const double * restrict uy = fluid->uy;
  const double * restrict p = fluid->p;
  double * restrict srcuya = fluid->srcuya;
  double * restrict srcuyb = fluid->srcuyb;
  double * restrict srcuyg = fluid->srcuyg;
  /* ! previous k-step source term of uy is copied ! 1 ! */
  memcpy(srcuyb, srcuya, SRCUYA_MEMSIZE);
  /* ! uy is computed from i=1 to itot ! 2 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      /* ! velocity-gradient tensor L_yx ! 2 ! */
      double duydx_xm = (-UY(i-1, j  )+UY(i  , j  ))/DXC(i  );
      double duydx_xp = (-UY(i  , j  )+UY(i+1, j  ))/DXC(i+1);
      /* ! velocity-gradient tensor L_yy ! 2 ! */
      double duydy_ym = (-UY(i  , j-1)+UY(i  , j  ))/dy;
      double duydy_yp = (-UY(i  , j  )+UY(i  , j+1))/dy;
      /* advection */
      /* ! y-momentum is transported by ux ! 10 ! */
      double adv1;
      {
        double c_xm = DXC(i  )/(2.*DXF(i));
        double c_xp = DXC(i+1)/(2.*DXF(i));
        double ux_xm = 0.5*UX(i  , j-1)+0.5*UX(i  , j  );
        double ux_xp = 0.5*UX(i+1, j-1)+0.5*UX(i+1, j  );
        adv1 =
          -c_xm*ux_xm*duydx_xm
          -c_xp*ux_xp*duydx_xp;
      }
      /* ! y-momentum is transported by uy ! 8 ! */
      double adv2;
      {
        double uy_ym = 0.5*UY(i  , j-1)+0.5*UY(i  , j  );
        double uy_yp = 0.5*UY(i  , j  )+0.5*UY(i  , j+1);
        adv2 =
          -0.5*uy_ym*duydy_ym
          -0.5*uy_yp*duydy_yp;
      }
      /* diffusion */
      /* ! y-momentum is diffused in x ! 5 ! */
      double dif1;
      {
        const double prefactor = 1./Re;
        dif1 = prefactor*(duydx_xp-duydx_xm)/DXF(i);
      }
      /* ! y-momentum is diffused in y ! 5 ! */
      double dif2;
      {
        const double prefactor = 1./Re;
        dif2 = prefactor*(duydy_yp-duydy_ym)/dy;
      }
      /* pressure gradient */
      /* ! pressure gradient in y ! 6 ! */
      double pre;
      {
        double p_ym = P(i  , j-1);
        double p_yp = P(i  , j  );
        pre = -(p_yp-p_ym)/dy;
      }
      /* summation */
      /* ! summation of uy explicit terms ! 4 ! */
      SRCUYA(i, j) =
        +(adv1+adv2)
        +(dif1+dif2)
        +ext;
      /* ! summation of uy implicit terms ! 1 ! */
      SRCUYG(i, j) = pre;
    }
  }
  return 0;
}

static int update_ux(const param_t *param, const parallel_t *parallel, const int rkstep, fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double alpha = param->rkcoefs[rkstep].alpha;
  const double beta  = param->rkcoefs[rkstep].beta;
  const double gamma = param->rkcoefs[rkstep].gamma;
  const double dt = param->dt;
  const double *srcuxa = fluid->srcuxa;
  const double *srcuxb = fluid->srcuxb;
  const double *srcuxg = fluid->srcuxg;
  double *ux = fluid->ux;
  /* ! compute increments of ux ! 8 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=2; i<=itot; i++){
      UX(i, j) +=
        +alpha*dt*SRCUXA(i, j)
        +beta *dt*SRCUXB(i, j)
        +gamma*dt*SRCUXG(i, j);
    }
  }
  fluid_update_boundaries_ux(param, parallel, ux);
  return 0;
}

static int update_uy(const param_t *param, const parallel_t *parallel, const int rkstep, fluid_t *fluid){
  const int mpisize = parallel->mpisize;
  const int mpirank = parallel->mpirank;
  const int itot = param->itot;
  const int jtot = param->jtot;
  const int jsize = parallel_get_size(jtot, mpisize, mpirank);
  const double alpha = param->rkcoefs[rkstep].alpha;
  const double beta  = param->rkcoefs[rkstep].beta;
  const double gamma = param->rkcoefs[rkstep].gamma;
  const double dt = param->dt;
  const double *srcuya = fluid->srcuya;
  const double *srcuyb = fluid->srcuyb;
  const double *srcuyg = fluid->srcuyg;
  double *uy = fluid->uy;
  /* ! compute increments of uy ! 8 ! */
  for(int j=1; j<=jsize; j++){
    for(int i=1; i<=itot; i++){
      UY(i, j) +=
        +alpha*dt*SRCUYA(i, j)
        +beta *dt*SRCUYB(i, j)
        +gamma*dt*SRCUYG(i, j);
    }
  }
  fluid_update_boundaries_uy(param, parallel, uy);
  return 0;
}

int fluid_update_velocity(const param_t *param, const parallel_t *parallel, const int rkstep, fluid_t *fluid){
  /* ! source terms of Runge-Kutta scheme are updated ! 2 ! */
  compute_src_ux(param, parallel, fluid);
  compute_src_uy(param, parallel, fluid);
  /* ! velocities are updated ! 2 ! */
  update_ux(param, parallel, rkstep, fluid);
  update_uy(param, parallel, rkstep, fluid);
  return 0;
}

