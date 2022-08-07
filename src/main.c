#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>
#include "param.h"
#include "parallel.h"
#include "fluid.h"
#include "suspensions.h"
#include "statistics.h"
#include "save.h"
#include "logging.h"


static int integrate(const param_t *param, const parallel_t *parallel, fluid_t *fluid, suspensions_t *suspensions){
  for(int rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
    suspensions_reset_particle_increments(suspensions);
    /* ! update boundary and halo values ! 3 ! */
    fluid_update_boundaries_ux(param, parallel, fluid->ux);
    fluid_update_boundaries_uy(param, parallel, fluid->uy);
    fluid_update_boundaries_p(param, parallel, fluid->p);
    // \int_{Vp^k} u_i^k d{Vp^k}
    suspensions_compute_inertia(param, parallel, 0, fluid, suspensions);
    // u_i^k -> u_i^*
    fluid_update_velocity(param, parallel, rkstep, fluid);
    // \alpha ( U_i - u_i^* ) / Delta t
    suspensions_exchange_momentum(param, parallel, fluid, suspensions);
    // u_i^* -> u_i^{**}
    suspensions_update_momentum_fleid(param, parallel, fluid, suspensions);
    /* ! compute scalar potential ! 1 ! */
    fluid_compute_potential(param, parallel, rkstep, fluid);
    /* ! correct velocity to be solenoidal ! 1 ! */
    fluid_correct_velocity(param, parallel, rkstep, fluid);
    /* ! update pressure ! 1 ! */
    fluid_update_pressure(param, parallel, fluid);
    /*** update particles iteratively ***/
    suspensions_compute_collision_force(param, parallel, 0, suspensions);
    for(int substep = 0; ; substep++){
      // \int_{Vp^{k+1}} u_i^{k+1} d{Vp^{k+1}}
      suspensions_compute_inertia(param, parallel, 1, fluid, suspensions);
      suspensions_compute_collision_force(param, parallel, 1, suspensions);
      const double residual_max = 1.e-8;
      double residual;
      suspensions_increment_particles(param, rkstep, suspensions, &residual);
      if(residual < residual_max){
        break;
      }
      const int substep_max = 100;
      if(substep > substep_max){
        break;
      }
    }
    // U_i^{k} -> U_i^{k+1}
    suspensions_update_particles(param, suspensions);
  }
  return 0;
}

int main(void){
  /* ! launch MPI, start timer ! 3 ! */
  MPI_Init(NULL, NULL);
  double wtimes[2] = {0.};
  wtimes[0] = parallel_get_wtime(MPI_MIN);
  /* ! initialise structures ! 5 ! */
  param_t       *param       = param_init();
  parallel_t    *parallel    = parallel_init();
  fluid_t       *fluid       = fluid_init(param, parallel);
  suspensions_t *suspensions = suspensions_init(param, parallel);
  statistics_t  *statistics  = statistics_init(param, parallel);
  /* main loop */
  for(;;){
    /* ! decide time step size ! 1 ! */
    param_decide_dt(param, parallel, fluid, suspensions);
    /* ! integrate mass, momentum, and motions of suspensions in time ! 1 ! */
    integrate(param, parallel, fluid, suspensions);
    /* ! step and time are incremented ! 2 ! */
    param->step += 1;
    param->time += param->dt;
    /* ! output log ! 4 ! */
    if(param->log.next < param->time){
      logging(param, parallel, fluid, suspensions);
      param->log.next += param->log.rate;
    }
    /* ! save flow fields ! 4 ! */
    if(param->save.next < param->time){
      save(param, parallel, fluid, suspensions);
      param->save.next += param->save.rate;
    }
    /* ! collect statistics ! 4 ! */
    if(param->stat.next < param->time){
      statistics_collect(param, parallel, fluid, suspensions, statistics);
      param->stat.next += param->stat.rate;
    }
    /* ! terminate when the simulation is finished ! 3 ! */
    if(param->time > param->timemax){
      break;
    }
    /* ! terminate when wall time limit is reached ! 4 ! */
    wtimes[1] = parallel_get_wtime(MPI_MAX);
    if(wtimes[1]-wtimes[0] > param->wtimemax){
      break;
    }
  }
  /* ! check duration ! 3 ! */
  if(parallel->mpirank == 0){
    printf("elapsed: %.2f [s]\n", wtimes[1]-wtimes[0]);
  }
  /* ! save restart file and statistics at last ! 2 ! */
  save(param, parallel, fluid, suspensions);
  statistics_output(param, parallel, statistics);
  /* ! finalise structures ! 5 ! */
  statistics_finalise(statistics);
  suspensions_finalise(suspensions);
  fluid_finalise(fluid);
  parallel_finalise(parallel);
  param_finalise(param);
  /* ! finalise MPI ! 1 ! */
  MPI_Finalize();
  return 0;
}

