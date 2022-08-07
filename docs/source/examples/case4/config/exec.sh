#!/bin/bash

## directory name containing flow fields to restart
# export dirname_restart=

## durations
# maximum duration (in free-fall time)
export timemax=1.0e+3
# maximum duration (in wall time [s])
export wtimemax=8.6e+4
# logging rate (in free-fall time)
export log_rate=1.0e-1
# logging after (in free-fall time)
export log_after=0.0e+0
# save rate (in free-fall time)
export save_rate=1.0e+2
# save after (in free-fall time)
export save_after=0.0e+0
# statistics collection rate (in free-fall time)
export stat_rate=1.0e-1
# statistics collection after (in free-fall time)
export stat_after=2.0e+3

## domain
# domain lengths
export ly=4.0e+0
# number of grids
export itot=256
export jtot=1024

## dt safe factors, adv and dif
export safefactor_adv=0.95
export safefactor_dif=0.95

## physical parameters
export Re=1.0e+3

## external forcing in y direction
export extfrcy=8.0e-3

mpirun -n 4 ./a.out
