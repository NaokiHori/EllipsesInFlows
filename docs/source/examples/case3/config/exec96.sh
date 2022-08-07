#!/bin/bash

## directory name containing flow fields to restart
# export dirname_restart="output/save/stepxxxxxxxxxx"

## durations
# maximum duration (in free-fall time)
export timemax=1.0e+1
# maximum duration (in wall time [s])
export wtimemax=6.0e+2
# logging rate (in free-fall time)
export log_rate=1.0e-2
# logging after (in free-fall time)
export log_after=0.0e+0
# save rate (in free-fall time)
export save_rate=5.0e+1
# save after (in free-fall time)
export save_after=0.0e+0
# statistics collection rate (in free-fall time)
export stat_rate=1.0e-1
# statistics collection after (in free-fall time)
export stat_after=2.0e+3

## domain
# domain lengths
export ly=2.0e+0
# number of grids
export itot=96
export jtot=192

## dt safe factors, adv and dif
export safefactor_adv=0.9
export safefactor_dif=0.9

## physical parameters
export Re=2.5e+1
# export Fr=1.e+0

## external forcing in y direction
export extfrcy=0.

mpirun -n 2 --oversubscribe ./a.out
