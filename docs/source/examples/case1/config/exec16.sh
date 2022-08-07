#!/bin/bash

## directory name containing flow fields to restart
# export dirname_restart="output/save/stepxxxxxxxxxx"

## durations
# maximum duration (in free-fall time)
export timemax=1.5e+2
# maximum duration (in wall time [s])
export wtimemax=6.0e+2
# logging rate (in free-fall time)
export log_rate=1.0e+0
# logging after (in free-fall time)
export log_after=0.0e+0
# save rate (in free-fall time)
export save_rate=1.0e+3
# save after (in free-fall time)
export save_after=0.0e+0
# statistics collection rate (in free-fall time)
export stat_rate=1.0e-1
# statistics collection after (in free-fall time)
export stat_after=2.0e+3

## domain
# domain lengths
export ly=1.0e+0
# number of grids
export itot=16
export jtot=16

## dt safe factors, adv and dif
export safefactor_adv=0.9
export safefactor_dif=0.9

## physical parameters
export Re=4.e+1
# export Fr=1.e+0

## external forcing in y direction
export extfrcy=0.

mpirun -n 1 --oversubscribe ./a.out
