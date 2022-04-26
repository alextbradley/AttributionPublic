#!/bin/bash

# clean run directory and link all required files
./prep_run.sh

# set job id
JOBNO=915

# record start times
TIMEQSTART="$(date +%s)"
echo Start-time `date` >> times


# submit the job
sbatch -J ATTR_$JOBNO \
       -A $HECACC \
       --export HECACC=$HECACC,TIMEQSTART=$TIMEQSTART,IMGNAME=$IMGNAME,JDEPOT=$JDEPOT,JOBNO=$JOBNO,TIMEQSTART=$TIMEQSTART \
       ./run_repeat.sh

