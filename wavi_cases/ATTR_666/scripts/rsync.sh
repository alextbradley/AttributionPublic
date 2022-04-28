#!/bin/bash

JOBNO=666

# move to run directory
cd ../run

# make target directory
W_HOMEDIR=$W_HOMEROOT/ATTR_${JOBNO}/run
ssh $W_HOMEHOST "mkdir -p $W_HOMEDIR"

rsync -avzL *.nc $W_HOMEHOST:$W_HOMEDIR

