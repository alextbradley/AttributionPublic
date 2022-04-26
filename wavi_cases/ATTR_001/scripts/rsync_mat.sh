#!/bin/bash

JOBNO=024

# move to run directory
cd ../run

# copy the julia script from utilities
cp ../../../utilities/zip_mat.jl .

#execute the zipping script
export SINGULARITYENV_JULIA_DEPOT_PATH="/opt/julia"
singularity exec -B $JDEPOT:/opt/julia $IMGNAME julia zip_jld2.jl

# make target directory
W_HOMEDIR=$W_HOMEROOT/MISMIP_${JOBNO}/run
ssh $W_HOMEHOST "mkdir -p $W_HOMEDIR"

rsync -avzL *.nc $W_HOMEHOST:$W_HOMEDIR

