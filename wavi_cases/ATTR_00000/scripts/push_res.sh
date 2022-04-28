#! /bin/bash
#
# Make an outfile of the current run directory and push the results. 
# 

JOBNO=00000

cd ../run
cp $W_ROOT/utilities/zip_mat.jl .
export SINGULARITYENV_JULIA_DEPOT_PATH="/opt/julia"
singularity exec -B ${JDEPOT}:/opt/julia,$(pwd) ${IMGNAME} julia zip_mat.jl

# make target directory
W_HOMEDIR=$W_HOMEROOT/ATTR_${JOBNO}/run
ssh $W_HOMEHOST "mkdir -p $W_HOMEDIR"

rsync -avzL *.nc $W_HOMEHOST:$W_HOMEDIR


