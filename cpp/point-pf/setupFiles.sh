#!/bin/sh

# Arguments:
# 1: point-pf folder
#        platt
# 2: timesteps
# 3: particles

POINT_PF_TYPE=$1
TIMESTEPS=$2
PARTICLES=$3

MPC_FILE_NAME="${POINT_PF_TYPE}MPC"

MPC_C_FILE="${POINT_PF_TYPE}/mpc/${MPC_FILE_NAME}${TIMESTEPS}_${PARTICLES}.c"
MPC_H_FILE="${POINT_PF_TYPE}/mpc/${MPC_FILE_NAME}${TIMESTEPS}_${PARTICLES}.h"

if [ ! -f $MPC_C_FILE ] || [ ! -f $MPC_H_FILE ]; then
    echo "mpc files do not exist. generating with matlab..."
    cd ${POINT_PF_TYPE}
    
    matlab -nodisplay -nosplash -nodesktop -r "${POINT_PF_TYPE}_mpc_gen(${TIMESTEPS},${PARTICLES}); exit"
    cd ..
fi

echo "Copying ${MPC_FILE_NAME}${TIMESTEPS}_${PARTICLES}.c/h to ${MPC_FILE_NAME}.c/h"
cp ${MPC_C_FILE} ${POINT_PF_TYPE}/${MPC_FILE_NAME}".c"
cp ${MPC_H_FILE} ${POINT_PF_TYPE}/${MPC_FILE_NAME}".h"

# modify point-pf.h
CPP_TIMESTEPS_DEF="#define TIMESTEPS"
CPP_PARTICLES_DEF="#define PARTICLES"

echo "replacing TIMESTEPS definition with new TIMESTEPS for point-pf.h"
sed -i "s/^${CPP_TIMESTEPS_DEF}.*/${CPP_TIMESTEPS_DEF} ${TIMESTEPS}/" point-pf.h

echo "replacing PARTICLES definition with new PARTICLES for point-pf.h"
sed -i "s/^${CPP_PARTICLES_DEF}.*/${CPP_PARTICLES_DEF} ${PARTICLES}/" point-pf.h
