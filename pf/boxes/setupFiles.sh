#!/bin/sh

# Arguments:
# 1: timesteps

TIMESTEPS=$1

MPC_FILE_NAME="boxesMPC"

MPC_C_FILE="mpc/${MPC_FILE_NAME}${TIMESTEPS}_${AGENTS}.c"
MPC_H_FILE="mpc/${MPC_FILE_NAME}${TIMESTEPS}_${AGENTS}.h"

echo ${MPC_C_FILE}
echo ${MPC_H_FILE}

if [ ! -f $MPC_C_FILE ] || [ ! -f $MPC_H_FILE ]; then
    echo "mpc files do not exist. generating with matlab..."
    
    matlab -nodisplay -nosplash -nodesktop -r "$boxes_mpc_gen(${TIMESTEPS}); exit"
fi

echo "Copying ${MPC_FILE_NAME}${TIMESTEPS}_${AGENTS}.c/h to ${MPC_FILE_NAME}.c/h"
cp ${MPC_C_FILE} ${MPC_FILE_NAME}".c"
cp ${MPC_H_FILE} ${MPC_FILE_NAME}".h"

CPP_TIMESTEPS_DEF="#define TIMESTEPS"

echo "replacing TIMESTEPS definition with new TIMESTEPS for explore.cpp"
sed -i "s/^${CPP_TIMESTEPS_DEF}.*/${CPP_TIMESTEPS_DEF} ${TIMESTEPS}/" explore.cpp
