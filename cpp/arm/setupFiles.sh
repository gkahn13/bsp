#!/bin/sh

# Arguments:
# 1: arm folder
#        belief, state
# 2: timesteps

ARM_TYPE=$1
TIMESTEPS=$2


BSP_DIR=$(sed s:/"bsp"/.*:/"bsp": <<< $(pwd))
ARM_DIR=" $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DIR=$ARM_DIR/$ARM_TYPE

MPC_DIR=$DIR/mpc
if [ $ARM_TYPE = "belief" ]; then
MPC_FILE_NAME="${ARM_TYPE}PenaltyMPC"
else
MPC_FILE_NAME="${ARM_TYPE}MPC"
fi

CPP_TIMESTEP_DEF="#define TIMESTEPS"

cd $ARM_DIR

# deal with FORCES files
MPC_C_FILE="${MPC_DIR}/${MPC_FILE_NAME}${TIMESTEPS}.c"
MPC_H_FILE="${MPC_DIR}/${MPC_FILE_NAME}${TIMESTEPS}.h"

if [ ! -f $MPC_C_FILE ] || [ ! -f $MPC_H_FILE ]; then
    echo "mpc files do not exist. generating with matlab..."
    cd $DIR
    if [ $ARM_TYPE = "belief" ]; then MATLAB_FILE_NAME="${ARM_TYPE}_penalty"
    else MATLAB_FILE_NAME="{$ARM_TYPE}"; fi
    
    matlab -nodisplay -nosplash -nodesktop -r "arm_${MATLAB_FILE_NAME}_mpc_gen(${TIMESTEPS}); exit"
    cd $ARM_DIR
fi

echo "Copying ${MPC_FILE_NAME}${TIMESTEPS}.c/h to ${MPC_FILE_NAME}.c/h"
cp ${MPC_C_FILE} ${DIR}/${MPC_FILE_NAME}".c"
cp ${MPC_H_FILE} ${DIR}/${MPC_FILE_NAME}".h"

# modify arm.h
echo "replacing TIMESTEPS definition with new TIMESTEPS for arm.h in ${ARM_DIR}"
H_WRITE="${ARM_DIR}/arm.h"
sed -i "s/^${CPP_TIMESTEP_DEF}.*/${CPP_TIMESTEP_DEF} ${TIMESTEPS}/" $H_WRITE
