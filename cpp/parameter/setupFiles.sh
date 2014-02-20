#!/bin/sh

# Arguments:
# 1: parameter folder
#        belief, state
# 2: timesteps

PARAM_TYPE=$1
TIMESTEPS=$2


BSP_DIR=$(sed s:/"bsp"/.*:/"bsp": <<< $(pwd))
PARAM_DIR=" $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DIR=$PARAM_DIR/$PARAM_TYPE

MPC_DIR=$DIR/mpc
if [ $PARAM_TYPE = "belief" ]; then
MPC_FILE_NAME="${PARAM_TYPE}PenaltyMPC"
else
MPC_FILE_NAME="${PARAM_TYPE}MPC"
fi

CPP_TIMESTEP_DEF="#define TIMESTEPS"

cd $PARAM_DIR

# deal with FORCES files
MPC_C_FILE="${MPC_DIR}/${MPC_FILE_NAME}${TIMESTEPS}.c"
MPC_H_FILE="${MPC_DIR}/${MPC_FILE_NAME}${TIMESTEPS}.h"

if [ ! -f $MPC_C_FILE ] || [ ! -f $MPC_H_FILE ]; then
    echo "mpc files do not exist. generating with matlab..."
    cd $DIR
    if [ $PARAM_TYPE = "belief" ]; then MATLAB_FILE_NAME="${PARAM_TYPE}_penalty"
    else MATLAB_FILE_NAME="{$PARAM_TYPE}"; fi
    
    matlab -nodisplay -nosplash -nodesktop -r "parameter_${MATLAB_FILE_NAME}_mpc_gen(${TIMESTEPS}); exit"
    cd $PARAM_DIR
fi

echo "Copying ${MPC_FILE_NAME}${TIMESTEPS}.c/h to ${MPC_FILE_NAME}.c/h"
cp ${MPC_C_FILE} ${DIR}/${MPC_FILE_NAME}".c"
cp ${MPC_H_FILE} ${DIR}/${MPC_FILE_NAME}".h"

# modify parameter.h
echo "replacing TIMESTEPS definition with new TIMESTEPS for parameter.h in ${PARAM_DIR}"
H_WRITE="${PARAM_DIR}/parameter.h"
sed -i .bk "s/^${CPP_TIMESTEP_DEF}.*/${CPP_TIMESTEP_DEF} ${TIMESTEPS}/" $H_WRITE
rm *.bk