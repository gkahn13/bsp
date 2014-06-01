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
MPC_FILE_NAME="${ARM_TYPE}PenaltyMPC"
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
sed -i .bk "s/^${CPP_TIMESTEP_DEF}.*/${CPP_TIMESTEP_DEF} ${TIMESTEPS}/" $H_WRITE
rm -f *.bk

# make casadi files
# move casadi files over if different
if [ $ARM_TYPE = "control" ] || [ $ARM_TYPE = "state" ]; then
    CASADI_ARM_DIR="${ARM_DIR}/casadi/"
    echo "Making casadi files"
    make -C ${CASADI_ARM_DIR} all T=${TIMESTEPS}

    CASADI_FILES=${CASADI_ARM_DIR}/*.c
    for CASADI_FILE in $CASADI_FILES
    do
    CASADI_FILE_NAME=$(basename $CASADI_FILE)
    ARM_CASADI_FILE=${ARM_DIR}/${ARM_TYPE}/${CASADI_FILE_NAME}

    if [ ! -f $ARM_CASADI_FILE ]; then
        echo "casadi file ${CASADI_FILE_NAME} does not exit, copying over"
        cp $CASADI_FILE $ARM_CASADI_FILE
    fi

    if [ $(diff $CASADI_FILE $ARM_CASADI_FILE | wc -w) -gt 0 ]; then
        echo "casadi file ${CASADI_FILE_NAME} differs, copying new one over"
        cp $CASADI_FILE $ARM_CASADI_FILE
    fi
   done
fi
