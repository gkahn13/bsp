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
MPC_FILE_NAME="${PARAM_TYPE}PenaltyMPC"
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
    else MATLAB_FILE_NAME="{$PARAM_TYPE}_penalty"; fi
    
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
rm -f *.bk


# make casadi files
# move casadi files over if different
if [ $PARAM_TYPE = "state" ]; then
    CASADI_PARAM_DIR="${BSP_DIR}/casadi/parameter"
    echo "Making casadi files"
    make -C ${CASADI_PARAM_DIR} all T=${TIMESTEPS}

    CASADI_FILES=${CASADI_PARAM_DIR}/*.c
    for CASADI_FILE in $CASADI_FILES
    do
    CASADI_FILE_NAME=$(basename $CASADI_FILE)
    PARAM_CASADI_FILE=${PARAM_DIR}/state/${CASADI_FILE_NAME}

    if [ ! -f $PARAM_CASADI_FILE ]; then
        echo "casadi file ${CASADI_FILE_NAME} does not exit, copying over"
        cp $CASADI_FILE $PARAM_CASADI_FILE
    fi

    if [ $(diff $CASADI_FILE $PARAM_CASADI_FILE | wc -w) -gt 0 ]; then
        echo "casadi file ${CASADI_FILE_NAME} differs, copying new one over"
        cp $CASADI_FILE $PARAM_CASADI_FILE
    fi
   done
fi

# make casadi files
# move casadi files over if different
if [ $PARAM_TYPE = "controls" ]; then
    CASADI_PARAM_DIR="${BSP_DIR}/casadi/parameter/controls"
    echo "Making casadi files"
    make -C ${CASADI_PARAM_DIR} all T=${TIMESTEPS}

    CASADI_FILES=${CASADI_PARAM_DIR}/*.c
    for CASADI_FILE in $CASADI_FILES
    do
    CASADI_FILE_NAME=$(basename $CASADI_FILE)
    PARAM_CASADI_FILE=${PARAM_DIR}/controls/${CASADI_FILE_NAME}

    if [ ! -f $PARAM_CASADI_FILE ]; then
        echo "casadi file ${CASADI_FILE_NAME} does not exit, copying over"
        cp $CASADI_FILE $PARAM_CASADI_FILE
    fi

    if [ $(diff $CASADI_FILE $PARAM_CASADI_FILE | wc -w) -gt 0 ]; then
        echo "casadi file ${CASADI_FILE_NAME} differs, copying new one over"
        cp $CASADI_FILE $PARAM_CASADI_FILE
    fi
   done
fi