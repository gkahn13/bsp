#!/bin/sh

# Arguments:
# 1: slam folder
#        belief, state, lp, control
# 2: timesteps
# 3: landmarks

SLAM_TYPE=$1
TIMESTEPS=$2
NUM_LANDMARKS=$3


BSP_DIR=$(sed s:/"bsp"/.*:/"bsp": <<< $(pwd))
SLAM_DIR=" $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DIR=$SLAM_DIR/$SLAM_TYPE

MPC_DIR=$DIR/mpc
if [ $SLAM_TYPE = "belief" ]; then
MPC_FILE_NAME="${SLAM_TYPE}PenaltyMPC"
else
MPC_FILE_NAME="${SLAM_TYPE}MPC"
fi

CPP_TIMESTEP_DEF="#define TIMESTEPS"
CPP_LANDMARK_DEF="#define NUM_LANDMARKS"

cd $SLAM_DIR

# deal with FORCES files
MPC_C_FILE="${MPC_DIR}/${MPC_FILE_NAME}${TIMESTEPS}.c"
MPC_H_FILE="${MPC_DIR}/${MPC_FILE_NAME}${TIMESTEPS}.h"

# must always make for belief because landmarks is part of state
if [ ! -f $MPC_C_FILE ] || [ ! -f $MPC_H_FILE ] || [ $SLAM_TYPE = "belief" ]; then
    echo "mpc files do not exist. generating with matlab..."
    cd $DIR
    if [ $SLAM_TYPE = "belief" ]; then MATLAB_FILE_NAME="${SLAM_TYPE}_penalty"
    else MATLAB_FILE_NAME="${SLAM_TYPE}"; fi
    
    if [ $SLAM_TYPE = "belief" ]; then matlab -nodisplay -nosplash -nodesktop -r "slam_${MATLAB_FILE_NAME}_mpc_gen(${TIMESTEPS},${NUM_LANDMARKS}); exit"
    else matlab -nodisplay -nosplash -nodesktop -r "slam_${MATLAB_FILE_NAME}_mpc_gen(${TIMESTEPS}); exit"; fi
    cd $SLAM_DIR
fi


if [ ! -f ${DIR}/${MPC_FILE_NAME}".c" ] || [ $(diff ${MPC_C_FILE} ${DIR}/${MPC_FILE_NAME}".c" | wc -w) -gt 0 ]; then
    echo "Copying ${MPC_FILE_NAME}${TIMESTEPS}.c to ${MPC_FILE_NAME}.c"
    cp ${MPC_C_FILE} ${DIR}/${MPC_FILE_NAME}".c"
fi
if [ ! -f ${DIR}/${MPC_FILE_NAME}".h" ] || [ $(diff ${MPC_H_FILE} ${DIR}/${MPC_FILE_NAME}".h" | wc -w) -gt 0 ]; then
    echo "Copying ${MPC_FILE_NAME}${TIMESTEPS}.h to ${MPC_FILE_NAME}.h"
    cp ${MPC_H_FILE} ${DIR}/${MPC_FILE_NAME}".h"
fi

# modify slam.h
echo "replacing TIMESTEPS definition with new TIMESTEPS for slam.h in ${SLAM_DIR}"
H_WRITE="${SLAM_DIR}/slam.h"
sed -i "s/^${CPP_TIMESTEP_DEF}.*/${CPP_TIMESTEP_DEF} ${TIMESTEPS}/" $H_WRITE
echo "replacing NUM_LANDMARKS definition with new NUM_LANDMARKS for slam.h in ${SLAM_DIR}"
sed -i "s/^${CPP_LANDMARK_DEF}.*/${CPP_LANDMARK_DEF} ${NUM_LANDMARKS}/" $H_WRITE

# make casadi files
# move casadi files over if different
if [ $SLAM_TYPE = "state" ] || [ $SLAM_TYPE = "control" ]; then
    CASADI_SLAM_DIR="${BSP_DIR}/casadi/slam"
    echo "Making casadi files"
    make -C ${CASADI_SLAM_DIR} "slam-${SLAM_TYPE}" T=${TIMESTEPS} NUM_LANDMARKS=${NUM_LANDMARKS}

    CASADI_FILES=${CASADI_SLAM_DIR}/*.c
    for CASADI_FILE in $CASADI_FILES
    do
	CASADI_FILE_NAME=$(basename $CASADI_FILE)
	SLAM_CASADI_FILE=${SLAM_DIR}/${SLAM_TYPE}/${CASADI_FILE_NAME}

	if [ ! -f $SLAM_CASADI_FILE ]; then
	    echo "casadi file ${CASADI_FILE_NAME} does not exit, copying over"
	    cp $CASADI_FILE $SLAM_CASADI_FILE
	fi

	if [ $(diff $CASADI_FILE $SLAM_CASADI_FILE | wc -w) -gt 0 ]; then
	    echo "casadi file ${CASADI_FILE_NAME} differs, copying new one over"
	    cp $CASADI_FILE $SLAM_CASADI_FILE
	fi
   done
fi