#!/bin/sh

# Arguments:
# 1: point folder
#        belief, state, lp, control
# 2: timesteps

POINT_TYPE=$1
TIMESTEPS=$2


BSP_DIR=$(sed s:/"bsp"/.*:/"bsp": <<< $(pwd))
POINT_DIR=" $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DIR=$POINT_DIR/$POINT_TYPE

MPC_DIR=$DIR/mpc
MPC_FILE_NAMES[0]="${POINT_TYPE}MPC"
if [ $POINT_TYPE = "belief" ]; then MPC_FILE_NAMES+=("${POINT_TYPE}PenaltyMPC"); fi

CPP_TIMESTEP_DEF="#define TIMESTEPS"

cd $POINT_DIR

for MPC_FILE_NAME in ${MPC_FILE_NAMES[@]}
do
    MPC_C_FILE="${MPC_DIR}/${MPC_FILE_NAME}${TIMESTEPS}.c"
    MPC_H_FILE="${MPC_DIR}/${MPC_FILE_NAME}${TIMESTEPS}.h"

    if [ ! -f $MPC_C_FILE ] || [ ! -f $MPC_H_FILE ]; then
	echo "mpc files do not exist. generating with matlab..."
	cd $DIR
	matlab -nodisplay -nosplash -nodesktop -r "point_${POINT_TYPE}_mpc_gen(${TIMESTEPS}); exit"
	cd $POINT_DIR
    fi

    echo "Copying ${MPC_FILE_NAME}${TIMESTEPS}.c/h to ${MPC_FILE_NAME}.c/h"
    cp ${MPC_C_FILE} ${DIR}/${MPC_FILE_NAME}".c"
    cp ${MPC_H_FILE} ${DIR}/${MPC_FILE_NAME}".h"
done

echo "replacing TIMESTEPS definition with new TIMESTEPS for point.h in ${POINT_DIR}"
H_WRITE="${POINT_DIR}/point.h"
sed -i "s/^${CPP_TIMESTEP_DEF}.*/${CPP_TIMESTEP_DEF} ${TIMESTEPS}/" $H_WRITE

if [ $POINT_TYPE = "state" ] || [ $POINT_TYPE = "lp" ] || [ $POINT_TYPE = "control" ]; then
    # since lp uses same sym-eval and mask as state
    if [ $POINT_TYPE = "state" ] || [ $POINT_TYPE = "lp" ]; then
	SYM_MASK_TYPE="state"
    else
	SYM_MASK_TYPE="control"
    fi

    DSTAR_PT_SYM_DIR="${BSP_DIR}/dstar/point/sym"
    SYM_FILE_TIMESTEP="${DSTAR_PT_SYM_DIR}/${SYM_MASK_TYPE}-symeval-${TIMESTEPS}.c"
    SYM_FILE="${POINT_DIR}/sym/${SYM_MASK_TYPE}-symeval.c"

    if [ ! -f $SYM_FILE_TIMESTEP ]; then
	echo "sym file ${SYM_FILE_TIMESTEP} does not exist"
	exit -1
    fi

    if [ $(diff $SYM_FILE_TIMESTEP $SYM_FILE | wc -w) -gt 0 ]; then
	echo "sym_eval file differs, copying over correct file..."
	cp $SYM_FILE_TIMESTEP $SYM_FILE
    else
	echo "sym_eval file is the same. not copying over"
    fi

fi
