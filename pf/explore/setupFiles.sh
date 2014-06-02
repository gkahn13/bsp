#!/bin/sh

# Arguments:
# 1: timesteps
# 2: agents

TIMESTEPS=$1
AGENTS=$2

MPC_FILE_NAME="exploreMPC"

MPC_C_FILE="mpc/${MPC_FILE_NAME}${TIMESTEPS}_${AGENTS}.c"
MPC_H_FILE="mpc/${MPC_FILE_NAME}${TIMESTEPS}_${AGENTS}.h"

echo ${MPC_C_FILE}
echo ${MPC_H_FILE}

if [ ! -f $MPC_C_FILE ] || [ ! -f $MPC_H_FILE ]; then
    echo "mpc files do not exist. generating with matlab..."
    
    matlab -nodisplay -nosplash -nodesktop -r "$explore_mpc_gen(${TIMESTEPS},${AGENTS}); exit"
fi

echo "Copying ${MPC_FILE_NAME}${TIMESTEPS}_${AGENTS}.c/h to ${MPC_FILE_NAME}.c/h"
cp ${MPC_C_FILE} ${MPC_FILE_NAME}".c"
cp ${MPC_H_FILE} ${MPC_FILE_NAME}".h"

# modify point-pf.h
CPP_TIMESTEPS_DEF="#define TIMESTEPS"
CPP_AGENTS_DEF="#define AGENTS"

echo "replacing TIMESTEPS definition with new TIMESTEPS for explore.cpp"
sed -i "s/^${CPP_TIMESTEPS_DEF}.*/${CPP_TIMESTEPS_DEF} ${TIMESTEPS}/" explore.cpp

echo "replacing AGENTS definition with new AGENTS for explore.cpp"
sed -i "s/^${CPP_AGENTS_DEF}.*/${CPP_AGENTS_DEF} ${AGENTS}/" explore.cpp
