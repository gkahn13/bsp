#!/bin/sh

# Arguments:
# 1: point-explore folder
#        explore
# 2: timesteps
# 3: particles
# 4: agents

POINT_EXPLORE_TYPE=$1
TIMESTEPS=$2
PARTICLES=$3
AGENTS=$4

MPC_FILE_NAME="${POINT_EXPLORE_TYPE}MPC"

MPC_C_FILE="${POINT_EXPLORE_TYPE}/mpc/${MPC_FILE_NAME}${TIMESTEPS}_${AGENTS}.c"
MPC_H_FILE="${POINT_EXPLORE_TYPE}/mpc/${MPC_FILE_NAME}${TIMESTEPS}_${AGENTS}.h"

echo ${MPC_C_FILE}
echo ${MPC_H_FILE}

if [ ! -f $MPC_C_FILE ] || [ ! -f $MPC_H_FILE ]; then
    echo "mpc files do not exist. generating with matlab..."
    cd ${POINT_EXPLORE_TYPE}
    
    matlab -nodisplay -nosplash -nodesktop -r "${POINT_EXPLORE_TYPE}_mpc_gen(${TIMESTEPS},${AGENTS}); exit"
    cd ..
fi

echo "Copying ${MPC_FILE_NAME}${TIMESTEPS}_${AGENTS}.c/h to ${MPC_FILE_NAME}.c/h"
cp ${MPC_C_FILE} ${POINT_EXPLORE_TYPE}/${MPC_FILE_NAME}".c"
cp ${MPC_H_FILE} ${POINT_EXPLORE_TYPE}/${MPC_FILE_NAME}".h"

# modify point-pf.h
CPP_TIMESTEPS_DEF="#define TIMESTEPS"
CPP_PARTICLES_DEF="#define PARTICLES"
CPP_AGENTS_DEF="#define AGENTS"

echo "replacing TIMESTEPS definition with new TIMESTEPS for point-explore.h"
sed -i "s/^${CPP_TIMESTEPS_DEF}.*/${CPP_TIMESTEPS_DEF} ${TIMESTEPS}/" point-explore-system.h

echo "replacing PARTICLES definition with new PARTICLES for point-explore.h"
sed -i "s/^${CPP_PARTICLES_DEF}.*/${CPP_PARTICLES_DEF} ${PARTICLES}/" point-explore-system.h

echo "replacing AGENTS definition with new AGENTS for point-explore.h"
sed -i "s/^${CPP_AGENTS_DEF}.*/${CPP_AGENTS_DEF} ${AGENTS}/" point-explore-system.h