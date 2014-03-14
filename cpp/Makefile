CXX = g++
CC = gcc
INTELCC = icc

OS := mac

LINKER_FLAGS.mac = -lm
LINKER_FLAGS.linux = -lrt
LINKER_FLAGS = $(LINKER_FLAGS.$(OS))

C_FLAGS = -std=gnu99
CPP_FLAGS = -std=c++0x -fpermissive -I.
PYTHON_FLAGS = -I/usr/include/python2.7 -lpython2.7
BOOST_FLAGS = -lboost_python -lboost_filesystem -lboost_system
ADOLC_FLAGS = -I/home/gkahn/adolc_base/include -L/home/gkahn/adolc_base/lib64 -ladolc -lColPack
#ADOLC_FLAGS = -I/usr/include/adolc -ladolc -lColPack
EIGEN_FLAGS = -I/usr/include/eigen3
CASADI_FLAGS = -I$(CASADIPATH) -L$(CASADIPATH)/build/lib
CASADI_LIBS = -lcasadi -ldl

# default to debug build
BUILD := debug

BFLAGS.debug = -g -Wall
BFLAGS.release = -O3
BFLAGS := $(BFLAGS.$(BUILD))

OBJ_DIR = build/$(BUILD)
BIN_DIR = bin/$(BUILD)

# default to 15 T
T = 15

###### DIRECTORY CREATION ###########
dir:
	mkdir -p build/debug
	mkdir -p build/release
	mkdir -p bin/debug
	mkdir -p bin/release

###### UTIL ############

$(OBJ_DIR)/logging.o : util/logging.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $^
	
###### ARM ############

# share arm functions
ARM_SETUP_FILE = arm/setupFiles.sh

ARM_HEADERS = arm/arm.h arm/matrix.h
ARM_BELIEF_DIR = arm/belief

# make arm-belief-joint
ARM_BELIEF_JOINT_FILES = arm-belief-joint-goal arm-belief-joint logging
ARM_BELIEF_JOINT_OBJS = $(ARM_BELIEF_JOINT_FILES:%=$(OBJ_DIR)/%.o)

arm-belief-joint: $(ARM_BELIEF_JOINT_OBJS)
	$(CXX) $(BFLAGS) $(ARM_BELIEF_JOINT_OBJS) -o $(BIN_DIR)/arm-belief-joint $(BOOST_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/arm-belief-joint-goal.o : $(ARM_BELIEF_DIR)/arm-belief-joint-goal-MPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/arm-belief-joint.o : $(ARM_BELIEF_DIR)/arm-belief-joint-goal.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $<

# make arm-belief-pos
ARM_BELIEF_POS_FILES = arm-belief-pos-goal arm-belief-pos logging
ARM_BELIEF_POS_OBJS = $(ARM_BELIEF_POS_FILES:%=$(OBJ_DIR)/%.o)

arm-belief-pos: $(ARM_BELIEF_POS_OBJS)
	$(CXX) $(BFLAGS) $(ARM_BELIEF_POS_OBJS) -o $(BIN_DIR)/arm-belief-pos $(BOOST_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/arm-belief-pos-goal.o : $(ARM_BELIEF_DIR)/arm-belief-pos-goal-MPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/arm-belief-pos.o : $(ARM_BELIEF_DIR)/arm-belief-pos-goal.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $<

# make arm-ilqg
ARM_ILQG_DIR = arm/ilqg
ARM_ILQG_FILES = arm-ilqg logging
ARM_ILQG_OBJS = $(ARM_ILQG_FILES:%=$(OBJ_DIR)/%.o)

arm-ilqg: | $(ARM_ILQG_OBJS)
	$(CXX) $(BFLAGS) $(ARM_ILQG_OBJS) -o $(BIN_DIR)/arm-ilqg $(BOOST_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/arm-ilqg.o : $(ARM_ILQG_DIR)/arm-ilqg-pos-goal.cpp $(ARM_ILQG_DIR)/ilqg.h
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $<
	
	
# make arm-state
ARM_STATE_DIR = arm/state
ARM_STATE_FILES = arm-state-cost arm-state-grad arm-statePenaltyMPC arm-state logging
ARM_STATE_OBJS = $(ARM_STATE_FILES:%=$(OBJ_DIR)/%.o)

#touch-arm-state:
#	touch $(ARM_STATE_DIR)/arm-state-pos-goal.cpp
arm-state-setup:
	bash $(ARM_SETUP_FILE) state $(T)

arm-state: | arm-state-setup $(ARM_STATE_OBJS)
	$(CXX) $(BFLAGS) $(ARM_STATE_OBJS) -o $(BIN_DIR)/arm-state $(BOOST_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/arm-state-cost.o : $(ARM_STATE_DIR)/arm-state-cost.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/arm-state-grad.o : $(ARM_STATE_DIR)/arm-state-grad.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/arm-statePenaltyMPC.o : $(ARM_STATE_DIR)/statePenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/arm-state.o : $(ARM_STATE_DIR)/arm-state-pos-goal.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $<

#make arm-control 
ARM_CONTROL_DIR = arm/control
ARM_CONTROL_FILES = arm-control-cost arm-control-grad arm-statePenaltyMPC arm-control logging
ARM_CONTROL_OBJS = $(ARM_CONTROL_FILES:%=$(OBJ_DIR)/%.o)

#touch-arm-CONTROL:
#	touch $(ARM_CONTROL_DIR)/arm-CONTROL-pos-goal.cpp
arm-control-setup:
	bash $(ARM_SETUP_FILE) control $(T)

arm-control: | arm-control-setup $(ARM_CONTROL_OBJS)
	$(CXX) $(BFLAGS) $(ARM_CONTROL_OBJS) -o $(BIN_DIR)/arm-control $(BOOST_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/arm-control-cost.o : $(ARM_CONTROL_DIR)/arm-control-cost.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/arm-control-grad.o : $(ARM_CONTROL_DIR)/arm-control-grad.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/arm-control-pos-goal.o : $(ARM_CONTROL_DIR)/controlPenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/arm-control.o : $(ARM_CONTROL_DIR)/arm-control-pos-goal.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $<
	


###### PARAMETER ###########

# share parameter functions
PARAM_SETUP_FILE = parameter/setupFiles.sh

PARAM_HEADERS = parameter/parameter.h

# make parameter-belief
PARAM_BELIEF_DIR = parameter/belief
PARAM_BELIEF_FILES = parameterBeliefPenaltyMPC parameter-belief logging
PARAM_BELIEF_OBJS = $(PARAM_BELIEF_FILES:%=$(OBJ_DIR)/%.o)

parameter-belief-setup:
	touch $(PARAM_BELIEF_DIR)/parameter-belief.cpp
	bash $(PARAM_SETUP_FILE) belief $(T)

parameter-belief: | parameter-belief-setup $(PARAM_BELIEF_OBJS)
	$(CXX) $(BFLAGS) $(PARAM_BELIEF_OBJS) -o $(BIN_DIR)/parameter-belief $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/parameterBeliefPenaltyMPC.o : $(PARAM_BELIEF_DIR)/beliefPenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/parameter-belief.o : $(PARAM_BELIEF_DIR)/parameter-belief.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	

# make parameter-controls
PARAM_CONTROLS_DIR = parameter/controls
PARAM_CONTROLS_FILES = parameter-controls logging parameter-controls-cost parameter-controls-grad controlsPenaltyMPC
PARAM_CONTROLS_OBJS = $(PARAM_CONTROLS_FILES:%=$(OBJ_DIR)/%.o)

parameter-controls-setup:
	bash $(PARAM_SETUP_FILE) controls $(T)

parameter-controls: | parameter-controls-setup $(PARAM_CONTROLS_OBJS)
	$(CXX) $(BFLAGS) $(PARAM_CONTROLS_OBJS) -o $(BIN_DIR)/parameter-controls $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/parameter-controls-cost.o : $(PARAM_CONTROLS_DIR)/parameter-controls-cost.c
	$(CC) $(C_FLAGS) -c -o $@ $^

$(OBJ_DIR)/parameter-controls-grad.o : $(PARAM_CONTROLS_DIR)/parameter-controls-grad.c
	$(CC) $(C_FLAGS) -c -o $@ $^

$(OBJ_DIR)/controlsPenaltyMPC.o : $(PARAM_CONTROLS_DIR)/controlsPenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/parameter-controls.o : $(PARAM_CONTROLS_DIR)/parameter-controls.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<


# make parameter-state
PARAM_STATE_DIR = parameter/state
PARAM_STATE_FILES = parameter-state logging parameter-state-cost parameter-state-grad statePenaltyMPC
PARAM_STATE_OBJS = $(PARAM_STATE_FILES:%=$(OBJ_DIR)/%.o)

parameter-state-setup:
	bash $(PARAM_SETUP_FILE) state $(T)

parameter-state: | parameter-state-setup $(PARAM_STATE_OBJS)
	$(CXX) $(BFLAGS) $(PARAM_STATE_OBJS) -o $(BIN_DIR)/parameter-state $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/parameter-state-cost.o : $(PARAM_STATE_DIR)/parameter-state-cost.c
	$(CC) $(C_FLAGS) -c -o $@ $^

$(OBJ_DIR)/parameter-state-grad.o : $(PARAM_STATE_DIR)/parameter-state-grad.c
	$(CC) $(C_FLAGS) -c -o $@ $^

$(OBJ_DIR)/statePenaltyMPC.o : $(PARAM_STATE_DIR)/statePenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/parameter-state.o : $(PARAM_STATE_DIR)/parameter-state.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<


# make parameter-ilqg
PARAM_ILQG_DIR = parameter/ilqg
PARAM_ILQG_FILES = parameter-ilqg
PARAM_ILQG_OBJS = $(PARAM_ILQG_FILES:%=$(OBJ_DIR)/%.o)

parameter-ilqg: $(PARAM_ILQG_OBJS)
	$(CXX) $(BFLAGS) $(PARAM_ILQG_OBJS) -o $(BIN_DIR)/parameter-ilqg $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/parameter-ilqg.o : $(PARAM_ILQG_DIR)/parameter-ilqg.cpp $(PARAM_ILQG_DIR)/ilqg.h
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) -c -o $@ $<

# make parameter-random
PARAM_RANDOM_DIR = parameter/random
PARAM_RANDOM_FILES = parameter-random logging
PARAM_RANDOM_OBJS = $(PARAM_RANDOM_FILES:%=$(OBJ_DIR)/%.o)

touch-parameter-random:
	touch $(PARAM_RANDOM_DIR)/parameter-random.cpp

parameter-random: | touch-parameter-random $(PARAM_RANDOM_OBJS)
	$(CXX) $(BFLAGS) $(PARAM_RANDOM_OBJS) -o $(BIN_DIR)/parameter-random $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/parameter-random.o : $(PARAM_RANDOM_DIR)/parameter-random.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) -c -o $@ $<




###### POINT ##########

# setupFiles bash
POINT_SETUP_FILE = point/setupFiles.sh

# share point functions
PT_HEADERS = point/point.h
PT_ADOLC_HEADERS = point/point.h point/point-adouble.h

# symeval	
PT_SYM_DIR = point/sym

$(OBJ_DIR)/state_symeval.o : $(PT_SYM_DIR)/state-symeval.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/control_symeval.o : $(PT_SYM_DIR)/control-symeval.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
# make point-belief
PT_BELIEF_DIR = point/belief
PT_BELIEF_FILES = pointBeliefMPC pointBeliefPenaltyMPC point-belief logging
PT_BELIEF_OBJS = $(PT_BELIEF_FILES:%=$(OBJ_DIR)/%.o)

point-belief-setup:
	bash $(POINT_SETUP_FILE) belief $(T)

point-belief: | point-belief-setup $(PT_BELIEF_OBJS)
	$(CXX) $(BFLAGS) $(PT_BELIEF_OBJS) -o $(BIN_DIR)/point-belief $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/pointBeliefMPC.o : $(PT_BELIEF_DIR)/beliefMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/pointBeliefPenaltyMPC.o : $(PT_BELIEF_DIR)/beliefPenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-belief.o : $(PT_BELIEF_DIR)/point-belief.cpp $(POINT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
	
# make point-state
PT_STATE_DIR = point/state
PT_STATE_FILES = pointStateMPC state_symeval point-state logging
PT_STATE_OBJS = $(PT_STATE_FILES:%=$(OBJ_DIR)/%.o)

point-state-setup:
	bash $(POINT_SETUP_FILE) state $(T)

point-state: | point-state-setup $(PT_STATE_OBJS)
	$(CXX) $(BFLAGS) $(PT_STATE_OBJS) -o $(BIN_DIR)/point-state $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/pointStateMPC.o : $(PT_STATE_DIR)/stateMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-state.o : $(PT_STATE_DIR)/point-state.cpp $(PT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
# make point-state-finitediff
PT_STATE_FD_DIR = $(PT_STATE_DIR)
PT_STATE_FD_FILES = pointStateMPC state_symeval point-state-finitediff logging
PT_STATE_FD_OBJS = $(PT_STATE_FD_FILES:%=$(OBJ_DIR)/%.o)

point-state-finitediff: | point-state-setup $(PT_STATE_FD_OBJS) 
	$(CXX) $(BFLAGS) $(PT_STATE_FD_OBJS) -o $(BIN_DIR)/point-state-finitediff $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/point-state-finitediff.o : $(PT_STATE_FD_DIR)/point-state-finitediff.cpp $(POINT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
# make point-state-adolc
PT_STATE_ADOLC_DIR = point/state
PT_STATE_ADOLC_FILES = pointStateMPC point-state-adolc logging
PT_STATE_ADOLC_OBJS = $(PT_STATE_ADOLC_FILES:%=$(OBJ_DIR)/%.o)

point-state-adolc: | $(PT_STATE_ADOLC_OBJS)
	$(CXX) $(BFLAGS) $(PT_STATE_ADOLC_OBJS) -o $(BIN_DIR)/point-state-adolc $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS) $(ADOLC_FLAGS)
	
$(OBJ_DIR)/point-state-adolc.o : $(PT_STATE_ADOLC_DIR)/point-state-adolc.cpp $(PT_ADOLC_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) $(ADOLC_FLAGS) -c -o $@ $<

# make point-lp
PT_LP_DIR = point/lp
PT_LP_FILES = pointLpMPC state_symeval point-lp logging
PT_LP_OBJS = $(PT_LP_FILES:%=$(OBJ_DIR)/%.o)

point-lp-setup:
	bash $(POINT_SETUP_FILE) lp $(T)

point-lp: | point-lp-setup $(PT_LP_OBJS)
	$(CXX) $(BFLAGS) $(PT_LP_OBJS) -o $(BIN_DIR)/point-lp $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/pointLpMPC.o : $(PT_LP_DIR)/lpMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-lp.o : $(PT_LP_DIR)/point-lp.cpp $(POINT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
# make point-lp-adolc
PT_LP_ADOLC_DIR = point/lp
PT_LP_ADOLC_FILES = pointLpMPC point-lp-adolc logging
PT_LP_ADOLC_OBJS = $(PT_LP_ADOLC_FILES:%=$(OBJ_DIR)/%.o)

point-lp-adolc: | $(PT_LP_ADOLC_OBJS)
	$(CXX) $(BFLAGS) $(PT_LP_ADOLC_OBJS) -o $(BIN_DIR)/point-lp-adolc $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS) $(ADOLC_FLAGS)
	
$(OBJ_DIR)/point-lp-adolc.o : $(PT_LP_ADOLC_DIR)/point-lp-adolc.cpp $(PT_ADOLC_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) $(ADOLC_FLAGS) -c -o $@ $<

# make point-control
PT_CONTROL_DIR = point/control
PT_CONTROL_FILES = pointControlMPC control_symeval point-control logging
PT_CONTROL_OBJS = $(PT_CONTROL_FILES:%=$(OBJ_DIR)/%.o)

point-control-setup:
	bash $(POINT_SETUP_FILE) control $(T)

point-control: | point-control-setup $(PT_CONTROL_OBJS) 
	$(CXX) $(BFLAGS) $(PT_CONTROL_OBJS) -o $(BIN_DIR)/point-control $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/pointControlMPC.o : $(PT_CONTROL_DIR)/controlMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-control.o : $(PT_CONTROL_DIR)/point-control.cpp $(PT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
# make point-ilqg
PT_ILQG_DIR = point/ilqg
PT_ILQG_FILES = point-ilqg
PT_ILQG_OBJS = $(PT_ILQG_FILES:%=$(OBJ_DIR)/%.o)

point-ilqg: $(PT_ILQG_OBJS)
	$(CXX) $(BFLAGS) $(PT_ILQG_OBJS) -o $(BIN_DIR)/point-ilqg $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/point-ilqg.o : $(PT_ILQG_DIR)/point-ilqg.cpp $(PT_ILQG_DIR)/ilqg.h
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) -c -o $@ $<
	
	
###### POINT-SLAM ##########

# share point functions
PT_SLAM_HEADERS = point-slam/point-slam.h

# make point-slam-test
PT_SLAM_TEST_DIR = point-slam/test
PT_SLAM_TEST_FILES = point-slam-test
PT_SLAM_TEST_OBJS = $(PT_SLAM_TEST_FILES:%=$(OBJ_DIR)/%.o)

touch-file:
	touch point-slam/test/point-slam-test.cpp

point-slam-test: | touch-file $(PT_SLAM_TEST_OBJS)
	$(CXX) $(BFLAGS) $(PT_SLAM_TEST_OBJS) -o $(BIN_DIR)/point-slam-test $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/point-slam-test.o : $(PT_SLAM_TEST_DIR)/point-slam-test.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) -c -o $@ $<

###### SLAM ##########

NUM_LANDMARKS = 3

# setupFiles bash
SLAM_SETUP_FILE = slam/setupFiles.sh

# share point functions
SLAM_HEADERS = slam/slam.h
SLAM_SRC = slam/slam.cpp

SLAM_CASADI_DIR = slam/casadi


$(OBJ_DIR)/slam.o : $(SLAM_SRC) $(SLAM_HEADERS)
	 $(CXX) $(CPP_FLAGS) $(BFLAGS)  $(PYTHON_FLAGS) $(BOOST_FLAGS) -c -o $@ $<

# make slam-test
SLAM_TEST_DIR = slam/test
SLAM_TEST_FILES = test-slam
SLAM_TEST_OBJS = $(SLAM_TEST_FILES:%=$(OBJ_DIR)/%.o)

touch-file-slam:
	touch slam/test/test-slam.cpp

test-slam: | touch-file-slam $(SLAM_TEST_OBJS)
	$(CXX) $(BFLAGS) $(SLAM_TEST_OBJS) -o $(BIN_DIR)/test-slam $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/test-slam.o : $(SLAM_TEST_DIR)/test-slam.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) -c -o $@ $<
	
# make slam-traj
SLAM_TRAJ_DIR = slam/traj
SLAM_TRAJ_FILES = trajMPC slam-traj logging slam
SLAM_TRAJ_OBJS = $(SLAM_TRAJ_FILES:%=$(OBJ_DIR)/%.o)

slam-traj-setup:
	bash $(SLAM_SETUP_FILE) traj $(T) $(NUM_LANDMARKS)

slam-traj: | slam-traj-setup $(SLAM_TRAJ_OBJS)
	$(CXX) $(SLAM_TRAJ_OBJS) $(BFLAGS) -o $(BIN_DIR)/slam-traj $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/trajMPC.o : $(SLAM_TRAJ_DIR)/trajMPC.c
	$(CC) $(C_FLAGS) -c -o $@ $^
	
$(OBJ_DIR)/slam-traj.o : $(SLAM_TRAJ_DIR)/slam-traj.cpp $(SLAM_TRAJ_DIR)/trajMPC.c
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
# make slam-traj-plan
SLAM_TRAJ_PLAN_DIR = slam/traj
SLAM_TRAJ_PLAN_FILES = trajMPC slam-traj slam-traj-plan logging slam
SLAM_TRAJ_PLAN_OBJS = $(SLAM_TRAJ_PLAN_FILES:%=$(OBJ_DIR)/%.o)

slam-traj-plan: | slam-traj-setup $(SLAM_TRAJ_PLAN_OBJS)
	$(CXX) $(SLAM_TRAJ_PLAN_OBJS) $(BFLAGS) -o $(BIN_DIR)/slam-traj-plan $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/slam-traj-plan.o : $(SLAM_TRAJ_PLAN_DIR)/slam-traj-plan.cpp
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
	
# make slam-smooth
SLAM_SMOOTH_DIR = slam/smooth
SLAM_SMOOTH_FILES = smoothMPC slam-smooth logging slam
SLAM_SMOOTH_OBJS = $(SLAM_SMOOTH_FILES:%=$(OBJ_DIR)/%.o)

slam-smooth-setup:
	bash $(SLAM_SETUP_FILE) smooth $(T) $(NUM_LANDMARKS)

slam-smooth: | $(SLAM_SMOOTH_OBJS)
	$(CXX) $(SLAM_SMOOTH_OBJS) $(BFLAGS) -o $(BIN_DIR)/slam-smooth $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/smoothMPC.o : $(SLAM_SMOOTH_DIR)/smoothMPC.c
	$(CC) $(C_FLAGS) -c -o $@ $^
	
$(OBJ_DIR)/slam-smooth.o : $(SLAM_SMOOTH_DIR)/slam-smooth.cpp $(SLAM_SMOOTH_DIR)/smoothMPC.c
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<

# make slam-belief
SLAM_BELIEF_DIR = slam/belief
SLAM_BELIEF_FILES = slamBeliefPenaltyMPC slam-belief logging slam slam-traj trajMPC
SLAM_BELIEF_OBJS = $(SLAM_BELIEF_FILES:%=$(OBJ_DIR)/%.o)

slam-belief-setup:
	bash $(SLAM_SETUP_FILE) belief $(T) $(NUM_LANDMARKS)

slam-belief: | slam-belief-setup slam-traj-setup $(SLAM_BELIEF_OBJS)
	$(CXX) $(BFLAGS) $(SLAM_BELIEF_OBJS) -o $(BIN_DIR)/slam-belief $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/slamBeliefPenaltyMPC.o : $(SLAM_BELIEF_DIR)/beliefPenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/slam-belief.o : $(SLAM_BELIEF_DIR)/slam-belief.cpp $(SLAM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
# make slam-state
SLAM_STATE_DIR = slam/state
# casadi-slam-state 
SLAM_STATE_FILES = slam-state-cost slam-state-grad slamStateMPC slam-state logging slam trajMPC slam-traj
SLAM_STATE_OBJS = $(SLAM_STATE_FILES:%=$(OBJ_DIR)/%.o)

slam-state-setup:
	bash $(SLAM_SETUP_FILE) state $(T) $(NUM_LANDMARKS)

slam-state: | slam-state-setup slam-traj-setup $(SLAM_STATE_OBJS)
	$(CXX) $(BFLAGS) $(SLAM_STATE_OBJS) -o $(BIN_DIR)/slam-state $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS) $(CASADI_FLAGS) $(CASADI_LIBS)
	
$(OBJ_DIR)/slam-state-cost.o : $(SLAM_STATE_DIR)/slam-state-cost.c
	$(CC) $(C_FLAGS) -c -o $@ $^

$(OBJ_DIR)/slam-state-grad.o : $(SLAM_STATE_DIR)/slam-state-grad.c
	$(CC) $(C_FLAGS) -c -o $@ $^
	
$(OBJ_DIR)/slamStateMPC.o : $(SLAM_STATE_DIR)/stateMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/casadi-slam-state.o : $(SLAM_CASADI_DIR)/slam-state.cpp
	$(CXX) $(CPP_FLAGS) -DNDEBUG -fopenmp $(BFLAGS) $(CASADI_FLAGS) -c -o $@ $< $(CASADI_LIBS) 
	
$(OBJ_DIR)/slam-state.o : $(SLAM_STATE_DIR)/slam-state.cpp $(SLAM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) $(CASADI_FLAGS) -c -o $@ $< $(CASADI_LIBS)
	
#casadi-slam-state:
#	$(CXX) -O3 -DNDEBUG -fopenmp -I$(CASADIPATH) -L$(CASADIPATH)/build/lib slam/casadi/slam-state.cpp -o casadi-slam-state -lcasadi -ldl
	
# make slam-state-mpc
SLAM_STATE_MPC_DIR = slam/state
SLAM_STATE_MPC_FILES = slam-state-cost slam-state-grad slamStateMPC slam-state-mpc logging slam trajMPC slam-traj
SLAM_STATE_MPC_OBJS = $(SLAM_STATE_MPC_FILES:%=$(OBJ_DIR)/%.o)

slam-state-mpc: | slam-state-setup slam-traj-setup $(SLAM_STATE_MPC_OBJS)
	$(CXX) $(BFLAGS) $(SLAM_STATE_MPC_OBJS) -o $(BIN_DIR)/slam-state-mpc $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/slam-state-mpc.o : $(SLAM_STATE_DIR)/slam-state-mpc.cpp $(SLAM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
# make slam-control
SLAM_CONTROL_DIR = slam/control
SLAM_CONTROL_FILES = slam-control-cost slam-control-grad slamControlMPC slam-control logging slam trajMPC slam-traj
SLAM_CONTROL_OBJS = $(SLAM_CONTROL_FILES:%=$(OBJ_DIR)/%.o)

slam-control-setup:
	bash $(SLAM_SETUP_FILE) control $(T) $(NUM_LANDMARKS)

slam-control: | slam-control-setup slam-traj-setup $(SLAM_CONTROL_OBJS)
	$(CXX) $(BFLAGS) $(SLAM_CONTROL_OBJS) -o $(BIN_DIR)/slam-control $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/slam-control-cost.o : $(SLAM_CONTROL_DIR)/slam-control-cost.c
	$(CC) $(C_FLAGS) -c -o $@ $^

$(OBJ_DIR)/slam-control-grad.o : $(SLAM_CONTROL_DIR)/slam-control-grad.c
	$(CC) $(C_FLAGS) -c -o $@ $^
	
$(OBJ_DIR)/slamControlMPC.o : $(SLAM_CONTROL_DIR)/controlMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/slam-control.o : $(SLAM_CONTROL_DIR)/slam-control.cpp $(SLAM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<

# make slam-ilqg
SLAM_ILQG_DIR = slam/ilqg
SLAM_ILQG_FILES = slam-ilqg logging slam
SLAM_ILQG_OBJS = $(SLAM_ILQG_FILES:%=$(OBJ_DIR)/%.o)

slam-ilqg: | $(SLAM_ILQG_OBJS)
	$(CXX) $(BFLAGS) $(SLAM_ILQG_OBJS) -o $(BIN_DIR)/slam-ilqg $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/slam-ilqg.o : $(SLAM_ILQG_DIR)/slam-ilqg.cpp $(SLAM_ILQG_DIR)/ilqg.h
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	

	



###### CLEAN ########

clean:
	rm -f -r bin/debug/*
	rm -f -r bin/release/*
	rm -f build/debug/*.o
	rm -f build/release/*.o
	
#### NOT FILES ######
.PHONY: point-belief-setup point-state-setup point-lp-setup point-control-setup arm-belief-setup
