CXX = g++
CC = gcc

C_FLAGS = -std=gnu99
CPP_FLAGS = -std=c++0x -fpermissive -I.
PYTHON_FLAGS = -I/usr/include/python2.7 -lpython2.7
BOOST_FLAGS = -lboost_python -lboost_filesystem -lboost_system
LINKER_FLAGS = -lrt
ADOLC_FLAGS = -I/home/gkahn/adolc_base/include -L/home/gkahn/adolc_base/lib64 -ladolc -lColPack
#ADOLC_FLAGS = -I/usr/include/adolc -ladolc -lColPack
EIGEN_FLAGS = -I/usr/include/eigen3

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
ARM_STATE_FILES = arm-state-pos-goal arm-state-cost arm-state-grad arm-state logging
ARM_STATE_OBJS = $(ARM_STATE_FILES:%=$(OBJ_DIR)/%.o)

arm-state: | $(ARM_STATE_OBJS)
	$(CXX) $(BFLAGS) $(ARM_STATE_OBJS) -o $(BIN_DIR)/arm-state $(BOOST_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/arm-state-cost.o : $(ARM_STATE_DIR)/arm-state-cost.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/arm-state-grad.o : $(ARM_STATE_DIR)/arm-state-grad.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/arm-state-pos-goal.o : $(ARM_STATE_DIR)/arm-state-pos-goal-MPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/arm-state.o : $(ARM_STATE_DIR)/arm-state-pos-goal.cpp $(ARM_HEADERS)
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
	bash $(PARAM_SETUP_FILE) belief $(T)

parameter-belief: | parameter-belief-setup $(PARAM_BELIEF_OBJS)
	$(CXX) $(BFLAGS) $(PARAM_BELIEF_OBJS) -o $(BIN_DIR)/parameter-belief $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/parameterBeliefPenaltyMPC.o : $(PARAM_BELIEF_DIR)/beliefPenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/parameter-belief.o : $(PARAM_BELIEF_DIR)/parameter-belief.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	

# make parameter-ilqg
PARAM_ILQG_DIR = parameter/ilqg
PARAM_ILQG_FILES = parameter-ilqg
PARAM_ILQG_OBJS = $(PARAM_ILQG_FILES:%=$(OBJ_DIR)/%.o)

parameter-ilqg: $(PARAM_ILQG_OBJS)
	$(CXX) $(BFLAGS) $(PARAM_ILQG_OBJS) -o $(BIN_DIR)/parameter-ilqg $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/parameter-ilqg.o : $(PARAM_ILQG_DIR)/parameter-ilqg.cpp $(PARAM_ILQG_DIR)/ilqg.h
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



###### CLEAN ########

clean:
	rm -f $(BIN_DIR)/*
	rm -f $(OBJ_DIR)/*.o
	
#### NOT FILES ######
.PHONY: point-belief-setup point-state-setup point-lp-setup point-control-setup arm-belief-setup
