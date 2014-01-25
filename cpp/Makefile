CXX = g++
CC = gcc

C_FLAGS = -std=gnu99
CPP_FLAGS = -std=c++0x -fpermissive -I.
PYTHON_FLAGS = -I/usr/include/python2.7 -lpython2.7
BOOST_FLAGS = -lboost_python -lboost_filesystem -lboost_system
LINKER_FLAGS = -lrt
ADOLC_FLAGS = -I/usr/include/adolc -ladolc
EIGEN_FLAGS = -I/usr/include/eigen3

# default to debug build
BUILD := debug

BFLAGS.debug = 
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
ARM_HEADERS = arm/arm.h

# make arm-state
ARM_STATE_DIR = arm/state
ARM_STATE_FILES = armStateMPC arm-state logging
ARM_STATE_OBJS = $(ARM_STATE_FILES:%=$(OBJ_DIR)/%.o)

arm-state: | $(ARM_STATE_OBJS)
	$(CXX) $(BFLAGS) $(ARM_STATE_OBJS) -o $(BIN_DIR)/arm-state $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS) $(ADOLC_FLAGS)

$(OBJ_DIR)/armStateMPC.o : $(ARM_STATE_DIR)/armStateMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/arm-state.o : $(ARM_STATE_DIR)/arm-state.cpp $(ARM_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) $(ADOLC_FLAGS) -c -o $@ $<



###### POINT ##########

# setupFiles bash
POINT_SETUP_FILE = point/setupFiles.sh

# share point functions
PT_HEADERS = point/point.h

# symeval	
PT_SYM_DIR = point/sym

$(OBJ_DIR)/state_symeval.o : $(PT_SYM_DIR)/state-symeval.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/control_symeval.o : $(PT_SYM_DIR)/control-symeval.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
# make point-belief
PT_BELIEF_DIR = point/belief
PT_BELIEF_FILES = beliefMPC beliefPenaltyMPC point-belief logging
PT_BELIEF_OBJS = $(PT_BELIEF_FILES:%=$(OBJ_DIR)/%.o)

belief-setup:
	bash $(POINT_SETUP_FILE) belief $(T)

point-belief: | belief-setup $(PT_BELIEF_OBJS)
	$(CXX) $(BFLAGS) $(PT_BELIEF_OBJS) -o $(BIN_DIR)/point-belief $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/beliefMPC.o : $(PT_BELIEF_DIR)/beliefMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/beliefPenaltyMPC.o : $(PT_BELIEF_DIR)/beliefPenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-belief.o : $(PT_BELIEF_DIR)/point-belief.cpp $(POINT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
	
# make point-state
PT_STATE_DIR = point/state
PT_STATE_FILES = stateMPC state_symeval point-state logging
PT_STATE_OBJS = $(PT_STATE_FILES:%=$(OBJ_DIR)/%.o)

state-setup:
	bash $(POINT_SETUP_FILE) state $(T)

point-state: | state-setup $(PT_STATE_OBJS)
	$(CXX) $(BFLAGS) $(PT_STATE_OBJS) -o $(BIN_DIR)/point-state $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/stateMPC.o : $(PT_STATE_DIR)/stateMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-state.o : $(PT_STATE_DIR)/point-state.cpp $(PT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
# make point-state-finitediff
PT_STATE_FD_DIR = $(PT_STATE_DIR)
PT_STATE_FD_FILES = stateMPC state_symeval point-state-finitediff logging
PT_STATE_FD_OBJS = $(PT_STATE_FD_FILES:%=$(OBJ_DIR)/%.o)

point-state-finitediff: | state-setup $(PT_STATE_FD_OBJS) 
	$(CXX) $(BFLAGS) $(PT_STATE_FD_OBJS) -o $(BIN_DIR)/point-state-finitediff $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/point-state-finitediff.o : $(PT_STATE_FD_DIR)/point-state-finitediff.cpp $(POINT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<
	
# make point-state-adolc
PT_ADOLC_HEADERS = point/point.h
PT_STATE_ADOLC_DIR = point/state
PT_STATE_ADOLC_FILES = stateMPC point-state-adolc logging
PT_STATE_ADOLC_OBJS = $(PT_STATE_ADOLC_FILES:%=$(OBJ_DIR)/%.o)

point-state-adolc: | $(PT_STATE_ADOLC_OBJS)
	$(CXX) $(BFLAGS) $(PT_STATE_ADOLC_OBJS) -o $(BIN_DIR)/point-state-adolc $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS) $(ADOLC_FLAGS)
	
$(OBJ_DIR)/point-state-adolc.o : $(PT_STATE_ADOLC_DIR)/point-state-adolc.cpp $(PT_ADOLC_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) $(ADOLC_FLAGS) -c -o $@ $<

# make point-lp
PT_LP_DIR = point/lp
PT_LP_FILES = lpMPC state_symeval point-lp logging
PT_LP_OBJS = $(PT_LP_FILES:%=$(OBJ_DIR)/%.o)

lp-setup:
	bash $(POINT_SETUP_FILE) lp $(T)

point-lp: | lp-setup $(PT_LP_OBJS)
	$(CXX) $(BFLAGS) $(PT_LP_OBJS) -o $(BIN_DIR)/point-lp $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/lpMPC.o : $(PT_LP_DIR)/lpMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-lp.o : $(PT_LP_DIR)/point-lp.cpp $(POINT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<

# make point-control
PT_CONTROL_DIR = point/control
PT_CONTROL_FILES = controlMPC control_symeval point-control logging
PT_CONTROL_OBJS = $(PT_CONTROL_FILES:%=$(OBJ_DIR)/%.o)

control-setup:
	bash $(POINT_SETUP_FILE) control $(T)

point-control: | control-setup $(PT_CONTROL_OBJS) 
	$(CXX) $(BFLAGS) $(PT_CONTROL_OBJS) -o $(BIN_DIR)/point-control $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/controlMPC.o : $(PT_CONTROL_DIR)/controlMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-control.o : $(PT_CONTROL_DIR)/point-control.cpp $(PT_HEADERS)
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $<

###### CLEAN ########

clean:
	rm -f $(BIN_DIR)/*
	rm -f $(OBJ_DIR)/*.o
	
#### NOT FILES ######
.PHONY: belief-setup state-setup lp-setup control-setup
