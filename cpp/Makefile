CXX = g++
CC = gcc

# default to debug build
BUILD := debug

BFLAGS.debug = -g -Wall
BFLAGS.release = -O3
BFLAGS := $(BFLAGS.$(BUILD))

OBJ_DIR = build/$(BUILD)
BIN_DIR = bin/$(BUILD)

C_FLAGS = -std=gnu99
CPP_FLAGS = -std=c++0x -I.
PYTHON_FLAGS = -I/usr/include/python2.7 -lpython2.7
BOOST_FLAGS = -lboost_python -lboost_filesystem -lboost_system
LINKER_FLAGS = -lrt

###### DIRECTORY CREATION ###########
dir:
	mkdir -p build/debug
	mkdir -p build/release
	mkdir -p bin/debug
	mkdir -p bin/release

###### UTIL ############

$(OBJ_DIR)/logging.o : util/logging.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $^

###### POINT ##########

# symeval	
PT_SYM_DIR = point/sym

$(OBJ_DIR)/pt_symeval.o : $(PT_SYM_DIR)/symeval.c
	gcc -c $(C_FLAGS) $(BFLAGS) -o $(OBJ_DIR)/pt_symeval.o $(PT_SYM_DIR)/symeval.c

# make point-belief
PT_BELIEF_DIR = point/belief
PT_BELIEF_FILES = beliefMPC beliefPenaltyMPC point-belief logging
PT_BELIEF_OBJS = $(PT_BELIEF_FILES:%=$(OBJ_DIR)/%.o)

belief: $(PT_BELIEF_OBJS)
	$(CXX) $(BFLAGS) $(PT_BELIEF_OBJS) -o $(BIN_DIR)/point-belief $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/beliefMPC.o : $(PT_BELIEF_DIR)/beliefMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/beliefPenaltyMPC.o : $(PT_BELIEF_DIR)/beliefPenaltyMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-belief.o : $(PT_BELIEF_DIR)/point-belief.cpp
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $^
	
	
# make point-state
PT_STATE_DIR = point/state
PT_STATE_FILES = stateMPC pt_symeval point-state logging
PT_STATE_OBJS = $(PT_STATE_FILES:%=$(OBJ_DIR)/%.o)

state: $(PT_STATE_OBJS)
	$(CXX) $(BFLAGS) $(PT_STATE_OBJS) -o $(BIN_DIR)/point-state $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/stateMPC.o : $(PT_STATE_DIR)/stateMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-state.o : $(PT_STATE_DIR)/point-state.cpp
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $^
	
# make point-state-finitediff
PT_STATE_FD_DIR = $(PT_STATE_DIR)
PT_STATE_FD_FILES = stateMPC pt_symeval point-state-finitediff logging
PT_STATE_FD_OBJS = $(PT_STATE_FD_FILES:%=$(OBJ_DIR)/%.o)

state-finitediff: $(PT_STATE_FD_OBJS)
	$(CXX) $(BFLAGS) $(PT_STATE_FD_OBJS) -o $(BIN_DIR)/point-state-finitediff $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/point-state-finitediff.o : $(PT_STATE_FD_DIR)/point-state-finitediff.cpp
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $^


# make point-lp
PT_LP_DIR = point/lp
PT_LP_FILES = lpMPC pt_symeval point-lp logging
PT_LP_OBJS = $(PT_LP_FILES:%=$(OBJ_DIR)/%.o)

lp: $(PT_LP_OBJS)
	$(CXX) $(BFLAGS) $(PT_LP_OBJS) -o $(BIN_DIR)/point-lp $(BOOST_FLAGS) $(PYTHON_FLAGS) $(LINKER_FLAGS)

$(OBJ_DIR)/lpMPC.o : $(PT_LP_DIR)/lpMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^
	
$(OBJ_DIR)/point-lp.o : $(PT_LP_DIR)/point-lp.cpp
	$(CXX) $(CPP_FLAGS) $(PYTHON_FLAGS) $(BFLAGS) -c -o $@ $^


###### CLEAN ########

clean:
	rm -f $(BIN_DIR)/*
	rm -f $(OBJ_DIR)/*.o
