CXX = g++
CC = gcc

OS := linux

LINKER_FLAGS.mac = -lm
LINKER_FLAGS.linux = -lrt
LINKER_FLAGS = $(LINKER_FLAGS.$(OS))

C_FLAGS = -std=gnu99
CPP_FLAGS = -std=c++0x -fpermissive -I.
PYTHON_FLAGS = -I/usr/include/python2.7 -lpython2.7
BOOST_FLAGS = -lboost_python -lboost_filesystem -lboost_system
CASADI_FLAGS = -I$(CASADIPATH) -L$(CASADIPATH)/build/lib
CASADI_LIBS = -lcasadi -ldl

# default to debug build
BUILD := debug

BFLAGS.debug = -g -Wall -Wno-sign-compare 
BFLAGS.release = -O3
BFLAGS := $(BFLAGS.$(BUILD))

OBJ_DIR = build/$(BUILD)
BIN_DIR = bin/$(BUILD)

# default TIMESTEPS
T = 10

###### DIRECTORY CREATION ###########
dir:
	mkdir -p build/debug
	mkdir -p build/release
	mkdir -p bin/debug
	mkdir -p bin/release

###### UTIL ############

$(OBJ_DIR)/logging.o : ../util/logging.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $^
	
####### test-point-explore #######

TEST_POINT_EXPLORE_FILES = test-point-explore logging
TEST_POINT_EXPLORE_OBJS = $(TEST_POINT_EXPLORE_FILES:%=$(OBJ_DIR)/%.o)

test-point-explore: $(TEST_POINT_EXPLORE_OBJS)
	$(CXX) $(BFLAGS) $(TEST_POINT_EXPLORE_OBJS) -o $(BIN_DIR)/test-point-explore $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS) 
	
$(OBJ_DIR)/test-point-explore.o : test-point-explore.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) -c -o $@ $^
	

#######  explore ########

POINT_EXPLORE_FILES = explore-cpp exploreMPC logging
POINT_EXPLORE_OBJS = $(POINT_EXPLORE_FILES:%=$(OBJ_DIR)/%.o)

explore-setup:
	bash setupFiles.sh explore $(T)

explore: explore-setup $(POINT_EXPLORE_OBJS) 
	$(CXX) $(BFLAGS) $(POINT_EXPLORE_OBJS) -o $(BIN_DIR)/explore $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS) $(CASADI_FLAGS) $(CASADI_LIBS)

$(OBJ_DIR)/exploreMPC.o : explore/exploreMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/explore-cpp.o : explore/explore.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)


###### CLEAN ########

clean:
	rm -f bin/debug/*
	rm -f bin/release/*
	rm -f build/debug/*
	rm -f build/release/*
	
.PHONY:
	platt-setup


