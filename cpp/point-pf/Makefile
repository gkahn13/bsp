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
T = 15
# default PARTICLES
M = 5

###### DIRECTORY CREATION ###########
dir:
	mkdir -p build/debug
	mkdir -p build/release
	mkdir -p bin/debug
	mkdir -p bin/release

###### UTIL ############

$(OBJ_DIR)/logging.o : ../util/logging.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $^
	
####### test-point-pf #######

TEST_POINT_PF_FILES = test-point-pf logging
TEST_POINT_PF_OBJS = $(TEST_POINT_PF_FILES:%=$(OBJ_DIR)/%.o)

test-point-pf: $(TEST_POINT_PF_OBJS)
	$(CXX) $(BFLAGS) $(TEST_POINT_PF_OBJS) -o $(BIN_DIR)/test-point-pf $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) $(LINKER_FLAGS) $(CASADI_LIBS)
	
$(OBJ_DIR)/test-point-pf.o : test-point-pf.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)
	
#######  point-platt ########

POINT_PLATT_FILES = platt-cpp plattMPC logging
POINT_PLATT_OBJS = $(POINT_PLATT_FILES:%=$(OBJ_DIR)/%.o)

platt-setup:
	bash setupFiles.sh platt $(T) $(M)

platt: platt-setup $(POINT_PLATT_OBJS) 
	$(CXX) $(BFLAGS) $(POINT_PLATT_OBJS) -o $(BIN_DIR)/point-platt $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) $(LINKER_FLAGS) $(CASADI_LIBS)

$(OBJ_DIR)/plattMPC.o : platt/plattMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/platt-cpp.o : platt/platt.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)


###### CLEAN ########

clean:
	rm -f bin/debug/*
	rm -f bin/release/*
	rm -f build/debug/*
	rm -f build/release/*
	
.PHONY:
	platt-setup


