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

$(OBJ_DIR)/logging.o : ../util/logging.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $^
	
####### test-point-pf #######

TEST_POINT_PF_FILES = test-point-pf logging
TEST_POINT_PF_OBJS = $(TEST_POINT_PF_FILES:%=$(OBJ_DIR)/%.o)

test-point-pf: $(TEST_POINT_PF_OBJS)
	$(CXX) $(BFLAGS) $(TEST_POINT_PF_OBJS) -o $(BIN_DIR)/test-point-pf $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS)
	
$(OBJ_DIR)/test-point-pf.o : test-point-pf.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) -c -o $@ $^
	

###### CLEAN ########

clean:
	rm -f build/debug/*.o
	rm -f build/release/*.o


