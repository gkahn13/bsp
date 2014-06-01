CXX = g++
CC = gcc

OS := linux

LINKER_FLAGS.mac = -lm
LINKER_FLAGS.linux = -lrt
LINKER_FLAGS = $(LINKER_FLAGS.$(OS)) -larmadillo

C_FLAGS = -std=gnu99
CPP_FLAGS = -std=c++0x -fpermissive -I.
PYTHON_FLAGS = -I/usr/include/python2.7 -lpython2.7
BOOST_FLAGS = -lboost_python -lboost_filesystem -lboost_system -lboost_program_options
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
# default AGENTS
N = 1

###### DIRECTORY CREATION ###########
dir:
	mkdir -p build/debug
	mkdir -p build/release
	mkdir -p bin/debug
	mkdir -p bin/release

###### util ############

$(OBJ_DIR)/logging.o : ../util/logging.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $^
	
##### system #####

$(OBJ_DIR)/system.o : system.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)
	
$(OBJ_DIR)/casadi-system.o : casadi-system.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)
	
#######  explore ########

EXPLORE_FILES = explore-cpp exploreMPC system casadi-system casadi-explore-system explore-system logging
EXPLORE_OBJS = $(EXPLORE_FILES:%=$(OBJ_DIR)/%.o)

explore-setup:
	bash explore/setupFiles.sh $(T) $(N)

explore: explore-setup $(EXPLORE_OBJS) 
	$(CXX) $(BFLAGS) $(EXPLORE_OBJS) -o $(BIN_DIR)/explore $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS) $(CASADI_FLAGS) $(CASADI_LIBS)

$(OBJ_DIR)/exploreMPC.o : explore/exploreMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/explore-system.o : explore/explore-system.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)

$(OBJ_DIR)/casadi-explore-system.o : explore/casadi-explore-system.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)
	
$(OBJ_DIR)/explore-cpp.o : explore/explore.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)

###### boxes ########

BOXES_FILES = boxes-cpp boxesMPC system casadi-system boxes-system casadi-boxes-system logging
BOXES_OBJS = $(BOXES_FILES:%=$(OBJ_DIR)/%.o)

boxes-setup:
	#bash boxes/setupFiles.sh $(T) $(M)

boxes: boxes-setup $(BOXES_OBJS) 
	$(CXX) $(BFLAGS) $(BOXES_OBJS) -o $(BIN_DIR)/boxes $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS) $(CASADI_FLAGS) $(CASADI_LIBS)

$(OBJ_DIR)/boxesMPC.o : boxes/boxesMPC.c
	$(CC) $(C_FLAGS) $(BFLAGS) -c -o $@ $^

$(OBJ_DIR)/boxes-system.o : boxes/boxes-system.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)

$(OBJ_DIR)/casadi-boxes-system.o : boxes/casadi-boxes-system.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)

$(OBJ_DIR)/boxes-cpp.o : boxes/boxes.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)



##### kdpee ##############

$(OBJ_DIR)/kdpee.o : kdpee/src/kdpee.c
	$(CXX) $(CPP_FLAGS) $(BFLAGS) -c -o $@ $^
	

	
####### test-point-explore #######

TEST_POINT_EXPLORE_FILES = test-point-explore point-explore-system casadi-point-explore-system logging
TEST_POINT_EXPLORE_OBJS = $(TEST_POINT_EXPLORE_FILES:%=$(OBJ_DIR)/%.o)

test-point-explore: $(TEST_POINT_EXPLORE_OBJS)
	$(CXX) $(BFLAGS) $(TEST_POINT_EXPLORE_OBJS) -o $(BIN_DIR)/test-point-explore $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS) $(CASADI_FLAGS) $(CASADI_LIBS) 
	
$(OBJ_DIR)/test-point-explore.o : tests/test-point-explore.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)
	
####### test-entropy #######

TEST_ENTROPY_FILES = test-entropy kdpee logging
TEST_ENTROPY_OBJS = $(TEST_ENTROPY_FILES:%=$(OBJ_DIR)/%.o)

test-entropy: $(TEST_ENTROPY_OBJS)
	$(CXX) $(BFLAGS) $(TEST_ENTROPY_OBJS) -o $(BIN_DIR)/test-entropy $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS) $(CASADI_FLAGS) $(CASADI_LIBS) 
	
$(OBJ_DIR)/test-entropy.o : tests/test-entropy.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)


##### test-casadi-point-explore

TEST_CASADI_POINT_EXPLORE_FILES = test-casadi-point-explore casadi-point-explore point-explore logging
TEST_CASADI_POINT_EXPLORE_OBJS = $(TEST_CASADI_POINT_EXPLORE_FILES:%=$(OBJ_DIR)/%.o)

test-casadi-point-explore: $(TEST_CASADI_POINT_EXPLORE_OBJS) 
	$(CXX) $(BFLAGS) $(TEST_CASADI_POINT_EXPLORE_OBJS) -o $(BIN_DIR)/test-casadi-point-explore $(PYTHON_FLAGS) $(BOOST_FLAGS) $(LINKER_FLAGS) $(CASADI_FLAGS) $(CASADI_LIBS)

$(OBJ_DIR)/test-casadi-point-explore.o : casadi/test-casadi-point-explore.cpp
	$(CXX) $(CPP_FLAGS) $(BFLAGS) $(PYTHON_FLAGS) $(BOOST_FLAGS) $(CASADI_FLAGS) -c -o $@ $^ $(CASADI_LIBS)
	

###### CLEAN ########

clean:
	rm -f bin/debug/*
	rm -f bin/release/*
	rm -f build/debug/*
	rm -f build/release/*
	
.PHONY:
	platt-setup


