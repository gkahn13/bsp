CXX = g++

CPP_FLAGS = -DNDEBUG 
CASADI_FLAGS = -I/Users/laskeymd/Documents/Research/BSP/casadi -L/Users/laskeymd/Documents/Research/BSP/casadi/build/lib
CASADI_LIBS = -lcasadi -ldl

OS := mac

LINKER_FLAGS.mac = -lm
LINKER_FLAGS.linux = -lrt
LINKER_FLAGS = $(LINKER_FLAGS.$(OS))

# default to debug build
BUILD := debug

BFLAGS.debug = -g
BFLAGS.release = -O3
BFLAGS := $(BFLAGS.$(BUILD))

# default to 15 T
T = 15

PARAM_STATE_FILES = parameter-controls 
PARAM_STATE_OBJS = $(PARAM_STATE_FILES:%= %.o)

parameter-controls: parameter-controls.cpp
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) -g parameter-controls.cpp -o parameter-controls $(CASADI_LIBS)


all: .FORCE
	./parameter-controls $(T)
	bash modifyCasadi.sh parameter-controls-cost.c Cost parameter
	bash modifyCasadi.sh parameter-controls-grad.c CostGrad parameter
	
clean:
	rm -f parameter-controls
	rm -f *.o
	
.FORCE:
	
.PHONY: .FORCE