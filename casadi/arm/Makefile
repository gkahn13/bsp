CXX = g++

CPP_FLAGS = -DNDEBUG 
CASADI_FLAGS = -I/Users/laskeymd/Documents/Research/BSP/casadi -L/Users/laskeymd/Documents/Research/BSP/casadi/build/lib
CASADI_LIBS = -lcasadi -ldl

OS := mac

LINKER_FLAGS.mac = -lm
LINKER_FLAGS.linux = -lrt
LINKER_FLAGS = $(LINKER_FLAGS.$(OS))

# default to debug build
BUILD := release

BFLAGS.debug = -g
BFLAGS.release = -O3
BFLAGS := $(BFLAGS.$(BUILD))

# default to 15 T
T = 15

PARAM_STATE_FILES = arm-control 
PARAM_STATE_OBJS = $(PARAM_STATE_FILES:%= %.o)

arm-control: arm-control.cpp
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) -g arm-control.cpp -o arm-control $(CASADI_LIBS)

arm-state: arm-state.cpp
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) -g arm-state.cpp -o arm-state $(CASADI_LIBS)

all: .FORCE
	./arm-control $(T)
	./arm-state $(T)
	bash modifyCasadi.sh arm-control-cost.c Cost arm control
	bash modifyCasadi.sh arm-control-grad.c CostGrad arm control
	bash modifyCasadi.sh arm-state-cost.c Cost arm state 
	bash modifyCasadi.sh arm-state-grad.c CostGrad arm state 
	
clean:
	rm -f parameter-jac
	rm -f arm-control
	rm -f *.o
	
.FORCE:
	
.PHONY: .FORCE