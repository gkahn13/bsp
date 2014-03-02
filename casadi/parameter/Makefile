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

PARAM_STATE_FILES = parameter-state parameter-dyndx parameter-obsdx
PARAM_STATE_OBJS = $(PARAM_STATE_FILES:%= %.o)


parameter-state: | $(PARAM_STATE_OBJS)
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) $(PARAM_STATE_OBJS) -o parameter-state $(CASADI_LIBS)  $(LINKER_FLAGS)
	
parameter-dyndx.o : parameter-dyndx.c
	$(CC) $(C_FLAGS) -c -o $@ $^

parameter-obsdx.o : parameter-obsdx.c
	$(CC) $(C_FLAGS) -c -o $@ $^

parameter-state.o : parameter-state.cpp 
	$(CC) $(C_FLAGS) $(CASADI_FLAGS) $(BFLAGS) -c -o $@ $<

#$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) $(PARAM_STATE_OBJS) -o parameter-state $(CASADI_LIBS)

parameter-jac: parameter-jac.cpp
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) -g parameter-jac.cpp -o parameter-jac $(CASADI_LIBS)

fix_jac: 
	./parameter-jac
	bash modifyCasadi.sh parameter-obsdx.c Obs parameter
	bash modifyCasadi.sh parameter-dyndx.c Dyn parameter

all: .FORCE
	./parameter-state $(T)
	bash modifyCasadi.sh parameter-state-cost.c Cost parameter
	bash modifyCasadi.sh parameter-state-grad.c CostGrad parameter
	
clean:
	rm -f parameter-jac
	rm -f parameter-state
	rm -f *.o
	
.FORCE:
	
.PHONY: .FORCE