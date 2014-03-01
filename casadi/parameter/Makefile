CXX = g++

CPP_FLAGS = -DNDEBUG 
CASADI_FLAGS = -I$(CASADIPATH) -L$(CASADIPATH)/build/lib
CASADI_LIBS = -lcasadi -ldl

# default to debug build
BUILD := release

BFLAGS.debug = -g
BFLAGS.release = -O3
BFLAGS := $(BFLAGS.$(BUILD))

# default to 15 T
T = 15

PARAM_STATE_FILES = parameter-state parameter-dyndx parameter-obsdx
PARAM_STATE_OBJS = $(PARAM_STATE_FILES:%= %.o)


parameter-state: parameter-state.cpp
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) -g parameter-state.cpp $(PARAM_STATE_OBJS) -o parameter-state $(CASADI_LIBS)

	parameter-dyndx.o : parameter-dyndx.c
	$(CC) $(C_FLAGS) -c -o $@ $^

	parameter-obsdx.o : parameter-obsdx.c
	$(CC) $(C_FLAGS) -c -o $@ 


parameter-jac: parameter-jac.cpp
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) -g parameter-jac.cpp -o parameter-jac $(CASADI_LIBS)


all: .FORCE
	./parameter-state $(T)
	bash modifyCasadi.sh parameter-state-cost.c Cost parameter
	bash modifyCasadi.sh parameter-state-grad.c CostGrad parameter
	bash modifyCasadi.sh parameter-obsdx.c Obs parameter
	bash modifyCasadi.sh parameter-dyndx.c Dyn parameter
clean:
	rm -f *.c
	rm -f -state
	
.FORCE:
	
.PHONY: .FORCE