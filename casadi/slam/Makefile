CXX = g++

CPP_FLAGS = -DNDEBUG -fopenmp
CASADI_FLAGS = -I$(CASADIPATH) -L$(CASADIPATH)/build/lib
CASADI_LIBS = -lcasadi -ldl

# default to debug build
BUILD := release

BFLAGS.debug = -g
BFLAGS.release = -O3
BFLAGS := $(BFLAGS.$(BUILD))

# default to 15 T
T = 15

slam-state: .FORCE
	rm -f *.c
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) slam-state.cpp -o slam-state $(CASADI_LIBS)
	./slam-state $(T)
	bash modifyCasadi.sh slam-state-cost.c Cost state
	bash modifyCasadi.sh slam-state-grad.c CostGrad state

slam-control: .FORCE
	rm -f *.c
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) slam-control.cpp -o slam-control $(CASADI_LIBS)
	./slam-control $(T)
	bash modifyCasadi.sh slam-control-cost.c Cost control
	bash modifyCasadi.sh slam-control-grad.c CostGrad control
	
clean:
	rm -f *.c
	rm -f slam-state
	
.FORCE:
	
.PHONY: .FORCE