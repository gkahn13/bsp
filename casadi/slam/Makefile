CXX = g++

CPP_FLAGS = -DNDEBUG -fopenmp
CASADI_FLAGS = -I/home/gkahn/source/casadi/ -L/home/gkahn/source/casadi/build/lib
CASADI_LIBS = -lcasadi -ldl

# default to debug build
BUILD := release

BFLAGS.debug = -g
BFLAGS.release = -O3
BFLAGS := $(BFLAGS.$(BUILD))

# default to 15 T
T = 15

all: slam-state.cpp | clean
	$(CXX) $(BFLAGS) $(CPP_FLAGS) $(CASADI_FLAGS) slam-state.cpp -o slam-state $(CASADI_LIBS)
	./slam-state $(T)
	bash modifyCasadi.sh slam-state-cost.c Cost slam
	bash modifyCasadi.sh slam-state-grad.c CostGrad slam
	
clean:
	rm -f *.c
	rm -f slam-state