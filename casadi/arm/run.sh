#!/bin/bash

rm -f f.c grad_f.c hess_f.c f grad hess
rm -f arm-state arm-dynobsjac 

g++ -O3 -DNDEBUG -fopenmp -I/home/sachin/Workspace/casadi/ -L/home/sachin/Workspace/casadi/build/lib arm-state.cpp -o arm-state -lcasadi -ldl
echo Finished compilation
./arm-state

#g++ -I/home/sachin/Workspace/casadi/ -L/home/sachin/Workspace/casadi/build/lib arm-dynobsjac.cpp -o arm-dynobsjac -lcasadi -ldl
#echo Finished compilation
#./arm-dynobsjac
