#!/bin/bash

rm -f f.c grad_f.c hess_f.c f grad hess
rm -f slam-state

g++ -g -DNDEBUG -fopenmp -I/home/gkahn/source/casadi/ -L/home/gkahn/source/casadi/build/lib slam-state.cpp -o slam-state -lcasadi -ldl
echo Finished compilation