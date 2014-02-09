#!/bin/bash

# $1 := source file
# $2 := casadi type (Cost, Grad, Hess)
# $3 := example name (i.e. slam)

SRC_C=$1
CASADI_TYPE=$2
EXAMPLE_NAME=$3

# replace #include<math.h> with our own header file
sed -i "s/math.h/${EXAMPLE_NAME}-casadi.h/" $SRC_C

# delete these two functions
sed -i '/sq(d x)/d' $SRC_C
sed -i '/sign(d x)/d' $SRC_C

# replace function names
sed -i "s/getSparsity/get${CASADI_TYPE}Sparsity/" $SRC_C
sed -i "s/evaluate/evaluate${CASADI_TYPE}/" $SRC_C

# delete main at end
START_TAIL_CUT=$(grep -n "stdio.h" $SRC_C | cut -f1 -d: | head -1)
sed -i ${START_TAIL_CUT}',$d' $SRC_C