#!/bin/bash

# $1 := source file
# $2 := casadi type (Cost, Grad, Hess,Dyn,Obs)
# $3 := example name (i.e. slam)

SRC_C=$1
CASADI_TYPE=$2
EXAMPLE_NAME=$3
EXAMPLE_TYPE=$4

# replace #include<math.h> with our own header file
sed -i "s/<math.h>/\"${EXAMPLE_NAME}-${EXAMPLE_TYPE}-casadi.h\"/" $SRC_C

# replace sparsity names so no conflicts between files
sed -i "s/s[0-9]/&${CASADI_TYPE}/" $SRC_C

# delete these two functions
sed -i '/sq(d x)/d' $SRC_C
sed -i '/sign(d x)/d' $SRC_C

# replace function names
sed -i "s/init/init${CASADI_TYPE}/" $SRC_C
sed -i "s/getSparsity/get${CASADI_TYPE}Sparsity/" $SRC_C
sed -i "s/evaluate/evaluate${CASADI_TYPE}/" $SRC_C

# delete main at end
START_TAIL_CUT=$(grep -n "stdio.h" $SRC_C | cut -f1 -d: | head -1)
sed -i ${START_TAIL_CUT}',$d' $SRC_C


rm -f *.bk

