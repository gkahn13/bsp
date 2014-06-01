#!/bin/bash

# $1 := source file
# $2 := casadi type (Cost, Grad, Hess,Dyn,Obs)
# $3 := example name (i.e. slam)

SRC_C=$1
CASADI_TYPE=$2
EXAMPLE_NAME=$3

# replace #include<math.h> with our own header file
sed -i .bk "s/<math.h>/\"${EXAMPLE_NAME}-state-casadi.h\"/" $SRC_C

# replace sparsity names so no conflicts between files
sed -i .bk "s/s[0-9]/&${CASADI_TYPE}/" $SRC_C

# delete these two functions
sed -i .bk '/sq(d x)/d' $SRC_C
sed -i .bk '/sign(d x)/d' $SRC_C

# replace function names
sed -i .bk "s/init/init${CASADI_TYPE}/" $SRC_C
sed -i .bk "s/getSparsity/get${CASADI_TYPE}Sparsity/" $SRC_C
sed -i .bk "s/evaluate/evaluate${CASADI_TYPE}/" $SRC_C

# delete main at end
START_TAIL_CUT=$(grep -n "stdio.h" $SRC_C | cut -f1 -d: | head -1)
sed -i .bk ${START_TAIL_CUT}',$d' $SRC_C


rm -f *.bk

