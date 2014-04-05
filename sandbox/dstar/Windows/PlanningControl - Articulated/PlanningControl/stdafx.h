// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#define _CRT_RAND_S

#define DIM 3
#define X_DIM 6
#define U_DIM 6
#define Z_DIM 4

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <vector>
#include <map>
#include <queue>
#include <time.h>
#include <float.h>

#include "callisto.h"
#include "callistoTypes.h"
#include "matrix.h"
#include "igammaf.h"

struct RRTNode {
  Matrix<X_DIM> x;
  Matrix<U_DIM> u;
  int bp;
  Matrix<DIM> p;
};

struct PathNode {
  Matrix<X_DIM> x;
  Matrix<U_DIM> u;
};



// TODO: reference additional headers your program requires here
