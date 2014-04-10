/* Copyright (C) 2012 Benjamin Charrow
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <sys/time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "ros_gsl/random.hpp"

static const gsl_rng_type *global_rng_T;
static gsl_rng *global_rng = NULL;

void __attribute__ ((constructor)) init_random(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);

  global_rng_T = gsl_rng_mt19937;
  global_rng = gsl_rng_alloc(global_rng_T);
  gsl_rng_set(global_rng, tv.tv_usec);
}

void __attribute__ ((destructor)) free_random(void) {
  gsl_rng_free(global_rng);
}

namespace gsl {
  double uniform(double a, double b) {
    return gsl_ran_flat(global_rng, a, b);
  }

  double beta(double a, double b) {
    return gsl_ran_beta(global_rng, a, b);
  }
  
  double normal(double sd) {
    return gsl_ran_gaussian(global_rng, sd);
  }

  gsl_rng *rng() {
    return global_rng;
  }
}
