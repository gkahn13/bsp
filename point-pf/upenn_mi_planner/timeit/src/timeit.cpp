/** Copyright (C) 2012 Benjamin Charrow
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

#include "timeit/timeit.hpp"

#include <sys/time.h>
#include <cstdio>

namespace timeit {
  static struct timeval tv1, tv2;

  void tic() {
    gettimeofday(&tv1, NULL);
  }

  void toc(const char *msg) {
    gettimeofday(&tv2, NULL);
    double t = (static_cast<double>((tv2.tv_sec - tv1.tv_sec)) +
                static_cast<double>(tv2.tv_usec - tv1.tv_usec) * 1e-6);
    fprintf(stderr, "%s = %10.4f msec\n", msg, 1e3 * t);
  }

  void toc(struct duration *d) {
    gettimeofday(&tv2, NULL);
    double t = (static_cast<double>((tv2.tv_sec - tv1.tv_sec)) +
                static_cast<double>(tv2.tv_usec - tv1.tv_usec) * 1e-6);
    d->start.tv_sec = tv1.tv_sec;
    d->start.tv_usec = tv1.tv_usec;
    d->stop.tv_sec = tv2.tv_sec;
    d->stop.tv_usec = tv2.tv_usec;
    d->msec = t * 1e3;
  }
}
