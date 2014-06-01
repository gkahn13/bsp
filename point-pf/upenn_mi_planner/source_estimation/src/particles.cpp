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

#include "particles.hpp"

namespace rf {
  std::ostream& operator<<(std::ostream &s, const struct StationaryParticle &sp) {
    s << "(" << sp.point.x << ", " << sp.point.y << ")" << " weight: " << sp.weight;
    return s;
  }

  //============================ ParticleBuilder ============================//
  template<>
  ParticleArray<StationaryParticle>*
  ParticleBuilder<StationaryParticle>::NewParticles(int num,
                                                    const Measurement &m) {
    ParticleArray<StationaryParticle> *pa =
      new ParticleArray<StationaryParticle>(num);

    for (int i = 0; i < num; ++i) {
      double theta = gsl::uniform(0, 2 * M_PI);
      double r = gsl::normal(rsd_) + m.dist();
      pa->ps[i].point.x = m.sourceLoc().x() + r * cos(theta);
      pa->ps[i].point.y = m.sourceLoc().y() + r * sin(theta);
    }
    return pa;
  }

}
