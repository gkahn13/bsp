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

#include "gtest/gtest.h"

#include "particles.hpp"

using namespace rf;

TEST(ParticleArray, FilterCDFRemoveNonec) {
  ParticleArray<> pa_orig;
  pa_orig.resize(2);
  pa_orig.ps[0].point.x = 2.29;
  pa_orig.ps[0].point.y = -8.12;
  pa_orig.ps[0].weight = 0.68;
  
  pa_orig.ps[1].point.x = 12.34;
  pa_orig.ps[1].point.y = 5.29;
  pa_orig.ps[1].weight = 0.32;


  ParticleArray<> *sub = pa_orig.FilterCDF(0.1);
  ASSERT_TRUE(sub != NULL);
  ASSERT_EQ(2, sub->size);
  EXPECT_DOUBLE_EQ(sub->ps[0].point.x, 12.34);
  EXPECT_DOUBLE_EQ(sub->ps[0].point.y, 5.29);
  EXPECT_DOUBLE_EQ(sub->ps[0].weight, 0.32);
  EXPECT_DOUBLE_EQ(sub->ps[1].point.x, 2.29);
  EXPECT_DOUBLE_EQ(sub->ps[1].point.y, -8.12);
  EXPECT_DOUBLE_EQ(sub->ps[1].weight, 0.68);
  
  delete sub;
}

TEST(ParticleArray, FilterCDFRemoveOne) {
  ParticleArray<> pa_orig;
  pa_orig.resize(3);
  pa_orig.ps[0].point.x = 1.0;
  pa_orig.ps[0].point.y = -2.10;
  pa_orig.ps[0].weight = 0.56;
  
  pa_orig.ps[1].point.x = -1.0;
  pa_orig.ps[1].point.y = 2.10;
  pa_orig.ps[1].weight = 0.01;

  pa_orig.ps[2].point.x = 1.0;
  pa_orig.ps[2].point.y = -2.09;
  pa_orig.ps[2].weight = 0.43;

  ParticleArray<> *sub = pa_orig.FilterCDF(0.1);
  ASSERT_TRUE(sub != NULL);
  ASSERT_EQ(2, sub->size);
  EXPECT_DOUBLE_EQ(sub->ps[0].point.x, 1.0);
  EXPECT_DOUBLE_EQ(sub->ps[0].point.y, -2.09);
  EXPECT_DOUBLE_EQ(sub->ps[0].weight, 0.43 / (0.43 + 0.56));
  EXPECT_DOUBLE_EQ(sub->ps[1].point.x, 1.0);
  EXPECT_DOUBLE_EQ(sub->ps[1].point.y, -2.10);
  EXPECT_DOUBLE_EQ(sub->ps[1].weight, 0.56 / (0.43 + 0.56));
  
  delete sub;
}

TEST(STAgg, Time1Bot2Msmts) {
  double ts1 = 12.0, ts3 = 17.1;
  Pose2D start(0.0, 0.0, 1.5), end(50.0, 50.0, 0.0);
  Measurement m1(ts1, 0, 1, start, 1.0);
  Measurement m2(ts3, 0, 1, end, 6.0);

  STAggregator st(5.0, 50);

  EXPECT_FALSE(st.addMeasurement(m1));
  EXPECT_TRUE(st.addMeasurement(m2));
  const Measurement *agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), 0);
  EXPECT_EQ(agg->sourceID(), 1);
  EXPECT_DOUBLE_EQ(m1.dist(), agg->dist());
}


TEST(STAgg, Time1Bot3Msmts) {
  double ts1 = 12.0, ts2 = 13.2, ts3 = 17.1;
  Pose2D start(0.0, 0.0, 1.5), end(50.0, 50.0, 0.0);
  Measurement m1(ts1, 0, 1, start, 1.0);
  Measurement m2(ts2, 0, 1, end, 2.0);
  Measurement m3(ts3, 0, 1, end, 3.0);

  STAggregator st(5.0, 1.0);

  EXPECT_FALSE(st.addMeasurement(m1));
  EXPECT_FALSE(st.addMeasurement(m2));
  EXPECT_TRUE(st.addMeasurement(m3));
  const Measurement *agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), 0);
  EXPECT_EQ(agg->sourceID(), 1);
  EXPECT_DOUBLE_EQ((m1.dist() + m2.dist()) / 2.0, agg->dist());

  const std::vector<Measurement*> &msmts = st.getOrig();
  ASSERT_EQ(2, msmts.size());
  EXPECT_EQ(m1.timestamp(), msmts[0]->timestamp());
  EXPECT_EQ(m2.timestamp(), msmts[1]->timestamp());
}

TEST(STAgg, Dist1Bot2Msmts) {
  Pose2D p1(1.0, 2.3, 0.5);
  Pose2D p2(1.0, 2.4, -3.0);
  
  Measurement m1(0.0, 1, 0, p1, 1.0);
  Measurement m2(0.001, 1, 0, p2, -5.0); 

  STAggregator st(0.0, 0.01);

  EXPECT_FALSE(st.addMeasurement(m1));
  EXPECT_TRUE(st.addMeasurement(m2));  
  const Measurement *agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), 1);
  EXPECT_EQ(agg->sourceID(), 0);
  EXPECT_DOUBLE_EQ(m1.dist(), agg->dist());
}

TEST(STAgg, Dist1Bot3Msmts) {
  Pose2D p1(1.0, 2.3, 0.5);
  Pose2D p2(1.0, 2.4, -3.0);
  Pose2D p3(12.0, 3.0, 3.0);
  
  Measurement m1(0.0, 1, 0, p1, 1.0);
  Measurement m2(1.0, 1, 0, p2, -5.0);
  Measurement m3(2.0, 1, 0, p3, 6.0);  

  STAggregator st(0.5, 10.0);

  EXPECT_FALSE(st.addMeasurement(m1));
  EXPECT_FALSE(st.addMeasurement(m2));  
  EXPECT_TRUE(st.addMeasurement(m3));
  const Measurement *agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), 1);
  EXPECT_EQ(agg->sourceID(), 0);
  EXPECT_DOUBLE_EQ((m1.dist() + m2.dist()) / 2.0, agg->dist());
}

TEST(STAgg, OutOfOrderMeasurement) {
  Pose2D p1(2.0, 2.3, 0.5);
  Pose2D p2(1.0, -50, 0.0);
  
  Measurement m1(10.0, 1, 0, p1, 1.0);
  Measurement m2(9.0, 1, 0, p2, -5.0);
  Measurement m3(15.0, 1, 0, p2, 6.0);  

  STAggregator st(0.5, 10.0);

  EXPECT_FALSE(st.addMeasurement(m1));
  EXPECT_FALSE(st.addMeasurement(m2));  
  EXPECT_TRUE(st.addMeasurement(m3));
  const Measurement *agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), 1);
  EXPECT_EQ(agg->sourceID(), 0);
  EXPECT_DOUBLE_EQ(m1.dist(), agg->dist());
}

TEST(STAgg, MultipleSources) {
  int src1 = 1;
  int src2 = 2;
  int recv = 3;
  
  Pose2D start(0, 0, 0), end(10, 10, 0);
  
  Measurement m1(10.0, src1, recv, start, 1.2);
  Measurement m2(12.0, src2, recv, end, -5.0);
  Measurement m3(15.0, src1, recv, end, 6.0);
  Measurement m4(18.0, src2, recv, start, 0.5);  

  STAggregator st(0.5, 10.0);

  EXPECT_FALSE(st.addMeasurement(m1));
  EXPECT_FALSE(st.addMeasurement(m2));
  
  EXPECT_TRUE(st.addMeasurement(m3));
  const Measurement *agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), src1);
  EXPECT_EQ(agg->sourceID(), recv);
  EXPECT_DOUBLE_EQ(m1.dist(), agg->dist());

  EXPECT_TRUE(st.addMeasurement(m4));
  agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), src2);
  EXPECT_EQ(agg->sourceID(), recv);
  EXPECT_DOUBLE_EQ(m2.dist(), agg->dist());
}

TEST(STAgg, MultipleReceivers) {
  int src = 1;
  int recv1 = 3;
  int recv2 = 5;
  
  Pose2D start(0, 0, 0), end(10, 10, 0);
  
  Measurement m1(10.0, src, recv1, start, 1.2);
  Measurement m2(12.0, src, recv2, end, -5.0);
  Measurement m3(14.0, src, recv1, end, 6.0);
  Measurement m5(16.0, src, recv1, end, 0.1);
  Measurement m4(18.0, src, recv2, start, 0.5);

  STAggregator st(5.0, 10.0);

  EXPECT_FALSE(st.addMeasurement(m1));
  EXPECT_FALSE(st.addMeasurement(m2));
  EXPECT_FALSE(st.addMeasurement(m3));
  
  EXPECT_TRUE(st.addMeasurement(m4));  
  const Measurement *agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), src);
  EXPECT_EQ(agg->sourceID(), recv2);
  EXPECT_DOUBLE_EQ(m2.dist(), agg->dist());

  EXPECT_FALSE(st.addMeasurement(m4));  
  EXPECT_TRUE(st.addMeasurement(m5));
  agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), src);
  EXPECT_EQ(agg->sourceID(), recv1);
  EXPECT_DOUBLE_EQ((m1.dist() + m3.dist()) / 2.0, agg->dist());
}

TEST(STAgg, RemoveOldEntries) {
  double ts1 = 12.0, ts2 = 73.0, ts3 = 76.0;
  double max_age = 60.0;
  Pose2D start(0.0, 0.0, 1.5), end(50.0, 50.0, 0.0);
  Measurement m1(ts1, 0, 1, start, 1.0);
  Measurement m2(ts2, 0, 1, end, 6.0);
  Measurement m3(ts3, 0, 1, start, 6.0);
  
  STAggregator st(2.0, 2.0, max_age);

  EXPECT_FALSE(st.addMeasurement(m1));
  EXPECT_FALSE(st.addMeasurement(m2));
  EXPECT_TRUE(st.addMeasurement(m3));  
  const Measurement *agg = st.getAggregate();
  ASSERT_TRUE(agg != NULL);
  EXPECT_EQ(agg->destID(), 0);
  EXPECT_EQ(agg->sourceID(), 1);
  EXPECT_DOUBLE_EQ(m2.dist(), agg->dist());
}


