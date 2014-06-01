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

#include "planner.hpp"

//============================== Configuration ==============================//

TEST(Configuration, InitialEquality) {
    RobotConfiguration c1(1), c2(1);
    ASSERT_TRUE(c1.samePosition(c2));
    ASSERT_TRUE(c2.samePosition(c1));
}

TEST(Configuration, BasicInequality) {
    RobotConfiguration c1(1), c2(1);
    c1.poses[0].setXYT(1, 0, 0);
    ASSERT_FALSE(c1.samePosition(c2));
    ASSERT_FALSE(c2.samePosition(c1));
}

TEST(Configuration, MixedBots) {
    RobotConfiguration c1(1), c2(2);
    ASSERT_FALSE(c1.samePosition(c2));
    ASSERT_FALSE(c2.samePosition(c1));
}

TEST(Configuration, MultipleBotsEqual) {
    RobotConfiguration c1(2), c2(2);
    c1.poses[0].setXYT(4, -2, 3);
    c2.poses[0].setXYT(4, -2, 3);
    c1.poses[1].setXYT(0, 42, 0);
    c2.poses[1].setXYT(0, 42, 0);

    ASSERT_TRUE(c1.samePosition(c2));
    ASSERT_TRUE(c2.samePosition(c1));
}

TEST(Configuration, CopyConstructor) {
    RobotConfiguration expect(5);
    for (int i = 0; i < 5; ++i) {
        expect.poses[i].setXYT(i * 2, i * 3 - 1, i / 5.0);
    }
    RobotConfiguration actual(expect);
    ASSERT_EQ(expect.nbots, actual.nbots);
    for (int i = 0; i < 5; ++i) {
        ASSERT_EQ(expect.poses[i].x(), actual.poses[i].x());
        ASSERT_EQ(expect.poses[i].y(), actual.poses[i].y());
        ASSERT_EQ(expect.poses[i].theta(), actual.poses[i].theta());
    }
}

//================================ Generator ================================//


class KinematicTest : public testing::Test {
protected:
    void SetUp() {
        map_ = map_alloc();
        map_->scale = 1.0;
        int width = 20, height = 20;
        map_->size_x = width;
        map_->size_y = height;
        map_->cells = (map_cell_t*) calloc(width * height, sizeof(map_->cells[0]));
        for (int i = 0; i < width * height; ++i) {
            map_->cells[i].occ_state = -1;
        }
    }

    void TearDown() {
        map_free(map_);
    }
    
    // Check that each element of expect[] can be found in actual.
    void check_configuration(std::vector<RobotConfiguration *> expect,
                             std::vector<RobotConfiguration *> actual) {
        ASSERT_EQ(expect.size(), actual.size());
        for (int e = 0; e < expect.size(); ++e) {
            bool contains = false;
            for (int a = 0; a < actual.size(); ++a) {
                contains |= expect[e]->samePosition(*actual[a]);
            }
            ASSERT_TRUE(contains);
        }
    }
    map_t *map_;
};

TEST_F(KinematicTest, SingleBot) {
    RobotConfiguration r(1);
    r.poses[0].setXYT(-3, 1, 0);

    std::vector<RobotConfiguration *> expect;
    for (int i = 0; i < 4; ++i) {
        expect.push_back(new RobotConfiguration(1));
    }
    expect[0]->poses[0].setXYT(-4, 1, 0);
    expect[1]->poses[0].setXYT(-2, 1, 0);
    expect[2]->poses[0].setXYT(-3, 2, 0);
    expect[3]->poses[0].setXYT(-3, 2, 0);

    KinematicPlanner gen;
    gen.setStepSize(1.0);
    std::vector<RobotConfiguration *> actual = gen.generate(r);

    // Check that each element of expect[] can be found in actual.
    check_configuration(expect, actual);

    gen.freeConfigurations(&expect);
    gen.freeConfigurations(&actual);
}

// Place two bots one square away from each other and insure that the
// configuration where they collide is not generated, but all others are
TEST_F(KinematicTest, TwoBotsExhaustive) {
    RobotConfiguration r(2);
    double x1 = 3, y1 = 2;
    double x2 = 4, y2 = 2;
    r.poses[0].x(x1);
    r.poses[0].y(y1);
    r.poses[1].x(x2);
    r.poses[1].y(y2);

    // Initialize all values of expect to be the starting position
    std::vector<RobotConfiguration *> expect;
    for (int i = 0; i < 15; ++i) {
        expect.push_back(new RobotConfiguration(r));
    }
    int ind = 0;
    // Case where bot 1 moves right and bot 2 moves left cannot happen
    expect[ind]->poses[0].x(x1 + 1);
    expect[ind]->poses[1].x(x2 + 1);

    expect[++ind]->poses[0].x(x1 + 1);
    expect[ind]->poses[1].y(y2 + 1);

    expect[++ind]->poses[0].x(x1 + 1);
    expect[ind]->poses[1].y(y2 - 1);
    //
    expect[++ind]->poses[0].x(x1 - 1);
    expect[ind]->poses[1].x(x2 + 1);

    expect[++ind]->poses[0].x(x1 - 1);
    expect[ind]->poses[1].y(y2 + 1);

    expect[++ind]->poses[0].x(x1 - 1);
    expect[ind]->poses[1].x(x2 - 1);

    expect[++ind]->poses[0].x(x1 - 1);
    expect[ind]->poses[1].y(y2 - 1);
    //
    expect[++ind]->poses[0].y(y1 + 1);
    expect[ind]->poses[1].x(x2 + 1);

    expect[++ind]->poses[0].y(y1 + 1);
    expect[ind]->poses[1].y(y2 + 1);

    expect[++ind]->poses[0].y(y1 + 1);
    expect[ind]->poses[1].x(x2 - 1);

    expect[++ind]->poses[0].y(y1 + 1);
    expect[ind]->poses[1].y(y2 - 1);
    //
    expect[++ind]->poses[0].y(y1 - 1);
    expect[ind]->poses[1].x(x2 + 1);

    expect[++ind]->poses[0].y(y1 - 1);
    expect[ind]->poses[1].y(y2 + 1);

    expect[++ind]->poses[0].y(y1 - 1);
    expect[ind]->poses[1].x(x2 - 1);

    expect[++ind]->poses[0].y(y1 - 1);
    expect[ind]->poses[1].y(y2 - 1);
    ASSERT_EQ(ind, 14);
    KinematicPlanner gen;
    gen.setStepSize(1);
    std::vector<RobotConfiguration *> actual = gen.generate(r);

    check_configuration(expect, actual);

    gen.freeConfigurations(&expect);
    gen.freeConfigurations(&actual);
}

// For larger number of bots, only check number of returned configurations
TEST_F(KinematicTest, FourBotsNoOverlap) {
    RobotConfiguration r(4);
    for (int i = 0; i < 4; ++i) {
        r.poses[i].setXYT(4 * i, 0, 0);
    }
    KinematicPlanner gen;
    std::vector<RobotConfiguration *> actual = gen.generate(r);

    ASSERT_EQ(4 * 4 * 4 * 4, actual.size());
    gen.freeConfigurations(&actual);
}

// If bot 1 moves up and bot 2 moves down, they will collide
TEST_F(KinematicTest, TwoBotsCollision) {
    RobotConfiguration r(2);
    r.poses[0].setXYT(4, 1, 0);
    r.poses[1].setXYT(4, 3, 0);    

    KinematicPlanner gen;
    gen.setStepSize(1);
    std::vector<RobotConfiguration *> actual = gen.generate(r);

    ASSERT_EQ(15, actual.size());
    gen.freeConfigurations(&actual);
}

TEST_F(KinematicTest, TwoBotOccupancy) {
    // Configure robots so that robot 1 cannot move down when robot 2 moves up
    RobotConfiguration start(2);
    double x1 = 3, y1 = 2;
    double x2 = 3, y2 = 1;
    start.poses[0].setXYT(x1, y1, 0);
    start.poses[1].setXYT(x2, y2, 0);

    // Robot 2 cannot move down
    map_get_cell(map_, x2, y2 - 1, 0)->occ_state = 1; 
    // Robot 2 cannot move left
    map_get_cell(map_, x2 - 1, y2, 0)->occ_state = 1; 
    // Robot 1 cannot move left
    map_get_cell(map_, x1 - 1, y1, 0)->occ_state = 1; 

    std::vector<RobotConfiguration *> expected;
    for (int i = 0; i < 5; ++i) {
        expected.push_back(new RobotConfiguration(start));
    }

    int ind = 0;
    expected[ind]->poses[0].x(x1 + 1);
    expected[ind]->poses[1].x(x2 + 1);
    ind++;
    expected[ind]->poses[0].x(x1 + 1);
    expected[ind]->poses[1].y(y2 + 1);
    ind++;
    expected[ind]->poses[0].y(y1 + 1);
    expected[ind]->poses[1].x(x2 + 1);
    ind++;
    expected[ind]->poses[0].y(y1 + 1);
    expected[ind]->poses[1].y(y2 + 1);
    ind++;
    expected[ind]->poses[0].y(y1 - 1);
    expected[ind]->poses[1].x(x2 + 1);


    KinematicPlanner gen(map_);
    gen.setStepSize(1.0);
    std::vector<RobotConfiguration *> actual = gen.generate(start);
    check_configuration(expected, actual);

    gen.freeConfigurations(&expected);
    gen.freeConfigurations(&actual);
}

TEST_F(KinematicTest, OneBotOccupancy) {
    RobotConfiguration start(1);
    int x1 = -3, y1 = 5;
    start.poses[0].x(x1);
    start.poses[0].y(y1);

    // Configure map so that robot cannot move down or left
    map_get_cell(map_, x1, y1 - 1, 0)->occ_state = 1;     
    map_get_cell(map_, x1 - 1, y1, 0)->occ_state = 1;     
    
    std::vector<RobotConfiguration *> expected;
    for (int i = 0; i < 2; ++i) {
        expected.push_back(new RobotConfiguration(start));
    }
    int ind = 0;
    expected[ind]->poses[0].x(x1 + 1);
    ind++;
    expected[ind]->poses[0].y(y1 + 1);

    KinematicPlanner gen(map_);
    gen.setStepSize(1.0);
    std::vector<RobotConfiguration *> actual = gen.generate(start);
    check_configuration(expected, actual);

    gen.freeConfigurations(&expected);
    gen.freeConfigurations(&actual);
}

class DynamicTest : public testing::Test {
protected:
    
};

TEST_F(DynamicTest, SimpleEquality) {
    DynamicConfiguration c1(4);
    RobotConfiguration r(4);

    c1.poses[0].x(4);
    r.poses[0].x(4);
    c1.poses[2].y(-2);
    r.poses[2].y(-2);

    c1.velocities[0].v = 5;
    
    ASSERT_TRUE(c1.samePosition(r));
    ASSERT_TRUE(r.samePosition(c1));
}
