#include "../include/back_sub.hpp"

#include <iostream>
#include <gtest/gtest.h>

TEST(BackSubTest, BasicAssertions) {
    std::vector<float> u = {1.0,5.0,0.2};
    std::vector<float> b = {1.0,8.0};
    std::vector<float> y = {6.0, -3.0, 7.8};
    std::vector<float> x = backward_substitution(u, b, y);


    EXPECT_NEAR(x[0], 69.0f, 1e-6);
    EXPECT_NEAR(x[1],-63.0f, 1e-6);
    EXPECT_NEAR(x[2], 39.0f, 1e-6);
}