#include "../include/forward_sub.hpp"

#include <iostream>
#include <gtest/gtest.h>

TEST(ForwardSubTest, BasicAssertions) {
    std::vector<float> l = {2.0, 0.6};
    std::vector<float> b = {6.0, 9.0, 6.0};
    std::vector<float> y = forward_substitution(l, b);


    EXPECT_NEAR(y[0], 6.0f, 1e-6);
    EXPECT_NEAR(y[1],-3.0f, 1e-6);
    EXPECT_NEAR(y[2], 7.8f, 1e-6);
}