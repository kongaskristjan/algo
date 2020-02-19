
#include <gtest/gtest.h>
#include "algo.h"

using namespace std;
using namespace algo;

TEST(arithmetic, most_significant_bit_test) {
    EXPECT_EQ(-1, most_significant_bit(0));
    EXPECT_EQ(0, most_significant_bit(1));
    EXPECT_EQ(1, most_significant_bit(2));
    EXPECT_EQ(1, most_significant_bit(3));
    EXPECT_EQ(2, most_significant_bit(4));
    EXPECT_EQ(2, most_significant_bit(5));
    EXPECT_EQ(2, most_significant_bit(6));
    EXPECT_EQ(2, most_significant_bit(7));
}

TEST(range_minimum_query_test, test) {
    vector<int> v({-2, 5, 1, 2, 9, 4, 3});
    range_minimum_query<int> r(v);
    for (int i = 0; i < v.size(); ++i) {
        int mn = 100;
        for (int j = i; j < v.size(); ++j) {
            mn = min(mn, v[j]);
            EXPECT_EQ(mn, r.query(i, j + 1));
        }
    }
}
