
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

TEST(sparse_table_test, test) {
    vector<int> v({-2, 5, 1, 2, 9, 4, 3});
    sparse_table<int> r(v);
    for (size_t i = 0; i < v.size(); ++i) {
        int mn = 100;
        for (size_t j = i; j < v.size(); ++j) {
            mn = min(mn, v[j]);
            EXPECT_EQ(mn, r.query(i, j + 1));
        }
    }
}

TEST(sparse_table_test, test2d) {
    vector<vector<int>> v({
        {-2, 5, 1, 2, 7, 2},
        {2, 9, 4, 4, 5, 1},
        {3, 0, 2, -3, -2, -1},
        {1, 2, -2, -1, 0, -1}
    });
    sparse_table<sparse_table<int>> r(v);
    for (int i = 0; i < v.size(); ++i)
        for (int j = 0; j < v[i].size(); ++j)
            for (int iMax = i; iMax < v.size(); ++iMax)
                for (int jMax = j; jMax < v.size(); ++jMax) {
                    int mn = 100;
                    for(int iIdx = i; iIdx <= iMax; ++iIdx)
                        for(int jIdx = j; jIdx <= jMax; ++jIdx)
                            mn = min(mn, v[iIdx][jIdx]);
                    EXPECT_EQ(mn, r.query(i, iMax + 1, j, jMax + 1));
                }
}
