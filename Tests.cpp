
#include <gtest/gtest.h>
#include "algo.h"

using namespace std;
using namespace algo;

// ARITHMETICS

TEST(arithmetic_test, most_significant_bit_test) {
    EXPECT_EQ(-1, most_significant_bit(0));
    EXPECT_EQ(0, most_significant_bit(1));
    EXPECT_EQ(1, most_significant_bit(2));
    EXPECT_EQ(1, most_significant_bit(3));
    EXPECT_EQ(2, most_significant_bit(4));
    EXPECT_EQ(2, most_significant_bit(5));
    EXPECT_EQ(2, most_significant_bit(6));
    EXPECT_EQ(2, most_significant_bit(7));
}

TEST(modular_test, basic_arithm) {
    using m10 = modular<int, 10>;
    using u10 = modular<unsigned, 10>;

    EXPECT_EQ(m10(-17), 3);
    EXPECT_EQ(m10(13), 3);
    EXPECT_EQ(m10(3) * 4, 2);
    EXPECT_EQ(m10(3) - 4, 9);

    EXPECT_EQ(m10(-17).get(), 3);
    EXPECT_EQ(m10(103), 3);
    EXPECT_NE(m10(5), 3);
    EXPECT_EQ(m10(4), -6);

    EXPECT_EQ(m10(5) * 5, 5);
    EXPECT_EQ(m10(2) + 2, 4);
    EXPECT_EQ(m10(2) + 8, 0);
    EXPECT_EQ(m10(2) - 2, 0);
    EXPECT_EQ(m10(2) - 8, 4);
    EXPECT_EQ(-m10(3), 7);
}

TEST(modular_test, division) {
    using m13 = modular<int, 13>;
    using u13 = modular<unsigned, 13>;

    for(int i = 1; i < 13; ++i) {
        EXPECT_EQ(m13(i).inverse() * i, 1);
        EXPECT_EQ(u13(i).inverse() * i, 1);
    }

    for(int i = 0; i < 13; ++i)
        for(int j = 1; j < 13; ++j) {
            EXPECT_EQ(m13(i) / j * j, i);
            EXPECT_EQ(u13(i) / j * j, i);
        }
}

TEST(modular_test, power) {
    using m13 = modular<int, 13>;
    using u13 = modular<unsigned, 13>;

    for(int i = 1; i < 13; ++i) {
        m13 m(1);
        u13 u(1);
        for (int j = 0; j <= 50; ++j) {
            ASSERT_EQ(m.get(), u.get());
            ASSERT_EQ(pow(m13(i), j), m);
            ASSERT_EQ(pow(u13(i), j), u);
            ASSERT_EQ(pow(m13(i), -j), m.inverse());
            ASSERT_EQ(pow(u13(i), -j), u.inverse());
            m *= i;
            u *= i;
        }
    }
}

TEST(modular_test, stream_test) {
    using m13 = modular<int, 13>;

    std::stringstream ss;
    m13 a(4), b;
    ss << a;
    ss >> b;
    ASSERT_EQ(a, b);
}

// TABULAR

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
