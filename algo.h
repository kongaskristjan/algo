#ifndef __ALGO_H__
#define __ALGO_H__

#include <cassert>
#include <algorithm>
#include <vector>

namespace algo {
    template <typename T>
    constexpr int most_significant_bit(T value) {
        assert(value >= 0);

        int msb = -1;
        while(value) {
            value >>= 1;
            ++msb;
        }
        return msb;
    }

    template <typename T>
    class range_minimum_query {
    private:
        std::vector<std::vector<T>> sparse_; // sparse_[log(query size)][start position]

    public:
        explicit range_minimum_query(const std::vector<T> &v) {
            assert(v.size() > 0);
            size_t levels = 1 + most_significant_bit(v.size() - 1);
            for(size_t i = v.size() - 1; i; i /= 2) ++levels;

            sparse_.resize(levels, std::vector<T>(v.size()));
            sparse_[0] = v;
            size_t power_i = 1;
            for(size_t i = 0; i < levels - 1; ++i) {
                for(size_t j = 0; j < sparse_[i].size(); ++j) {
                    if(j + power_i < sparse_[i].size()) {
                        size_t second_index = std::min(sparse_[i].size() - 1, j + power_i);
                        sparse_[i + 1][j] = std::min(sparse_[i][j], sparse_[i][second_index]);
                    }
                    else
                        sparse_[i + 1][j] = sparse_[i][j];
                }
                power_i *= 2;
            }
        }

        T query(size_t begin, size_t end) {
            assert(end - begin > 0);
            if(end - begin == 1) return sparse_[0][begin];

            size_t level = most_significant_bit(end - begin - 1);
            size_t second_index = std::min(sparse_[level].size() - 1, end - (size_t(1) << level));

            return std::min(sparse_[level][begin], sparse_[level][second_index]);
        }
    };
};

#endif
