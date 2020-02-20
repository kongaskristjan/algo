#ifndef __ALGO_H__
#define __ALGO_H__

#include <cassert>
#include <algorithm>
#include <vector>
#include <array>

// DECLARATIONS

namespace algo {
    template<typename T>
    constexpr int most_significant_bit(T value);

    template<typename T>
    class range_minimum_query {
    private:
        std::vector<std::vector<T>> sparse_; // sparse_[log(query size)][start position]

        static auto query_if_appropriate_(const T &sub_rmq);

        template<typename... Args>
        static auto query_if_appropriate_(const T &sub_rmq, Args... sub_range);

    public:
        range_minimum_query() = default;

        template<typename VecT>
        explicit range_minimum_query(const VecT &v);

        template<typename... Args>
        auto query(size_t begin, size_t end, Args... sub_ranges) const;

        size_t size() const { return sparse_[0].size(); }

        static range_minimum_query<T> min(const range_minimum_query<T> &a, const range_minimum_query<T> &b);
    };
};

template<typename T>
algo::range_minimum_query<T> std::min(const algo::range_minimum_query<T> &a, const algo::range_minimum_query<T> &b);

// DEFINITIONS

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
    auto range_minimum_query<T>::query_if_appropriate_(const T &sub_rmq) {
        return sub_rmq;
    }

    template <typename T> template <typename... Args>
    auto range_minimum_query<T>::query_if_appropriate_(const T &sub_rmq, Args... sub_range) {
        return sub_rmq.query(sub_range...);
    }

    template <typename T> template <typename VecT>
    range_minimum_query<T>::range_minimum_query(const VecT &v) {
        assert(v.size() > 0);
        size_t levels = 1 + most_significant_bit(v.size() - 1);
        for(size_t i = v.size() - 1; i; i /= 2) ++levels;

        sparse_.resize(levels, std::vector<T>(v.size()));

        // If T == range_minimum_query: convert vector<vector<>> v to range_minimum_query<range_minimum_query<>>
        // If T == numeric type: copy
        for(size_t i = 0; i < v.size(); ++i)
            sparse_[0][i] = T(v[i]);

        size_t power_i = 1;
        for(size_t i = 0; i < levels - 1; ++i) {
            for(size_t j = 0; j < sparse_[i].size(); ++j) {
                size_t second_index = j + power_i;
                if(second_index < sparse_[i].size())
                    sparse_[i + 1][j] = std::min(sparse_[i][j], sparse_[i][second_index]);
                else
                    sparse_[i + 1][j] = sparse_[i][j];
            }
            power_i *= 2;
        }
    }

    template <typename T> template <typename... Args>
    auto range_minimum_query<T>::query(size_t begin, size_t end, Args... sub_range) const {
        assert(end - begin > 0);

        if(end - begin == 1) return query_if_appropriate_(sparse_[0][begin], sub_range...);

        size_t level = most_significant_bit(end - begin - 1);
        size_t second_index = std::min(sparse_[level].size() - 1, end - (size_t(1) << level));

        auto query1 = query_if_appropriate_(sparse_[level][begin], sub_range...);
        auto query2 = query_if_appropriate_(sparse_[level][second_index], sub_range...);
        return std::min(query1, query2);
    }

    template <typename T>
    range_minimum_query<T> range_minimum_query<T>::min(
            const range_minimum_query<T> &a, const range_minimum_query<T> &b) {
        assert(a.size() == b.size());
        range_minimum_query<T> ret;
        ret.sparse_.resize(a.sparse_.size(), std::vector<T>(a.sparse_[0].size()));
        for(size_t i = 0; i < ret.sparse_.size(); ++i)
            for (size_t j = 0; j < ret.sparse_[i].size(); ++j)
                ret.sparse_[i][j] = std::min(a.sparse_[i][j], b.sparse_[i][j]);

        return ret;
    }
};

template<typename T>
algo::range_minimum_query<T> std::min(const algo::range_minimum_query<T> &a, const algo::range_minimum_query<T> &b) {
    return algo::range_minimum_query<T>::min(a, b);
}

#endif
