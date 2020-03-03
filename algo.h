#ifndef __ALGO_H__
#define __ALGO_H__

#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>
#include <array>
#include <iostream>
#include <type_traits>

// DECLARATIONS

namespace algo {
    // ARITHMETIC

    template<typename T>
    constexpr int most_significant_bit(T value);

    template<typename T>
    T gcd(T a, T b, typename std::make_signed<T>::type &x, typename std::make_signed<T>::type &y);

    // std::make_unsigned<T>::type must support at least Mod^2 sized integers
    // Sign of T only matters for the interface, computations are mostly done with unsigned integers anyway
    template<typename T, typename std::make_unsigned<T>::type Mod>
    class modular {
        using S = typename std::make_signed<T>::type;
        using U = typename std::make_unsigned<T>::type;
    private:
        U value_ = 0; // 0 <= value_ < Mod

    public:
        modular& operator=(T value__);
        modular() = default;
        modular(T value__) { *this = value__; }

        bool operator==(const modular<T, Mod> &other) const { return value_ == other.value_; }
        bool operator!=(const modular<T, Mod> &other) const { return value_ != other.value_; }

        modular<T, Mod> inverse() const;
        modular<T, Mod>& operator+=(const modular<T, Mod> &other);
        modular<T, Mod>& operator-=(const modular<T, Mod> &other);
        modular<T, Mod>& operator*=(const modular<T, Mod> &other);
        modular<T, Mod>& operator/=(const modular<T, Mod> &other);
        modular<T, Mod> operator+(const modular<T, Mod> &rhs) const { modular<T, Mod> ret(*this); return ret += rhs; }
        modular<T, Mod> operator-(const modular<T, Mod> &rhs) const { modular<T, Mod> ret(*this); return ret -= rhs; }
        modular<T, Mod> operator*(const modular<T, Mod> &rhs) const { modular<T, Mod> ret(*this); return ret *= rhs; }
        modular<T, Mod> operator/(const modular<T, Mod> &rhs) const { modular<T, Mod> ret(*this); return ret /= rhs; }
        modular<T, Mod> operator-() const;

        T get() const { return value_; }
    };

    template<typename T, typename std::make_unsigned<T>::type Mod>
    std::ostream& operator<<(std::ostream &os, const modular<T, Mod> &md) { return os << md.get(); }
    template<typename T, typename std::make_unsigned<T>::type Mod>
    std::istream& operator>>(std::istream &is, modular<T, Mod> &md) { T x; is >> x; md = x; return is; }

    // TABULAR

    template<typename T>
    class sparse_table {
    private:
        std::vector<std::vector<T>> sparse_; // sparse_[log(query size)][start position]

        static auto query_if_appropriate_(const T &sub_rmq);

        template<typename... Args>
        static auto query_if_appropriate_(const T &sub_rmq, Args... sub_range);

    public:
        sparse_table() = default;

        template<typename VecT>
        explicit sparse_table(const VecT &v);

        template<typename... Args>
        auto query(size_t begin, size_t end, Args... sub_ranges) const;

        size_t size() const { return sparse_[0].size(); }

        static sparse_table<T> min(const sparse_table<T> &a, const sparse_table<T> &b);
    };

    template<typename T>
    class segment_tree {
    private:
        // linear representation of binary tree from index 1, children of N are 2N, 2N+1
        // elements are in leaf nodes at SIZE, SIZE+1, ..., 2*SIZE-1
        std::vector<T> tree_;
        size_t size_;

    public:
        segment_tree() = default;
        explicit segment_tree(const std::vector<T> &v);

        T query(size_t begin, size_t end) const;
        void update(size_t idx, T value);
        size_t size() const { return size_; }
    };
};

template<typename T, typename std::make_unsigned<T>::type Mod>
algo::modular<T, Mod> pow(algo::modular<T, Mod> x, typename std::make_signed<T>::type n);

template<typename T>
algo::sparse_table<T> std::min(const algo::sparse_table<T> &a, const algo::sparse_table<T> &b);

// DEFINITIONS

namespace algo {
    // ARITHMETIC

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

    // a >= b
    template<typename T>
    T gcd(T a, T b, typename std::make_signed<T>::type &x, typename std::make_signed<T>::type &y) {
        if(b == 0) {
            x = 1;
            y = 0;
            return a;
        }
        typename std::make_signed<T>::type x1, y1;
        T g = gcd(b, a % b, x1, y1);
        x = y1;
        y = x1 - a / b * y1;
        return g;
    }

    template<typename T, typename std::make_unsigned<T>::type Mod>
    modular<T, Mod>& modular<T, Mod>::operator=(T value__) {
        value__ = value__ % (T) Mod;
        if(value__ < 0) value__ += Mod;
        value_ = value__;
        return *this;
    }

    template<typename T, typename std::make_unsigned<T>::type Mod>
    modular<T, Mod> modular<T, Mod>::inverse() const {
        assert(value_ != 0);
        S inv, _;
        gcd(Mod, value_, _, inv);
        return modular<T, Mod>(inv + Mod);
    }

    template<typename T, typename std::make_unsigned<T>::type Mod>
    modular<T, Mod>& modular<T, Mod>::operator+=(const modular<T, Mod> &other) {
        value_ += other.value_;
        if(value_ >= Mod)
            value_ -= Mod;
        return *this;
    }

    template<typename T, typename std::make_unsigned<T>::type Mod>
    modular<T, Mod>& modular<T, Mod>::operator-=(const modular<T, Mod> &other) {
        value_ -= other.value_;
        if(value_ >= Mod) // overflow
            value_ += Mod;
        return *this;
    }

    template<typename T, typename std::make_unsigned<T>::type Mod>
    modular<T, Mod>& modular<T, Mod>::operator*=(const modular<T, Mod> &other) {
        value_ *= other.value_;
        value_ %= Mod;
        return *this;
    }

    template<typename T, typename std::make_unsigned<T>::type Mod>
    modular<T, Mod>& modular<T, Mod>::operator/=(const modular<T, Mod> &other) {
        modular<T, Mod> inv = other.inverse();
        *this *= inv;
        return *this;
    }

    template<typename T, typename std::make_unsigned<T>::type Mod>
    modular<T, Mod> modular<T, Mod>::operator-() const {
        modular<T, Mod> ret;
        ret.value_ = Mod - value_;
        return ret;
    }

    // TABULAR

    template <typename T>
    auto sparse_table<T>::query_if_appropriate_(const T &sub_rmq) {
        return sub_rmq;
    }

    template <typename T> template <typename... Args>
    auto sparse_table<T>::query_if_appropriate_(const T &sub_rmq, Args... sub_range) {
        return sub_rmq.query(sub_range...);
    }

    template <typename T> template <typename VecT>
    sparse_table<T>::sparse_table(const VecT &v) {
        assert(v.size() > 0);
        size_t levels = 1 + most_significant_bit(v.size() - 1);
        for(size_t i = v.size() - 1; i; i /= 2) ++levels;

        sparse_.resize(levels, std::vector<T>(v.size()));

        // If T == sparse_table: convert vector<vector<>> v to sparse_table<sparse_table<>>
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
    auto sparse_table<T>::query(size_t begin, size_t end, Args... sub_range) const {
        assert(begin < end);

        if(end - begin == 1) return query_if_appropriate_(sparse_[0][begin], sub_range...);

        size_t level = most_significant_bit(end - begin - 1);
        size_t second_index = std::min(sparse_[level].size() - 1, end - (size_t(1) << level));

        auto query1 = query_if_appropriate_(sparse_[level][begin], sub_range...);
        auto query2 = query_if_appropriate_(sparse_[level][second_index], sub_range...);
        return std::min(query1, query2);
    }

    template <typename T>
    sparse_table<T> sparse_table<T>::min(
            const sparse_table<T> &a, const sparse_table<T> &b) {
        assert(a.size() == b.size());
        sparse_table<T> ret;
        ret.sparse_.resize(a.sparse_.size(), std::vector<T>(a.sparse_[0].size()));
        for(size_t i = 0; i < ret.sparse_.size(); ++i)
            for (size_t j = 0; j < ret.sparse_[i].size(); ++j)
                ret.sparse_[i][j] = std::min(a.sparse_[i][j], b.sparse_[i][j]);

        return ret;
    }


    template <typename T>
    segment_tree<T>::segment_tree(const std::vector<T> &v): tree_(v.size() * 2), size_(v.size()) {
        std::copy(v.begin(), v.end(), tree_.begin() + v.size());
        for(size_t i = v.size() - 1; i > 0; --i)
            tree_[i] = std::min(tree_[2 * i], tree_[2 * i + 1]);
    }

    template <typename T>
    T segment_tree<T>::query(size_t begin, size_t end) const {
        assert(begin < end);

        size_t first = size_ + begin, last = size_ + end - 1;
        T value = tree_[first];
        while(first <= last) {
            value = std::min(value, tree_[first]);
            value = std::min(value, tree_[last]);
            first = (first + 1) / 2;
            last = (last - 1) / 2;
        }
        return value;
    }

    template <typename T>
    void segment_tree<T>::update(size_t idx, T value) {
        assert(idx < size_);

        size_t i = size_ + idx;
        tree_[i] = value;
        for(i /= 2; i > 0; i /= 2)
            tree_[i] = std::min(tree_[2 * i], tree_[2 * i + 1]);
    }
};

template<typename T, typename std::make_unsigned<T>::type Mod>
algo::modular<T, Mod> pow(algo::modular<T, Mod> x, typename std::make_signed<T>::type n) {
    if(n < 0) {
        n = -n;
        x = x.inverse();
    }
    algo::modular<T, Mod> ret = 1;
    algo::modular<T, Mod> mul = x;
    while(n) {
        if(n & 1) ret *= mul;
        mul *= mul;
        n >>= 1;
    }
    return ret;
}

template<typename T>
algo::sparse_table<T> std::min(const algo::sparse_table<T> &a, const algo::sparse_table<T> &b) {
    return algo::sparse_table<T>::min(a, b);
}

#endif
