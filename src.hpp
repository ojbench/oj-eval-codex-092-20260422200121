#ifndef SRC_HPP
#define SRC_HPP

#include <vector>
#include <string>
#include <iostream>
#include <exception>
#include <algorithm>
#include <cmath>
#include "fraction.hpp"

class resistive_network {
private:
    int n, m;
    std::vector<int> u, v;
    std::vector<fraction> R;
    std::vector<std::vector<fraction>> L;

    static std::vector<fraction> solve_linear(std::vector<std::vector<fraction>> A,
                                              const std::vector<fraction>& b) {
        int N = (int)A.size();
        if (N == 0) throw matrix_error();
        for (auto &row : A) if ((int)row.size() != N) throw matrix_error();
        if ((int)b.size() != N) throw matrix_error();

        std::vector<std::vector<fraction>> aug(N, std::vector<fraction>(N + 1, fraction(0)));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) aug[i][j] = A[i][j];
            aug[i][N] = b[i];
        }

        int row = 0;
        for (int col = 0; col < N && row < N; ++col) {
            int pivot = -1;
            for (int r = row; r < N; ++r) {
                if (!(aug[r][col] == fraction(0))) { pivot = r; break; }
            }
            if (pivot == -1) continue;
            if (pivot != row) std::swap(aug[pivot], aug[row]);

            fraction piv = aug[row][col];
            for (int j = col; j <= N; ++j) aug[row][j] = aug[row][j] / piv;

            for (int r = 0; r < N; ++r) {
                if (r == row) continue;
                fraction factor = aug[r][col];
                if (factor == fraction(0)) continue;
                for (int j = col; j <= N; ++j) {
                    aug[r][j] = aug[r][j] - factor * aug[row][j];
                }
            }
            ++row;
        }

        std::vector<fraction> x(N, fraction(0));
        for (int i = 0; i < N; ++i) x[i] = aug[i][N];
        return x;
    }

    std::vector<std::vector<fraction>> reduced_L(int removed) const {
        int N = n - 1;
        std::vector<std::vector<fraction>> M(N, std::vector<fraction>(N, fraction(0)));
        int ii = 0;
        for (int i = 0; i < n; ++i) {
            if (i == removed) continue;
            int jj = 0;
            for (int j = 0; j < n; ++j) {
                if (j == removed) continue;
                M[ii][jj] = L[i][j];
                ++jj;
            }
            ++ii;
        }
        return M;
    }

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        n = interface_size_;
        m = connection_size_;
        u.resize(m); v.resize(m); R.resize(m);
        for (int i = 0; i < m; ++i) {
            u[i] = from[i];
            v[i] = to[i];
            R[i] = resistance[i];
        }
        L.assign(n, std::vector<fraction>(n, fraction(0)));
        for (int i = 0; i < m; ++i) {
            int a = u[i] - 1;
            int b = v[i] - 1;
            fraction g = fraction(1) / R[i];
            L[a][a] = L[a][a] + g;
            L[b][b] = L[b][b] + g;
            L[a][b] = L[a][b] - g;
            L[b][a] = L[b][a] - g;
        }
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int a = interface_id1;
        int b = interface_id2;
        if (a == b) return fraction(0);
        auto A = reduced_L(b - 1);
        std::vector<fraction> rhs(n - 1, fraction(0));
        int ai = (a < b) ? (a - 1) : (a - 2);
        rhs[ai] = fraction(1);
        auto x = solve_linear(A, rhs);
        return x[ai];
    }

    fraction get_voltage(int id, fraction current[]) {
        auto A = reduced_L(n - 1);
        std::vector<fraction> rhs(n - 1, fraction(0));
        for (int i = 0; i < n - 1; ++i) rhs[i] = current[i];
        auto x = solve_linear(A, rhs);
        return x[id - 1];
    }

    fraction get_power(fraction voltage[]) {
        fraction total(0);
        for (int i = 0; i < m; ++i) {
            int a = u[i] - 1;
            int b = v[i] - 1;
            fraction delta = voltage[a] - voltage[b];
            fraction term = (delta * delta) / R[i];
            total = total + term;
        }
        return total;
    }
};

#endif // SRC_HPP

