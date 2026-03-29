#ifndef MATH__HPP
#define MATH__HPP

#include "mipp.h"
#include <vector>
#include <stdexcept>

template<typename T>
class Matrix {
public:
    std::vector<T> data;
    int rows;
    int cols;

    Matrix(int rows, int cols);
    ~Matrix() = default;

    T& operator()(int i, int j);
    const T& operator()(int i, int j) const;

    Matrix<T> operator+(const Matrix<T>& other) const;
    Matrix<T> operator-(const Matrix<T>& other) const;
    Matrix<T> operator*(const Matrix<T>& other) const;
};

std::vector<int> Generate_Poisson_Knuth_SIMD(int n, int lambda);

#endif