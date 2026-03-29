#include "MATH__.hpp"
#include <cmath>    
#include <cstdlib>   

template<typename T>
Matrix<T>::Matrix(int rows, int cols) : rows(rows), cols(cols) {
    data.resize(rows * cols, T(0));
}

template<typename T>
T& Matrix<T>::operator()(int i, int j) {
    return data[i * cols + j];
}

template<typename T>
const T& Matrix<T>::operator()(int i, int j) const {
    return data[i * cols + j];
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {
    if (rows != other.rows || cols != other.cols)
        throw std::invalid_argument("Matrix addition: dimension mismatch");

    Matrix<T> result(rows, cols);
    constexpr int N = mipp::N<T>();
    int vec_end = (cols / N) * N;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < vec_end; j += N) {
            mipp::Reg<T> rA = mipp::loadu<T>(&(*this)(i, j));
            mipp::Reg<T> rB = mipp::loadu<T>(&other(i, j));
            mipp::Reg<T> rC = rA + rB;
            mipp::storeu(&result(i, j), rC);
        }
        for (int j = vec_end; j < cols; ++j) {
            result(i, j) = (*this)(i, j) + other(i, j);
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const {
    if (rows != other.rows || cols != other.cols)
        throw std::invalid_argument("Matrix subtraction: dimension mismatch");

    Matrix<T> result(rows, cols);
    constexpr int N = mipp::N<T>();
    int vec_end = (cols / N) * N;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < vec_end; j += N) {
            mipp::Reg<T> rA = mipp::loadu<T>(&(*this)(i, j));
            mipp::Reg<T> rB = mipp::loadu<T>(&other(i, j));
            mipp::Reg<T> rC = rA - rB;
            mipp::storeu(&result(i, j), rC);
        }
        for (int j = vec_end; j < cols; ++j) {
            result(i, j) = (*this)(i, j) - other(i, j);
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {
    if (cols != other.rows)
        throw std::invalid_argument("Matrix multiplication: incompatible dimensions");

    Matrix<T> result(rows, other.cols);
    constexpr int N = mipp::N<T>();
    int vec_end = (other.cols / N) * N;

    for (int i = 0; i < rows; ++i) {
        for (int k = 0; k < cols; ++k) {
            T aik = (*this)(i, k);
            mipp::Reg<T> aikk_reg(aik);
            for (int j = 0; j < vec_end; j += N) {
                mipp::Reg<T> rResult = mipp::loadu<T>(&result(i, j));
                mipp::Reg<T> rOther = mipp::loadu<T>(&other(k, j));
                mipp::Reg<T> rProduct = aikk_reg * rOther;
                mipp::Reg<T> rNewResult = rResult + rProduct;
                mipp::storeu(&result(i, j), rNewResult);
            }
            for (int j = vec_end; j < other.cols; ++j) {
                result(i, j) += aik * other(k, j);
            }
        }
    }
    return result;
}

std::vector<int> Generate_Poisson_Knuth_SIMD(int n, int lambda) {
    constexpr int N = mipp::N<float>();
    const float L = std::exp(-lambda);
    mipp::Reg<float> rL(L);

    std::vector<int> result(n);

    int i = 0;
    for (; i + N <= n; i += N) {
        mipp::Reg<float> p(1.0f);
        mipp::Reg<int>   k(0);
        mipp::Msk<N> active(true);

        while (!mipp::testz(active)) {
            float u_arr[N];
            for (int j = 0; j < N; ++j)
                u_arr[j] = static_cast<float>(std::rand()) / RAND_MAX;
            mipp::Reg<float> u(u_arr);

            mipp::Msk<N> cond = p > rL;
            p = mipp::blend(p * u, p, cond);
            k = mipp::blend(k + 1, k, cond);
            active = cond;
        }

        mipp::Reg<int> res = k - 1;
        res.store(&result[i]);
    }

    for (; i < n; i++) {
        float p = 1.0f;
        int k = 0;
        while (p > L) {
            p *= static_cast<float>(std::rand()) / RAND_MAX;
            ++k;
        }
        result[i] = k - 1;
    }

    return result;
}
template class Matrix<int>;