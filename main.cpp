#include "MATH__.hpp"
#include <iostream>

int main() {
    srand(static_cast<unsigned>(time(nullptr)));

    Matrix<int> A(3, 3);
    Matrix<int> B(3, 3);

    A(0,0) = 1; A(0,1) = 0; A(0,2) = 0;
    A(1,0) = 0; A(1,1) = 2; A(1,2) = 0;
    A(2,0) = 0; A(2,1) = 0; A(2,2) = 3;

    B(0,0) = 6; B(0,1) = 7; B(0,2) =12;
    B(1,0) = 8; B(1,1) = 9; B(1,2) = 13;
    B(2,0) = 10; B(2,1) = 11; B(2,2) =14;

    Matrix<int> C = A + B; 
    Matrix<int> D = A - B; 
    Matrix<int> E = A * B;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            std::cout << "C[" << i << "][" << j << "] = " << C(i,j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            std::cout << "D[" << i << "][" << j << "] = " << D(i,j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            std::cout << "E[" << i << "][" << j << "] = " << E(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std :: cout << std::endl;

    std :: vector <int> poisson_samples = Generate_Poisson_Knuth_SIMD(10, 1);
    for (size_t i = 0; i < poisson_samples.size(); ++i) {
        std::cout << "Poisson sample " << i << ": " << poisson_samples[i] << std::endl;
    }
    return 0;
}