#pragma once
#include "Matrix.hpp"
#include <vector>

struct LUDecomposition {
    Matrix L;
    Matrix U;
    std::vector<int> P;   // P[i] — исходный номер строки, которая оказалась на позиции i
};

LUDecomposition luDecompose(const Matrix& A, std::ostream& log = std::cout);

std::vector<double> solveLU(const LUDecomposition& lu, const std::vector<double>& b);

double computeDeterminant(const LUDecomposition& lu);

Matrix computeInverse(const LUDecomposition& lu);

void run_1_1(const std::string& inputFile);