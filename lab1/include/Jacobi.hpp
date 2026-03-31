#pragma once
#include "Matrix.hpp"
#include <vector>
#include <ostream>
#include <utility>
#include <string>

struct EigenResult {
    std::vector<double> eigenvalues;
    Matrix eigenvectors; // столбцы — собственные векторы
};

struct EigenInput {
    Matrix A;
    double eps;
};

EigenInput loadEigenInputFromFile(const std::string& filename);

EigenResult jacobiEigen(const Matrix& A, double eps, std::ostream& log);

void run_1_4(const std::string& inputFile);