#pragma once
#include "Matrix.hpp"
#include <vector>
#include <string>

const double EPS_ZERO = 1e-12;

struct LinearSystem {
    Matrix A;
    std::vector<double> b;
};

LinearSystem loadSystemFromFile(const std::string& filename);

double vectorNorm2(const std::vector<double>& x);

double vectorNormC(const std::vector<double>& x);

double matrixNormC(const Matrix& A);

std::vector<double> subtractVectors(const std::vector<double>& a,
                                    const std::vector<double>& b);

std::vector<double> multiplyMatrixVector(const Matrix& A,
                                         const std::vector<double>& x);

std::vector<double> addVectors(const std::vector<double>& a,
                               const std::vector<double>& b);

bool isSymmetric(const Matrix& A, double tol = EPS_ZERO);

Matrix transpose(const Matrix& A);

Matrix multiply(const Matrix& A, const Matrix& B);

double matrixDiffNorm1(const Matrix& A, const Matrix& B);

double signNonZero(double x);