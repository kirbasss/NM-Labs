#pragma once
#include "Matrix.hpp"
#include <vector>
#include <complex>
#include <string>
#include <ostream>

struct QRDecomposition {
    Matrix Q;
    Matrix R;
};

struct QREigenResult {
    std::vector<std::complex<double>> eigenvalues;
    Matrix finalA;      // финальная почти квазитреугольная матрица
    int iterations = 0;
};

// Критерий для вещественного собственного значения
// sqrt(sum_{m=l+1}^{n-1} a_{m,l}^2)
double columnSubdiagNorm(const Matrix& A, size_t l);

std::pair<std::complex<double>, std::complex<double>>
eigenvalues2x2(double a, double b, double c, double d);

double strictLowerNorm(const Matrix& A);

QRDecomposition qrDecomposeHouseholder(const Matrix& A, std::ostream& log);

Matrix reduceToHessenberg(const Matrix& A, std::ostream& log);

QREigenResult qrAlgorithm(const Matrix& A, double eps, std::ostream& log);

void run_1_5(const std::string& inputFile);