#pragma once
#include "Matrix.hpp"
#include "utils.hpp"
#include <vector>
#include <ostream>
#include <utility>

Matrix buildAlpha(const Matrix& A);

std::vector<double> buildBeta(const Matrix& A, const std::vector<double>& b);

void ensureNonZeroDiagonal(Matrix& A, std::vector<double>& b, std::ostream& log);

double residualNorm1(const Matrix& A, const std::vector<double>& b, const std::vector<double>& x);

// Метод простых итераций (Якоби)
std::pair<std::vector<double>, int> simpleIteration(
    const Matrix& A, const std::vector<double>& b, double eps,
    std::ostream& log);

// Метод Зейделя (Гаусс-Зейдель)
std::pair<std::vector<double>, int> gaussSeidel(
    const Matrix& A, const std::vector<double>& b, double eps,
    std::ostream& log);

void run_1_3(const std::string& inputFile);