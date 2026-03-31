#include "utils.hpp"
#include <fstream>
#include <stdexcept>
#include <cmath>

extern const double EPS_ZERO;

LinearSystem loadSystemFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Не удалось открыть файл " + filename);

    size_t n;
    file >> n;
    Matrix A(n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            file >> A(i, j);

    std::vector<double> b(n);
    for (size_t i = 0; i < n; ++i)
        file >> b[i];

    return {A, b};
}

double vectorNorm2(const std::vector<double>& x) {
    double s = 0.0;
    for (double v : x) s += v * v;
    return std::sqrt(s);
}

double vectorNormC(const std::vector<double>& x) {
    double norm = 0.0;
    for (double v : x) {
        norm = std::max(norm, std::fabs(v));
    }
    return norm;
}

// ||A||_c = max_i sum_j |a_ij|
double matrixNormC(const Matrix& A) {
    size_t n = A.size();
    double norm = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            row_sum += std::fabs(A(i, j));
        }
        norm = std::max(norm, row_sum);
    }
    return norm;
}

std::vector<double> subtractVectors(const std::vector<double>& a,
                                    const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("subtractVectors: incompatible sizes");
    }

    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<double> multiplyMatrixVector(const Matrix& A,
                                         const std::vector<double>& x) {
    size_t n = A.size();
    if (x.size() != n) {
        throw std::invalid_argument("multiplyMatrixVector: incompatible sizes");
    }

    std::vector<double> result(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result[i] += A(i, j) * x[j];
        }
    }
    return result;
}

std::vector<double> addVectors(const std::vector<double>& a,
                               const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("addVectors: incompatible sizes");
    }

    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

bool isSymmetric(const Matrix& A, double tol) {
    size_t n = A.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (std::fabs(A(i, j) - A(j, i)) > tol) {
                return false;
            }
        }
    }
    return true;
}

Matrix transpose(const Matrix& A) {
    size_t n = A.size();
    Matrix T(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            T(i, j) = A(j, i);
        }
    }
    return T;
}

Matrix multiply(const Matrix& A, const Matrix& B) {
    size_t n = A.size();
    Matrix C(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (size_t k = 0; k < n; ++k) {
                s += A(i, k) * B(k, j);
            }
            C(i, j) = s;
        }
    }
    return C;
}

// Поэлементная 1-норма разности двух матриц
double matrixDiffNorm1(const Matrix& A, const Matrix& B) {
    size_t n = A.size();
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            sum += std::fabs(A(i, j) - B(i, j));
        }
    }
    return sum;
}

double signNonZero(double x) {
    return (x >= 0.0 ? 1.0 : -1.0);
}