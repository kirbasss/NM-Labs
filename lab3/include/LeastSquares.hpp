#pragma once
#include <vector>
#include <string>
#include <ostream>

struct LeastSquaresInput {
    std::vector<double> x;
    std::vector<double> y;
};

struct PolynomialFit {
    std::vector<double> coeffs;   // a0, a1, ..., am
    double squaredError = 0.0;
    std::vector<double> valuesAtNodes;
};

LeastSquaresInput loadLeastSquaresInputFromFile(const std::string& filename);

double evaluatePolynomial(const std::vector<double>& coeffs, double x);

PolynomialFit leastSquaresFit(const std::vector<double>& x,
                              const std::vector<double>& y,
                              int degree,
                              std::ostream& log);

std::string formatPolynomial(const std::vector<double>& coeffs,
                             const std::string& var = "x",
                             int precision = 10);

void run_3_3(const std::string& inputFile);