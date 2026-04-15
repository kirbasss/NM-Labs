#pragma once
#include <vector>
#include <string>
#include <ostream>

struct DifferentiationInput {
    std::vector<double> x;
    std::vector<double> y;
    double xStar = 0.0;
};

struct DerivativeResult {
    double firstDerivative = 0.0;
    double secondDerivative = 0.0;

    size_t leftIndex = 0;   // используем точки x[i], x[i+1], x[i+2]
};

DifferentiationInput loadDifferentiationInputFromFile(const std::string& filename);

DerivativeResult differentiateAtPoint(const std::vector<double>& x,
                                      const std::vector<double>& y,
                                      double xStar,
                                      std::ostream& log);

void run_3_4(const std::string& inputFile);