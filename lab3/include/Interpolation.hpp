#pragma once
#include <vector>
#include <string>
#include <ostream>

struct PointTable {
    std::vector<double> x;
    std::vector<double> y;
};

struct InterpolationInput {
    PointTable table;
    double xStar;
};

struct InterpolationErrorEstimates {
    double actualErrorLagrange = 0.0;
    double actualErrorNewton = 0.0;

    double aPrioriBound = 0.0;          // по формуле с M_{n+1}
    double aPosterioriEstimate = 0.0;   // по первому отброшенному члену Ньютона
};

struct InterpolationResult {
    double lagrangeValue = 0.0;
    double newtonValue = 0.0;
    double exactValue = 0.0;

    std::vector<std::vector<double>> dividedDifferences;
    InterpolationErrorEstimates errors;
};

InterpolationInput loadInterpolationInputFromFile(const std::string& filename);

double lagrangeInterpolate(const std::vector<double>& x,
                           const std::vector<double>& y,
                           double xStar,
                           std::ostream& log);

std::vector<std::vector<double>> buildDividedDifferences(
    const std::vector<double>& x,
    const std::vector<double>& y,
    std::ostream& log);

double newtonInterpolate(const std::vector<double>& x,
                         const std::vector<std::vector<double>>& dd,
                         double xStar,
                         std::ostream& log);

double omegaProduct(const std::vector<double>& x, double xStar);

double computeInterpolationAPrioriBound(const std::vector<double>& x,
                                        double xStar,
                                        double maxAbsDerivative);

void run_3_1(const std::string& inputFile);