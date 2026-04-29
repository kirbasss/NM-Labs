#pragma once

#include <ostream>
#include <string>
#include <vector>

struct BvpInput {
    double a = 0.0;
    double b = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double leftValue = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    double rightValue = 0.0;
    double h1 = 0.0;
    double h2 = 0.0;
    double eps = 0.0;
    double guess0 = 0.0;
    double guess1 = 0.0;
};

struct BvpSolution {
    std::vector<double> x;
    std::vector<double> y;
    int iterations = 0;
    double residual = 0.0;
};

BvpInput loadBvpInputFromFile(const std::string& filename);

BvpSolution solveShooting(const BvpInput& input, double h, std::ostream& log);

BvpSolution solveFiniteDifference(const BvpInput& input, double h, std::ostream& log);

void run_4_2(const std::string& inputFile);
