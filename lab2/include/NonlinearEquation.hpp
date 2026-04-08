#pragma once
#include <vector>
#include <string>
#include <ostream>

struct RootResult {
    double root;
    int iterations;
};

RootResult solveNewton(double a, double b, double eps, std::ostream& log);

RootResult solveSimpleIteration(double a, double b, double eps, std::ostream& log);

void run_2_1(const std::string& inputFile);