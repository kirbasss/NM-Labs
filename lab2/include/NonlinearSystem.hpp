#pragma once

#include "Matrix.hpp"
#include <functional>
#include <ostream>
#include <string>
#include <vector>

struct SystemResult {
    std::vector<double> x;
    int iterations;
};

using VectorFunction   = std::function<std::vector<double>(const std::vector<double>&)>;
using JacobianFunction = std::function<Matrix(const std::vector<double>&)>;

// Простая итерация для конкретной системы 2x2 из задания 2.2
SystemResult solveSimpleIterationSystem(double a1, double b1,
                                        double a2, double b2,
                                        const std::vector<double>& x0,
                                        double eps,
                                        std::ostream& log);

// Обобщённый метод Ньютона для системы произвольного размера
SystemResult solveNewtonSystem(const std::vector<double>& x0,
                               double eps,
                               const VectorFunction& F,
                               const JacobianFunction& J,
                               std::ostream& log);

void run_2_2(const std::string& inputFile);