#pragma once

#include <ostream>
#include <string>
#include <vector>

struct CauchyInput {
    double x0 = 0.0;
    double xk = 0.0;
    double y0 = 0.0;
    double z0 = 0.0;
    double h1 = 0.0;
    double h2 = 0.0;
};

struct ODESolution {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

CauchyInput loadCauchyInputFromFile(const std::string& filename);

ODESolution solveEuler(double x0, double xk, double y0, double z0, double h, std::ostream& log);

ODESolution solveRungeKutta4(double x0, double xk, double y0, double z0, double h, std::ostream& log);

ODESolution solveAdams4(double x0, double xk, double y0, double z0, double h, std::ostream& log);

void run_4_1(const std::string& inputFile);
