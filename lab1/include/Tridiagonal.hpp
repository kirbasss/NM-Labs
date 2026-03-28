#pragma once
#include <vector>
#include <ostream>
#include <string>

struct TridiagonalSystem {
    std::vector<double> a; // нижняя диагональ, размер n-1
    std::vector<double> b; // главная диагональ, размер n
    std::vector<double> c; // верхняя диагональ, размер n-1
    std::vector<double> d; // правая часть, размер n
};

TridiagonalSystem loadTridiagonalSystemFromFile(const std::string& filename);

std::vector<double> tridiagSolve(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c,
                                const std::vector<double>& d,
                                std::ostream& log);

double computeTridiagonalResidual(const std::vector<double>& a,
                                  const std::vector<double>& b,
                                  const std::vector<double>& c,
                                  const std::vector<double>& d,
                                  const std::vector<double>& x);

void run_1_2(const std::string& inputFile);