#pragma once
#include <vector>
#include <ostream>
#include <string>

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