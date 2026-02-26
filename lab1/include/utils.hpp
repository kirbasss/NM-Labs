#pragma once
#include "Matrix.hpp"
#include <vector>
#include <string>

struct LinearSystem {
    Matrix A;
    std::vector<double> b;
};

LinearSystem loadSystemFromFile(const std::string& filename);