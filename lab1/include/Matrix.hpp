#pragma once
#include <vector>
#include <string>
#include <iostream>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    size_t n = 0;

public:
    Matrix(size_t size = 0);
    Matrix(const std::vector<std::vector<double>>& d);

    size_t size() const { return n; }

    double& operator()(size_t i, size_t j);
    const double& operator()(size_t i, size_t j) const;

    void swapRows(size_t r1, size_t r2);

    void print(std::ostream& os = std::cout, int prec = 8) const;
    void saveToFile(const std::string& filename, int prec = 8) const;
};