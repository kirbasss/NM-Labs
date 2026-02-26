#include "Matrix.hpp"
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include "utils.hpp"

Matrix::Matrix(size_t size) : n(size), data(size, std::vector<double>(size, 0.0)) {}

Matrix::Matrix(const std::vector<std::vector<double>>& d) 
    : n(d.size()), data(d) {}

double& Matrix::operator()(size_t i, size_t j) { return data[i][j]; }
const double& Matrix::operator()(size_t i, size_t j) const { return data[i][j]; }

void Matrix::swapRows(size_t r1, size_t r2) {
    if (r1 != r2) std::swap(data[r1], data[r2]);
}

void Matrix::print(std::ostream& os, int prec) const {
    os << std::fixed << std::setprecision(prec);
    for (const auto& row : data) {
        for (double v : row) os << std::setw(12) << v << " ";
        os << '\n';
    }
}

void Matrix::saveToFile(const std::string& filename, int prec) const {
    std::ofstream file(filename);
    if (!file) return;
    print(file, prec);
}
