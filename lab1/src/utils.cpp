#include "utils.hpp"
#include <fstream>
#include <stdexcept>

LinearSystem loadSystemFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Не удалось открыть файл " + filename);

    size_t n;
    file >> n;
    Matrix A(n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            file >> A(i, j);

    std::vector<double> b(n);
    for (size_t i = 0; i < n; ++i)
        file >> b[i];

    return {A, b};
}