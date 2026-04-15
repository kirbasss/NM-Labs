#include "Differentiation.hpp"

#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace {
    const double EPS = 1e-12;

    void validateInput(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Размеры x и y не совпадают");
        }
        if (x.size() < 3) {
            throw std::invalid_argument("Для численного дифференцирования нужно хотя бы 3 точки");
        }
        for (size_t i = 1; i < x.size(); ++i) {
            if (!(x[i] > x[i - 1])) {
                throw std::invalid_argument("Узлы x должны быть строго возрастающими");
            }
        }
    }

    bool almostEqual(double a, double b, double eps = EPS) {
        return std::fabs(a - b) < eps;
    }

    size_t findLeftIndex(const std::vector<double>& x, double xStar) {
        const size_t n = x.size();

        for (size_t i = 0; i < n; ++i) {
            if (almostEqual(x[i], xStar)) {
                if (i == 0) {
                    return 0;
                }
                if (i == n - 1) {
                    return n - 3;
                }
                return i - 1;
            }
        }

        for (size_t i = 0; i + 1 < n; ++i) {
            if (xStar >= x[i] - EPS && xStar <= x[i + 1] + EPS) {
                if (i + 2 < n) {
                    return i;
                }
                return n - 3;
            }
        }

        throw std::out_of_range("Точка x* находится вне диапазона таблицы");
    }
}

DifferentiationInput loadDifferentiationInputFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filename);
    }

    size_t n;
    file >> n;
    if (!file) {
        throw std::runtime_error("Ошибка чтения числа точек из файла " + filename);
    }

    DifferentiationInput input;
    input.x.resize(n);
    input.y.resize(n);

    for (size_t i = 0; i < n; ++i) {
        file >> input.x[i];
    }
    for (size_t i = 0; i < n; ++i) {
        file >> input.y[i];
    }
    file >> input.xStar;

    if (!file) {
        throw std::runtime_error("Ошибка чтения данных дифференцирования из файла " + filename);
    }

    validateInput(input.x, input.y);
    return input;
}

DerivativeResult differentiateAtPoint(const std::vector<double>& x,
                                      const std::vector<double>& y,
                                      double xStar,
                                      std::ostream& log) {
    validateInput(x, y);

    size_t i = findLeftIndex(x, xStar);

    const double x0 = x[i];
    const double x1 = x[i + 1];
    const double x2 = x[i + 2];

    const double y0 = y[i];
    const double y1 = y[i + 1];
    const double y2 = y[i + 2];

    const double h01 = x1 - x0;
    const double h12 = x2 - x1;
    const double h02 = x2 - x0;

    if (std::fabs(h01) < EPS || std::fabs(h12) < EPS || std::fabs(h02) < EPS) {
        throw std::runtime_error("Некорректные расстояния между узлами");
    }

    log << "=== Численное дифференцирование ===\n";
    log << "Используются точки:\n";
    log << "x0 = " << x0 << ", y0 = " << y0 << "\n";
    log << "x1 = " << x1 << ", y1 = " << y1 << "\n";
    log << "x2 = " << x2 << ", y2 = " << y2 << "\n";
    log << "x* = " << xStar << "\n\n";

    log << "Формула для первой производной:\n";
    log << "f'(x) ≈ (y1-y0)/(x1-x0)"
        << " + [ ((y2-y1)/(x2-x1) - (y1-y0)/(x1-x0)) / (x2-x0) ] * (2x - x0 - x1)\n\n";

    const double slope01 = (y1 - y0) / h01;
    const double slope12 = (y2 - y1) / h12;
    const double correction = (slope12 - slope01) / h02;

    const double firstDerivative =
        slope01 + correction * (2.0 * xStar - x0 - x1);

    log << "(y1 - y0)/(x1 - x0) = (" << y1 << " - " << y0 << ")/(" << x1 << " - " << x0
        << ") = " << slope01 << "\n";
    log << "(y2 - y1)/(x2 - x1) = (" << y2 << " - " << y1 << ")/(" << x2 << " - " << x1
        << ") = " << slope12 << "\n";
    log << "correction = (" << slope12 << " - " << slope01 << ")/(" << x2 << " - " << x0
        << ") = " << correction << "\n";
    log << "f'(x*) = " << slope01 << " + " << correction << " * (2*" << xStar
        << " - " << x0 << " - " << x1 << ") = " << firstDerivative << "\n\n";

    log << "Формула для второй производной:\n";
    log << "f''(x) ≈ 2 * [ ((y2-y1)/(x2-x1) - (y1-y0)/(x1-x0)) / (x2-x0) ]\n\n";

    const double secondDerivative = 2.0 * correction;

    log << "f''(x*) = 2 * " << correction << " = " << secondDerivative << "\n\n";

    DerivativeResult result;
    result.firstDerivative = firstDerivative;
    result.secondDerivative = secondDerivative;
    result.leftIndex = i;
    return result;
}

void run_3_4(const std::string& inputFile) {
    namespace fs = std::filesystem;

    std::string base    = fs::path(inputFile).stem().string();
    std::string outFile = "output/" + base + ".txt";
    std::string logFile = "output/" + base + "_log.txt";

    fs::create_directories("output");

    std::ofstream log(logFile);
    if (!log) {
        std::cerr << "Не удалось создать " << logFile << std::endl;
        return;
    }

    try {
        auto input = loadDifferentiationInputFromFile(inputFile);

        log << std::fixed << std::setprecision(10);
        log << "Задача 3.4. Численное дифференцирование\n";
        log << "Табличные данные:\n";
        for (size_t i = 0; i < input.x.size(); ++i) {
            log << "i = " << i
                << ": x = " << input.x[i]
                << ", y = " << input.y[i] << "\n";
        }
        log << "\nx* = " << input.xStar << "\n\n";

        DerivativeResult result = differentiateAtPoint(input.x, input.y, input.xStar, log);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile);
        }

        out << std::fixed << std::setprecision(10);

        out << "=== Задача 3.4. Численное дифференцирование ===\n\n";

        out << "Исходные данные:\n";
        for (size_t i = 0; i < input.x.size(); ++i) {
            out << "i = " << i
                << ", x_i = " << input.x[i]
                << ", y_i = " << input.y[i] << "\n";
        }

        out << "\nТочка x* = " << input.xStar << "\n";
        out << "Использованы узлы:\n";
        out << "x[" << result.leftIndex << "] = " << input.x[result.leftIndex] << "\n";
        out << "x[" << (result.leftIndex + 1) << "] = " << input.x[result.leftIndex + 1] << "\n";
        out << "x[" << (result.leftIndex + 2) << "] = " << input.x[result.leftIndex + 2] << "\n";

        out << "\n=== Результат ===\n";
        out << "f'(x*)  = " << result.firstDerivative << "\n";
        out << "f''(x*) = " << result.secondDerivative << "\n";

        std::cout << "Алгоритм 3.4 (численное дифференцирование) завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << "\n";
        std::cerr << e.what() << std::endl;
    }
}