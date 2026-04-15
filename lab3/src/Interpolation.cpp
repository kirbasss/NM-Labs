#include "Interpolation.hpp"
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace {
    const double EPS = 1e-12;
    const double PI = 3.14159265358979323846;

    double factorial(int n) {
        if (n < 0) {
            throw std::invalid_argument("factorial: n must be non-negative");
        }
        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= static_cast<double>(i);
        }
        return result;
    }

    void validateNodes(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.empty()) {
            throw std::invalid_argument("Список узлов пуст");
        }
        if (x.size() != y.size()) {
            throw std::invalid_argument("Размеры x и y не совпадают");
        }
        for (size_t i = 0; i < x.size(); ++i) {
            for (size_t j = i + 1; j < x.size(); ++j) {
                if (std::fabs(x[i] - x[j]) < EPS) {
                    throw std::invalid_argument("Обнаружены совпадающие узлы интерполяции");
                }
            }
        }
    }
}

InterpolationInput loadInterpolationInputFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filename);
    }

    size_t n;
    file >> n;
    if (!file) {
        throw std::runtime_error("Ошибка чтения числа узлов из файла " + filename);
    }

    InterpolationInput input;
    input.table.x.resize(n);
    input.table.y.resize(n);

    for (size_t i = 0; i < n; ++i) file >> input.table.x[i];
    for (size_t i = 0; i < n; ++i) file >> input.table.y[i];
    file >> input.xStar;

    if (!file) {
        throw std::runtime_error("Ошибка чтения данных интерполяции из файла " + filename);
    }

    validateNodes(input.table.x, input.table.y);
    return input;
}

double lagrangeInterpolate(const std::vector<double>& x,
                           const std::vector<double>& y,
                           double xStar,
                           std::ostream& log) {
    validateNodes(x, y);

    const size_t n = x.size();
    double result = 0.0;

    log << "=== Интерполяционный многочлен Лагранжа ===\n";
    log << "L(x) = sum_{i=0}^{n-1} y[i] * l[i](x)\n\n";

    for (size_t i = 0; i < n; ++i) {
        double basis = 1.0;

        log << "i = " << i << ":\n";
        log << "  l[" << i << "](" << xStar << ") = ";

        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;

            double denom = x[i] - x[j];
            if (std::fabs(denom) < EPS) {
                throw std::runtime_error("Нулевой знаменатель в Лагранже");
            }

            double factor = (xStar - x[j]) / denom;
            basis *= factor;

            log << "((" << xStar << " - " << x[j] << ")/(" << x[i] << " - " << x[j] << "))";
            if (j + 1 < n) log << " * ";
        }
        log << " = " << basis << "\n";

        double term = y[i] * basis;
        result += term;

        log << "  term = y[" << i << "] * l[" << i << "] = "
            << y[i] << " * " << basis << " = " << term << "\n\n";
    }

    log << "L(" << xStar << ") = " << result << "\n\n";
    return result;
}

std::vector<std::vector<double>> buildDividedDifferences(
    const std::vector<double>& x,
    const std::vector<double>& y,
    std::ostream& log) {
    validateNodes(x, y);

    const size_t n = x.size();
    std::vector<std::vector<double>> dd(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        dd[i][0] = y[i];
    }

    log << "=== Таблица разделённых разностей ===\n";
    log << "Порядок 0:\n";
    for (size_t i = 0; i < n; ++i) {
        log << "f[x" << i << "] = " << dd[i][0] << "\n";
    }
    log << "\n";

    for (size_t order = 1; order < n; ++order) {
        log << "Порядок " << order << ":\n";
        for (size_t i = 0; i + order < n; ++i) {
            double denom = x[i + order] - x[i];
            if (std::fabs(denom) < EPS) {
                throw std::runtime_error("Нулевой знаменатель в таблице разделённых разностей");
            }

            dd[i][order] = (dd[i + 1][order - 1] - dd[i][order - 1]) / denom;

            log << "f[x" << i;
            for (size_t k = 1; k <= order; ++k) {
                log << ",x" << (i + k);
            }
            log << "] = (" << dd[i + 1][order - 1]
                << " - " << dd[i][order - 1]
                << ") / (" << x[i + order] << " - " << x[i] << ") = "
                << dd[i][order] << "\n";
        }
        log << "\n";
    }

    return dd;
}

double newtonInterpolate(const std::vector<double>& x,
                         const std::vector<std::vector<double>>& dd,
                         double xStar,
                         std::ostream& log) {
    const size_t n = x.size();
    if (dd.size() != n) {
        throw std::invalid_argument("Некорректная таблица разделённых разностей");
    }

    log << "=== Интерполяционный многочлен Ньютона ===\n";

    double result = dd[0][0];
    log << "P(x) = f[x0]";
    for (size_t order = 1; order < n; ++order) {
        log << " + f[x0";
        for (size_t j = 1; j <= order; ++j) {
            log << ",x" << j;
        }
        log << "]";
        for (size_t j = 0; j < order; ++j) {
            log << "*(x-x" << j << ")";
        }
    }
    log << "\n\n";

    log << "P(" << xStar << ") = " << dd[0][0];
    for (size_t order = 1; order < n; ++order) {
        double product = 1.0;
        for (size_t j = 0; j < order; ++j) {
            product *= (xStar - x[j]);
        }

        double term = dd[0][order] * product;
        result += term;

        log << " + (" << dd[0][order] << ")";
        for (size_t j = 0; j < order; ++j) {
            log << "*(" << xStar << " - " << x[j] << ")";
        }
        log << " = term " << term << "\n";
    }

    log << "\nP(" << xStar << ") = " << result << "\n\n";
    return result;
}

double omegaProduct(const std::vector<double>& x, double xStar) {
    double w = 1.0;
    for (double xi : x) {
        w *= (xStar - xi);
    }
    return w;
}

double computeInterpolationAPrioriBound(const std::vector<double>& x,
                                        double xStar,
                                        double maxAbsDerivative) {
    const int n = static_cast<int>(x.size()) - 1;
    double w = std::fabs(omegaProduct(x, xStar));
    return maxAbsDerivative * w / factorial(n + 1);
}

void run_3_1(const std::string& inputFile) {
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
        auto input = loadInterpolationInputFromFile(inputFile);
        const auto& x = input.table.x;
        const auto& y = input.table.y;
        const double xStar = input.xStar;

        log << std::fixed << std::setprecision(10);
        log << "Задача 3.1. Интерполяция\n";
        log << "Число узлов: " << x.size() << "\n";
        log << "x* = " << xStar << "\n\n";

        log << "Таблица узлов:\n";
        for (size_t i = 0; i < x.size(); ++i) {
            log << "i = " << i << ": x = " << x[i] << ", y = " << y[i] << "\n";
        }
        log << "\n";

        double lagrangeValue = lagrangeInterpolate(x, y, xStar, log);
        auto dd = buildDividedDifferences(x, y, log);
        double newtonValue = newtonInterpolate(x, dd, xStar, log);

        double exactValue = std::sin(xStar);

        // Для sin(x): |f^(k)(x)| <= 1 для любого порядка
        double maxAbsDerivative = 1.0;
        double aPrioriBound = computeInterpolationAPrioriBound(x, xStar, maxAbsDerivative);

        double actualErrorLagrange = std::fabs(exactValue - lagrangeValue);
        double actualErrorNewton   = std::fabs(exactValue - newtonValue);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile);
        }

        out << std::fixed << std::setprecision(10);

        out << "=== Задача 3.1. Интерполяция ===\n\n";

        out << "Узлы интерполяции:\n";
        for (size_t i = 0; i < x.size(); ++i) {
            out << "i = " << i
                << ", x_i = " << x[i]
                << ", y_i = " << y[i] << "\n";
        }

        out << "\nТочка x* = " << xStar << "\n";
        out << "Точное значение f(x*) = sin(x*) = " << exactValue << "\n";

        out << "\n=== Значение интерполяционного многочлена Лагранжа ===\n";
        out << lagrangeValue << "\n";

        out << "\n=== Значение интерполяционного многочлена Ньютона ===\n";
        out << newtonValue << "\n";

        out << "\n=== Фактическая абсолютная погрешность ===\n";
        out << "|f(x*) - L(x*)| = " << actualErrorLagrange << "\n";
        out << "|f(x*) - P(x*)| = " << actualErrorNewton << "\n";

        out << "\n=== Априорная оценка погрешности ===\n";
        out << "eps(x*) <= M_(n+1) / (n+1)! * |omega(x*)|\n";
        out << "Для sin(x): max |f^(n+1)(x)| <= 1\n";
        out << "Оценка = " << aPrioriBound << "\n";

        std::cout << "Алгоритм 3.1 (интерполяция) завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << "\n";
        std::cerr << e.what() << std::endl;
    }
}