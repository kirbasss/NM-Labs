#include "LeastSquares.hpp"
#include "LU.hpp"
#include "Matrix.hpp"

#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace {
    const double EPS = 1e-12;

    void validateInput(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Размеры x и y не совпадают");
        }
        if (x.empty()) {
            throw std::invalid_argument("Пустая таблица данных");
        }
    }

    std::vector<double> computePowerSums(const std::vector<double>& x, int maxPower) {
        std::vector<double> sums(maxPower + 1, 0.0);
        for (double xi : x) {
            double p = 1.0;
            for (int k = 0; k <= maxPower; ++k) {
                sums[k] += p;
                p *= xi;
            }
        }
        return sums;
    }

    std::vector<double> computeWeightedPowerSums(const std::vector<double>& x,
                                                 const std::vector<double>& y,
                                                 int maxPower) {
        std::vector<double> sums(maxPower + 1, 0.0);
        for (size_t i = 0; i < x.size(); ++i) {
            double p = 1.0;
            for (int k = 0; k <= maxPower; ++k) {
                sums[k] += y[i] * p;
                p *= x[i];
            }
        }
        return sums;
    }
}

LeastSquaresInput loadLeastSquaresInputFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filename);
    }

    size_t n;
    file >> n;
    if (!file) {
        throw std::runtime_error("Ошибка чтения числа точек из файла " + filename);
    }

    LeastSquaresInput input;
    input.x.resize(n);
    input.y.resize(n);

    for (size_t i = 0; i < n; ++i) {
        file >> input.x[i];
    }
    for (size_t i = 0; i < n; ++i) {
        file >> input.y[i];
    }

    if (!file) {
        throw std::runtime_error("Ошибка чтения табличных данных из файла " + filename);
    }

    validateInput(input.x, input.y);
    return input;
}

double evaluatePolynomial(const std::vector<double>& coeffs, double x) {
    double value = 0.0;
    double p = 1.0;
    for (double a : coeffs) {
        value += a * p;
        p *= x;
    }
    return value;
}

PolynomialFit leastSquaresFit(const std::vector<double>& x,
                              const std::vector<double>& y,
                              int degree,
                              std::ostream& log) {
    validateInput(x, y);

    if (degree < 0) {
        throw std::invalid_argument("Степень полинома должна быть неотрицательной");
    }

    const int m = degree + 1;
    if (static_cast<size_t>(m) > x.size()) {
        throw std::invalid_argument("Слишком высокая степень полинома для данного числа точек");
    }

    log << "=== Метод наименьших квадратов, степень " << degree << " ===\n";

    std::vector<double> powerSums = computePowerSums(x, 2 * degree);
    std::vector<double> weightedSums = computeWeightedPowerSums(x, y, degree);

    log << "Суммы степеней x:\n";
    for (int k = 0; k <= 2 * degree; ++k) {
        log << "sum x^" << k << " = " << powerSums[k] << "\n";
    }
    log << "\n";

    log << "Суммы y*x^k:\n";
    for (int k = 0; k <= degree; ++k) {
        log << "sum y*x^" << k << " = " << weightedSums[k] << "\n";
    }
    log << "\n";

    Matrix A(m);
    std::vector<double> b(m, 0.0);

    for (int row = 0; row < m; ++row) {
        for (int col = 0; col < m; ++col) {
            A(row, col) = powerSums[row + col];
        }
        b[row] = weightedSums[row];
    }

    log << "Нормальная система:\n";
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            log << std::setw(14) << A(i, j) << " ";
        }
        log << " | " << b[i] << "\n";
    }
    log << "\n";

    auto lu = luDecompose(A, log);
    std::vector<double> coeffs = solveLU(lu, b);

    log << "Коэффициенты полинома:\n";
    for (size_t i = 0; i < coeffs.size(); ++i) {
        log << "a[" << i << "] = " << coeffs[i] << "\n";
    }
    log << "Полином: " << formatPolynomial(coeffs) << "\n\n";

    PolynomialFit fit;
    fit.coeffs = coeffs;
    fit.valuesAtNodes.resize(x.size(), 0.0);
    fit.squaredError = 0.0;

    log << "Значения полинома в узлах и вклад в ошибку:\n";
    for (size_t i = 0; i < x.size(); ++i) {
        fit.valuesAtNodes[i] = evaluatePolynomial(coeffs, x[i]);
        double diff = y[i] - fit.valuesAtNodes[i];
        double sq = diff * diff;
        fit.squaredError += sq;

        log << "i = " << i
            << ", x = " << x[i]
            << ", y = " << y[i]
            << ", F(x) = " << fit.valuesAtNodes[i]
            << ", (y - F)^2 = " << sq << "\n";
    }

    log << "\nСумма квадратов ошибок Phi = " << fit.squaredError << "\n\n";
    return fit;
}

std::string formatPolynomial(const std::vector<double>& coeffs,
                             const std::string& var,
                             int precision) {
    if (coeffs.empty()) {
        return "0";
    }

    std::ostringstream out;
    out << std::fixed << std::setprecision(precision);

    bool firstTerm = true;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        double a = coeffs[i];

        if (std::fabs(a) < EPS) {
            continue;
        }

        if (firstTerm) {
            if (a < 0.0) {
                out << "-";
            }
        } else {
            out << (a >= 0.0 ? " + " : " - ");
        }

        double absA = std::fabs(a);

        if (i == 0) {
            out << absA;
        } else if (i == 1) {
            out << absA << "*" << var;
        } else {
            out << absA << "*" << var << "^" << i;
        }

        firstTerm = false;
    }

    if (firstTerm) {
        return "0";
    }

    return out.str();
}

void run_3_3(const std::string& inputFile) {
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
        auto input = loadLeastSquaresInputFromFile(inputFile);

        log << std::fixed << std::setprecision(10);
        log << "Задача 3.3. Метод наименьших квадратов\n";
        log << "Табличные данные:\n";
        for (size_t i = 0; i < input.x.size(); ++i) {
            log << "i = " << i
                << ": x = " << input.x[i]
                << ", y = " << input.y[i] << "\n";
        }
        log << "\n";

        PolynomialFit fit1 = leastSquaresFit(input.x, input.y, 1, log);
        PolynomialFit fit2 = leastSquaresFit(input.x, input.y, 2, log);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile);
        }

        out << std::fixed << std::setprecision(10);

        out << "=== Задача 3.3. Метод наименьших квадратов ===\n\n";

        out << "Исходные данные:\n";
        for (size_t i = 0; i < input.x.size(); ++i) {
            out << "i = " << i
                << ", x_i = " << input.x[i]
                << ", y_i = " << input.y[i] << "\n";
        }

        out << "\n=== Приближающий многочлен 1-й степени ===\n";
        out << "F1(x) = " << formatPolynomial(fit1.coeffs) << "\n";
        out << "Коэффициенты:\n";
        for (size_t i = 0; i < fit1.coeffs.size(); ++i) {
            out << "a[" << i << "] = " << fit1.coeffs[i] << "\n";
        }
        out << "Сумма квадратов ошибок Phi1 = " << fit1.squaredError << "\n";

        out << "\nЗначения F1(x_i):\n";
        for (size_t i = 0; i < input.x.size(); ++i) {
            out << "x = " << input.x[i] << ", F1(x) = " << fit1.valuesAtNodes[i] << "\n";
        }

        out << "\n=== Приближающий многочлен 2-й степени ===\n";
        out << "F2(x) = " << formatPolynomial(fit2.coeffs) << "\n";
        out << "Коэффициенты:\n";
        for (size_t i = 0; i < fit2.coeffs.size(); ++i) {
            out << "a[" << i << "] = " << fit2.coeffs[i] << "\n";
        }
        out << "Сумма квадратов ошибок Phi2 = " << fit2.squaredError << "\n";

        out << "\nЗначения F2(x_i):\n";
        for (size_t i = 0; i < input.x.size(); ++i) {
            out << "x = " << input.x[i] << ", F2(x) = " << fit2.valuesAtNodes[i] << "\n";
        }

        out << "\n=== Сравнение ошибок ===\n";
        if (fit2.squaredError < fit1.squaredError) {
            out << "Полином 2-й степени даёт лучшее приближение по МНК.\n";
        } else if (fit2.squaredError > fit1.squaredError) {
            out << "Полином 1-й степени даёт лучшее приближение по МНК.\n";
        } else {
            out << "Суммы квадратов ошибок совпадают.\n";
        }

        // Машиночитаемый блок для Python
        out << "\n=== PLOT_DATA ===\n";
        out << "N " << input.x.size() << "\n";

        out << "X ";
        for (double v : input.x) out << v << " ";
        out << "\n";

        out << "Y ";
        for (double v : input.y) out << v << " ";
        out << "\n";

        out << "C1 ";
        for (double v : fit1.coeffs) out << v << " ";
        out << "\n";

        out << "C2 ";
        for (double v : fit2.coeffs) out << v << " ";
        out << "\n";

        std::cout << "Алгоритм 3.3 (МНК) завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << "\n";
        std::cerr << e.what() << std::endl;
    }
}