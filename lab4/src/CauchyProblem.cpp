#include "CauchyProblem.hpp"
#include "OdeHelpers.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
double rhs(double x, double y) {
    return std::sin(3.0 * x) - y;
}

double exactY(double x) {
    return std::cos(x) + 11.0 * std::sin(x) / 8.0 - std::sin(3.0 * x) / 8.0;
}

double exactZ(double x) {
    return -std::sin(x) + 11.0 * std::cos(x) / 8.0 - 3.0 * std::cos(3.0 * x) / 8.0;
}

double maxAbsErrorToExact(const ODESolution& solution) {
    double error = 0.0;
    for (size_t i = 0; i < solution.x.size(); ++i) {
        error = std::max(error, std::fabs(solution.y[i] - exactY(solution.x[i])));
    }
    return error;
}

void writeMethodTable(std::ostream& out,
                      const std::string& title,
                      const ODESolution& solution,
                      bool includeDerivative = false) {
    out << title << "\n";
    out << "x\t\ty\t\ty_exact\t\t|y - y_exact|";
    if (includeDerivative) {
        out << "\t\tz = y'\t\tz_exact\t\t|z - z_exact|";
    }
    out << "\n";
    for (size_t i = 0; i < solution.x.size(); ++i) {
        double yExact = exactY(solution.x[i]);
        double zExact = exactZ(solution.x[i]);
        out << solution.x[i] << "\t"
            << solution.y[i] << "\t"
            << yExact << "\t"
            << std::fabs(solution.y[i] - yExact);
        if (includeDerivative) {
            out << "\t" << solution.z[i]
                << "\t" << zExact
                << "\t" << std::fabs(solution.z[i] - zExact);
        }
        out << "\n";
    }
    out << "\n";
}

void writeSummary(std::ostream& out,
                  const std::string& title,
                  const ODESolution& coarse,
                  const ODESolution& fine,
                  int order) {
    out << title << "\n";
    out << "max |y_h - y_exact|     = " << maxAbsErrorToExact(coarse) << "\n";
    out << "max |y_h/2 - y_exact|   = " << maxAbsErrorToExact(fine) << "\n";
    out << "Runge-Romberg estimate  = "
        << maxAbsErrorRungeRomberg(coarse.x, coarse.y, fine.x, fine.y, order) << "\n\n";
}
}

CauchyInput loadCauchyInputFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filename);
    }

    CauchyInput input;
    file >> input.x0 >> input.xk;
    file >> input.y0 >> input.z0;
    file >> input.h1 >> input.h2;

    if (!file) {
        throw std::runtime_error("Ошибка чтения параметров задачи Коши из файла " + filename);
    }

    return input;
}

ODESolution solveEuler(double x0, double xk, double y0, double z0, double h, std::ostream& log) {
    int n = computeStepCount(x0, xk, h);

    ODESolution solution;
    solution.x.resize(n + 1);
    solution.y.resize(n + 1);
    solution.z.resize(n + 1);

    solution.x[0] = x0;
    solution.y[0] = y0;
    solution.z[0] = z0;

    log << "=== Метод Эйлера (явный), h = " << h << " ===\n";
    log << "k = 0, x = " << solution.x[0]
        << ", y = " << solution.y[0]
        << ", z = " << solution.z[0] << "\n";

    for (int i = 0; i < n; ++i) {
        double x = solution.x[i];
        double y = solution.y[i];
        double z = solution.z[i];

        solution.x[i + 1] = x + h;
        solution.y[i + 1] = y + h * z;
        solution.z[i + 1] = z + h * rhs(x, y);

        log << "k = " << (i + 1)
            << ", x = " << solution.x[i + 1]
            << ", y = " << solution.y[i + 1]
            << ", z = " << solution.z[i + 1]
            << ", |dy| = " << std::fabs(solution.y[i + 1] - y)
            << "\n";
    }

    log << "\n";
    return solution;
}

ODESolution solveRungeKutta4(double x0, double xk, double y0, double z0, double h, std::ostream& log) {
    int n = computeStepCount(x0, xk, h);

    ODESolution solution;
    solution.x.resize(n + 1);
    solution.y.resize(n + 1);
    solution.z.resize(n + 1);

    solution.x[0] = x0;
    solution.y[0] = y0;
    solution.z[0] = z0;

    log << "=== Метод Рунге-Кутты 4-го порядка, h = " << h << " ===\n";
    log << "k = 0, x = " << solution.x[0]
        << ", y = " << solution.y[0]
        << ", z = " << solution.z[0] << "\n";

    for (int i = 0; i < n; ++i) {
        double x = solution.x[i];
        double y = solution.y[i];
        double z = solution.z[i];

        double K1 = h * z;
        double L1 = h * rhs(x, y);

        double K2 = h * (z + L1 / 2.0);
        double L2 = h * rhs(x + h / 2.0, y + K1 / 2.0);

        double K3 = h * (z + L2 / 2.0);
        double L3 = h * rhs(x + h / 2.0, y + K2 / 2.0);

        double K4 = h * (z + L3);
        double L4 = h * rhs(x + h, y + K3);

        solution.x[i + 1] = x + h;
        solution.y[i + 1] = y + (K1 + 2.0 * K2 + 2.0 * K3 + K4) / 6.0;
        solution.z[i + 1] = z + (L1 + 2.0 * L2 + 2.0 * L3 + L4) / 6.0;

        log << "k = " << (i + 1)
            << ", x = " << solution.x[i + 1]
            << ", y = " << solution.y[i + 1]
            << ", z = " << solution.z[i + 1]
            << ", K = (" << K1 << ", " << K2 << ", " << K3 << ", " << K4 << ")"
            << ", L = (" << L1 << ", " << L2 << ", " << L3 << ", " << L4 << ")"
            << "\n";
    }

    log << "\n";
    return solution;
}

ODESolution solveAdams4(double x0, double xk, double y0, double z0, double h, std::ostream& log) {
    int n = computeStepCount(x0, xk, h);
    if (n < 4) {
        throw std::invalid_argument("Для метода Адамса 4-го порядка требуется не менее 4 шагов");
    }

    log << "=== Метод Адамса 4-го порядка, h = " << h << " ===\n";
    log << "Первые четыре точки получаем методом Рунге-Кутты 4-го порядка\n";

    ODESolution solution = solveRungeKutta4(x0, x0 + 3.0 * h, y0, z0, h, log);
    solution.x.resize(n + 1);
    solution.y.resize(n + 1);
    solution.z.resize(n + 1);
    for (int i = 4; i <= n; ++i) {
        solution.x[i] = x0 + i * h;
    }

    for (int i = 3; i < n; ++i) {
        double f0 = rhs(solution.x[i], solution.y[i]);
        double f1 = rhs(solution.x[i - 1], solution.y[i - 1]);
        double f2 = rhs(solution.x[i - 2], solution.y[i - 2]);
        double f3 = rhs(solution.x[i - 3], solution.y[i - 3]);

        solution.y[i + 1] = solution.y[i] + h * (
            55.0 * solution.z[i]
            - 59.0 * solution.z[i - 1]
            + 37.0 * solution.z[i - 2]
            - 9.0 * solution.z[i - 3]
        ) / 24.0;

        solution.z[i + 1] = solution.z[i] + h * (
            55.0 * f0
            - 59.0 * f1
            + 37.0 * f2
            - 9.0 * f3
        ) / 24.0;

        log << "k = " << (i + 1)
            << ", x = " << solution.x[i + 1]
            << ", y = " << solution.y[i + 1]
            << ", z = " << solution.z[i + 1]
            << ", f = (" << f3 << ", " << f2 << ", " << f1 << ", " << f0 << ")"
            << "\n";
    }

    log << "\n";
    return solution;
}

void run_4_1(const std::string& inputFile) {
    namespace fs = std::filesystem;

    const fs::path resolvedInput = resolveInputPath(inputFile);
    const fs::path labRoot = resolvedInput.parent_path().parent_path();
    const std::string base = resolvedInput.stem().string();
    const fs::path outFile = labRoot / "output" / (base + ".txt");
    const fs::path logFile = labRoot / "output" / (base + "_log.txt");

    fs::create_directories(outFile.parent_path());

    std::ofstream log(logFile);
    if (!log) {
        std::cerr << "Не удалось создать " << logFile << "\n";
        return;
    }

    try {
        CauchyInput input = loadCauchyInputFromFile(resolvedInput.string());

        log << std::fixed << std::setprecision(10);
        log << "Задача 4.1. Решение задачи Коши\n";
        log << "y'' + y - sin(3x) = 0\n";
        log << "y(" << input.x0 << ") = " << input.y0
            << ", y'(" << input.x0 << ") = " << input.z0 << "\n";
        log << "Интервал: [" << input.x0 << ", " << input.xk << "]\n";
        log << "h1 = " << input.h1 << ", h2 = " << input.h2 << "\n\n";

        ODESolution euler1 = solveEuler(input.x0, input.xk, input.y0, input.z0, input.h1, log);
        ODESolution euler2 = solveEuler(input.x0, input.xk, input.y0, input.z0, input.h2, log);

        ODESolution rk1 = solveRungeKutta4(input.x0, input.xk, input.y0, input.z0, input.h1, log);
        ODESolution rk2 = solveRungeKutta4(input.x0, input.xk, input.y0, input.z0, input.h2, log);

        ODESolution adams1 = solveAdams4(input.x0, input.xk, input.y0, input.z0, input.h1, log);
        ODESolution adams2 = solveAdams4(input.x0, input.xk, input.y0, input.z0, input.h2, log);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile.string());
        }

        out << std::fixed << std::setprecision(10);
        out << "=== Задача 4.1. Решение задачи Коши ===\n\n";
        out << "Точное решение: y = cos(x) + 11/8 sin(x) - 1/8 sin(3x)\n";
        out << "Интервал: [" << input.x0 << ", " << input.xk << "]\n";
        out << "h1 = " << input.h1 << ", h2 = " << input.h2 << "\n\n";
        out << "В обозначениях метода Рунге-Кутты: K_i -- приращения для y, L_i -- приращения для z = y'.\n\n";

        writeMethodTable(out, "=== Метод Эйлера (h1) ===", euler1, true);
        writeMethodTable(out, "=== Метод Эйлера (h2) ===", euler2, true);
        writeMethodTable(out, "=== Метод Рунге-Кутты 4-го порядка (h1) ===", rk1, true);
        writeMethodTable(out, "=== Метод Рунге-Кутты 4-го порядка (h2) ===", rk2, true);
        writeMethodTable(out, "=== Метод Адамса 4-го порядка (h1) ===", adams1, true);
        writeMethodTable(out, "=== Метод Адамса 4-го порядка (h2) ===", adams2, true);

        out << "=== Оценка погрешности ===\n\n";
        writeSummary(out, "Метод Эйлера", euler1, euler2, 1);
        writeSummary(out, "Метод Рунге-Кутты 4-го порядка", rk1, rk2, 4);
        writeSummary(out, "Метод Адамса 4-го порядка", adams1, adams2, 4);

        std::cout << "Алгоритм 4.1 (задача Коши для ОДУ) завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";
    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << "\n";
        std::cerr << e.what() << "\n";
    }
}
