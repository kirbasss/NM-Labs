#include "BoundaryValueProblem.hpp"
#include "OdeHelpers.hpp"

#include "Tridiagonal.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
const int MAX_ITER = 1000;

double p(double x) {
    return 2.0 / x;
}

double q(double x) {
    return -1.0;
}

double f(double x) {
    (void)x;
    return 0.0;
}

double exactY(double x) {
    return std::exp(x) / x;
}

double exactDy(double x) {
    return std::exp(x) * (x - 1.0) / (x * x);
}

double rhs(double x, double y, double z) {
    return f(x) - p(x) * z - q(x) * y;
}

struct Trajectory {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

Trajectory solveSystemRungeKutta(double a,
                                 double b,
                                 double y0,
                                 double z0,
                                 double h) {
    int n = computeStepCount(a, b, h);

    Trajectory result;
    result.x.resize(n + 1);
    result.y.resize(n + 1);
    result.z.resize(n + 1);

    result.x[0] = a;
    result.y[0] = y0;
    result.z[0] = z0;

    for (int i = 0; i < n; ++i) {
        double x = result.x[i];
        double y = result.y[i];
        double z = result.z[i];

        double K1 = h * z;
        double L1 = h * rhs(x, y, z);

        double K2 = h * (z + L1 / 2.0);
        double L2 = h * rhs(x + h / 2.0, y + K1 / 2.0, z + L1 / 2.0);

        double K3 = h * (z + L2 / 2.0);
        double L3 = h * rhs(x + h / 2.0, y + K2 / 2.0, z + L2 / 2.0);

        double K4 = h * (z + L3);
        double L4 = h * rhs(x + h, y + K3, z + L3);

        result.x[i + 1] = x + h;
        result.y[i + 1] = y + (K1 + 2.0 * K2 + 2.0 * K3 + K4) / 6.0;
        result.z[i + 1] = z + (L1 + 2.0 * L2 + 2.0 * L3 + L4) / 6.0;
    }

    return result;
}

double maxAbsErrorToExact(const std::vector<double>& x, const std::vector<double>& y) {
    double error = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        error = std::max(error, std::fabs(y[i] - exactY(x[i])));
    }
    return error;
}

void writeTable(std::ostream& out, const std::string& title, const BvpSolution& solution) {
    out << title << "\n";
    out << "x\t\ty\t\ty_exact\t\t|y - y_exact|\n";
    for (size_t i = 0; i < solution.x.size(); ++i) {
        double yExact = exactY(solution.x[i]);
        out << solution.x[i] << "\t"
            << solution.y[i] << "\t"
            << yExact << "\t"
            << std::fabs(solution.y[i] - yExact) << "\n";
    }
    out << "\n";
}
}

BvpInput loadBvpInputFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filename);
    }

    BvpInput input;
    file >> input.a >> input.b;
    file >> input.alpha >> input.beta >> input.leftValue;
    file >> input.delta >> input.gamma >> input.rightValue;
    file >> input.h1 >> input.h2;
    file >> input.eps;
    file >> input.guess0 >> input.guess1;

    if (!file) {
        throw std::runtime_error("Ошибка чтения параметров краевой задачи из файла " + filename);
    }

    return input;
}

BvpSolution solveShooting(const BvpInput& input, double h, std::ostream& log) {
    if (std::fabs(input.beta) < 1e-12) {
        throw std::invalid_argument("Для метода стрельбы по методичке требуется beta != 0");
    }

    log << "=== Метод стрельбы, h = " << h << " ===\n";
    log << "Решаем задачу Коши с параметром eta = y(a).\n";
    log << "Начальная производная выражается из левого граничного условия:\n";
    log << "y'(a) = (A - alpha * eta) / beta.\n";
    log << "Строим Phi(eta) = delta * y(b, eta) + gamma * y'(b, eta) - B\n";
    log << "и уточняем eta методом секущих.\n";

    double eta0 = input.guess0;
    double eta1 = input.guess1;

    double zEta0 = (input.leftValue - input.alpha * eta0) / input.beta;
    Trajectory tr0 = solveSystemRungeKutta(input.a, input.b, eta0, zEta0, h);
    double phi0 = input.delta * tr0.y.back() + input.gamma * tr0.z.back() - input.rightValue;

    double zEta1 = (input.leftValue - input.alpha * eta1) / input.beta;
    Trajectory tr1 = solveSystemRungeKutta(input.a, input.b, eta1, zEta1, h);
    double phi1 = input.delta * tr1.y.back() + input.gamma * tr1.z.back() - input.rightValue;

    log << "eta0 = " << eta0 << ", Phi(eta0) = " << phi0 << "\n";
    log << "eta1 = " << eta1 << ", Phi(eta1) = " << phi1 << "\n";

    int iter = 0;
    while (iter < MAX_ITER) {
        if (std::fabs(phi1 - phi0) < 1e-12) {
            throw std::runtime_error("В методе стрельбы возникло деление на нуль в формуле секущих");
        }

        const double eta2 = eta1 - (eta1 - eta0) * phi1 / (phi1 - phi0);
        const double zEta2 = (input.leftValue - input.alpha * eta2) / input.beta;
        Trajectory tr2 = solveSystemRungeKutta(input.a, input.b, eta2, zEta2, h);
        const double phi2 = input.delta * tr2.y.back() + input.gamma * tr2.z.back() - input.rightValue;

        log << "k = " << iter
            << ", eta_k+1 = " << eta2
            << ", y(b) = " << tr2.y.back()
            << ", z(b) = " << tr2.z.back()
            << ", Phi(eta_k+1) = " << phi2
            << ", |eta_k+1 - eta_k| = " << std::fabs(eta2 - eta1)
            << "\n";

        if (std::fabs(phi2) < input.eps || std::fabs(eta2 - eta1) < input.eps) {
            return {tr2.x, tr2.y, iter + 1, std::fabs(phi2)};
        }

        eta0 = eta1;
        phi0 = phi1;
        eta1 = eta2;
        phi1 = phi2;
        ++iter;
    }

    throw std::runtime_error("Метод стрельбы не сошёлся");
}

BvpSolution solveFiniteDifference(const BvpInput& input, double h, std::ostream& log) {
    int n = computeStepCount(input.a, input.b, h);
    int size = n + 1;

    std::vector<double> x(size);
    for (int i = 0; i <= n; ++i) {
        x[i] = input.a + i * h;
    }

    std::vector<double> aDiag(size - 1, 0.0);
    std::vector<double> bDiag(size, 0.0);
    std::vector<double> cDiag(size - 1, 0.0);
    std::vector<double> dDiag(size, 0.0);

    log << "=== Конечно-разностный метод, h = " << h << " ===\n";

    bDiag[0] = input.alpha * h - input.beta;
    cDiag[0] = input.beta;
    dDiag[0] = input.leftValue * h;
    log << "Граничное условие слева: (" << bDiag[0] << ") * y0 + ("
        << cDiag[0] << ") * y1 = " << dDiag[0] << "\n";

    for (int i = 1; i < n; ++i) {
        double xi = x[i];
        aDiag[i - 1] = 1.0 - p(xi) * h / 2.0;
        bDiag[i] = q(xi) * h * h - 2.0;
        cDiag[i] = 1.0 + p(xi) * h / 2.0;
        dDiag[i] = f(xi) * h * h;

        log << "i = " << i
            << ": A = " << aDiag[i - 1]
            << ", B = " << bDiag[i]
            << ", C = " << cDiag[i]
            << ", D = " << dDiag[i]
            << "\n";
    }

    aDiag[n - 1] = -input.gamma;
    bDiag[n] = input.delta * h + input.gamma;
    dDiag[n] = input.rightValue * h;
    log << "Граничное условие справа: (" << aDiag[n - 1] << ") * y[n-1] + ("
        << bDiag[n] << ") * y[n] = " << dDiag[n] << "\n";

    std::vector<double> y = tridiagSolve(aDiag, bDiag, cDiag, dDiag, log);
    double residual = std::fabs(input.alpha * y.front() + input.beta * (y[1] - y[0]) / h - input.leftValue)
        + std::fabs(input.delta * y.back() + input.gamma * (y.back() - y[n - 1]) / h - input.rightValue);

    return {x, y, 1, residual};
}

void run_4_2(const std::string& inputFile) {
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
        BvpInput input = loadBvpInputFromFile(resolvedInput.string());

        log << std::fixed << std::setprecision(10);
        log << "Задача 4.2. Решение краевой задачи\n";
        log << "x y'' + 2 y' - x y = 0\n";
        log << "Интервал: [" << input.a << ", " << input.b << "]\n";
        log << "eps = " << input.eps << ", h1 = " << input.h1 << ", h2 = " << input.h2 << "\n\n";

        BvpSolution shooting1 = solveShooting(input, input.h1, log);
        BvpSolution shooting2 = solveShooting(input, input.h2, log);

        BvpSolution diff1 = solveFiniteDifference(input, input.h1, log);
        BvpSolution diff2 = solveFiniteDifference(input, input.h2, log);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile.string());
        }

        out << std::fixed << std::setprecision(10);
        out << "=== Задача 4.2. Решение краевой задачи ===\n\n";
        out << "Точное решение: y = e^x / x\n";
        out << "Интервал: [" << input.a << ", " << input.b << "]\n";
        out << "h1 = " << input.h1 << ", h2 = " << input.h2 << "\n\n";

        writeTable(out, "=== Метод стрельбы (h1) ===", shooting1);
        out << "iter = " << shooting1.iterations << "\n";
        out << "Невязка граничного условия = " << shooting1.residual << "\n\n";

        writeTable(out, "=== Метод стрельбы (h2) ===", shooting2);
        out << "iter = " << shooting2.iterations << "\n";
        out << "Невязка граничного условия = " << shooting2.residual << "\n\n";

        writeTable(out, "=== Конечно-разностный метод (h1) ===", diff1);
        out << "Невязка граничных условий = " << diff1.residual << "\n\n";

        writeTable(out, "=== Конечно-разностный метод (h2) ===", diff2);
        out << "Невязка граничных условий = " << diff2.residual << "\n\n";

        out << "=== Оценка погрешности ===\n\n";
        out << "Метод стрельбы:\n";
        out << "max |y_h - y_exact|     = " << maxAbsErrorToExact(shooting1.x, shooting1.y) << "\n";
        out << "max |y_h/2 - y_exact|   = " << maxAbsErrorToExact(shooting2.x, shooting2.y) << "\n";
        out << "Runge-Romberg estimate  = "
            << maxAbsErrorRungeRomberg(shooting1.x, shooting1.y, shooting2.x, shooting2.y, 4) << "\n\n";

        out << "Конечно-разностный метод:\n";
        out << "max |y_h - y_exact|     = " << maxAbsErrorToExact(diff1.x, diff1.y) << "\n";
        out << "max |y_h/2 - y_exact|   = " << maxAbsErrorToExact(diff2.x, diff2.y) << "\n";
        out << "Runge-Romberg estimate  = "
            << maxAbsErrorRungeRomberg(diff1.x, diff1.y, diff2.x, diff2.y, 1) << "\n";

        std::cout << "Алгоритм 4.2 (краевая задача для ОДУ) завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";
    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << "\n";
        std::cerr << e.what() << "\n";
    }
}
