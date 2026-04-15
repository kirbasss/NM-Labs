#include "Integration.hpp"

#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace {
    const double EPS = 1e-12;

    bool isIntegerWithinTolerance(double value, double tol = 1e-9) {
        return std::fabs(value - std::round(value)) < tol;
    }

    int computeStepCount(double a, double b, double h) {
        if (h <= 0.0) {
            throw std::invalid_argument("Шаг h должен быть положительным");
        }
        if (b <= a) {
            throw std::invalid_argument("Требуется xk > x0");
        }

        double n = (b - a) / h;
        if (!isIntegerWithinTolerance(n)) {
            throw std::invalid_argument("Длина интервала должна делиться на шаг без остатка");
        }

        int steps = static_cast<int>(std::round(n));
        if (steps <= 0) {
            throw std::invalid_argument("Число шагов должно быть положительным");
        }
        return steps;
    }

    double exactIntegral(double a, double b) {
        // ∫ x/(2x+5) dx = x/2 - (5/4) ln|2x+5| + C
        auto F = [](double x) {
            return x / 2.0 - 1.25 * std::log(std::fabs(2.0 * x + 5.0));
        };
        return F(b) - F(a);
    }
}

IntegrationInput loadIntegrationInputFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filename);
    }

    IntegrationInput input;
    file >> input.x0 >> input.xk;
    file >> input.h1 >> input.h2;

    if (!file) {
        throw std::runtime_error("Ошибка чтения параметров интегрирования из файла " + filename);
    }

    if (input.h1 <= 0.0 || input.h2 <= 0.0) {
        throw std::invalid_argument("Шаги должны быть положительными");
    }
    if (input.xk <= input.x0) {
        throw std::invalid_argument("Правая граница должна быть больше левой");
    }

    return input;
}

double integrand(double x) {
    return x / (2.0 * x + 5.0);
}

double integrateRectangles(const std::function<double(double)>& f,
                           double a,
                           double b,
                           double h,
                           std::ostream& log) {
    int n = computeStepCount(a, b, h);

    log << "=== Метод прямоугольников, h = " << h << " ===\n";

    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double xl = a + i * h;
        double xr = xl + h;
        double xm = 0.5 * (xl + xr);
        double fxm = f(xm);
        double term = h * fxm;
        sum += term;

        log << "i = " << i
            << ", [" << xl << ", " << xr << "]"
            << ", mid = " << xm
            << ", f(mid) = " << fxm
            << ", term = " << term << "\n";
    }

    log << "Итог: I_rect = " << sum << "\n\n";
    return sum;
}

double integrateTrapezoids(const std::function<double(double)>& f,
                           double a,
                           double b,
                           double h,
                           std::ostream& log) {
    int n = computeStepCount(a, b, h);

    log << "=== Метод трапеций, h = " << h << " ===\n";

    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        double xi = a + i * h;
        double fi = f(xi);
        double coeff = (i == 0 || i == n) ? 1.0 : 2.0;
        sum += coeff * fi;

        log << "i = " << i
            << ", x = " << xi
            << ", f(x) = " << fi
            << ", coeff = " << coeff << "\n";
    }

    double result = 0.5 * h * sum;
    log << "Итог: I_trap = h/2 * sum = " << result << "\n\n";
    return result;
}

double integrateSimpson(const std::function<double(double)>& f,
                        double a,
                        double b,
                        double h,
                        std::ostream& log) {
    int n = computeStepCount(a, b, h);
    if (n % 2 != 0) {
        throw std::invalid_argument("Для метода Симпсона число интервалов должно быть чётным");
    }

    log << "=== Метод Симпсона, h = " << h << " ===\n";

    log << "i = 0, x = " << a << ", f(x) = " << f(a) << ", coeff = 1\n";

    double sum = 0.0;
    int pairCount = n / 2;

    for (int k = 0; k < pairCount; ++k) {
        double x0 = a + (2 * k) * h;
        double x1 = x0 + h;
        double x2 = x0 + 2 * h;

        double f0 = f(x0);
        double f1 = f(x1);
        double f2 = f(x2);

        double term = f0 + 4.0 * f1 + f2;
        sum += term;

        log << "k = " << k
            << ": [" << x0 << ", " << x2 << "], midpoint = " << x1 << "\n";
        log << "term = (" << h << "/3) * ("
            << f0 << " + 4*" << f1 << " + " << f2
            << ") = " << term << "\n";
    }

    double result = (h / 3.0) * sum;
    log << "Итог: I_simp = h/3 * sum = " << result << "\n\n";
    return result;
}

double rungeRomberg(double FhCoarse,
                    double FhFine,
                    double hCoarse,
                    double hFine,
                    int p) {
    if (hCoarse <= 0.0 || hFine <= 0.0) {
        throw std::invalid_argument("Шаги должны быть положительными");
    }
    if (p <= 0) {
        throw std::invalid_argument("Порядок точности p должен быть положительным");
    }

    double k = hCoarse / hFine;
    if (k <= 1.0 + EPS) {
        throw std::invalid_argument("Для Рунге-Ромберга требуется hCoarse > hFine");
    }

    return FhFine + (FhFine - FhCoarse) / (std::pow(k, p) - 1.0);
}

double estimateRectError(double M2, double a, double b, double h) {
    return M2 * (b - a) * h * h / 24.0;
}

double estimateTrapError(double M2, double a, double b, double h) {
    return M2 * (b - a) * h * h / 12.0;
}

double estimateSimpsonError(double M4, double a, double b, double h) {
    return M4 * (b - a) * std::pow(h, 4) / 180.0;
}

void run_3_5(const std::string& inputFile) {
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
        auto input = loadIntegrationInputFromFile(inputFile);

        log << std::fixed << std::setprecision(10);
        log << "Задача 3.5. Численное интегрирование\n";
        log << "f(x) = x / (2x + 5)\n";
        log << "x0 = " << input.x0 << ", xk = " << input.xk << "\n";
        log << "h1 = " << input.h1 << ", h2 = " << input.h2 << "\n\n";

        auto f = [](double x) { return integrand(x); };

        double rect_h1 = integrateRectangles(f, input.x0, input.xk, input.h1, log);
        double trap_h1 = integrateTrapezoids(f, input.x0, input.xk, input.h1, log);
        double simp_h1 = integrateSimpson(f, input.x0, input.xk, input.h1, log);

        double rect_h2 = integrateRectangles(f, input.x0, input.xk, input.h2, log);
        double trap_h2 = integrateTrapezoids(f, input.x0, input.xk, input.h2, log);
        double simp_h2 = integrateSimpson(f, input.x0, input.xk, input.h2, log);

        double rect_rr = rungeRomberg(rect_h1, rect_h2, input.h1, input.h2, 2);
        double trap_rr = rungeRomberg(trap_h1, trap_h2, input.h1, input.h2, 2);
        double simp_rr = rungeRomberg(simp_h1, simp_h2, input.h1, input.h2, 4);

        double exact = exactIntegral(input.x0, input.xk);

        double err_rect = std::fabs(exact - rect_rr);
        double err_trap = std::fabs(exact - trap_rr);
        double err_simp = std::fabs(exact - simp_rr);

        // Теоретические оценки
        double M2 = 20.0 / 27.0;
        double M4 = 960.0 / 243.0;

        double rect_err_est_h1 = estimateRectError(M2, input.x0, input.xk, input.h1);
        double trap_err_est_h1 = estimateTrapError(M2, input.x0, input.xk, input.h1);
        double simp_err_est_h1 = estimateSimpsonError(M4, input.x0, input.xk, input.h1);

        double rect_err_est_h2 = estimateRectError(M2, input.x0, input.xk, input.h2);
        double trap_err_est_h2 = estimateTrapError(M2, input.x0, input.xk, input.h2);
        double simp_err_est_h2 = estimateSimpsonError(M4, input.x0, input.xk, input.h2);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile);
        }

        out << std::fixed << std::setprecision(10);

        out << "=== Задача 3.5. Численное интегрирование ===\n\n";
        out << "f(x) = x / (2x + 5)\n";
        out << "Интервал: [" << input.x0 << ", " << input.xk << "]\n";
        out << "h1 = " << input.h1 << "\n";
        out << "h2 = " << input.h2 << "\n\n";
        
        out << "=== Результаты для h1 ===\n";
        out << "Метод прямоугольников: " << rect_h1 << "\n";
        out << "Метод трапеций:       " << trap_h1 << "\n";
        out << "Метод Симпсона:       " << simp_h1 << "\n\n";

        out << "=== Результаты для h2 ===\n";
        out << "Метод прямоугольников: " << rect_h2 << "\n";
        out << "Метод трапеций:       " << trap_h2 << "\n";
        out << "Метод Симпсона:       " << simp_h2 << "\n\n";
        
        out << "=== Теоретические оценки погрешностей ===\n";
        out << "Метод прямоугольников: " << "для h1 " << rect_err_est_h1 << ", для h2 " << rect_err_est_h2 << '\n';
        out << "Метод трапеций: " << "для h1 " << rect_err_est_h1 << ", для h2 " << rect_err_est_h2 << '\n';
        out << "Метод Симпсона: " << "для h1 " << rect_err_est_h1 << ", для h2 " << rect_err_est_h2 << "\n\n";

        out << "=== Уточнение по Рунге-Ромбергу ===\n";
        out << "Прямоугольники (p=2): " << rect_rr << "\n";
        out << "Трапеции       (p=2): " << trap_rr << "\n";
        out << "Симпсон        (p=4): " << simp_rr << "\n\n";

        out << "=== Точное значение интеграла ===\n";
        out << exact << "\n\n";

        out << "=== Абсолютная погрешность уточнённых значений ===\n";
        out << "Прямоугольники: " << err_rect << "\n";
        out << "Трапеции:       " << err_trap << "\n";
        out << "Симпсон:        " << err_simp << "\n";

        std::cout << "Алгоритм 3.5 (численное интегрирование) завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << "\n";
        std::cerr << e.what() << std::endl;
    }
}