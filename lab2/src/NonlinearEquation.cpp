#include "NonlinearEquation.hpp"
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <stdexcept>
#include <iostream>

const int MAX_ITER = 10000;
const double EPS_ZERO = 1e-12;

// Исходная функция
double f(double x) {
    return std::pow(2.0, x) - x * x - 0.5;
}

double df(double x) {
    return std::pow(2.0, x) * std::log(2.0) - 2.0 * x;
}

double ddf(double x) {
    return std::pow(2.0, x) * std::log(2.0) * std::log(2.0) - 2.0;
}

// phi(x) для простой итерации
enum class PhiType {
    NegativeRoot,
    PositiveMiddleRoot,
    PositiveLargeRoot
};

PhiType choosePhi(double a, double b) {
    if (b <= 0.0) return PhiType::NegativeRoot;
    if (a >= 3.0) return PhiType::PositiveLargeRoot;
    return PhiType::PositiveMiddleRoot;
}

double phi(double x, PhiType type) {
    if (type == PhiType::NegativeRoot) {
        double val = std::pow(2.0, x) - 0.5;
        if (val < 0.0) throw std::runtime_error("phi not defined");
        return -std::sqrt(val);
    }

    if (type == PhiType::PositiveMiddleRoot) {
        double val = std::pow(2.0, x) - 0.5;
        if (val < 0.0) throw std::runtime_error("phi not defined");
        return std::sqrt(val);
    }

    // PhiType::PositiveLargeRoot
    return std::log2(x * x + 0.5);
}


double dphi(double x, PhiType type) {
    if (type == PhiType::NegativeRoot) {
        double val = std::pow(2.0, x) - 0.5;
        if (val <= 0.0) throw std::runtime_error("dphi not defined");
        return -std::pow(2.0, x) * std::log(2.0) / (2.0 * std::sqrt(val));
    }

    if (type == PhiType::PositiveMiddleRoot) {
        double val = std::pow(2.0, x) - 0.5;
        if (val <= 0.0) throw std::runtime_error("dphi not defined");
        return std::pow(2.0, x) * std::log(2.0) / (2.0 * std::sqrt(val));
    }

    // PhiType::PositiveLargeRoot
    return (2.0 * x) / ((x * x + 0.5) * std::log(2.0));
}


bool phiValue(double x, PhiType type, double& y) {
    try {
        y = phi(x, type);
        return true;
    } catch (...) {
        return false;
    }
}

bool dphiValue(double x, PhiType type, double& y) {
    try {
        y = dphi(x, type);
        return true;
    } catch (...) {
        return false;
    }
}

bool mapsIntervalToItself(double a, double b, PhiType type) {
    const int N = 1000;
    for (int i = 0; i <= N; ++i) {
        double x = a + (b - a) * i / N;
        double y;
        if (!phiValue(x, type, y)) {
            return false;
        }
        if (y < a || y > b) {
            return false;
        }
    }
    return true;
}

// Оценка q на интервале
bool estimateQ(double a, double b, double& q, PhiType type) {
    const int N = 1000;
    q = 0.0;

    for (int i = 0; i <= N; ++i) {
        double x = a + (b - a) * i / N;
        double val;
        if (!dphiValue(x, type, val)) {
            return false;
        }
        val = std::fabs(val);
        if (val > q) q = val;
    }
    return true;
}

// МН
RootResult solveNewton(double a, double b, double eps, std::ostream& log) {
    log << "--- Метод Ньютона ---\n";

    double x;

    if (f(a) * ddf(a) > 0) x = a;
    else if (f(b) * ddf(b) > 0) x = b;
    else {
        log << "Не удалось выбрать начальное приближение строго по условию\n";
        x = (a + b) / 2.0;
    }

    log << "Начальное приближение x0 = " << x << "\n";

    int iter = 0;

    while (iter < MAX_ITER) {
        double fx = f(x);
        double dfx = df(x);

        if (std::fabs(dfx) < EPS_ZERO)
            throw std::runtime_error("df(x) ~ 0");

        double x_new = x - fx / dfx;

        double diff = std::fabs(x_new - x);

        log << "k=" << iter
            << " x=" << x_new
            << " f(x)=" << f(x_new)
            << " |dx|=" << diff << "\n";

        if (diff < eps)
            return {x_new, iter + 1};

        x = x_new;
        iter++;
    }

    log << "Не сошлось\n";
    return {x, MAX_ITER};
}

// МПИ
RootResult solveSimpleIteration(double a, double b, double eps, std::ostream& log) {
    log << "--- Метод простой итерации ---\n";

    PhiType type = choosePhi(a, b);

    log << "Выбранная phi: ";
    if (type == PhiType::NegativeRoot) {
        log << "-sqrt(2^x - 0.5)\n";
    } else if (type == PhiType::PositiveMiddleRoot) {
        log << "sqrt(2^x - 0.5)\n";
    } else {
        log << "log2(x^2 + 0.5)\n";
    }

    double q = 0.0;
    bool q_ok = estimateQ(a, b, q, type);
    bool invariant = mapsIntervalToItself(a, b, type);

    if (!q_ok) {
        log << "phi или phi' не определены на всём интервале\n";
        return {std::numeric_limits<double>::quiet_NaN(), 0};
    }

    log << "Оценка q = " << q << "\n";
    log << "phi([a,b]) subset [a,b]: " << (invariant ? "да" : "нет") << "\n";

    if (!invariant) {
        log << "ВНИМАНИЕ: phi не отображает интервал в себя\n";
        return {std::numeric_limits<double>::quiet_NaN(), 0};
    }

    if (q >= 1.0) {
        log << "ВНИМАНИЕ: условие сходимости |phi'(x)| < 1 не выполнено\n";
        return {std::numeric_limits<double>::quiet_NaN(), 0};
    }

    double x = (a + b) / 2.0;
    log << "Начальное приближение x0 = " << x << "\n";

    int iter = 0;
    while (iter < MAX_ITER) {
        double x_new;
        if (!phiValue(x, type, x_new)) {
            log << "Итерация вышла в область, где phi не определена\n";
            return {std::numeric_limits<double>::quiet_NaN(), iter};
        }

        double diff = std::fabs(x_new - x);
        double est = q / (1.0 - q) * diff;

        log << "k=" << iter
            << " x=" << x_new
            << " f(x)=" << f(x_new)
            << " |dx|=" << diff
            << " est=" << est << "\n";

        if (est < eps) {
            return {x_new, iter + 1};
        }

        x = x_new;
        ++iter;
    }

    log << "Не сошлось\n";
    return {x, MAX_ITER};
}

void run_2_1(const std::string& inputFile) {
    namespace fs = std::filesystem;

    std::string base = fs::path(inputFile).stem().string();
    std::string outFile = "output/" + base + ".txt";
    std::string logFile = "output/" + base + "_log.txt";

    fs::create_directories("output");

    std::ifstream in(inputFile);
    std::ofstream log(logFile);

    if (!in) throw std::runtime_error("Не открыть input");
    if (!log) throw std::runtime_error("Не создать log");

    double eps;
    int m;

    in >> eps >> m;

    std::ofstream out(outFile);
    out << std::fixed << std::setprecision(8);

    for (int i = 0; i < m; ++i) {
        double a, b;
        in >> a >> b;

        log << "\n=== Интервал [" << a << ", " << b << "] ===\n";

        if (f(a) * f(b) >= 0) {
            log << "Нет смены знака, пропуск\n";
            continue;
        }

        out << "=== Интервал [" << a << ", " << b << "] ===\n";

        try {
            auto res_iter = solveSimpleIteration(a, b, eps, log);
            double residual = std::fabs(f(res_iter.root));

            out << "Простая итерация:\n";
            out << "root = " << res_iter.root
                << ", iter = " << res_iter.iterations << "\n";
            out << "|f(root)| = " << residual << "\n";
        } catch (const std::exception& e) {
            log << "Ошибка в методе простой итерации: " << e.what() << "\n";
            out << "Простая итерация:\n";
            out << "ошибка: " << e.what() << "\n";
        }

        try {
            auto res_newt = solveNewton(a, b, eps, log);
            double residual = std::fabs(f(res_newt.root));

            out << "Ньютон:\n";
            out << "root = " << res_newt.root
                << ", iter = " << res_newt.iterations << "\n";
            out << "|f(root)| = " << residual << "\n\n";
        } catch (const std::exception& e) {
            log << "Ошибка в методе Ньютона: " << e.what() << "\n";
            out << "Ньютон:\n";
            out << "ошибка: " << e.what() << "\n\n";
        }
    }

    std::cout << "Алгоритм 2.1 (решение нелинейного уравнения) завершён. Результаты в " << outFile << '\n';
    std::cout << "Лог: " << logFile << "\n";
}