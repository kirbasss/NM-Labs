#include "NonlinearSystem.hpp"
#include "LU.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

const int MAX_ITER = 10000;
extern const double EPS_ZERO;

// Система:
// (x1^2 + 4) x2 - 8 = 0
// (x1 - 1)^2 + (x2 - 1)^2 - 4 = 0

double f1(double x1, double x2) {
    return (x1 * x1 + 4.0) * x2 - 8.0;
}

double f2(double x1, double x2) {
    return (x1 - 1.0) * (x1 - 1.0) + (x2 - 1.0) * (x2 - 1.0) - 4.0;
}

std::vector<double> F_system(const std::vector<double>& x) {
    if (x.size() != 2) {
        throw std::invalid_argument("F_system: x must have size 2");
    }
    return {
        f1(x[0], x[1]),
        f2(x[0], x[1])
    };
}

Matrix J_system(const std::vector<double>& x) {
    if (x.size() != 2) {
        throw std::invalid_argument("J_system: x must have size 2");
    }

    Matrix J(2);
    J(0, 0) = 2.0 * x[0] * x[1];
    J(0, 1) = x[0] * x[0] + 4.0;
    J(1, 0) = 2.0 * (x[0] - 1.0);
    J(1, 1) = 2.0 * (x[1] - 1.0);
    return J;
}

double residualNormC(const VectorFunction& F, const std::vector<double>& x) {
    return vectorNormC(F(x));
}

// МПИ
// Ищем положительное решение, поэтому берём ветвь:
// x1 = 1 + sqrt(4 - (x2 - 1)^2)
// x2 = 8 / (x1^2 + 4)

double phi1(double x1, double x2) {
    (void)x1;
    double inside = 4.0 - (x2 - 1.0) * (x2 - 1.0);
    if (inside < 0.0) {
        throw std::runtime_error("phi1: подкоренное выражение отрицательно");
    }
    return 1.0 + std::sqrt(inside);
}

double phi2(double x1, double x2) {
    return 8.0 / (x1 * x1 + 4.0);
}

// Частные производные phi
double dphi1_dx1(double x1, double x2) {
    (void)x1;
    (void)x2;
    return 0.0;
}

double dphi1_dx2(double x1, double x2) {
    (void)x1;
    double inside = 4.0 - (x2 - 1.0) * (x2 - 1.0);
    if (inside <= 0.0) {
        throw std::runtime_error("dphi1_dx2: подкоренное выражение неположительно");
    }
    return -(x2 - 1.0) / std::sqrt(inside);
}

double dphi2_dx1(double x1, double x2) {
    (void)x2;
    double denom = x1 * x1 + 4.0;
    return -16.0 * x1 / (denom * denom);
}

double dphi2_dx2(double x1, double x2) {
    (void)x1;
    (void)x2;
    return 0.0;
}

// max-row-sum норма матрицы Якоби phi на прямоугольнике
double estimateQ(double a1, double b1, double a2, double b2) {
    const int N = 200;
    double q = 0.0;

    for (int i = 0; i <= N; ++i) {
        double x1 = a1 + (b1 - a1) * i / N;
        for (int j = 0; j <= N; ++j) {
            double x2 = a2 + (b2 - a2) * j / N;

            double row1 = 0.0;
            double row2 = 0.0;

            try {
                row1 = std::fabs(dphi1_dx1(x1, x2)) + std::fabs(dphi1_dx2(x1, x2));
                row2 = std::fabs(dphi2_dx1(x1, x2)) + std::fabs(dphi2_dx2(x1, x2));
            } catch (...) {
                return std::numeric_limits<double>::infinity();
            }

            q = std::max(q, std::max(row1, row2));
        }
    }

    return q;
}

bool mapsRectangleToItself(double a1, double b1, double a2, double b2) {
    const int N = 200;

    for (int i = 0; i <= N; ++i) {
        double x1 = a1 + (b1 - a1) * i / N;
        for (int j = 0; j <= N; ++j) {
            double x2 = a2 + (b2 - a2) * j / N;

            double y1, y2;
            try {
                y1 = phi1(x1, x2);
                y2 = phi2(x1, x2);
            } catch (...) {
                return false;
            }

            if (y1 < a1 || y1 > b1 || y2 < a2 || y2 > b2) {
                return false;
            }
        }
    }

    return true;
}

// МПИ
SystemResult solveSimpleIterationSystem(double a1, double b1,
                                        double a2, double b2,
                                        const std::vector<double>& x0,
                                        double eps,
                                        std::ostream& log) {
    if (x0.size() != 2) {
        throw std::invalid_argument("solveSimpleIterationSystem: x0 must have size 2");
    }

    log << "--- Метод простой итерации для системы ---\n";

    double q = estimateQ(a1, b1, a2, b2);
    bool invariant = mapsRectangleToItself(a1, b1, a2, b2);

    log << "Оценка q = " << q << "\n";
    log << "phi(G) subset G: " << (invariant ? "да" : "нет") << "\n";

    if (!std::isfinite(q)) {
        log << "phi или производные phi не определены на всей области\n";
        return {{std::numeric_limits<double>::quiet_NaN(),
                 std::numeric_limits<double>::quiet_NaN()}, 0};
    }

    if (!invariant) {
        log << "ВНИМАНИЕ: phi не отображает область в себя\n";
        return {{std::numeric_limits<double>::quiet_NaN(),
                 std::numeric_limits<double>::quiet_NaN()}, 0};
    }

    if (q >= 1.0) {
        log << "ВНИМАНИЕ: достаточное условие сходимости q < 1 не выполнено\n";
        return {{std::numeric_limits<double>::quiet_NaN(),
                 std::numeric_limits<double>::quiet_NaN()}, 0};
    }

    std::vector<double> x = x0;
    log << "Начальное приближение x^(0) = (" << x[0] << ", " << x[1] << ")\n";

    int iter = 0;
    while (iter < MAX_ITER) {
        std::vector<double> x_new(2);
        x_new[0] = phi1(x[0], x[1]);
        x_new[1] = phi2(x[0], x[1]);

        auto diff = subtractVectors(x_new, x);
        double diffNorm = vectorNormC(diff);
        double est = q / (1.0 - q) * diffNorm;

        log << "k=" << iter
            << " x=(" << x_new[0] << ", " << x_new[1] << ")"
            << " F=(" << f1(x_new[0], x_new[1]) << ", " << f2(x_new[0], x_new[1]) << ")"
            << " ||dx||_c=" << diffNorm
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

// МН
// J(x^(k)) * delta^(k) = -F(x^(k))
// x^(k+1) = x^(k) + delta^(k)
SystemResult solveNewtonSystem(const std::vector<double>& x0,
                               double eps,
                               const VectorFunction& F,
                               const JacobianFunction& J,
                               std::ostream& log) {
    if (x0.empty()) {
        throw std::invalid_argument("solveNewtonSystem: x0 must be non-empty");
    }

    std::vector<double> x = x0;
    log << "--- Метод Ньютона для системы ---\n";
    log << "Начальное приближение x^(0) = (";
    for (size_t i = 0; i < x.size(); ++i) {
        log << x[i] << (i + 1 == x.size() ? "" : ", ");
    }
    log << ")\n";

    int iter = 0;
    while (iter < MAX_ITER) {
        auto Fx = F(x);
        Matrix Jx = J(x);

        if (Fx.size() != x.size() || Jx.size() != x.size()) {
            throw std::runtime_error("solveNewtonSystem: inconsistent sizes of F or J");
        }

        std::vector<double> rhs = Fx;
        for (double& v : rhs) {
            v = -v;
        }

        auto lu = luDecompose(Jx, log);
        auto delta = solveLU(lu, rhs);
        auto x_new = addVectors(x, delta);

        double diffNorm = vectorNormC(delta);

        log << "k=" << iter << " x=(";
        for (size_t i = 0; i < x_new.size(); ++i) {
            log << x_new[i] << (i + 1 == x_new.size() ? "" : ", ");
        }
        log << ") F=(";
        auto Fnew = F(x_new);
        for (size_t i = 0; i < Fnew.size(); ++i) {
            log << Fnew[i] << (i + 1 == Fnew.size() ? "" : ", ");
        }
        log << ") ||delta||_c=" << diffNorm << "\n";

        if (diffNorm < eps) {
            return {x_new, iter + 1};
        }

        x = x_new;
        ++iter;
    }

    log << "Не сошлось\n";
    return {x, MAX_ITER};
}

// Формат входного файла:
// eps
// a1 b1
// a2 b2
// x10 x20

void run_2_2(const std::string& inputFile) {
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
    in >> eps;

    double a1, b1, a2, b2;
    in >> a1 >> b1;
    in >> a2 >> b2;

    std::vector<double> x0(2);
    in >> x0[0] >> x0[1];

    std::ofstream out(outFile);
    out << std::fixed << std::setprecision(8);

    log << "Прямоугольник G: [" << a1 << ", " << b1 << "] x ["
        << a2 << ", " << b2 << "]\n";
    log << "eps = " << eps << "\n\n";

    try {
        auto res_iter = solveSimpleIterationSystem(a1, b1, a2, b2, x0, eps, log);

        out << "=== Метод простой итерации ===\n";
        out << "x1 = " << res_iter.x[0] << "\n";
        out << "x2 = " << res_iter.x[1] << "\n";
        out << "iter = " << res_iter.iterations << "\n";

        if (std::isfinite(res_iter.x[0]) && std::isfinite(res_iter.x[1])) {
            out << "||F(x)||_c = "
                << std::max(std::fabs(f1(res_iter.x[0], res_iter.x[1])),
                            std::fabs(f2(res_iter.x[0], res_iter.x[1])))
                << "\n";
        } else {
            out << "||F(x)||_c = nan\n";
        }

        out << "\n";
    } catch (const std::exception& e) {
        log << "Ошибка в методе простой итерации: " << e.what() << "\n";
        out << "=== Метод простой итерации ===\n";
        out << "ошибка: " << e.what() << "\n\n";
    }

    try {
        auto res_newton = solveNewtonSystem(x0, eps, F_system, J_system, log);

        out << "=== Метод Ньютона ===\n";
        out << "x1 = " << res_newton.x[0] << "\n";
        out << "x2 = " << res_newton.x[1] << "\n";
        out << "iter = " << res_newton.iterations << "\n";
        out << "||F(x)||_c = " << residualNormC(F_system, res_newton.x) << "\n";
    } catch (const std::exception& e) {
        log << "Ошибка в методе Ньютона: " << e.what() << "\n";
        out << "=== Метод Ньютона ===\n";
        out << "ошибка: " << e.what() << "\n";
    }

    std::cout << "Алгоритм 2.2 (решение системы нелинейных уравнений) завершён. Результаты в "
              << outFile << "\n";
    std::cout << "Лог: " << logFile << "\n";
}