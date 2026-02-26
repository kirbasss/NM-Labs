#include "Tridiagonal.hpp"
#include "utils.hpp"
#include "Matrix.hpp"
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <stdexcept>

const double EPS = 1e-12;

std::vector<double> tridiagSolve(const std::vector<double>& a,   // нижняя диагональ (n-1)
                                const std::vector<double>& b,   // главная диагональ (n)
                                const std::vector<double>& c,   // верхняя диагональ (n-1)
                                const std::vector<double>& d,   // правая часть (n)
                                std::ostream& log) {
    size_t n = b.size();
    if (n == 0) return {};

    std::vector<double> alpha(n), beta(n), x(n);

    log << "=== Метод прогонки ===\n";
    log << "Прямой ход (вычисление коэффициентов alpha и beta):\n";

    // Специальный случай n == 1
    if (n == 1) {
        if (std::fabs(b[0]) < EPS) throw std::runtime_error("b[0] == 0");
        beta[0] = d[0] / b[0];
        alpha[0] = 0.0;
        log << "  i=0: alpha[0] = 0, beta[0] = " << beta[0] << '\n';
        x[0] = beta[0];
        log << "Обратный ход не нужен (n=1)\n";
        return x;
    }

    // n >= 2
    if (std::fabs(b[0]) < EPS) throw std::runtime_error("b[0] == 0");
    alpha[0] = -c[0] / b[0];
    beta[0]  =  d[0] / b[0];
    log << "  i=0: alpha[0] = " << alpha[0] << ", beta[0] = " << beta[0] << '\n';

    for (size_t i = 1; i < n; ++i) {
        double denom = b[i] + a[i-1] * alpha[i-1];
        if (std::fabs(denom) < EPS)
            throw std::runtime_error("Нулевой знаменатель на шаге i=" + std::to_string(i));

        beta[i] = (d[i] - a[i-1] * beta[i-1]) / denom;

        if (i < n - 1) {
            alpha[i] = -c[i] / denom;
        } else {
            alpha[i] = 0.0; // не используется
        }

        log << "  i=" << i << ": alpha[" << i << "] = " << alpha[i]
            << ", beta[" << i << "] = " << beta[i] << '\n';
    }

    log << "\nОбратный ход (нахождение решения x):\n";
    x[n-1] = beta[n-1];
    log << "  x[" << n-1 << "] = " << x[n-1] << '\n';

    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i+1] + beta[i];
        log << "  x[" << i << "] = " << x[i] << '\n';
    }

    log << "Метод прогонки завершён успешно\n";
    return x;
}

double computeTridiagonalResidual(const std::vector<double>& a,
                                  const std::vector<double>& b,
                                  const std::vector<double>& c,
                                  const std::vector<double>& d,
                                  const std::vector<double>& x) {
    size_t n = x.size();
    double res = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double sum = b[i] * x[i];
        if (i > 0)   sum += a[i-1] * x[i-1];
        if (i < n-1) sum += c[i]   * x[i+1];
        res += std::fabs(sum - d[i]);
    }
    return res;
}

void run_1_2(const std::string& inputFile) {
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
        auto sys = loadSystemFromFile(inputFile);
        size_t n = sys.A.size();

        log << "Загружена трёхдиагональная система n = " << n 
            << " из " << inputFile << "\n\n";

        // === Извлекаем три диагонали ===
        std::vector<double> sub(n > 0 ? n-1 : 0);
        std::vector<double> diag(n);
        std::vector<double> super(n > 0 ? n-1 : 0);
        for (size_t i = 0; i < n; ++i) {
            diag[i] = sys.A(i, i);
            if (i > 0)   sub[i-1] = sys.A(i, i-1);
            if (i < n-1) super[i] = sys.A(i, i+1);
        }

        auto x = tridiagSolve(sub, diag, super, sys.b, log);

        double residual = computeTridiagonalResidual(sub, diag, super, sys.b, x);

        std::ofstream out(outFile);
        out << std::fixed << std::setprecision(8);
        out << "=== Решение системы x ===\n";
        for (double v : x) out << v << "\n";

        out << "\n=== Проверка ||A x - b|| ===\n";
        out << residual << " (должно быть близко к 0)\n";

        std::cout << "Алгоритм 1.2 (метод прогонки) завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << std::endl;
        std::cerr << e.what() << std::endl;
    }
}