#include "Tridiagonal.hpp"
#include "utils.hpp"
#include "Matrix.hpp"
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <stdexcept>

const double EPS = 1e-12;

TridiagonalSystem loadTridiagonalSystemFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filename);
    }

    size_t n;
    file >> n;
    if (!file) {
        throw std::runtime_error("Ошибка чтения размера системы из файла " + filename);
    }

    TridiagonalSystem sys;
    sys.a.resize(n > 0 ? n - 1 : 0);
    sys.b.resize(n);
    sys.c.resize(n > 0 ? n - 1 : 0);
    sys.d.resize(n);

    for (size_t i = 0; i + 1 < n; ++i) {
        file >> sys.a[i];
    }
    for (size_t i = 0; i < n; ++i) {
        file >> sys.b[i];
    }
    for (size_t i = 0; i + 1 < n; ++i) {
        file >> sys.c[i];
    }
    for (size_t i = 0; i < n; ++i) {
        file >> sys.d[i];
    }

    if (!file) {
        throw std::runtime_error("Ошибка чтения данных трёхдиагональной системы из файла " + filename);
    }

    return sys;
}

std::vector<double> tridiagSolve(const std::vector<double>& a,
                                 const std::vector<double>& b,
                                 const std::vector<double>& c,
                                 const std::vector<double>& d,
                                 std::ostream& log) {
    const size_t n = b.size();
    if (n == 0) {
        return {};
    }

    if (d.size() != n) {
        throw std::invalid_argument("Размер d должен совпадать с размером b");
    }
    if (a.size() != (n > 0 ? n - 1 : 0)) {
        throw std::invalid_argument("Размер a должен быть n-1");
    }
    if (c.size() != (n > 0 ? n - 1 : 0)) {
        throw std::invalid_argument("Размер c должен быть n-1");
    }

    std::vector<double> p(n, 0.0);
    std::vector<double> q(n, 0.0);
    std::vector<double> x(n, 0.0);

    log << "=== Метод прогонки ===\n";
    log << "Решение ищется в виде x[i] = p[i] * x[i+1] + q[i]\n\n";

    if (n == 1) {
        if (std::fabs(b[0]) < EPS) {
            throw std::runtime_error("Нулевой элемент на главной диагонали: b[0] = 0");
        }

        p[0] = 0.0;
        q[0] = d[0] / b[0];
        x[0] = q[0];

        log << "n = 1:\n";
        log << "p[0] = 0\n";
        log << "q[0] = d[0] / b[0] = " << q[0] << "\n";
        log << "x[0] = q[0] = " << x[0] << "\n";

        return x;
    }

    log << "Прямой ход:\n";

    if (std::fabs(b[0]) < EPS) {
        throw std::runtime_error("Нулевой знаменатель на первом шаге: b[0] = 0");
    }

    p[0] = -c[0] / b[0];
    q[0] =  d[0] / b[0];

    log << "i = 0: "
        << "p[0] = " << p[0]
        << ", q[0] = " << q[0] << "\n";

    for (size_t i = 1; i < n; ++i) {
        double denom = b[i] + a[i - 1] * p[i - 1];
        if (std::fabs(denom) < EPS) {
            throw std::runtime_error("Нулевой знаменатель на шаге i = " + std::to_string(i));
        }

        if (i < n - 1) {
            p[i] = -c[i] / denom;
        } else {
            p[i] = 0.0; // p[n-1] = 0
        }

        q[i] = (d[i] - a[i - 1] * q[i - 1]) / denom;

        log << "i = " << i << ": "
            << "denom = b[" << i << "] + a[" << (i - 1) << "] * p[" << (i - 1) << "] = " << denom
            << ", p[" << i << "] = " << p[i]
            << ", q[" << i << "] = " << q[i] << "\n";
    }

    log << "\nОбратный ход:\n";

    x[n - 1] = q[n - 1];
    log << "x[" << (n - 1) << "] = q[" << (n - 1) << "] = " << x[n - 1] << "\n";

    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
        log << "x[" << i << "] = p[" << i << "] * x[" << (i + 1) << "] + q[" << i << "] = " << x[i] << "\n";
    }

    log << "\nМетод прогонки завершён успешно\n";
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
        auto sys = loadTridiagonalSystemFromFile(inputFile);
        const size_t n = sys.b.size();

        log << "Загружена трёхдиагональная система n = " << n 
            << " из " << inputFile << "\n\n";

        auto x = tridiagSolve(sys.a, sys.b, sys.c, sys.d, log);
        double residual = computeTridiagonalResidual(sys.a, sys.b, sys.c, sys.d, x);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile);
        }

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