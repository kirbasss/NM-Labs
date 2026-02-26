#include "LU.hpp"
#include "Matrix.hpp"
#include <fstream>
#include <cmath>
#include <filesystem>
#include "utils.hpp"

const double EPS = 1e-12;

LUDecomposition luDecompose(const Matrix& A_orig, std::ostream& log) {
    int n = A_orig.size();
    Matrix A = A_orig;                     // рабочая копия
    Matrix L(n), U(n);
    std::vector<int> P(n);

    for (int i = 0; i < n; ++i) {
        P[i] = i;
        L(i, i) = 1.0;
    }

    log << "=== LU-разложение с выбором главного элемента ===\n";

    for (int k = 0; k < n; ++k) {
        // поиск главного элемента
        int pivot = k;
        double max_val = std::fabs(A(k, k));
        for (int i = k + 1; i < n; ++i) {
            double val = std::fabs(A(i, k));
            if (val > max_val) {
                max_val = val;
                pivot = i;
            }
        }

        log << "Шаг k=" << k << ": максимум в строке " << pivot 
            << " (значение " << max_val << ")\n";

        if (max_val < EPS) {
            log << "ОШИБКА: матрица вырожденная!\n";
            throw std::runtime_error("Вырожденная матрица");
        }

        // перестановка
        if (pivot != k) {
            A.swapRows(k, pivot);
            std::swap(P[k], P[pivot]);
            log << "   - переставили строки " << k << " и " << pivot << "\n";
            for (int j = 0; j < k; ++j) {          // только столбцы 0..k-1
                std::swap(L(k, j), L(pivot, j));
            }
            if (k > 0) {
                log << "   - также поменяли множители в L (строки " 
                    << k << " <-> " << pivot << ")\n";
            }
        }

        // зануление
        for (int i = k + 1; i < n; ++i) {
            double factor = A(i, k) / A(k, k);
            L(i, k) = factor;
            A(i, k) = 0.0;
            for (int j = k + 1; j < n; ++j) {
                A(i, j) -= factor * A(k, j);
            }
            log << "   L[" << i << "][" << k << "] = " << factor << "\n";
        }
    }

    // заполняем U
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j)
            U(i, j) = A(i, j);

    log << "LU-разложение завершено успешно\n";
    return {L, U, P};
}

std::vector<double> solveLU(const LUDecomposition& lu, const std::vector<double>& b) {
    int n = lu.L.size();
    std::vector<double> bp(n);
    for (int i = 0; i < n; ++i) bp[i] = b[lu.P[i]];   // Pb

    // прямой ход L y = Pb
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) sum += lu.L(i, j) * y[j];
        y[i] = bp[i] - sum;
    }

    // обратный ход U x = y
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += lu.U(i, j) * x[j];
        x[i] = (y[i] - sum) / lu.U(i, i);
    }
    return x;
}

double computeDeterminant(const LUDecomposition& lu) {
    int n = lu.U.size();
    double det = 1.0;
    for (int i = 0; i < n; ++i) det *= lu.U(i, i);

    // подсчёт перестановок (чётность)
    int swaps = 0;
    std::vector<bool> visited(n, false);
    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            int j = i;
            while (!visited[j]) {
                visited[j] = true;
                j = lu.P[j];
                if (j != i) ++swaps;
            }
        }
    }
    if (swaps % 2 == 1) det = -det;
    return det;
}

Matrix computeInverse(const LUDecomposition& lu) {
    int n = lu.L.size();
    Matrix inv(n);
    std::vector<double> e(n, 0.0);

    for (int j = 0; j < n; ++j) {
        std::fill(e.begin(), e.end(), 0.0);
        e[j] = 1.0;
        auto col = solveLU(lu, e);
        for (int i = 0; i < n; ++i) inv(i, j) = col[i];
    }
    return inv;
}

void run_1_1(const std::string& inputFile) {
    namespace fs = std::filesystem;
    std::string base = fs::path(inputFile).stem().string();   // "1.1"
    std::string outFile  = "output/" + base + ".txt";
    std::string logFile  = "output/" + base + "_log.txt";

    fs::create_directories("output");

    std::ofstream log(logFile);
    if (!log) {
        std::cerr << "Не удалось создать лог " << logFile << std::endl;
        return;
    }

    try {
        auto sys = loadSystemFromFile(inputFile);
        log << "Загружена матрица " << sys.A.size() << "×" << sys.A.size() 
            << " из " << inputFile << "\n\n";

        auto lu = luDecompose(sys.A, log);

        auto x = solveLU(lu, sys.b);
        double det = computeDeterminant(lu);
        Matrix inv = computeInverse(lu);

        // Запись результата
        std::ofstream out(outFile);
        out << std::fixed << std::setprecision(8);
        out << "=== Решение системы (x) ===\n";
        for (double v : x) out << v << '\n';

        out << "\n=== Вектор перестановок P ===\n";
        for (size_t i : lu.P) out << i << '\n';

        out << "\n=== Матрица L ===\n";
        lu.L.print(out);

        out << "\n=== Матрица U ===\n";
        lu.U.print(out);

        out << "\n=== Определитель det(A) = " << det << " ===\n";

        out << "\n=== Обратная матрица A^-1 ===\n";
        inv.print(out);

        out << "\n=== Проверка: ||A x - b|| ===\n";
        double residual = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < x.size(); ++j) sum += sys.A(i, j) * x[j];
            residual += std::fabs(sum - sys.b[i]);
        }
        out << residual << " (должно быть близко к 0)\n";

        std::cout << "Алгоритм 1.1 завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << std::endl;
        std::cerr << e.what() << std::endl;
    }
}
