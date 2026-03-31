#include "Jacobi.hpp"
#include "Utils.hpp"
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <vector>

const int MAX_ITER = 10000;
extern const double EPS_ZERO;

#ifndef M_PI
constexpr float M_PI = 3.14159265358979323846f;
#endif

EigenInput loadEigenInputFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filename);
    }

    size_t n;
    file >> n;
    if (!file) {
        throw std::runtime_error("Ошибка чтения размера матрицы");
    }

    Matrix A(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            file >> A(i, j);
        }
    }

    double eps;
    file >> eps;
    if (!file) {
        throw std::runtime_error("Ошибка чтения eps");
    }
    if (eps <= 0.0) {
        throw std::runtime_error("eps должно быть положительным");
    }

    return {A, eps};
}

EigenResult jacobiEigen(const Matrix& A_orig, double eps, std::ostream& log) {
    size_t n = A_orig.size();

    // Проверка симметричности
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (std::fabs(A_orig(i, j) - A_orig(j, i)) > EPS_ZERO) {
                throw std::runtime_error("Метод вращений Якоби применим только к симметрическим матрицам");
            }
        }
    }

    Matrix A = A_orig;
    Matrix U(n); // накопленная матрица собственных векторов
    for (size_t i = 0; i < n; ++i) { U(i, i) = 1.0; }

    log << "=== Метод вращений (Якоби для собственных значений) ===\n";
    log << "Точность eps = " << eps << ", max итераций = " << MAX_ITER << "\n";

    std::vector<std::pair<double, int>> error_history;

    int iter = 0;
    while (iter < MAX_ITER) {
        // Критерий остановки:
        // t(A^(k)) = sqrt(sum_{l<m} (a_lm^(k))^2)
        double t = 0.0;
        for (size_t l = 0; l < n; ++l) {
            for (size_t m = l + 1; m < n; ++m) {
                t += A(l, m) * A(l, m);
            }
        }
        t = std::sqrt(t);

        log << "Итерация " << iter << ": t(A) = " << t << "\n";

        if (t < eps) {
            log << "Сошлось за " << iter << " итераций\n";
            break;
        }

        // Выбор максимального по модулю внедиагонального элемента a_lm
        size_t l = 0, m = 1;
        double max_val = std::fabs(A(0, 1));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                double val = std::fabs(A(i, j));
                if (val > max_val) {
                    max_val = val;
                    l = i;
                    m = j;
                }
            }
        }

        const double all = A(l, l);
        const double amm = A(m, m);
        const double alm = A(l, m);

        // Угол вращения
        double phi;
        if (std::fabs(all - amm) < EPS_ZERO) {
            phi = M_PI / 4.0;
        } else {
            phi = 0.5 * std::atan((2.0 * alm) / (all - amm));
        }

        const double c = std::cos(phi);
        const double s = std::sin(phi);

        log << "  max |a_lm| = " << max_val
            << " при l = " << l << ", m = " << m << "\n";
        log << "  phi = " << phi << "\n";
        
        // Обновление матрицы A
        for (size_t k = 0; k < n; ++k) {
            if (k == l || k == m) continue;

            const double akl = A(k, l);
            const double akm = A(k, m);

            const double new_kl = c * akl + s * akm;
            const double new_km = -s * akl + c * akm;

            A(k, l) = new_kl;
            A(l, k) = new_kl;
            A(k, m) = new_km;
            A(m, k) = new_km;
        }

        A(l, l) = c * c * all + 2.0 * c * s * alm + s * s * amm;
        A(m, m) = s * s * all - 2.0 * c * s * alm + c * c * amm;
        A(l, m) = 0.0;
        A(m, l) = 0.0;

        // Обновление матрицы собственных векторов
        for (size_t k = 0; k < n; ++k) {
            const double ukl = U(k, l);
            const double ukm = U(k, m);

            U(k, l) = c * ukl + s * ukm;
            U(k, m) = -s * ukl + c * ukm;
        }

        ++iter;
    }

    if (iter == MAX_ITER) {
        log << "ВНИМАНИЕ: не сошлось за " << MAX_ITER << " итераций!\n";
    }

    // Собственные значения — диагональ A
    std::vector<double> eigenvalues(n);
    for (size_t i = 0; i < n; ++i) {
        eigenvalues[i] = A(i, i);
    }

    return {eigenvalues, U};
}

void run_1_4(const std::string& inputFile) {
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
        auto input = loadEigenInputFromFile(inputFile);
        const Matrix& A = input.A;
        const double eps = input.eps;
        const size_t n = A.size();

        log << "Загружена матрица " << n << "x" << n
            << ", eps = " << eps << " из " << inputFile << "\n\n";

        auto result = jacobiEigen(A, eps, log);

        Matrix Lambda(n);
        for (size_t i = 0; i < n; ++i) {
            Lambda(i, i) = result.eigenvalues[i];
        }

        Matrix VT = transpose(result.eigenvectors);
        Matrix recon = multiply(multiply(result.eigenvectors, Lambda), VT);
        double reconError = matrixDiffNorm1(A, recon);

        std::ofstream out(outFile);
        out << std::fixed << std::setprecision(8);

        out << "=== Собственные значения ===\n";
        for (double lambda : result.eigenvalues) out << lambda << '\n';

        out << "\n=== Собственные векторы (столбцы) ===\n";
        result.eigenvectors.print(out);

        out << "\n=== Проверка реконструкции ===\n";
        out << "||A - U * Lambda * U^T||_1 = " << reconError << "\n";

        std::cout << "Алгоритм 1.4 (метод вращений) завершён. Результаты в " << outFile << '\n';
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << std::endl;
        std::cerr << e.what() << std::endl;
    }
}