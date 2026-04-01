#include "Iterative.hpp"
#include "utils.hpp"
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <limits>

const int MAX_ITER = 10000;
const double EPS = 1e-12;

// alpha_ii = 0, alpha_ij = -a_ij / a_ii
Matrix buildAlpha(const Matrix& A) {
    size_t n = A.size();
    Matrix alpha(n);

    for (size_t i = 0; i < n; ++i) {
        if (std::fabs(A(i, i)) < EPS) {
            throw std::runtime_error("Нулевой диагональный элемент при построении alpha");
        }

        for (size_t j = 0; j < n; ++j) {
            if (i == j) {
                alpha(i, j) = 0.0;
            } else {
                alpha(i, j) = -A(i, j) / A(i, i);
            }
        }
    }

    return alpha;
}

// beta_i = b_i / a_ii
std::vector<double> buildBeta(const Matrix& A, const std::vector<double>& b) {
    size_t n = A.size();
    std::vector<double> beta(n);

    for (size_t i = 0; i < n; ++i) {
        if (std::fabs(A(i, i)) < EPS) {
            throw std::runtime_error("Нулевой диагональный элемент при построении beta");
        }
        beta[i] = b[i] / A(i, i);
    }

    return beta;
}

// Перестановка строк, чтобы на диагонали не было нулей
void ensureNonZeroDiagonal(Matrix& A, std::vector<double>& b, std::ostream& log) {
    size_t n = A.size();

    for (size_t i = 0; i < n; ++i) {
        if (std::fabs(A(i, i)) >= EPS) {
            continue;
        }

        size_t pivot = i;
        double best = 0.0;

        for (size_t k = i + 1; k < n; ++k) {
            if (std::fabs(A(k, i)) > best) {
                best = std::fabs(A(k, i));
                pivot = k;
            }
        }

        if (pivot == i || best < EPS) {
            throw std::runtime_error(
                "Не удалось перестановкой строк получить ненулевой диагональный элемент на позиции " +
                std::to_string(i)
            );
        }

        A.swapRows(i, pivot);
        std::swap(b[i], b[pivot]);

        log << "Переставлены строки " << i << " и " << pivot
            << " для устранения нуля на диагонали\n";
    }
}

double residualNorm1(const Matrix& A, const std::vector<double>& b, const std::vector<double>& x) {
    size_t n = A.size();
    double res = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            sum += A(i, j) * x[j];
        }
        res += std::fabs(sum - b[i]);
    }

    return res;
}

std::pair<std::vector<double>, int> simpleIteration(
    const Matrix& A, const std::vector<double>& b, double eps,
    std::ostream& log) {

    size_t n = A.size();
    if (b.size() != n) {
        throw std::invalid_argument("simpleIteration: размер b не совпадает с размером A");
    }

    Matrix alpha = buildAlpha(A);
    std::vector<double> beta = buildBeta(A, b);
    std::vector<double> x = beta;  // x^(0) = beta
    std::vector<double> x_new(n);

    double alphaNorm = matrixNormC(alpha);

    log << "=== Метод простых итераций (Якоби) ===\n";
    log << "Начальное приближение: x^(0) = beta\n";
    log << "||alpha||_c = " << alphaNorm << "\n";
    log << "Точность eps = " << eps << ", max итераций = " << MAX_ITER << '\n';

    if (alphaNorm < 1.0) {
        log << "Достаточное условие сходимости ||alpha|| < 1 выполнено\n";
    } else {
        log << "ВНИМАНИЕ: ||alpha|| >= 1, достаточное условие сходимости не выполнено\n";
        log << "Используется эвристический критерий остановки: ||x^(k)-x^(k-1)||_c < eps\n";
    }

    int iter = 0;
    while (iter < MAX_ITER) {
        // x^(k+1) = beta + alpha * x^(k)
        std::vector<double> alphaX = multiplyMatrixVector(alpha, x);
        x_new = addVectors(beta, alphaX);

        std::vector<double> diff = subtractVectors(x_new, x);
        double diffNorm = vectorNormC(diff);
        ++iter;

        if (alphaNorm < 1.0) {
            double est = alphaNorm / (1.0 - alphaNorm) * diffNorm;
            log << "Итерация " << iter
                << ": ||x^(k)-x^(k-1)||_c = " << diffNorm
                << ", оценка погрешности = " << est << '\n';

            x = x_new;
            if (est <= eps) {
                log << "Сходимость достигнута по достаточному условию\n";
                return {x, iter};
            }
        } else {
            log << "Итерация " << iter
                << ": ||x^(k)-x^(k-1)||_c = " << diffNorm << '\n';

            x = x_new;
            if (diffNorm <= eps) {
                log << "Сходимость достигнута по эвристическому критерию\n";
                return {x, iter};
            }
        }
    }

    log << "ВНИМАНИЕ: метод простых итераций не сошёлся за " << MAX_ITER << " итераций\n";
    return {x, MAX_ITER};
}

std::pair<std::vector<double>, int> gaussSeidel(
    const Matrix& A, const std::vector<double>& b, double eps,
    std::ostream& log) {

    size_t n = A.size();
    if (b.size() != n) {
        throw std::invalid_argument("gaussSeidel: размер b не совпадает с размером A");
    }

    Matrix alpha = buildAlpha(A);
    std::vector<double> beta = buildBeta(A, b);
    std::vector<double> x = beta;
    double alphaNorm = matrixNormC(alpha);

    log << "=== Метод Зейделя ===\n";
    log << "Начальное приближение: x^(0) = beta\n";
    log << "||alpha||_c (для эквивалентной системы Якоби) = " << alphaNorm << "\n";
    log << "Точность eps = " << eps << ", max итераций = " << MAX_ITER << '\n';

    if (alphaNorm < 1.0) {
        log << "Достаточное условие сходимости выполнено\n";
    } else {
        log << "ВНИМАНИЕ: ||alpha|| >= 1, гарантии по достаточному условию нет\n";
        log << "Используется эвристический критерий остановки: ||x^(k)-x^(k-1)||_c < eps\n";
    }
        
    int iter = 0;
    while (iter < MAX_ITER) {
        std::vector<double> x_old = x;

        // x_i^(k+1) = beta_i + sum_{j<i} alpha_ij x_j^(k+1) + sum_{j>i} alpha_ij x_j^(k)
        for (size_t i = 0; i < n; ++i) {
            double sum = beta[i];

            for (size_t j = 0; j < i; ++j) {
                sum += alpha(i, j) * x[j];
            }
            for (size_t j = i + 1; j < n; ++j) {
                sum += alpha(i, j) * x_old[j];
            }

            x[i] = sum;
        }

        std::vector<double> diff = subtractVectors(x, x_old);
        double diffNorm = vectorNormC(diff);

        ++iter;

        if (alphaNorm < 1.0) {
            double est = alphaNorm / (1.0 - alphaNorm) * diffNorm; // можно брать норму верхнетреуг. матрицы C, тогда будет круче
            log << "Итерация " << iter
                << ": ||x^(k)-x^(k-1)||_c = " << diffNorm
                << ", оценка погрешности = " << est << '\n';

            if (est <= eps) {
                log << "Сходимость достигнута по достаточному условию\n";
                return {x, iter};
            }
        } else {
            log << "Итерация " << iter
                << ": ||x^(k)-x^(k-1)||_c = " << diffNorm << '\n';

            if (diffNorm <= eps) {
                log << "Сходимость достигнута по эвристическому критерию\n";
                return {x, iter};
            }
        }
    }

    log << "ВНИМАНИЕ: метод Зейделя не сошёлся за " << MAX_ITER << " итераций\n";
    return {x, MAX_ITER};
}

void run_1_3(const std::string& inputFile) {
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
        std::ifstream infile(inputFile);
        if (!infile) {
            throw std::runtime_error("Не открыт " + inputFile);
        }

        size_t n;
        infile >> n;

        Matrix A(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                infile >> A(i, j);
            }
        }

        std::vector<double> b(n);
        for (size_t i = 0; i < n; ++i) {
            infile >> b[i];
        }

        double eps;
        infile >> eps;
        if (eps <= 0) {
            eps = 1e-6;
        }

        log << "Загружена система " << n << "x" << n
            << ", eps = " << eps << " из " << inputFile << "\n\n";

        ensureNonZeroDiagonal(A, b, log);
        log << '\n';

        auto [x_jac, it_jac] = simpleIteration(A, b, eps, log);
        log << '\n';
        auto [x_sei, it_sei] = gaussSeidel(A, b, eps, log);
        log << '\n';

        double res_jac = residualNorm1(A, b, x_jac);
        double res_sei = residualNorm1(A, b, x_sei);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile);
        }

        out << std::fixed << std::setprecision(8);

        out << "=== Метод простых итераций ===\n";
        out << "Итераций: " << it_jac << "\n";
        out << "Решение x:\n";
        for (double v : x_jac) {
            out << v << "\n";
        }
        out << "||A x - b||_1 = " << res_jac << "\n\n";

        out << "=== Метод Зейделя ===\n";
        out << "Итераций: " << it_sei << "\n";
        out << "Решение x:\n";
        for (double v : x_sei) {
            out << v << "\n";
        }
        out << "||A x - b||_1 = " << res_sei << "\n\n";

        out << "=== Анализ ===\n";
        if (it_sei < it_jac) {
            out << "Метод Зейделя сошёлся быстрее\n";
        } else if (it_sei > it_jac) {
            out << "Метод простых итераций сошёлся быстрее\n";
        } else {
            out << "Оба метода сошлись за одинаковое число итераций\n";
        }

        std::cout << "Алгоритм 1.3 (метод простых итераций и метод Зейделя) завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << std::endl;
        std::cerr << e.what() << std::endl;
    }
}