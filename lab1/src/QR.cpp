#include "QR.hpp"
#include "Utils.hpp"
#include <fstream>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <stdexcept>

const int MAX_ITER = 10000;
extern const double EPS_ZERO;

double columnSubdiagNorm(const Matrix& A, size_t l) {
    size_t n = A.size();
    double s = 0.0;
    for (size_t m = l + 1; m < n; ++m) {
        s += A(m, l) * A(m, l);
    }
    return std::sqrt(s);
}

// Собственные значения блока 2x2:
// [ a b ]
// [ c d ]
std::pair<std::complex<double>, std::complex<double>>
eigenvalues2x2(double a, double b, double c, double d) {
    std::complex<double> tr = a + d;
    std::complex<double> det = a * d - b * c;
    std::complex<double> disc = tr * tr - std::complex<double>(4.0) * det;
    std::complex<double> root = std::sqrt(disc);
    return { (tr + root) / 2.0, (tr - root) / 2.0 };
}

// Норма поддиагональной части
double strictLowerNorm(const Matrix& A) {
    size_t n = A.size();
    double s = 0.0;
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            s += A(i, j) * A(i, j);
        }
    }
    return std::sqrt(s);
}

// Формируем отражение Хаусхолдера H = I - 2 vv^T/(v^T v),
// где v имеет размер n, но ненулевые компоненты только начиная с offset.
Matrix buildHouseholderFromTail(const std::vector<double>& x, size_t offset) {
    const size_t n = x.size();
    Matrix H = Matrix::identityMatrix(n);

    double normx = 0.0;
    for (size_t i = offset; i < n; ++i) {
        normx += x[i] * x[i];
    }
    normx = std::sqrt(normx);

    if (normx < EPS_ZERO) {
        return H;
    }

    std::vector<double> v = x;
    v[offset] += signNonZero(x[offset]) * normx;

    double vtv = 0.0;
    for (size_t i = offset; i < n; ++i) {
        vtv += v[i] * v[i];
    }

    if (vtv < EPS_ZERO) {
        return H;
    }

    for (size_t i = offset; i < n; ++i) {
        for (size_t j = offset; j < n; ++j) {
            H(i, j) -= 2.0 * v[i] * v[j] / vtv;
        }
    }

    return H;
}

QRDecomposition qrDecomposeHouseholder(const Matrix& A_orig, std::ostream& log) {
    const size_t n = A_orig.size();

    Matrix R = A_orig;
    Matrix Q = Matrix::identityMatrix(n);

    log << "=== QR-разложение (Хаусхолдер) ===\n";

    for (size_t k = 0; k + 1 < n; ++k) {
        std::vector<double> x(n, 0.0);
        for (size_t i = k; i < n; ++i) {
            x[i] = R(i, k);
        }

        Matrix H = buildHouseholderFromTail(x, k);

        R = multiply(H, R);
        Q = multiply(Q, H); // H^T = H

        log << "Шаг k = " << k << ": обнуление поддиагональных элементов столбца " << k << "\n";
    }

    return {Q, R};
}

QREigenResult qrAlgorithm(const Matrix& A_orig, double eps, std::ostream& log) {
    const size_t n = A_orig.size();
    if (eps <= 0.0) {
        throw std::invalid_argument("qrAlgorithm: eps must be positive");
    }

    log << "=== QR-алгоритм ===\n";

    Matrix A = A_orig;

    std::vector<std::complex<double>> prevComplex(n, std::complex<double>(0.0, 0.0));
    std::vector<bool> hasPrevComplex(n, false);

    int iter = 0;
    while (iter < MAX_ITER) {
        log << "\nИтерация " << iter
            << ": ||strict lower||_2 = " << strictLowerNorm(A) << "\n";

        // 1) Сначала попробуем выделить завершившиеся блоки
        bool all_done = true;
        size_t i = 0;
        while (i < n) {
            // последний 1x1 блок
            if (i == n - 1) {
                ++i;
                continue;
            }

            // вещественное собственное значение, критерий:
            // столбец i имеет малую норму под диагональю
            double col_norm = columnSubdiagNorm(A, i);
            if (col_norm <= eps) {
                ++i;
                continue;
            }

            // Иначе рассматриваем кандидат на 2x2 блок
            auto cur = eigenvalues2x2(
                A(i, i), A(i, i + 1),
                A(i + 1, i), A(i + 1, i + 1)
            );

            if (iter > 0 && hasPrevComplex[i]) {
                // Учитываем возможную перестановку корней между итерациями
                double directDiff =
                    std::abs(cur.first  - prevComplex[i]) +
                    std::abs(cur.second - prevComplex[i + 1]);

                double swappedDiff =
                    std::abs(cur.first  - prevComplex[i + 1]) +
                    std::abs(cur.second - prevComplex[i]);

                if (std::min(directDiff, swappedDiff) <= eps) {
                    i += 2;
                    continue;
                }
            }

            prevComplex[i]     = cur.first;
            prevComplex[i + 1] = cur.second;
            hasPrevComplex[i] = true;

            all_done = false;
            i += 2;
        }

        if (all_done) {
            log << "Критерий остановки выполнен\n";
            break;
        }

        // 2) QR-разложение и шаг A_{k+1} = R_k Q_k
        auto qr = qrDecomposeHouseholder(A, log);
        A = multiply(qr.R, qr.Q);

        ++iter;
    }

    if (iter == MAX_ITER) {
        log << "ВНИМАНИЕ: достигнут предел по числу итераций\n";
    }

    // Извлечение собственных значений из почти квазитреугольной матрицы
    std::vector<std::complex<double>> eigenvalues;
    size_t i = 0;
    while (i < n) {
        if (i == n - 1 || columnSubdiagNorm(A, i) <= eps) {
            eigenvalues.emplace_back(A(i, i), 0.0);
            ++i;
        } else {
            auto ev = eigenvalues2x2(
                A(i, i), A(i, i + 1),
                A(i + 1, i), A(i + 1, i + 1)
            );
            eigenvalues.push_back(ev.first);
            eigenvalues.push_back(ev.second);
            i += 2;
        }
    }

    return {eigenvalues, A, iter};
}

void run_1_5(const std::string& inputFile) {
    namespace fs = std::filesystem;
    std::string base = fs::path(inputFile).stem().string();
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

        double eps;
        infile >> eps;
        if (!infile || eps <= 0.0) {
            throw std::runtime_error("Некорректное значение eps");
        }

        log << "Загружена матрица " << n << "×" << n 
            << ", eps = " << eps << " из " << inputFile << "\n\n";

        auto result = qrAlgorithm(A, eps, log);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile);
        }
        out << std::fixed << std::setprecision(8);
        out << "=== Собственные значения ===\n";
        for (const auto& lambda : result.eigenvalues) {
            if (std::fabs(lambda.imag()) < EPS_ZERO) {
                out << lambda.real() << "\n";
            } else {
                out << lambda.real()
                    << (lambda.imag() >= 0.0 ? " + " : " - ")
                    << std::fabs(lambda.imag()) << "i\n";
            }
        }

        out << "\n=== Финальная матрица A_k ===\n";
        result.finalA.print(out);

        out << "\n=== Число итераций ===\n";
        out << result.iterations << "\n";

        std::cout << "Алгоритм 1.5 (QR-разложение) завершён. Результаты в " << outFile << '\n';
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << std::endl;
        std::cerr << e.what() << std::endl;
    }
}