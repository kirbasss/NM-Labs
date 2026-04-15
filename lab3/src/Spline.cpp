#include "Spline.hpp"
#include "Tridiagonal.hpp"

#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace {
    const double EPS = 1e-12;

    void validateSplineInput(const std::vector<double>& x,
                             const std::vector<double>& y) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Размеры x и y не совпадают");
        }
        if (x.size() < 2) {
            throw std::invalid_argument("Для построения сплайна нужно хотя бы 2 точки");
        }

        for (size_t i = 1; i < x.size(); ++i) {
            if (!(x[i] > x[i - 1])) {
                throw std::invalid_argument("Узлы x должны быть строго возрастающими");
            }
        }
    }

    size_t findSplineSegment(const CubicSpline& spline, double xStar) {
        if (spline.segments.empty()) {
            throw std::runtime_error("Сплайн не содержит сегментов");
        }

        for (size_t i = 0; i < spline.segments.size(); ++i) {
            const auto& seg = spline.segments[i];
            if (xStar >= seg.xLeft - EPS && xStar <= seg.xRight + EPS) {
                return i;
            }
        }

        throw std::out_of_range("Точка x* находится вне диапазона сплайна");
    }
}

SplineInput loadSplineInputFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filename);
    }

    size_t n;
    file >> n;
    if (!file) {
        throw std::runtime_error("Ошибка чтения числа узлов из файла " + filename);
    }

    SplineInput input;
    input.x.resize(n);
    input.y.resize(n);

    for (size_t i = 0; i < n; ++i) {
        file >> input.x[i];
    }
    for (size_t i = 0; i < n; ++i) {
        file >> input.y[i];
    }
    file >> input.xStar;

    if (!file) {
        throw std::runtime_error("Ошибка чтения данных сплайна из файла " + filename);
    }

    validateSplineInput(input.x, input.y);
    return input;
}

CubicSpline buildNaturalCubicSpline(const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    std::ostream& log) {
    validateSplineInput(x, y);

    const size_t n = x.size() - 1; // число интервалов
    std::vector<double> h(n);

    log << "=== Построение кубического сплайна ===\n";
    log << "Граничные условия: S''(x0)=0, S''(xn)=0\n\n";

    for (size_t i = 0; i < n; ++i) {
        h[i] = x[i + 1] - x[i];
        log << "h[" << i << "] = x[" << (i + 1) << "] - x[" << i << "] = " << h[i] << "\n";
    }
    log << "\n";

    std::vector<double> cFull(n + 1, 0.0);

    // Если только 2 точки, сплайн вырождается в линейную функцию на одном отрезке
    if (n == 1) {
        CubicSpline spline;
        CubicSplineSegment seg;
        seg.a = y[0];
        seg.b = (y[1] - y[0]) / h[0];
        seg.c = 0.0;
        seg.d = 0.0;
        seg.xLeft = x[0];
        seg.xRight = x[1];
        spline.segments.push_back(seg);

        log << "Только один интервал: сплайн линейный.\n";
        return spline;
    }

    // Система для внутренних c[1..n-1]
    // a[i] - нижняя диагональ, b[i] - главная, c[i] - верхняя, d[i] - правая часть
    const size_t m = n - 1;
    std::vector<double> lower(m > 0 ? m - 1 : 0, 0.0);
    std::vector<double> diag(m, 0.0);
    std::vector<double> upper(m > 0 ? m - 1 : 0, 0.0);
    std::vector<double> rhs(m, 0.0);

    log << "Составляем трёхдиагональную систему для внутренних коэффициентов c_i:\n\n";

    for (size_t row = 0; row < m; ++row) {
        size_t i = row + 1; // внутренний индекс c_i, i=1..n-1

        if (row > 0) {
            lower[row - 1] = h[i - 1];
        }

        diag[row] = 2.0 * (h[i - 1] + h[i]);

        if (row + 1 < m) {
            upper[row] = h[i];
        }

        rhs[row] = 3.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);

        log << "Уравнение для c[" << i << "]: ";
        if (i > 1) {
            log << h[i - 1] << " * c[" << (i - 1) << "] + ";
        }
        log << 2.0 * (h[i - 1] + h[i]) << " * c[" << i << "]";
        if (i < n - 1) {
            log << " + " << h[i] << " * c[" << (i + 1) << "]";
        }
        log << " = " << rhs[row] << "\n";
    }
    log << "\n";

    auto cInner = tridiagSolve(lower, diag, upper, rhs, log);

    for (size_t i = 1; i <= n - 1; ++i) {
        cFull[i] = cInner[i - 1];
    }
    cFull[0] = 0.0;
    cFull[n] = 0.0;

    log << "\nНайденные коэффициенты c_i:\n";
    for (size_t i = 0; i <= n; ++i) {
        log << "c[" << i << "] = " << cFull[i] << "\n";
    }
    log << "\n";

    CubicSpline spline;
    spline.segments.resize(n);

    for (size_t i = 0; i < n; ++i) {
        CubicSplineSegment seg;
        seg.a = y[i];
        seg.b = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * cFull[i] + cFull[i + 1]) / 3.0;
        seg.c = cFull[i];
        seg.d = (cFull[i + 1] - cFull[i]) / (3.0 * h[i]);
        seg.xLeft = x[i];
        seg.xRight = x[i + 1];

        spline.segments[i] = seg;

        log << "Интервал [" << seg.xLeft << ", " << seg.xRight << "]\n";
        log << "S_" << i << "(x) = a + b(x-x_i) + c(x-x_i)^2 + d(x-x_i)^3\n";
        log << "a = " << seg.a << "\n";
        log << "b = " << seg.b << "\n";
        log << "c = " << seg.c << "\n";
        log << "d = " << seg.d << "\n\n";
    }

    return spline;
}

double evaluateSpline(const CubicSpline& spline,
                      double xStar,
                      std::ostream& log) {
    size_t idx = findSplineSegment(spline, xStar);
    const auto& seg = spline.segments[idx];

    double t = xStar - seg.xLeft;
    double value = seg.a
                 + seg.b * t
                 + seg.c * t * t
                 + seg.d * t * t * t;

    log << "=== Вычисление значения сплайна ===\n";
    log << "x* = " << xStar << " принадлежит интервалу ["
        << seg.xLeft << ", " << seg.xRight << "]\n";
    log << "t = x* - xLeft = " << t << "\n";
    log << "S(x*) = " << seg.a
        << " + " << seg.b << " * " << t
        << " + " << seg.c << " * " << t << "^2"
        << " + " << seg.d << " * " << t << "^3"
        << " = " << value << "\n\n";

    return value;
}

void run_3_2(const std::string& inputFile) {
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
        auto input = loadSplineInputFromFile(inputFile);

        log << std::fixed << std::setprecision(10);

        log << "Задача 3.2. Кубический сплайн\n";
        log << "Узлы:\n";
        for (size_t i = 0; i < input.x.size(); ++i) {
            log << "i = " << i
                << ": x = " << input.x[i]
                << ", y = " << input.y[i] << "\n";
        }
        log << "x* = " << input.xStar << "\n\n";

        CubicSpline spline = buildNaturalCubicSpline(input.x, input.y, log);
        double splineValue = evaluateSpline(spline, input.xStar, log);

        std::ofstream out(outFile);
        if (!out) {
            throw std::runtime_error("Не удалось создать файл результата " + outFile);
        }

        out << std::fixed << std::setprecision(10);

        out << "=== Задача 3.2. Кубический сплайн ===\n\n";

        out << "Узлы интерполяции:\n";
        for (size_t i = 0; i < input.x.size(); ++i) {
            out << "i = " << i
                << ", x_i = " << input.x[i]
                << ", f_i = " << input.y[i] << "\n";
        }

        out << "\n=== Коэффициенты сплайна по интервалам ===\n";
        for (size_t i = 0; i < spline.segments.size(); ++i) {
            const auto& seg = spline.segments[i];
            out << "Интервал [" << seg.xLeft << ", " << seg.xRight << "]\n";
            out << "a = " << seg.a << "\n";
            out << "b = " << seg.b << "\n";
            out << "c = " << seg.c << "\n";
            out << "d = " << seg.d << "\n\n";
        }

        out << "=== Значение в точке x* ===\n";
        out << "x* = " << input.xStar << "\n";
        out << "S(x*) = " << splineValue << "\n";

        std::cout << "Алгоритм 3.2 (кубический сплайн) завершён. Результаты в " << outFile << "\n";
        std::cout << "Лог: " << logFile << "\n";

    } catch (const std::exception& e) {
        log << "ОШИБКА: " << e.what() << "\n";
        std::cerr << e.what() << std::endl;
    }
}