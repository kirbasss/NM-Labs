#pragma once

#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

inline bool isIntegerWithinTolerance(double value, double tol = 1e-9) {
    return std::fabs(value - std::round(value)) < tol;
}

inline int computeStepCount(double a, double b, double h) {
    if (h <= 0.0) {
        throw std::invalid_argument("Шаг h должен быть положительным");
    }
    if (b < a) {
        throw std::invalid_argument("Правая граница должна быть не меньше левой");
    }

    const double n = (b - a) / h;
    if (!isIntegerWithinTolerance(n)) {
        throw std::invalid_argument("Длина интервала должна делиться на шаг без остатка");
    }

    return static_cast<int>(std::round(n));
}

inline double maxAbsErrorRungeRomberg(const std::vector<double>& coarseX,
                                      const std::vector<double>& coarseY,
                                      const std::vector<double>& fineX,
                                      const std::vector<double>& fineY,
                                      int order) {
    if (coarseX.empty() || fineX.empty()) {
        return 0.0;
    }
    if (coarseX.size() != coarseY.size() || fineX.size() != fineY.size()) {
        throw std::invalid_argument("Размеры сетки и решения не совпадают");
    }
    if (coarseX.size() < 2 || fineX.size() < 2) {
        return 0.0;
    }

    const size_t coarseSteps = coarseX.size() - 1;
    const size_t fineSteps = fineX.size() - 1;
    if (fineSteps % coarseSteps != 0) {
        throw std::invalid_argument("Сетки не согласованы для метода Рунге-Ромберга");
    }

    const size_t ratio = fineSteps / coarseSteps;
    const double denom = std::pow(static_cast<double>(ratio), order) - 1.0;
    if (std::fabs(denom) < 1e-12) {
        throw std::invalid_argument("Некорректный порядок метода или отношение шагов");
    }

    double error = 0.0;
    for (size_t i = 0; i < coarseX.size(); ++i) {
        const size_t j = i * ratio;
        if (j >= fineY.size() || std::fabs(fineX[j] - coarseX[i]) > 1e-8) {
            throw std::invalid_argument("Сетки не согласованы по узлам");
        }
        error = std::max(error, std::fabs(fineY[j] - coarseY[i]) / denom);
    }

    return error;
}

inline std::filesystem::path resolveInputPath(const std::string& inputFile) {
    namespace fs = std::filesystem;

    const fs::path raw(inputFile);
    if (raw.is_absolute() && fs::exists(raw)) {
        return raw;
    }

    const fs::path current = fs::current_path();
    const std::vector<fs::path> candidates = {
        current / raw,
        current.parent_path() / raw,
        current.parent_path().parent_path() / raw
    };

    for (const auto& candidate : candidates) {
        if (fs::exists(candidate)) {
            return fs::weakly_canonical(candidate);
        }
    }

    return raw;
}
