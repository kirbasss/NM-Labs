#pragma once
#include <vector>
#include <string>
#include <ostream>

struct SplineInput {
    std::vector<double> x;
    std::vector<double> y;
    double xStar;
};

struct CubicSplineSegment {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double d = 0.0;
    double xLeft = 0.0;
    double xRight = 0.0;
};

struct CubicSpline {
    std::vector<CubicSplineSegment> segments;
};

SplineInput loadSplineInputFromFile(const std::string& filename);

CubicSpline buildNaturalCubicSpline(const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    std::ostream& log);

double evaluateSpline(const CubicSpline& spline,
                      double xStar,
                      std::ostream& log);

void run_3_2(const std::string& inputFile);