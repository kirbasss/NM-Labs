#pragma once
#include <string>
#include <ostream>
#include <functional>

struct IntegrationInput {
    double x0 = 0.0;
    double xk = 0.0;
    double h1 = 0.0;
    double h2 = 0.0;
};

struct IntegrationResult {
    double rectangles = 0.0;
    double trapezoids = 0.0;
    double simpson = 0.0;
};

IntegrationInput loadIntegrationInputFromFile(const std::string& filename);

double integrand(double x);

double integrateRectangles(const std::function<double(double)>& f,
                           double a,
                           double b,
                           double h,
                           std::ostream& log);

double integrateTrapezoids(const std::function<double(double)>& f,
                           double a,
                           double b,
                           double h,
                           std::ostream& log);

double integrateSimpson(const std::function<double(double)>& f,
                        double a,
                        double b,
                        double h,
                        std::ostream& log);

double rungeRomberg(double FhCoarse,
                    double FhFine,
                    double hCoarse,
                    double hFine,
                    int p);

void run_3_5(const std::string& inputFile);