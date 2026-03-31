#include "LU.hpp"
#include "Tridiagonal.hpp"
#include "Iterative.hpp"
#include "Jacobi.hpp"
#include "QR.hpp"

int main() {
    run_1_1("input/1.1.txt");
    run_1_2("input/1.2.txt");
    run_1_3("input/1.3.txt");
    run_1_4("input/1.4.txt");
    run_1_5("input/1.5.txt");

    return 0;
}