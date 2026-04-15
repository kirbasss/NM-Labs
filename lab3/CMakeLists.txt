CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

PROJECT(lab3)

INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(../lab1/include)

ADD_EXECUTABLE(main main.cpp src/Interpolation.cpp src/Spline.cpp ../lab1/src/Tridiagonal.cpp 
../lab1/src/Matrix.cpp ../lab1/src/LU.cpp src/LeastSquares.cpp ../lab1/src/utils.cpp src/Differentiation.cpp src/Integration.cpp)
TARGET_COMPILE_FEATURES(main PRIVATE cxx_std_23)