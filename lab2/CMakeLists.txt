CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

PROJECT(lab2)

INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(../lab1/include)

ADD_EXECUTABLE(main main.cpp src/NonlinearEquation.cpp src/NonlinearSystem.cpp ../lab1/src/Matrix.cpp ../lab1/src/LU.cpp ../lab1/src/utils.cpp)
TARGET_COMPILE_FEATURES(main PRIVATE cxx_std_23)