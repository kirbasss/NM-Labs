CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

PROJECT(lab4)

INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(../lab1/include)

ADD_EXECUTABLE(main
    main.cpp
    src/CauchyProblem.cpp
    src/BoundaryValueProblem.cpp
    ../lab1/src/Tridiagonal.cpp
)
TARGET_COMPILE_FEATURES(main PRIVATE cxx_std_23)
