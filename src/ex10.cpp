#include <iostream>
#include <cmath>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix;
    using Eigen::Vector;

    Matrix<double, 2, 3> A;

    A << 1, 2, 
         3, 4, 
         5, 6;

    utils::print("A", A);

    Matrix<double, 3, 2> B = A.reshaped(3, 2);

    utils::print("B", B);

    Vector<double, 6> v;

    v << 1, 2, 3, 4, 5, 6;

    utils::print("v", v);

    Matrix<double, 2, 3> C = v.reshaped(2, 3);

    utils::print("C", C);
}
