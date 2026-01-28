#include <iostream>
#include <cmath>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix4d;
    using Eigen::MatrixXd;
    using Eigen::Matrix2d;

    Matrix4d A;

    A << 1, 2, 3, 4, 
         5, 6, 7, 8, 
         9, 10, 11, 12, 
         13, 14, 15, 16;

    utils::print("A = ", A);

    MatrixXd B = A.block(0, 0, 2, 2);

    utils::print("B = ", B);

    MatrixXd C = A.block(2, 2, 2, 2);

    utils::print("C = ", C);

    Matrix2d D;

    D << 42, 42, 
         42, 42;

    utils::print("D = ", D);

    A.block(2, 2, 2, 2) = D;

    utils::print("A = ", A);
}
