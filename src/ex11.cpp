#include <iostream>
#include <cmath>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix2d;
    using Eigen::MatrixXd;

    Matrix2d A;

    A << 1, 2, 
         3, 4;

    utils::print("A", A);

    Matrix2d B;

    B << 5, 6, 
         7, 8;

    utils::print("B", B);

    MatrixXd C(A.rows(), A.cols() + B.cols());

    C << A, B;

    utils::print("C", C);

    MatrixXd D(A.rows() + B.rows(), A.cols());

    D << A, B;

    utils::print("D", D);
}
