#include <iostream>
#include <print>
#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix3d;
    using Eigen::Matrix;
    using std::cout;

    Matrix3d A;
    A.setZero();

    // A = Matrix3d::Zero();

    utils::print("A", A);

    Matrix<double, 2, 4> B;

    B << 1.0, 2.0, 3.0, 4.0,
         5.0, 6.0, 7.0, 8.0;

    utils::print("B", B);
}
