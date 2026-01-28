#include <iostream>
#include <cmath>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix;
    using Eigen::VectorXd;
    using std::cout;

    Matrix<double, 3, 4> A;

    A.setRandom();

    utils::print("A", A);

    utils::print("A.row(1)", A.row(1));
    utils::print("A.col(2)", A.col(2));

    A.row(1) << 1, 2, 3, 4;

    utils::print("A", A);

    VectorXd temp = A.col(2);
    A.col(2) = A.col(3);
    A.col(3) = temp;

    utils::print("A", A);

    A.col(2).swap(A.col(3));

    utils::print("A", A);
}
