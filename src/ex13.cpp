#include <iostream>
#include <cmath>
#include <vector>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    MatrixXd A(10, 10);

    A << MatrixXd::Random(10, 10);

    utils::print("A", A);

    VectorXd b(10);

    b << VectorXd::Random(10);

    utils::print("b", b);

    VectorXd x = A.fullPivLu().solve(b);

    utils::print("The solution is:", x);

    std::println("The relative error is: {}", (A * x - b).norm() / b.norm());

    std::println("A * x = {}", A * x);
}
