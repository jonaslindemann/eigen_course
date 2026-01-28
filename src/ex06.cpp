#include <iostream>
#include <cmath>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix3d;
    using Eigen::Vector3d;
    using std::cout;

    Matrix3d A;
    A.setRandom();

    utils::print("A", A);

    Vector3d b; // ColumnVector [3 x 1]
    b.setRandom();

    utils::print("b", b);
    utils::print("A * b", A * b);
}
