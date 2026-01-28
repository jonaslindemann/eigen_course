#include <iostream>
#include <cmath>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix3d;
    using std::cout;

    Matrix3d A;
    A.setRandom();

    utils::print("A", A);

    Matrix3d B;
    B.setRandom();

    utils::print("B", B);

    utils::print("A + B", A + B);
    utils::print("A - B", A - B);
    utils::print("A.array() * B.array()", A.array() * B.array());
    utils::print("A * B", A * B);
}
