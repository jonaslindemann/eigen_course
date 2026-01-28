#include <iostream>
#include <cmath>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::RowVector;
    using std::cout;

    RowVector<double, 10> v1;
    v1.setRandom();

    utils::print("v1", v1);

    RowVector<double, 10> v2;
    v2.setRandom();

    utils::print("v2", v2);

    utils::print("v1 + v2", v1 + v2);
    utils::print("v1 - v2", v1 - v2);
    utils::print("v1.cwiseProduct(v2)", v1.cwiseProduct(v2));
    std::println("v1.dot(v2) = {}", v1.dot(v2));
}
