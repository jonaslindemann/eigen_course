#include <iostream>
#include <cmath>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::RowVector;

    // -5 -4 -3 -2 -1 0 1 2 3 4 5
    RowVector<double, 11> v1;

    v1.setLinSpaced(11, -5.0, 5.0);
    // v1 = Vector<double, 11>::LinSpaced(11, -5.0, 5.0);

    utils::print("v1", v1);

    auto v2 = v1.cwiseAbs().cwiseSqrt();
    auto v2a = v1.array().abs().sqrt().matrix();

    utils::print("v2", v2);
    utils::print("v2a", v2a);

    auto v3 = v1;

    for (auto i = 0; i < v3.size(); ++i)
        v3[i] = exp(v3[i]);

    utils::print("v3", v3);

    std::println("v1.sum() = {}", v1.sum());
    std::println("v1.prod() = {}", v1.prod());
    std::println("v1.mean() = {}", v1.mean());
    std::println("v1.maxCoeff() = {}", v1.maxCoeff());
    std::println("v1.minCoeff() = {}", v1.minCoeff());
}
