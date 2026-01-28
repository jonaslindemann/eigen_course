#include <iostream>
#include <Eigen/Dense>
#include <print>

#include "utils_print.h"

int main()
{
    using Eigen::RowVector;
    using Eigen::VectorXi;
    
    RowVector<double, 5> v1{1.0, 2.0, 3.0, 4.0, 5.0};
    VectorXi v2{{1, 2, 3, 4, 5}};

    utils::print("v1", v1);
    utils::print("v2", v2);

    std::printf("v1[0] = %.1f\n", v1[0]);
    std::printf("v2[0] = %d\n", v2[0]);
}
