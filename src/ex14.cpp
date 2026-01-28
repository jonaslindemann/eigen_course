#include <iostream>
#include <cmath>
#include <vector>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix3d;
    using Eigen::SelfAdjointEigenSolver;

    Matrix3d A;

    A << 1, 2, 7, 
         2, 5, 6, 
         7, 6, 10;

    utils::print("A", A);

    SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);

    if (es.info() != Eigen::Success)
    {
        std::println("Eigen solver failed.");
        return 1;
    }

    utils::print("The eigenvalues of A are:", es.eigenvalues());
    utils::print("The matrix of eigenvectors, V, is:", es.eigenvectors());

    auto v1 = es.eigenvectors().col(0);
    auto v2 = es.eigenvectors().col(1);
    auto v3 = es.eigenvectors().col(2);

    auto lambda1 = es.eigenvalues()(0);
    auto lambda2 = es.eigenvalues()(1);
    auto lambda3 = es.eigenvalues()(2);

    utils::print("The first eigenvector is:", v1);
    utils::print("The second eigenvector is:", v2);
    utils::print("The third eigenvector is:", v3);

    std::println();

    std::println("The first eigenvalue is: {}", lambda1);
    std::println("The second eigenvalue is: {}", lambda2);
    std::println("The third eigenvalue is: {}\n", lambda3);

    Matrix3d I = Matrix3d::Identity();

    std::println("(A-lambda1 * I) * v1 = {}", (A - lambda1 * I) * v1);
    std::println("(A-lambda2 * I) * v2 = {}", (A - lambda2 * I) * v2);
    std::println("(A-lambda3 * I) * v3 = {}", (A - lambda3 * I) * v3);
}
