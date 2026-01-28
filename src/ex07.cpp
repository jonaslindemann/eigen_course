#include <iostream>
#include <cmath>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix4d;
    using Eigen::Matrix3d;
    using Eigen::MatrixXd;
    using Eigen::FullPivLU;
    using std::cout;

    Matrix4d A;
    A.setRandom();

    utils::print("A", A);

    std::println("det(A) = {}", A.determinant());

    utils::print("A^(-1)", A.inverse());
    utils::print("A^T", A.transpose());
    utils::print("A * A^(-1)", A * A.inverse());

    if ((A * A.inverse()).isIdentity())
        std::println("A * A^(-1) is identity");
    else
        std::println("A * A^(-1) is not identity");

    Matrix3d B;

    B << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    auto det = B.determinant();
    bool isSingular = std::abs(det) < 1e-10;

    if (isSingular)
        std::println("B is singular");
    else
        std::println("B is not singular");

    FullPivLU<MatrixXd> lu(B);

    bool isSingular_lu = !lu.isInvertible();

    if (isSingular_lu)
        std::println("B is singular (LU)");
    else
        std::println("B is not singular (LU)");
}
