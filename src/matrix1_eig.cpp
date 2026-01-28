#include <Eigen/Dense>
#include <print>

#include "utils_print.h"

int main()
{
    using Eigen::Matrix;
    using Eigen::Vector;
    using Eigen::Matrix3d;
    using Eigen::Matrix3i;
    using Eigen::Matrix3f;
    using Eigen::Vector3d;
    using Eigen::RowVector3d;
    using Eigen::Dynamic;

    Matrix3d A;

    utils::print("A uninitialised", A);

    A.setOnes();

    utils::print("A after setOnes", A);

    Matrix3i B = Matrix3i::Zero();

    utils::print("B after setZero", B);

    Matrix3f C;
    C.setIdentity();

    utils::print("C after setIdentity", C);

    Matrix3d D;

    D << 1, 2, 3,
        4, 5, 6,
        7, 8, 9;

    utils::print("D after initialization", D);

    D.col(0) << 1, 2, 3;
    D.col(1) << 4, 5, 6;
    D.col(2) << 7, 8, 9;

    utils::print("D after column initialization", D);

    Matrix<double, 6, 3> E;

    E.col(0) << 1, 2, 3, 4, 5, 6;
    E.col(1) << 7, 8, 9, 10, 11, 12;
    E.col(2) << 13, 14, 15, 16, 17, 18;

    utils::print("E after column initialization", E);

    Matrix<int, 3, 6> F;

    F.row(0) << 1, 2, 3, 4, 5, 6;
    F.row(1) << 7, 8, 9, 10, 11, 12;
    F.row(2) << 13, 14, 15, 16, 17, 18;

    utils::print("F after row initialization", F);

    Vector<double, 3> v;

    v << 1, 2, 3;

    utils::print("v after initialization", v);

    Vector3d w(1, 2, 3);

    utils::print("w after initialization", w);

    Matrix<double, Dynamic, Dynamic> A_dyn(3, 3); // Shorthand for MatrixXd

    A_dyn << 1, 2, 3,
        4, 5, 6,
        7, 8, 9;

    utils::print("A_dyn after initialization", A_dyn);

    A_dyn.resize(1, 9); // No reallocation

    utils::print("A_dyn after resizing", A_dyn);

    A_dyn.resize(6, 6); // Reallocation.

    utils::print("A_dyn after resizing", A_dyn);

    A_dyn.setZero();

    utils::print("A_dyn after setting to zero", A_dyn);

    Matrix<int, Dynamic, 3> B_dyn;

    B_dyn.resize(10, 3);

    B_dyn.setRandom();

    utils::print("B_dyn after initialization", B_dyn);

    RowVector3d v_row(1.0, 2.0, 3.0);

    utils::print("v_row after initialization", v_row);

    Vector3d v_col(1.0, 2.0, 3.0);

    utils::print("v_col after initialization", v_col);
}
