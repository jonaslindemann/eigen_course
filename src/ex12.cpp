#include <iostream>
#include <cmath>
#include <vector>
#include <print>

#include <Eigen/Dense>
#include "utils_print.h"

int main()
{
    using Eigen::Matrix;
    using Eigen::Matrix2d;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::Index;

    Matrix<double, 5, 5> A;

    A <<  1,  2,  3,  4,  5, 
          6,  7,  8,  9, 10, 
         11, 12, 13, 14, 15, 
         16, 17, 18, 19, 20, 
         21, 22, 23, 24, 25;

    utils::print("A", A);

    std::vector<Index> rows = {1, 3, 4};
    std::vector<Index> cols = {0, 2};

    MatrixXd B = A(rows, cols);

    utils::print("B", B);

    utils::print("diagonal of A", A.diagonal());

    auto n = A.rows();

    VectorXd antiDiagonal(n);

    for (int i = 0; i < n; i++)
        antiDiagonal(i) = A(i, n - i - 1);

    utils::print("anti-diagonal of A", antiDiagonal);

    auto tiledRows = 4;
    auto tiledCols = 3;

    Matrix2d C;

    C << 1, 2, 3, 4;

    MatrixXd D(tiledRows * C.rows(), tiledCols * C.cols());

    for (auto i = 0; i < tiledRows; i++)
        for (auto j = 0; j < tiledCols; j++)
            D.block(i * C.rows(), j * C.cols(), C.rows(), C.cols()) = C;

    utils::print("D", D);

    // Demonstrate different formatting options
    std::print("--- Different Format Examples ---\n\n");
    
    Matrix2d E;
    E << 1.23456789, 2.3456789,
         3.456789,   4.56789;
    
    utils::print("E (Default)", E);
    utils::print("E (Compact)", E, PrintFormat::Compact);
    utils::print("E (Python)", E, PrintFormat::Python);
    utils::print("E (Clean)", E, PrintFormat::Clean);
    utils::print("E (Octave)", E, PrintFormat::Octave);
    utils::print("E (CSV)", E, PrintFormat::CSV);

    // Custom format with specific precision
    Eigen::IOFormat CustomFmt(2, 0, ", ", "\n", "[", "]");
    utils::print("E (Custom 2 decimals)", E, CustomFmt);
}
