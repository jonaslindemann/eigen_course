#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    Matrix<double, 5, 5> A;

    A << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25;

    cout << "A = " << endl << A << endl;

    vector<int> rows = {1, 3, 4};
    vector<int> cols = {0, 2};

    MatrixXd B = A(rows, cols);

    cout << "B = " << endl << B << endl;

    cout << "diagonal of A" << endl << A.diagonal() << endl;

    auto n = A.rows();

    VectorXd antiDiagonal(n);

    for (int i = 0; i < n; i++)
        antiDiagonal(i) = A(i, n - i - 1);

    cout << "anti-diagonal of A" << endl << antiDiagonal << endl;

    auto tiledRows = 4;
    auto tiledCols = 3;

    Matrix2d C;

    C << 1, 2, 3, 4;

    MatrixXd D(tiledRows * C.rows(), tiledCols * C.cols());

    for (auto i = 0; i < tiledRows; i++)
        for (auto j = 0; j < tiledCols; j++)
            D.block(i * C.rows(), j * C.cols(), C.rows(), C.cols()) = C;

    cout << "D = " << endl << D << endl;

    return 0;
}
