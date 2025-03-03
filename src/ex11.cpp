#include <iostream>
#include <cmath>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    Matrix2d A;
    A << 1, 2, 3, 4;

    cout << "A = " << endl << A << endl;

    Matrix2d B;

    B << 5, 6, 7, 8;

    cout << "B = " << endl << B << endl;

    MatrixXd C(A.rows(), A.cols() + B.cols());

    C << A, B;

    cout << "C = " << endl << C << endl;

    MatrixXd D(A.rows() + B.rows(), A.cols());

    D << A, B;

    cout << "D = " << endl << D << endl;

    return 0;
}
