#include <iostream>
#include <cmath>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    Matrix4d A;

    A << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

    cout << "A = " << endl << A << endl;

    MatrixXd B = A.block(0, 0, 2, 2);

    cout << "B = " << endl << B << endl;

    MatrixXd C = A.block(2, 2, 2, 2);

    cout << "C = " << endl << C << endl;

    Matrix2d D;

    D << 42, 42, 42, 42;

    cout << "D = " << endl << D << endl;

    A.block(2, 2, 2, 2) = D;

    cout << "A = " << endl << A << endl;

    return 0;
}
