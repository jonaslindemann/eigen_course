#include <iostream>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    Matrix3d A;
    A.setZero();

    // A = Matrix3d::Zero();

    cout << "A = " << endl << A << endl;

    Matrix<double, 2, 4> B;

    B << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0;

    cout << "B = " << endl << B << endl;

    return 0;
}
