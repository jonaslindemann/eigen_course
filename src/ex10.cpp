#include <iostream>
#include <cmath>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    Matrix<double, 2, 3> A;

    A << 1, 2, 3, 4, 5, 6;

    cout << "A = " << endl << A << endl;

    Matrix<double, 3, 2> B = A.reshaped(3, 2);

    cout << "B = " << endl << B << endl;

    Vector<double, 6> v;

    v << 1, 2, 3, 4, 5, 6;

    cout << "v = " << endl << v << endl;

    Matrix<double, 2, 3> C = v.reshaped(2, 3);

    cout << "C = " << endl << C << endl;

    return 0;
}
