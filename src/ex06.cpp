#include <iostream>
#include <cmath>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    Matrix3d A;
    A.setRandom();

    cout << "A = " << endl << A << endl;

    Vector3d b; // ColumnVector [3 x 1]
    b.setRandom();

    cout << "b = " << endl << b << endl;

    cout << "A * b = " << endl << A * b << endl; // [3x3][3x1]

    return 0;
}
