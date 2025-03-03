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

    Vector3d b; // ColumnVector
    b.setRandom();

    cout << "b = " << endl << b << endl;

    cout << "A * b = " << endl << A * b << endl;

    return 0;
}
