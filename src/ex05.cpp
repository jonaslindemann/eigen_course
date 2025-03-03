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

    Matrix3d B;
    B.setRandom();

    cout << "B = " << endl << B << endl;

    cout << "A + B = " << endl << A + B << endl;
    cout << "A - B = " << endl << A - B << endl;
    cout << "A.array() * B.array() = " << endl << A.array() * B.array() << endl;
    cout << "A * B = " << endl << A * B << endl;

    return 0;
}
