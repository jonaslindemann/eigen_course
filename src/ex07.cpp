#include <iostream>
#include <cmath>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    Matrix4d A;
    A.setRandom();

    cout << "A = " << endl << A << endl;

    cout << "det(A) = " << A.determinant() << endl;
    cout << "A^(-1) = " << endl << A.inverse() << endl;
    cout << "A^T = " << endl << A.transpose() << endl;

    cout << "A * A^(-1) = " << endl << A * A.inverse() << endl;

    if ((A * A.inverse()).isIdentity())
        cout << "A * A^(-1) is identity" << endl;
    else
        cout << "A * A^(-1) is not identity" << endl;

    Matrix3d B;

    B << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    auto det = B.determinant();
    bool isSingular = std::abs(det) < 1e-10;

    if (isSingular)
    {
        cout << "B is singular" << endl;
    }
    else
    {
        cout << "B is not singular" << endl;
    }

    FullPivLU<MatrixXd> lu(B);
    bool isSingular_lu = !lu.isInvertible();

    if (isSingular_lu)
    {
        cout << "B is singular (LU)" << endl;
    }
    else
    {
        cout << "B is not singular (LU)" << endl;
    }

    return 0;
}
