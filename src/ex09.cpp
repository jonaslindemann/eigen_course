#include <iostream>
#include <cmath>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    Matrix<double, 3, 4> A;

    A.setRandom();

    cout << "A = " << endl << A << endl;

    cout << "A.row(1) = " << endl << A.row(1) << endl;
    cout << "A.col(2) = " << endl << A.col(2) << endl;

    A.row(1) << 1, 2, 3, 4;

    cout << "A = " << endl << A << endl;

    VectorXd temp = A.col(2);
    A.col(2) = A.col(3);
    A.col(3) = temp;

    cout << "A = " << endl << A << endl;

    A.col(2).swap(A.col(3));

    cout << "A = " << endl << A << endl;

    return 0;
}
