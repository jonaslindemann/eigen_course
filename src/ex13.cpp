#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    MatrixXd A(10, 10);

    A << MatrixXd::Random(10, 10);

    cout << "A = " << endl << A << endl;

    VectorXd b(10);

    b << VectorXd::Random(10);

    cout << "b = " << endl << b << endl;

    VectorXd x = A.fullPivLu().solve(b);

    cout << "The solution is:\n" << x << endl;

    cout << "The relative error is:\n" << (A * x - b).norm() / b.norm() << endl;

    cout << "A * x = " << endl << A * x << endl;

    return 0;
}
