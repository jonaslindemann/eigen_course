#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main()
{
    MatrixXd A(10, 10);
    A.setRandom();

    VectorXd b1(10);
    b1.setRandom();

    VectorXd b2(10);
    b2.setRandom();

    FullPivLU<MatrixXd> FullLU(A);

    auto x1 = FullLU.solve(b1);
    auto x2 = FullLU.solve(b2);

    cout << "x1: " << x1.transpose() << endl;
    cout << "x2: " << x2.transpose() << endl;

    MatrixXd B(10, 2);
    B.setRandom();

    auto X = FullLU.solve(B);

    cout << "X = " << endl << X << endl;

    return 0;
}
