#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    Matrix3d A;

    A << 1, 2, 7, 2, 5, 6, 7, 6, 10;

    cout << "A = " << endl << A << endl;

    SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);

    if (es.info() != Eigen::Success)
    {
        cout << "Eigen solver failed!" << endl;
        return 1;
    }

    cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
    cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl;

    auto v1 = es.eigenvectors().col(0);
    auto v2 = es.eigenvectors().col(1);
    auto v3 = es.eigenvectors().col(2);

    auto lambda1 = es.eigenvalues()(0);
    auto lambda2 = es.eigenvalues()(1);
    auto lambda3 = es.eigenvalues()(2);

    cout << "The first eigenvector is:" << endl << v1 << endl;
    cout << "The second eigenvector is:" << endl << v2 << endl;
    cout << "The third eigenvector is:" << endl << v3 << endl;

    cout << "The first eigenvalue is:" << endl << lambda1 << endl;
    cout << "The second eigenvalue is:" << endl << lambda2 << endl;
    cout << "The third eigenvalue is:" << endl << lambda3 << endl;

    Matrix3d I = Matrix3d::Identity();

    cout << "(A-lambda1 * I) * v1 = " << endl << (A - lambda1 * I) * v1 << endl;
    cout << "(A-lambda2 * I) * v2 = " << endl << (A - lambda2 * I) * v2 << endl;
    cout << "(A-lambda3 * I) * v3 = " << endl << (A - lambda3 * I) * v3 << endl;

    return 0;
}
