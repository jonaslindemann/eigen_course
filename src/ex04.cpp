#include <iostream>
#include <cmath>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    RowVector<double, 10> v1;
    v1.setRandom();

    cout << "v1 = " << v1 << "\n\n";

    RowVector<double, 10> v2;
    v2.setRandom();

    cout << "v2 = " << v2 << "\n\n";

    cout << "v1 + v2 = " << v1 + v2 << "\n\n";
    cout << "v1 - v2 = " << v1 + v2 << "\n\n";
    cout << "v1.cwiseProduct(v2) = " << v1.cwiseProduct(v2) << "\n\n";
    cout << "v1.dot(v2) = " << v1.dot(v2) << "\n\n";

    return 0;
}
