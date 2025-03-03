#include <iostream>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    RowVector<double, 5> v1{1.0, 2.0, 3.0, 4.0, 5.0};
    VectorXi v2{{1, 2, 3, 4, 5}};

    cout << "v1 = " << v1 << "\n\n";
    cout << "v2 = " << v2 << "\n\n";

    cout << "v1[0] = " << v1[0] << "\n";
    cout << "v2[0] = " << v2[0] << "\n";

    return 0;
}
