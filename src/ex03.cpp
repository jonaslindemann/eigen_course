#include <iostream>
#include <cmath>

#include <Eigen/Dense>

// Just for clarity of the examples.

using namespace Eigen;
using namespace std;

int main()
{
    // -5 -4 -3 -2 -1 0 1 2 3 4 5
    RowVector<double, 11> v1;

    v1.setLinSpaced(11, -5.0, 5.0);
    // v1 = Vector<double, 11>::LinSpaced(11, -5.0, 5.0);

    cout << "v1 = " << v1 << "\n\n";

    auto v2 = v1.cwiseAbs().cwiseSqrt();

    cout << "v2 = " << v2 << "\n\n";

    auto v3 = v1;

    for (int i = 0; i < v3.size(); ++i)
        v3[i] = exp(v3[i]);

    cout << "v3 = " << v3 << "\n\n";

    cout << "v1.sum() = " << v1.sum() << "\n";
    cout << "v1.maxCoeff() = " << v1.maxCoeff() << "\n";
    cout << "v1.minCoeff() = " << v1.minCoeff() << "\n";

    return 0;
}
