#include <Eigen/Dense>
#include <iostream>
#include <memory>

using namespace std;
using namespace Eigen;

void foo(double **data, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
            cout << data[i][j] << " ";
        cout << endl;
    }
}

int main()
{
    MatrixXd A(10, 10);
    A.setRandom();

    auto data = A.data();
    
    auto data2D = make_unique<double*[]>(A.rows());

    for (int i = 0; i < A.rows(); i++)
        data2D[i] = data + i * A.cols();

    auto rows = static_cast<int>(A.rows());
    auto cols = static_cast<int>(A.cols());

    foo(data2D.get(), rows, cols);
}
