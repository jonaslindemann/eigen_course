Solutions
=========

Basic Matrix and Vector operations
----------------------------------

**1. Vector Creation and Access**

.. code:: cpp

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


.. button-link:: https://godbolt.org/z/d1Pjo1cGe
    :color: primary
    :outline:

    Open solution in Compiler Explorer

**2. Matrix Initialization**

.. code:: cpp

    #include <iostream>

    #include <Eigen/Dense>

    // Just for clarity of the examples.

    using namespace Eigen;
    using namespace std;

    int main()
    {
        Matrix3d A;
        A.setZero();

        // A = Matrix3d::Zero();

        cout << "A = " << endl << A << endl;

        Matrix<double, 2, 4> B;

        B << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0;

        cout << "B = " << endl << B << endl;

        return 0;
    }


.. button-link:: https://godbolt.org/z/sPxaTs5hc
    :color: primary
    :outline:

    Open solution in Compiler Explorer
 
**3. Basic Math Functions**

.. code:: cpp

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

        auto v2a = v1.array().abs().sqrt();

        cout << "v2 = " << v2 << "\n\n";

        cout << "v2a = " << v2a << "\n\n";

        auto v3 = v1;

        for (int i = 0; i < v3.size(); ++i)
            v3[i] = exp(v3[i]);

        cout << "v3 = " << v3 << "\n\n";

        cout << "v1.sum() = " << v1.sum() << "\n";
        cout << "v1.maxCoeff() = " << v1.maxCoeff() << "\n";
        cout << "v1.minCoeff() = " << v1.minCoeff() << "\n";

        return 0;
    }

.. button-link:: https://godbolt.org/z/Kzh94ME6o
    :color: primary
    :outline:

**4. Coefficient-wise Operations**

.. code:: cpp

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

.. button-link:: https://godbolt.org/z/dv7e6ffxr
    :color: primary
    :outline:
 
**5. Matrix Operations**

.. code:: cpp

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

.. button-link:: https://godbolt.org/z/WxWMPjxc5
    :color: primary
    :outline:

**6. Vector-Matrix Operations**

.. code:: cpp
    
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

        Vector3d b; // ColumnVector [3 x 1]
        b.setRandom();

        cout << "b = " << endl << b << endl;

        cout << "A * b = " << endl << A * b << endl; // [3x3][3x1]

        return 0;
    }

.. button-link:: https://godbolt.org/z/YEzafd8Wq
    :color: primary
    :outline:

    Open solution in Compiler Explorer

**7. Special Matrix Operations**

.. code:: cpp

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



.. button-link:: https://godbolt.org/z/vM1nqo6Tf
    :color: primary
    :outline:

    Open solution in Compiler Explorer


Advanced Matrix operations
--------------------------

**8. Block Operations**

.. code:: cpp

    #include <iostream>
    #include <cmath>

    #include <Eigen/Dense>

    // Just for clarity of the examples.

    using namespace Eigen;
    using namespace std;

    int main()
    {
        Matrix4d A;

        A << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

        cout << "A = " << endl << A << endl;

        MatrixXd B = A.block(0, 0, 2, 2);

        cout << "B = " << endl << B << endl;

        MatrixXd C = A.block(2, 2, 2, 2);

        cout << "C = " << endl << C << endl;

        Matrix2d D;

        D << 42, 42, 42, 42;

        cout << "D = " << endl << D << endl;

        A.block(2, 2, 2, 2) = D;

        cout << "A = " << endl << A << endl;

        return 0;
    }

.. button-link:: https://godbolt.org/z/1bovzj7Gd
    :color: primary
    :outline:

    Open solution in Compiler Explorer

**9. Row and Column Operations**

.. code:: cpp

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



.. button-link:: https://godbolt.org/z/16cMqbY9Y
    :color: primary
    :outline:

    Open solution in Compiler Explorer


**10. Resizing and Reshaping**

.. code:: cpp

    #include <iostream>
    #include <cmath>

    #include <Eigen/Dense>

    // Just for clarity of the examples.

    using namespace Eigen;
    using namespace std;

    int main()
    {
        Matrix<double, 2, 3> A;

        A << 1, 2, 3, 4, 5, 6;

        cout << "A = " << endl << A << endl;

        Matrix<double, 3, 2> B = A.reshaped(3, 2);

        cout << "B = " << endl << B << endl;

        Vector<double, 6> v;

        v << 1, 2, 3, 4, 5, 6;

        cout << "v = " << endl << v << endl;

        Matrix<double, 2, 3> C = v.reshaped(2, 3);

        cout << "C = " << endl << C << endl;

        return 0;
    }

.. button-link:: https://godbolt.org/z/v44hhbM4n
    :color: primary
    :outline:

    Open solution in Compiler Explorer


**11. Concatenation**

.. code:: cpp

    #include <iostream>
    #include <cmath>

    #include <Eigen/Dense>

    // Just for clarity of the examples.

    using namespace Eigen;
    using namespace std;

    int main()
    {
        Matrix2d A;
        A << 1, 2, 3, 4;

        cout << "A = " << endl << A << endl;

        Matrix2d B;

        B << 5, 6, 7, 8;

        cout << "B = " << endl << B << endl;

        MatrixXd C(A.rows(), A.cols() + B.cols());

        C << A, B;

        cout << "C = " << endl << C << endl;

        MatrixXd D(A.rows() + B.rows(), A.cols());

        D << A, B;

        cout << "D = " << endl << D << endl;

        return 0;
    }


.. button-link:: https://godbolt.org/z/ber68r51E
    :color: primary
    :outline:

    Open solution in Compiler Explorer

**12. Advanced Slicing**

.. code:: cpp

    #include <iostream>
    #include <cmath>
    #include <vector>

    #include <Eigen/Dense>

    // Just for clarity of the examples.

    using namespace Eigen;
    using namespace std;

    int main()
    {
        Matrix<double, 5, 5> A;

        A << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25;

        cout << "A = " << endl << A << endl;

        vector<Index> rows = {1, 3, 4};
        vector<Index> cols = {0, 2};

        MatrixXd B = A(rows, cols);

        cout << "B = " << endl << B << endl;

        cout << "diagonal of A" << endl << A.diagonal() << endl;

        auto n = A.rows();

        VectorXd antiDiagonal(n);

        for (int i = 0; i < n; i++)
            antiDiagonal(i) = A(i, n - i - 1);

        cout << "anti-diagonal of A" << endl << antiDiagonal << endl;

        auto tiledRows = 4;
        auto tiledCols = 3;

        Matrix2d C;

        C << 1, 2, 3, 4;

        MatrixXd D(tiledRows * C.rows(), tiledCols * C.cols());

        for (auto i = 0; i < tiledRows; i++)
            for (auto j = 0; j < tiledCols; j++)
                D.block(i * C.rows(), j * C.cols(), C.rows(), C.cols()) = C;

        cout << "D = " << endl << D << endl;

        return 0;
    }

.. button-link:: https://godbolt.org/z/hsGWeYMWo
    :color: primary
    :outline:

    Open solution in Compiler Explorer

**13. Linear Algebra Operations**

.. code:: cpp

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


.. button-link:: https://godbolt.org/z/qdabzhdqf
    :color: primary
    :outline:

    Open solution in Compiler Explorer

**14. Eigenvalues and Eigenvectors**

.. code:: cpp

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


.. button-link:: https://godbolt.org/z/qx3xWffo6
    :color: primary
    :outline:

    Open solution in Compiler Explorer

Practical Applications
----------------------

**15. Image Processing Basics**

.. code:: cpp

    #include <iostream>
    #include <cmath>
    #include <vector>
    #include <numbers>

    #include <Eigen/Dense>

    #include <img/Image.h>

    // Just for clarity of the examples.

    using namespace Eigen;
    using namespace std;
    using namespace img;

    Eigen::MatrixXd adjustBrightnessContrast(const Eigen::MatrixXd &image, double alpha, double beta)
    {
        // Convert to array for element-wise operations
        Eigen::ArrayXXd result = image.array();

        // Apply contrast (multiply)
        result = alpha * result;

        // Apply brightness (add)
        result = result + beta;

        // Clamp values to valid pixel range [0, 1]
        result = result.min(1.0f).max(0.0f);

        return result.matrix();
    }

    Eigen::MatrixXd createGaussianKernel(int size, double sigma)
    {
        if (size % 2 == 0)
        {
            size++; // Make sure kernel size is odd
        }

        Eigen::MatrixXd kernel(size, size);
        int center = size / 2;
        double sum = 0.0;

        // Fill the kernel with Gaussian values
        for (int y = 0; y < size; y++)
        {
            for (int x = 0; x < size; x++)
            {
                int dx = x - center;
                int dy = y - center;
                float exponent = -(dx * dx + dy * dy) / (2 * sigma * sigma);
                kernel(y, x) = std::exp(exponent) / (2 * std::numbers::pi * sigma * sigma);
                sum += kernel(y, x);
            }
        }

        // Normalize the kernel so its sum equals 1
        kernel /= sum;

        return kernel;
    }

    Eigen::MatrixXd applyGaussianFilter(const Eigen::MatrixXd &image, int kernelSize, double sigma)
    {
        Eigen::MatrixXd result = image;
        Eigen::MatrixXd kernel = createGaussianKernel(kernelSize, sigma);

        int kernelRadius = kernelSize / 2;

        for (int y = kernelRadius; y < image.rows() - kernelRadius; y++)
        {
            for (int x = kernelRadius; x < image.cols() - kernelRadius; x++)
            {
                double sum = 0.0;
                for (int ky = 0; ky < kernelSize; ky++)
                {
                    for (int kx = 0; kx < kernelSize; kx++)
                    {
                        int imgX = x + kx - kernelRadius;
                        int imgY = y + ky - kernelRadius;
                        sum += image(imgY, imgX) * kernel(ky, kx);
                    }
                }
                result(y, x) = sum;
            }
        }
        return result;
    }

    ImageGd copyToImage(const Eigen::MatrixXd &mat)
    {
        ImageGd image;
        image.resize(mat.rows(), mat.cols());
        for (int i = 0; i < mat.rows(); i++)
            for (int j = 0; j < mat.cols(); j++)
                image(i, j) = mat(i, j);

        return image;
    }

    int main()
    {
        ImageGd image;

        cout << "Loading image..." << endl;

        if (!load("images/half-moon-986269.png", image))
        {
            cout << "Error loading image" << endl;
            return 1;
        }

        // Create a matrix with the image data

        auto mat = image.as_matrix();

        cout << "Adjusting brightness contrast..." << endl;

        MatrixXd brighterImageMat = adjustBrightnessContrast(mat, 1.0, 0.5);

        // Create a new image with the modified data

        cout << "Converting to image..." << endl;

        auto brighterImage = copyToImage(brighterImageMat);

        // Save the new image

        cout << "Saving image..." << endl;

        if (!save("images/brighter-half-moon-986269.png", brighterImage))
        {
            cout << "Error saving image" << endl;
            return 1;
        }

        cout << "Applying Gaussian filter..." << endl;

        MatrixXd blurredImageMat = applyGaussianFilter(mat, 10, 1.0);

        // Create a new image with the modified data

        cout << "Converting to image..." << endl;

        auto blurredImage = copyToImage(blurredImageMat);

        // Save the new image

        cout << "Saving image..." << endl;

        if (!save("images/blurred-half-moon-986269.png", blurredImage))
        {
            cout << "Error saving image" << endl;
            return 1;
        }
    }

.. button-link:: https://godbolt.org/z/qx3xWffo6
    :color: primary
    :outline:

    Open solution in Compiler Explorer