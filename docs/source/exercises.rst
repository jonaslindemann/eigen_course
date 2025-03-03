Exercises
=========

Basic Matrix and Vector operations
----------------------------------

**1. Vector Creation and Access**

- Create a 5-element vector using Eigen and initialize it with values 1 through 5.
- Print the vector's elements and demonstrate how to access individual elements.

**2. Matrix Initialization**

- Create a 3×3 matrix with all elements set to zero.
- Initialize a 2×4 matrix with specific values of your choice.
- Print both matrices to verify their contents.
 
**3. Basic Math Functions**

- Create a vector with values ranging from -5 to 5.
- Apply the absolute value, square root (for positive elements), and exponential functions.
- Calculate the sum, mean, minimum, and maximum values of the vector.

.. hint:: 

    If you want to apply sqrt to all components in a vector you have to use the **.cwiseSqrt()** method. The same applies to other functions like **.cwiseAbs()** and **.cwiseExp()**. If you want to apply a function to all values in a matrix, you can use the **.array()** method to convert the matrix to an array and then apply the function. **A.array().sqrt()** will apply the square root to all elements in the matrix.

    It is also possible to use the **.unaryExpr()** method to apply a function to all elements in a matrix. For example, **A.unaryExpr([](double x) { return foo(x); });** will apply the function **foo()** to all elements in the matrix A.

**4. Coefficient-wise Operations**

- Create two vectors of the same length.
- Perform element-wise multiplication, division, and power operations.
- Calculate the dot product and compare with element-wise multiplication.
 
**5. Matrix Operations**

- Create two 3×3 matrices with different values.
- Perform and print the results of addition, subtraction, and element-wise multiplication.
- Calculate and print the matrix product of the two matrices.

**6. Vector-Matrix Operations**

- Create a 3×3 matrix and a 3-element vector.
- Multiply the matrix by the vector and explain the result.
- Experiment with different ways to perform the multiplication.

**7. Special Matrix Operations**

- Create a 4×4 matrix and calculate its determinant, inverse, and transpose.
- Verify that A × A⁻¹ equals the identity matrix (within numerical precision).
- Explore how to check if a matrix is singular and handle potential errors.

.. hint:: 

    In addition to using the **.determinant()** to determine if a matrix is singular it is also possible to use the **.fullPivLu().isInvertible()** method. This method returns a boolean value indicating if the matrix is invertible. The **.fullPivLu()** method can also be used to calculate the inverse of a matrix. 

Advanced Matrix operations
--------------------------

**8. Block Operations**

- Create a 4×4 matrix with sequential numbers (1-16).
- Extract the 2×2 block in the top-left corner.
- Extract the 2×2 block in the bottom-right corner.
- Replace the bottom-left 2×2 block with a new 2×2 matrix.

.. hint:: 

    The **.block()** method can be used to extract a block from a matrix. The method takes four arguments: the row and column index of the top-left corner of the block, and the number of rows and columns in the block. For example, **A.block(0, 0, 2, 2)** will extract a 2×2 block from the top-left corner of the matrix A. It is also possible to assign a matrix to using the **.block()**-method. **A.block(0, 0, 2, 2) = B;** will replace the top-left 2×2 block in A with the matrix B.

**9. Row and Column Operations**

- Create a 3×4 matrix with random values.
- Extract the second row and third column.
- Replace the first row with new values.
- Swap two columns of your choice.

.. hint:: 

    The **.swap()** method can be used to swap two columns in a matrix.

**10. Resizing and Reshaping**

- Create a 2×3 matrix and then resize it to a 3×2 matrix.
- Create a 6-element vector and reshape it into a 2×3 matrix.
- Discuss what happens to the elements during these operations.

.. hint:: 

    Use the **.reshaped()** method to reshape a matrix. 

**11. Concatenation**

- Create two 2×2 matrices.
- Horizontally concatenate them to form a 2×4 matrix.
- Vertically concatenate them to form a 4×2 matrix.
- Create and demonstrate a function that can concatenate matrices of compatible dimensions.

.. hint:: 

    The **\<\<** operator can be used to concatenate matrices. **A \<\< B, C;** will concatenate the matrices B and C horizontally or vertically depending on the matrix shapes.

**12. Advanced Slicing**

- Create a 5×5 matrix with sequential values.
- Extract non-contiguous rows and columns (e.g., rows 1, 3, 4 and columns 0, 2).
- Extract a diagonal or anti-diagonal of the matrix.
- Create a tiled pattern by repeating a smaller matrix.

.. hint:: 

    **std::vector<int> indices = {1, 3, 4};** can be used as indices to extract non-contiguous areas of a matrix. 

**13. Linear Algebra Operations**

- Create a system of linear equations represented as Ax = b.
- Solve the system using Eigen's solver capabilities.
- Verify your solution by substituting it back into the original equations.

.. hint:: 

    Use **.fullPivLu()** to solve a system of linear equations.

**14. Eigenvalues and Eigenvectors**

- Create a symmetric 3×3 matrix.
- Calculate its eigenvalues and eigenvectors.
- Verify that Av = λv for each eigenvalue-eigenvector pair.

.. hint:: 

    Use the class **SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);** to compute eigen values and eiven vectors.  Use the **es.eigenvalues()** and **es.eigenvectors()** methods to calculate the eigenvalues and eigenvectors of a matrix.

Practical Applications
----------------------

**15. Image Processing Basics**

- Represent a grayscale image as an Eigen matrix.
- Implement basic operations like brightness adjustment and contrast enhancement.
- Apply simple filters (e.g., blur) using matrix operations.

.. hint::

    Use the Img library to load images that can converted to Eigen matrices. The library can be found here:

    `img: single-header-only C++ image library with Eigen interoperability <https://github.com/ThibaultLejemble/img>`_

    An image can be loaded and converted to a grayscale image (range [0.0..1.0]) using the following code snippet:

    .. code-block:: cpp

        Img::ImageGd image;

        std::cout << "Loading image..." << std::endl;

        if (!Img::load("images/half-moon-986269.png", image))
        {
            std::cout << "Error loading image" << std::endl;
            return 1;
        }

    Conversion from an Image class to an Eigen matrix can be done using the **.as_matrix()** method.

    .. code-block:: cpp

        MatrixXd matrix = image.as_matrix();

    Converting back to an Image class can be done using the following function:

    .. code-block:: cpp

        Img::ImageGd copyToImage(const Eigen::MatrixXd &mat)
        {
            ImageGd image;
            image.resize(mat.rows(), mat.cols());
            for (int i = 0; i < mat.rows(); i++)
                for (int j = 0; j < mat.cols(); j++)
                    image(i, j) = mat(i, j);

            return image;
        }

    Saving the image to a file can be done using the following code snippet:

    .. code-block:: cpp

        if (!Img::save("output.png", copyToImage(matrix)))
        {
            std::cout << "Error saving image" << std::endl;
            return 1;
        }

    The gaussian kernel can be generated using the following function:

    .. code-block:: cpp

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

