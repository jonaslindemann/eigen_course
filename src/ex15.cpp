#include <cmath>
#include <vector>
#include <numbers>
#include <print>

#include <Eigen/Dense>

#include "utils_print.h"

#include <img/Image.h>

Eigen::MatrixXd adjustBrightnessContrast(const Eigen::MatrixXd &image, double alpha, double beta)
{
    // Convert to array for element-wise operations
    Eigen::ArrayXXd result = image.array();

    // Apply contrast (multiply)
    result = alpha * result;

    // Apply brightness (add)
    result = result + beta;

    // Clamp values to valid pixel range [0, 1]
    result = result.max(0.0).min(1.0);

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
            double exponent = -(dx * dx + dy * dy) / (2 * sigma * sigma);
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

img::ImageGd copyToImage(const Eigen::MatrixXd &mat)
{
    img::ImageGd image;
    const int rows = static_cast<int>(mat.rows());
    const int cols = static_cast<int>(mat.cols());
    
    image.resize(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            image(i, j) = mat(i, j);

    return image;
}

int main()
{
    using std::cout;
    using Eigen::MatrixXd;
    using img::ImageGd;

    ImageGd image;

    std::println("Loading image...");

    if (!load("images/half-moon-986269.png", image))
    {
        std::println("Error loading image");
        return 1;
    }

    // Create a matrix with the image data

    auto mat = image.as_matrix();

    std::println("Adjusting brightness contrast...");

    MatrixXd brighterImageMat = adjustBrightnessContrast(mat, 1.0, 0.5);

    // Create a new image with the modified data

    std::println("Converting to image...");

    auto brighterImage = copyToImage(brighterImageMat);

    // Save the new image

    std::println("Saving image...");

    if (!save("images/brighter-half-moon-986269.png", brighterImage))
    {
        std::println("Error saving image");
        return 1;
    }

    std::println("Applying Gaussian filter...");

    MatrixXd blurredImageMat = applyGaussianFilter(mat, 10, 1.0);

    // Create a new image with the modified data

    std::println("Converting to image...");

    auto blurredImage = copyToImage(blurredImageMat);

    // Save the new image

    std::println("Saving image...");

    if (!save("images/blurred-half-moon-986269.png", blurredImage))
    {
        std::println("Error saving image");
        return 1;
    }
}
