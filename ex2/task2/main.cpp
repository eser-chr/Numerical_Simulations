// Task 2.2 (Group05)

// sudo apt-get install build-essential libeigen3-dev libopenblas-dev
// g++ main.cpp -std=c++17 -O3 -lopenblas -DEIGEN_DONT_PARALLELIZE -march=native -ffast-math -o mmm
// ./mmm CUSTOM 1024
// ./mmm EIGEN 1024
// export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 && ./mmm BLAS 1024

#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <chrono> // for time measurement
#include <cblas.h>
#include <eigen3/Eigen/Dense>


// namespace program_options
namespace program_options {
    std::tuple<std::string, std::size_t> parse(int argc, char *argv[]) {
        if (argc != 3)
            throw std::runtime_error("unexpected number of arguments");
        auto mode = std::string(argv[1]);
        const std::vector<std::string> modes = {"CUSTOM", "BLAS", "EIGEN"};
        if (std::find(modes.begin(), modes.end(), mode) == modes.end())
            throw std::runtime_error("unsupported mode");
        std::size_t size;
        if (std::sscanf(argv[2], "%zu", &size) != 1 && size != 0)
            throw std::runtime_error("invalid size");
        std::printf("arguments: impl=%s size=%zu\n", mode.data(), size);
        return {mode, size};
    }
} 


// Function to perform MMM using CUSTOM implementation
Eigen::MatrixXd custom_MMM(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B) {
    std::size_t rowsA = A.rows();
    std::size_t colsA = A.cols();
    std::size_t colsB = B.cols();

    // Check if matrix dimensions are compatible for matrix multiplication
    if (colsA != B.rows()) {
        throw std::invalid_argument("Matrix dimensions are not compatible for multiplication");
    }

    Eigen::MatrixXd result(rowsA, colsB);

    for (std::size_t i = 0; i < rowsA; ++i) {
        for (std::size_t j = 0; j < colsB; ++j) {
            result(i, j) = 0.0;
            for (std::size_t k = 0; k < colsA; ++k) {
                result(i, j) += A(i, k) * B(k, j);
            }
        }
    }

    return result;
}


// Function to perform MMM using BLAS library
Eigen::MatrixXd blas_MMM(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B) {
    Eigen::MatrixXd result(A.rows(), B.cols());
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A.rows(), B.cols(), A.cols(), 1.0, A.data(), A.rows(), B.data(), B.rows(), 0.0, result.data(), result.rows());
    return result;
}


// Function to perform MMM using EIGEN library
Eigen::MatrixXd eigen_MMM(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B) {
    return A * B;
}


int main(int argc, char *argv[]) try {

  auto [mode, size] = program_options::parse(argc, argv);

  std::size_t N = size;
  std::string impl = mode;

    
  // Creating and initializing the matrix M
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N, N);
  for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < N; ++j) {
          M(i, j) = static_cast<double>(i + j * N); 
      }
  }

  // Computing the transpose of M
  Eigen::MatrixXd M_T = M.transpose();


  // Start time measurement
  auto start_time = std::chrono::high_resolution_clock::now();

  // Perform MMM based on the selected implementation
  Eigen::MatrixXd result;
  if (impl == "CUSTOM") {
      result = custom_MMM(M, M_T);
  } else if (impl == "BLAS") {
      result = blas_MMM(M, M_T);
  } else if (impl == "EIGEN") {
      result = eigen_MMM(M, M_T);
  } else {
      throw std::runtime_error("Unsupported mode");
  }

  // End time measurement
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);


  // Printing the result of the runtime
  std::cout << "Runtime: " << duration.count() << " microseconds\n";


// Symmetry check for the product matrix
bool isSymmetric = true;

for (std::size_t i = 0; i < result.rows(); ++i) {
    for (std::size_t j = 0; j < result.cols(); ++j) {
        if (std::abs(result(i, j) - result(j, i)) > 1e-16) {
            isSymmetric = false;
            break;
        }
    }
    if (!isSymmetric) {
        break;
    }
}

// Printing the symmetry check result
std::cout << "Symmetry Check: " << (isSymmetric ? "Symmetric" : "Not Symmetric") << "\n";


  return EXIT_SUCCESS;
} catch (std::exception &e) {
  std::cout << e.what() << std::endl;
  return EXIT_FAILURE;
}
