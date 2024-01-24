// To compile this code : make
// To execute the solver without source : ./solver without_source 10 10000 0.0 1.0
// To execute the solver with source : ./solver with_source 50 100000 0.0 0.0 0.5 0.5 0.1

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <chrono>

namespace program_options {

    struct Options {
    std::string name;
    size_t N;
    size_t iters;
    double fix_west;
    double fix_east;
    bool has_source;
    double source_x;
    double source_y;
    double source_sigma;
    void print() const {
        std::printf("name: %s\n", name.c_str());
        std::printf("N: %zu\n", N);
        std::printf("iters: %zu\n", iters);
        std::printf("fix_west: %lf\n", fix_west);
        std::printf("fix_east: %lf\n", fix_east);
        std::printf("has_source: %s\n", has_source ? "true" : "false");
        std::printf("source_x: %lf\n", source_x);
        std::printf("source_y: %lf\n", source_y);
        std::printf("source_sigma: %lf\n", source_sigma);
    }
    };

    auto parse(int argc, char *argv[]) {
    if (argc != 9 && argc != 6)
        throw std::runtime_error("unexpected number of arguments");
    Options opts;
    opts.name = argv[1];
    if (std::sscanf(argv[2], "%zu", &opts.N) != 1 && opts.N >= 2)
        throw std::runtime_error("invalid parameter for N");
    if (std::sscanf(argv[3], "%zu", &opts.iters) != 1 && opts.iters != 0)
        throw std::runtime_error("invalid parameter for iters");
    if (std::sscanf(argv[4], "%lf", &opts.fix_west) != 1)
        throw std::runtime_error("invalid value for fix_west");
    if (std::sscanf(argv[5], "%lf", &opts.fix_east) != 1)
        throw std::runtime_error("invalid value for fix_east");
    if (argc == 6) {
        opts.has_source = false;
        opts.source_x = NAN;
        opts.source_y = NAN;
        opts.source_sigma = NAN;
        return opts;
    }
    if (std::sscanf(argv[6], "%lf", &opts.source_x) != 1)
        throw std::runtime_error("invalid value for source_x");
    if (std::sscanf(argv[7], "%lf", &opts.source_y) != 1)
        throw std::runtime_error("invalid value for source_y");
    if (std::sscanf(argv[8], "%lf", &opts.source_sigma) != 1)
        throw std::runtime_error("invalid value for source_sigma");
    opts.has_source = true;
    return opts;
    }

} // namespace program_options

void solver(std::vector<double>& input, const std::vector<double>& source,const program_options::Options& opts);
double gaussian_source(double x, double y, double sigma, double x_src, double y_src);
void print_vec(const std::vector<double>& vec, int n);
void compute_source_map(std::vector<double>& vec, const program_options::Options& opts);
std::vector<double> mul_mat(const std::vector<double>& A, const std::vector<double>& u, int N);
void compute_residuals(const std::vector<double>& u, const std::vector<double>& source, const program_options::Options& opts);
double norm_2(const std::vector<double>& res);
double norm_inf(const std::vector<double>& res);

int main(int argc, char *argv[]) try {

    auto opts = program_options::parse(argc, argv);

    opts.print();

    // write csv
    auto write = [ N = opts.N, name = opts.name ](const auto &x) -> auto {
        std::ofstream csv;
        csv.open(name + ".csv");
        for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < N - 1; ++i) {
            csv << x[i + j * N] << " ";
        }
        csv << x[(N - 1) + j * N];
        csv << "\n";
        }
        csv.close();
    };

    // Creating the gaussian source map
    std::vector<double> source(opts.N * opts.N, 0.0);
    compute_source_map(source, opts);

    // Creating the output vector
    std::vector<double> demo(opts.N * opts.N, 0.0);
    
    // Solving the equation on the vector and printing the runtime
    auto start = std::chrono::high_resolution_clock::now();
    solver(demo, source, opts);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "\nNumber of iterations : " << opts.iters << std::endl;
    std::cout << "Total runtime : " << duration.count() / 1e6 << " s\n" << std::endl;

    // Writing the solution in the file
    write(demo);

    // Printing the residual norms of the solution
    compute_residuals(demo, source, opts);


    return EXIT_SUCCESS;
    } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
}

double gaussian_source(double x, double y, double sigma, double x_src, double y_src) {
  double dx = pow(x - x_src, 2);
  double dy = pow(y - y_src, 2);
  double exponent = - (dx + dy) / (2 * pow(sigma, 2));
  double coefficient = 1 / (2 * M_PI * pow(sigma, 2));
  return coefficient * exp(exponent);
}

void solver(std::vector<double>& input, const std::vector<double>& source, const program_options::Options& opts) {

  // Matrix-free Jacobi solver

  /* Initalization of the buffer vector and of variables */
  size_t dim_1d = opts.N;
  size_t dim_2d = opts.N * opts.N;
  double h = 1 / ((double) dim_1d - 1);
  std::vector<double> buffer(dim_2d, 0.0);
  std::vector<double> gaussian_map(dim_2d, 0.0);

  /**** 1st step : applying boundary conditions ****/
  for (int i = 0; i < dim_2d; i++) {
    if ((i % dim_1d) == 0) { input[i] = opts.fix_west, buffer[i] = opts.fix_west; } // The point is on the west border
    if (i % dim_1d == (dim_1d - 1)) { input[i] = opts.fix_east, buffer[i] = opts.fix_east; } // The point is on the east border
  }

  /**** 2nd step : iterative solver ****/

  int w, e, n, s;
  bool valid = 0;
  double x2, y2;

  for (int k = 0; k < opts.iters; k++) {

    for (int i = 0; i < dim_2d; i++) {
      // Computing the directions
      w = i - 1;
      e = i + 1;
      n = i + dim_1d;
      s = i - dim_1d;

      // Checking if the point can be computed
      x2 = (i % dim_1d) * h;
      y2 = ((int) i / (int) dim_1d) * h;

      if (x2 != 0 and x2 != 1) {
        if (y2 != 0 and y2 != 1) { valid = 1; }
      }

      // Computing the point if possible (i.e. the point is not on an edge)
      if (valid) { buffer[i] = 0.25 * (pow(h,2) * source[i]  + input[n] + input[s] + input[w] + input[e]); }
      valid = 0;

    } 
    // Copying the buffer into the initial array
    input = buffer;
  }
}

void print_vec(const std::vector<double>& vec, int n) {

  for (int i = 0; i < n * n; i++) {
      std::cout << vec[i] << " ";

      if (i % n == n - 1) {
        std::cout << std::endl;
      }
  }

  std::cout << std::endl;

}

std::vector<double> mul_mat(const std::vector<double>& A, const std::vector<double>& u, int N) {
    // Useful function to compute Ah*uh
    std::vector<double> result(N, 0.0);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            result[i] += A[i * N + j] * u[j];
        }
    }

    return result;
}

void compute_residuals(const std::vector<double>& u, const std::vector<double>& source, const program_options::Options& opts) {

  // Compute and print in the console the residuals with the 2 different norms

  int N2 = u.size();
  int N = sqrt(N2);
  double h = 1 / (double) (N - 1);
  double coef = 1 / (double) pow(h, 2);

  int sw = N + 1, s = N + 2, se = N + 3;
  int w = 2 * N + 1, c = 2 * N + 2, e = 2 * N + 3;
  int nw = 3 * N + 1, n = 3 * N + 2, ne = 3 * N + 3;

  std::vector<double> uh;
  std::vector<double> bh(9, 0.0);
  std::vector<double> rh(9, 0.0);

  // Creating Ah
  std::vector<double> Ah = {
    4, -1, 0, -1, 0, 0, 0, 0, 0,
    -1, 4, -1, 0, -1, 0, 0, 0, 0,
    0, -1, 4, 0, 0, -1, 0, 0, 0,
    -1, 0, 0, 4, -1, 0, -1, 0, 0,
    0, -1, 0, -1, 4, -1, 0, -1, 0,
    0, 0, -1, 0, -1, 4, 0, 0, -1,
    0, 0, 0, -1, 0, 0, 4, -1, 0,
    0, 0, 0, 0, -1, 0, -1, 4, -1,
    0, 0, 0, 0, 0,  -1, 0, -1, 4
  };

  for (int i = 0; i < Ah.size(); i++) { Ah[i] *= coef; }

  // Creating uh 9x1 vector
  if (u.size() >= 25) {

    uh.push_back(u[sw]);
    uh.push_back(u[s]);
    uh.push_back(u[se]);
    uh.push_back(u[w]);
    uh.push_back(u[c]);
    uh.push_back(u[e]);
    uh.push_back(u[nw]);
    uh.push_back(u[n]);
    uh.push_back(u[ne]);
    
  }

  else {
    std::cout << "N is too small" << std::endl;
  }

  // Computing Ah * uh
  std::vector<double> AhUh = mul_mat(Ah, uh, 9);

  // Computing bh
  bh[0] = source[sw] + coef * (u[sw-N] + u[sw-1]);
  bh[1] = source[s] + coef * u[s - N];
  bh[2] = source[se] + coef * (u[se-N] + u[se+1]);
  bh[3] = source[w] + coef * u[w - 1];
  bh[4] = source[c];
  bh[5] = source[e] + coef * u[e + 1];
  bh[6] = source[nw] + coef * (u[nw - 1] + u[nw + N]);
  bh[7] = source[n] + coef * u[n + N];
  bh[8] = source[ne] + coef * (u[ne + 1] + u[ne + N]);

  // Computing Ah*uh - bh
  for (int i = 0; i < 9; i++) { rh[i] = AhUh[i] - bh[i];}

  std::cout << "Euclidean error : " << norm_2(rh) << std::endl;
  std::cout << "Infinite error : " << norm_inf(rh) << std::endl;

}

void compute_source_map(std::vector<double>& vec, const program_options::Options& opts) {
  
  int x = 0;
  int y = 0;
  double h = 1 / ((double) opts.N - 1);

  if (opts.has_source) {

    // Computing the gaussian source map
    for (int i = 0; i < opts.N * opts.N; i++) {

      x = i % opts.N;
      y = ((int) i / (int) opts.N);

      vec[i] = gaussian_source(x * h, y * h, opts.source_sigma, opts.source_x, opts.source_y);

    }
  }
}

double norm_2(const std::vector<double>& res) {
  // Return the euclidean norm of a vector
  double sum = 0.0;
  for (const auto& element : res) {
      sum += pow(element, 2);
  }
  return sqrt(sum);
}

double norm_inf(const std::vector<double>& res) {
  // Return the maximum element in a vector
  double maximum = 0.0;

    for (const auto& element : res) {

      double abs = std::abs(element);
      if (abs > maximum) maximum = abs;

    }

  return maximum;
}