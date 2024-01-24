#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>


Eigen::SparseMatrix<double> reader(const std::string &filename, int &N);

int main(int argc, char *argv[]) {
    if (argc != 3){
        std::cout << "Invalid number of arguments" << '\n';
        return 1;
    }

    int max_iterations = std::stoi(argv[2]);
    int N; // We need to use it later
    double tol = 1e-14; // Low tolerance with some buffer for preicsion

    Eigen::SparseMatrix<double> A = reader(argv[1], N);
    Eigen::VectorXd x(N), b, xC, xD;  // ---Initialize
    x = Eigen::VectorXd::Ones(N); // ---Initialize
    b = A.selfadjointView<Eigen::Lower>() * x;  // ---Initialize

    //-------------------------------------------------------------
    
    Eigen::ConjugateGradient<
    Eigen::SparseMatrix<double>, 
    Eigen::Lower, 
    Eigen::DiagonalPreconditioner<double>
    > cgD;
    cgD.compute(A);
    cgD.setTolerance(tol);
    cgD.setMaxIterations(max_iterations);
    xD = cgD.solve(b);
    
    //-------------------------------------------------------------
    
    Eigen::ConjugateGradient<
    Eigen::SparseMatrix<double>, 
    Eigen::Lower, 
    Eigen::IncompleteCholesky<double>
    > cgC;
    
    cgC.compute(A);
    cgC.setTolerance(tol);
    cgC.setMaxIterations(max_iterations);
    xC = cgC.solve(b); 

    //-------------------------------------------------------------
    
    std::cout<< cgC.iterations()<<","<<cgC.error()<<","<<cgD.iterations()<<","<<cgD.error()<<"\n";
    return 0;
}

Eigen::SparseMatrix<double> reader(const std::string &filename, int &N) {
    //Guard
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()){throw std::runtime_error("ERROR: Unable to open file ");}

    // std::string line;
    // while (getline(inputFile, line) && line[0] == '%') {} // skip comments. For some reason, it gives errors

    std::string header;
    std::getline(inputFile, header);
    int num_cols, num_rows, entries;
    inputFile >> num_cols >> num_rows >> entries;
    N = num_cols;
    
    if (num_cols != num_rows) {// SQUARE ??
        throw std::runtime_error("ERROR: Matrix is not square ");
    }

    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(entries);

    int row, col;
    double val;

    while (inputFile >> row){
        inputFile >> col >> val;
        trips.emplace_back(row - 1, col - 1, val);
    }

    Eigen::SparseMatrix<double> mat(N, N);
    mat.setFromTriplets(trips.begin(), trips.end());
    return mat;
}