#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>

using Vector = std::vector<double>;

double operator*(Vector v1, Vector v2) { // Dangerous move, but for our purposes its fine!
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

class SymCCS {
    private:
        Vector values; 
        std::vector<int> row_indices; 
        std::vector<int> col_ptrs; 
        int N; 
    
    public:
        int get_N() const { return N;}
       
        SymCCS(const std::string& filename) { // Constructor that populates the structure by reading the .mtx file
            // Guard
            std::ifstream inputFile(filename);
            if (!inputFile.is_open()) {
                throw std::invalid_argument("ERROR: Unable to open file ");
                            
            }

            std::string line;
            int num_rows, num_cols, entries;

            while (getline(inputFile, line) && line[0] == '%') {} // skip comments
            
            std::istringstream headerStream(line);               // Read header
            headerStream >> num_rows >> num_cols >> entries;
            
            if (num_cols != num_rows) {                         // Check if the matrix is square
                std::cerr << "ERROR: Matrix is not square." << std::endl;
                return;
            }

            values.resize(entries);
            row_indices.resize(entries);
            col_ptrs.resize(num_cols + 1, 0);

            // MAIN READING
            for (int i = 0; i < entries; ++i) {
                int row, col;
                double value;
                inputFile >> row >> col >> value;
                
                if(row >= col) {//Ensures to store only lower part
                    values[i] = value;
                    row_indices[i] = row - 1;
                    col_ptrs[col]++;
                }
            }            
            std::partial_sum(col_ptrs.begin(), col_ptrs.end(), col_ptrs.begin()); //CUMSUM            
            N = col_ptrs.size()-1;
        }

        Vector SMVM(const Vector& vector) const {     // Symmetric Matrix Vector Multiplication (SMVM)   
            if (vector.size() != static_cast<size_t>(N)) { // Guard
                std::cerr << "Invalid vector size for multiplication." << std::endl;
                return Vector();
            }

            Vector result(N, 0.0);

            for (int col = 0; col < N; ++col) {
                for (int i = col_ptrs[col]; i < col_ptrs[col + 1]; ++i) {
                    int row = row_indices[i];
                    result[row] += values[i] * vector[col];
                    if (row != col) { // Symmetry
                        result[col] += values[i] * vector[row];
                    }
                }
            }
            return result;
        }

        void printError (const Vector & x, const Vector & r, const Vector & b, int maxIterations){
            double rK_r0 = sqrt(r*r /(b*b));
            Vector xStar(N, 1.0);
            Vector ek(N, 0.0);
            for (size_t i = 0; i < N; ++i)
                    ek[i] = xStar[i] - x[i];
            double error_A_norm = sqrt(ek*SMVM(ek));
            std::cout << maxIterations <<  "," << rK_r0 << "," << error_A_norm << std::endl; // So later it would be stored in a csv file
        }

        // CG algorithm -- Acoording to prof. Faustmann it should have been a 5 liner. Now I am dissappointed.
        void CG(const Vector& b, Vector& x, int maxIterations, double tol=1e-20) {
            Vector r = b; // Initial residual
            Vector p = r; // Initial serach direction
            Vector Ap(N, 0.0);
            double alpha, beta, rDotRnew;
            double rDotR = r*r;
            
            // MAIN
            for (int k = 0; k < maxIterations && std::sqrt(rDotR) > tol; ++k) { // Check also for convergence                
                Ap = SMVM(p);
                alpha = rDotR / (p*Ap);

                for (size_t i = 0; i < N; ++i) 
                    x[i] += alpha * p[i];                

                for (size_t i = 0; i < N; ++i)
                    r[i] -= alpha * Ap[i];            

                rDotRnew = r*r;
                beta = rDotRnew / rDotR;

                for (size_t i = 0; i < N; ++i) 
                    p[i] = r[i] + beta * p[i];

                rDotR = rDotRnew;
            }
            printError(x, r, b, maxIterations);            
        }
};   

int main(int argc, char* argv[]) {
    if (argc != 3) {std::cerr << "Invalid num of arguments!\n"; return 1;}
    std::string filepath = argv[1];
    int iterations = std::stoi(argv[2]);

    std::string inv_iter = "Iterations should be positive and non zero\n";
    if(iterations<=0){std::cerr<< inv_iter; return 2;}
    
    SymCCS matrix(filepath);
    std::vector<double> x0(matrix.get_N(), 0.0);
    std::vector<double> xStar(matrix.get_N(), 1.0);
    std::vector<double> b = matrix.SMVM(xStar);
    
    
    matrix.CG(b, x0, iterations); 
    
    return 0;   
}