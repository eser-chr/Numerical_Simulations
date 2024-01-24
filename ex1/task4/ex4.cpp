#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>

using DATA_TYPE = double;  // Numerical
using Vector = std::vector<DATA_TYPE>;
using Matrix = std::vector<Vector>;


DATA_TYPE norm(const Matrix& matrix){
    int dim = matrix.size();
    DATA_TYPE max_col;

    for(int i = 0; i<dim; i++){
        DATA_TYPE sum=0;
        for(int j=0; j<dim; j++)
            sum+= std::abs(matrix[j][i]);

        if (sum>max_col)
            max_col=sum;
    }
    return max_col;
}

Matrix create_matrix(int N){
    Matrix matrix(N);
    for(int i=0; i< N; i++){
        Vector row(N);
        for(int j=0; j<N; j++)
            row[j]=(i+1.0)*j;     

        matrix[i]=row;
    }
    return matrix;
}

void write_val(std::string path, std::string name, DATA_TYPE val){
    std::ofstream file;
    file.open(path, std::ios_base::app);
    file << name;
    file<<",";
    file<<val;
    file<<"\n";
    file.close();
}



int main(int argc, char *argv[]){
    
    int N = 1000;
    int ITERATIONS = 8;
    std::string output_path = "./output.txt";
    
    Matrix matrix = create_matrix(N);

   
    // Could be better designed with functors or function pointers.
    auto start = std::chrono::steady_clock::now();
    double temp;
    for(int i =0 ; i<ITERATIONS; i++){
        temp = norm(matrix);
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;    
    auto avg = duration.count()/ITERATIONS;

    write_val(output_path, (std::string) argv[0], avg);

return 0;

}






















