#include <string>
#include <iostream>
#include <vector>
#include <fstream>

using DATA_TYPE = double;  // Numerical
using Vector = std::vector<DATA_TYPE>;
using Matrix = std::vector<Vector>;


int get_dimension(std::string path){
    int lines=0;
    std::ifstream file(path);
    std::string line;
    while(std::getline(file, line)){
        lines++;
    }
    file.close();
    return lines;
}


Matrix read_from_txt (std::string path){
    int dim = get_dimension(path);
    std::ifstream file(path);

    Matrix matrix(dim);
    DATA_TYPE temp;

    for(int i=0; i<dim; ++i){
        Vector row(dim);
        for(int j=0; j<dim; j++){
            file>>temp;
            row[j]=temp;
        }
        matrix[i]=row;
    }

    file.close();

    return matrix;

}

void print_matrix(const Matrix& matrix){
    int dim = matrix.size();
    for(int i=0; i<dim; i++){
        for(int j =0; j<dim; j++){
            std::cout<<matrix[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}

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

DATA_TYPE trace(const Matrix& matrix){
    int dim = matrix.size();
    DATA_TYPE sum=0;
    for(int i=0; i<dim; i++)
        sum += matrix[i][i];    
    return sum;
}






int main(){
    std::cout<<"Hello \n";
    Matrix matrix = read_from_txt("./MatrixA.txt");
    print_matrix(matrix);
    std::cout<<"The norm is: "<<norm(matrix)<<"\n";
    std::cout<<"The trace is: "<<trace(matrix)<<"\n";

    return 0;
}

