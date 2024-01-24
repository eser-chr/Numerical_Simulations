//main.cpp
#include<iostream>
#include<fstream>
#include "Benchmark.h"

int main() {
    //For adjusting vector size based on the system's cache
    const int vectorSize[] = {128, 256, 512, 1024, 2048, 4096, 8192, 16384, 20000};
    const int numSize = sizeof(vectorSize) / sizeof(vectorSize[0]);
    std::ofstream outputFile("benchmark_results.csv");
    outputFile << "Vector Size FLOPS/SEC." <<std::endl;

    for (int i = 0; i <numSize; ++i) {
        Benchmark Benchmark(vectorSize[i]);
        Benchmark.run();
        double flopsPerSecond = Benchmark.getFlopsPerSecond();

        std::cout << "vector Size: " << vectorSize[i] << ", FLOPS/SEC.: " << flopsPerSecond << std::endl;
        outputFile << vectorSize[i] << "," << flopsPerSecond << std::endl;
    }

    outputFile.close();

    return 0;
}