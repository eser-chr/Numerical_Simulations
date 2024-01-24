#include "Benchmark.h"

Benchmark::Benchmark(int vectorSize) : vectorSize(vectorSize), A(vectorSize), B(vectorSize), C(vectorSize), D(vectorSize), elapsedTime(0.0) {
}

void Benchmark::run() {
    auto start = std::chrono::high_resolution_clock::now();
    //Vector Traid Calculations
    for (int i = 0; i < vectorSize; ++i) {
        A[i] = B[i] + C[i] * D[i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e+06;
}

double Benchmark::getFlopsPerSecond() {
    // Copmuting FLOPS/SEC.
    return (static_cast<double>(vectorSize) / elapsedTime) / 1e+09;
}