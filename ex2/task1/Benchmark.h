#pragma once

#include <vector>
#include <chrono>

class Benchmark {
    public:
        Benchmark(int vectorSize);
        void run();
        double getFlopsPerSecond();

    private:
        int vectorSize;
        std::vector<double> A, B, C, D;
        double elapsedTime;
};