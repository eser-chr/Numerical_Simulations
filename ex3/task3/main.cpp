#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <random>

// Function to integrate
double function(const char* func, double x) {

    if (std::strcmp(func, "SINX") == 0) { return std::sin(x); }
    else if (std::strcmp(func, "COS2XINV") == 0) { return pow(std::cos(1 / x), 2); }
    else if (std::strcmp(func, "X4M5") == 0) { return 5.0 * std::pow(x, 4); }
    else { return 0; }
}

// Monte Carlo integration
double monteCarloIntegration(const char* func, int numSamples, double lowerBound, double upperBound, std::mt19937& exclusiveRNG) {

    double sum = 0.0;

    for (int i = 0; i < numSamples; ++i) {
        // Generate a random x within the bounds - With uniform distribution between lowerBound and upperBound
        double randomX = lowerBound + ((double)exclusiveRNG() / exclusiveRNG.max()) * (upperBound - lowerBound);
        // Evaluate the function at the random x and add it to the sum
        sum += function(func, randomX);
    }

    // Calculate the average value and multiply by the range
    double average = sum / numSamples;
    return average * (upperBound - lowerBound);
}

int main(int argc, char* argv[]) {

    // Getting arguments
    if (argc != 5) {
        std::cerr << "Not the right amount of argument." << std::endl;
        return 1;
    }

    const char* func = argv[1];
    double xmin = std::stod(argv[2]);
    double xmax = std::stod(argv[3]);
    int numSamples = std::stoi(argv[4]);

    std::cout << "\n" << std::endl;
    std::cout << func << " / " << xmin << " / " << xmax << " / " << numSamples << std::endl;;

    // Specify the number of threads
    int numThreads = omp_get_max_threads();

    // Start measuring time
    auto start_time = std::chrono::high_resolution_clock::now();

    // Perform Monte Carlo integration in parallel
    double totalResult = 0.0;

    // Single shared RNG to generate the seeds used for the exclusive RNG of each thread
    std::mt19937 sharedRNG(1234); // Use a fixed seed for the shared RNG: 1234
    std::uniform_int_distribution<int> seed_dist(0,1000);
    // Parallelize the computation of monte-carlo algorithm 
    #pragma omp parallel
    {

        int thread_Seed = seed_dist(sharedRNG)
        // Exclusive RNG for each thread with a seed generated by the sharedRNG
        std::mt19937 exclusiveRNG(thread_Seed);

        // Computing the result of the thread()
        double threadResult = monteCarloIntegration(func, numSamples / numThreads, xmin, xmax, exclusiveRNG); 

        // Get the threadID to be able to print it
        int threadID = omp_get_thread_num();
        
        // Critical on the next line so different threads do not update "totalResult" and print at the same time
        #pragma omp critical
        {
            // Output the approximation of the integral computed by the thread
            std::cout << "Thread " << threadID << ": " << threadResult << std::endl;
            totalResult += threadResult;
        }
    }

    // Stop measuring time
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    // Output the final approximation of the integral, runtime, number of threads, and number of total samples
    std::cout << "Final Approximation of the Integral: " << (double) (totalResult / numThreads) << std::endl;
    std::cout << "Runtime of the Integration " << std::fixed << std::setprecision(5) << duration.count() << " microseconds" << std::endl;
    std::cout << "Number of Threads used: " << numThreads << std::endl;
    std::cout << "Number of Total Samples: " << numSamples << std::endl;

    return 0;
}