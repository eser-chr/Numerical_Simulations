import numpy as np
import matplotlib.pyplot as plt

# size of matrices
N = [64, 128, 256, 500, 512, 1000, 1024, 1500]

# Peak Performance
peak_perf = 4.6 * 10**9 * 4 * 2

# runtimes in microseconds obtain by running the executable mmm for different command line argument "impl" and "size"
run_time_CUSTOM = [408,3618, 57639, 297939, 429750, 3556548, 3963387, 18469496]
run_time_BLAS = [129, 279, 2074, 6811, 9828, 51388, 66887, 155115]
run_time_EIGEN = [113, 647, 2247, 17007, 14941, 115108, 135694, 308144]

# theoretical runtime for peak performance of the CPU (in microseconds)
run_time_peak_performance = [((2*n**3) / (peak_perf))*10**6 for n in N]

# Plotting
plt.figure(figsize=(10, 6))
plt.semilogy(N, run_time_CUSTOM, label='Custom Implementation', linestyle='--', marker='o')
plt.semilogy(N, run_time_EIGEN, label='EIGEN Library', linestyle='--', marker='o')
plt.semilogy(N, run_time_BLAS, label='BLAS Library', linestyle='--', marker='o')
plt.semilogy(N, run_time_peak_performance, label='Theoretical runtime', linestyle='--', marker='o')

plt.xlabel('Matrix Size [N]')
plt.ylabel('Runtime [$\mu s$] in log-scale')
plt.title('Benchmark test (Runtime Comparison)')
plt.legend()
plt.grid(True)

# Save the plot as a PNG file
plt.savefig('runtime_comparison.png')

# Display the plot
plt.show()

