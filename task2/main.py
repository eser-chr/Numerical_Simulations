import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("input.txt", delimiter=",")

years = data[:,0]
cpu_count = data[:,1]
perf_number = data[:,2]

# Computing the log of the cpu count
L_cpu_count = np.log(cpu_count)

# Computing the linear interpolation 
slope, intercept = np.polyfit(years, L_cpu_count, 1)

L_cpu_count_pred = slope * years + intercept

sst = np.sum((L_cpu_count - np.mean(L_cpu_count))**2)

# Calculating the sum of residual square
srs = np.sum((L_cpu_count - L_cpu_count_pred)**2)

# Calculating the coefficient of determination R^2
r2 = 1 - (srs / sst)

# The coefficient of determination is close to 1, so it justifies the use of the log scale
print("Coefficient of determination : r2 = ", r2)
print("As the coefficient of determination is close to 1, the use of the log-scale is justified")

# Ploting the graph for the number of transistors
plt.bar(data[:,0], data[:,1])
plt.xlabel('Year')
plt.ylabel('Log of the number of transitors')
plt.title("Evolution of the number of transistors of CPU as a function of the year (LOG scale)")
plt.yscale('log')
plt.show()

# Ploting the graph for the performance number
plt.bar(data[:,0], data[:,2])
plt.xlabel('Year')
plt.ylabel('Log of the number of transitors')
plt.title("Evolution of the performance of CPU as a function of the year (LOG scale)")
plt.yscale('log')
plt.show()