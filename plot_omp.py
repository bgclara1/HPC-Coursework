import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Runtime Data
cores = np.array([1, 2, 4, 8, 12, 16, 32, 34, 36, 38, 40, 42, 46, 48])
runtimes = np.array([631.391, 576.238, 349.671, 211.066, 150.298,
                     120.99, 69.6578, 65.7052, 62.2505, 79.3753,
                     72.8861, 72.6114, 72.5133, 82.1919])

# Improved Rational Model with upward tail
def rational_tail_model(x, a, b, c, d):
    return (a + b * x + c * x**2) / (1 + d * x)

# Fit
params, _ = curve_fit(rational_tail_model, cores, runtimes, maxfev=10000)
a, b, c, d = params

# Smooth curve for plotting
x_smooth = np.linspace(1, max(cores), 500)
y_smooth = rational_tail_model(x_smooth, a, b, c, d)

# Plot
plt.figure(figsize=(12,6))
plt.plot(cores, runtimes, 'rx', markersize=8, label='Measured Runtime Data')
plt.plot(x_smooth, y_smooth, 'b-', linewidth=2.5,
         label=f'Rational Fit: (a + b·x + c·x²) / (1 + d·x)\n'
               f'a={a:.2f}, b={b:.2f}, c={c:.4f}, d={d:.4f}')

plt.title("OMP Runtime vs Number of Cores ", fontsize=14)
plt.xlabel("Number of Cores", fontsize=12)
plt.ylabel("Runtime (seconds)", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(cores)
plt.legend(fontsize=10)
plt.tight_layout()
plt.show()
