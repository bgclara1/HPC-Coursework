import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

with open("positions.txt", "r") as f:
    lines = [line for line in f if not line.startswith("runtime")]

data = np.loadtxt(StringIO("".join(lines)))
num_particles = (data.shape[1] - 1) // 2

plt.figure(figsize=(8, 6))
for i in range(num_particles):
    x = data[:, 1 + 2 * i]
    y = data[:, 1 + 2 * i + 1]
    plt.plot(x, y, label=f'Particle {i}')
plt.xlabel("x position")
plt.ylabel("y position")
plt.title("Trajectories in XY")
plt.legend()
plt.grid(True)
plt.show()
