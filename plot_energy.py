import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

with open("energy.txt", "r") as f:
    lines = [line for line in f if not line.startswith("runtime")]

data = np.loadtxt(StringIO("".join(lines)))

time = data[:, 0]
E_history = data[:, 1:]  

plt.figure(figsize=(8, 6))
num_particles = E_history.shape[1]
for i in range(num_particles):
    plt.plot(time, E_history[:, i], label=f'Particle {i}')

plt.xlabel("Time")
plt.ylabel("Kinetic Energy")
plt.title("Kinetic Energy vs Time for Each Particle")
plt.legend()
plt.grid(True)
plt.show()
