from io import StringIO
import numpy as np
import matplotlib.pyplot as plt

with open("energy.txt", "r") as f:
    lines = [line for line in f if not line.startswith("runtime")]

data = np.loadtxt(StringIO("".join(lines)))
time = data[:, 0]
E_history = data[:, 1:]

plt.figure(figsize=(8,6))
plt.plot(time, E_history[:, 0], label='Particle 0')
plt.xlabel("Time")
plt.ylabel("Kinetic Energy")
plt.title("Kinetic Energy vs Time for Each Particle")
plt.legend()
plt.grid(True)
plt.show()
