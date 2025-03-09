import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("energy.txt", skiprows=1)
time = data[:, 0]
E_history = data[:, 1:]

plt.figure(figsize=(8,6))
for i in range(E_history.shape[1]):
    plt.plot(time, E_history[:, i], label=f'Particle {i}')

plt.xlabel("Time")
plt.ylabel("Kinetic Energy")
plt.title("Kinetic Energy vs Time for Each Particle")
plt.legend()
plt.grid(True)
plt.show()
