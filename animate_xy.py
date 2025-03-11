
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

data = np.loadtxt("positions.txt", skiprows=1)
num_particles = (data.shape[1] - 1) // 2

fig, ax = plt.subplots(figsize=(8,6))

# Create markers for each particle (without a trail)
markers = []
for i in range(num_particles):
    (marker,) = ax.plot([], [], 'o', markersize=8, label=f'Particle {i}')
    markers.append(marker)

# Set axis limits based on your data
all_x = data[:, 1::2]
all_y = data[:, 2::2]
ax.set_xlim(np.min(all_x)-1, np.max(all_x)+1)
ax.set_ylim(np.min(all_y)-1, np.max(all_y)+1)
ax.set_xlabel("x position")
ax.set_ylabel("y position")
ax.set_title("Animated Particles in XY")
ax.legend()
ax.grid(True)

def init():
    for marker in markers:
        marker.set_data([], [])
    return markers

def update(frame):
    # Update only the current positions; no trailing history.
    for i in range(num_particles):
        x = data[frame, 1 + 2*i]
        y = data[frame, 1 + 2*i + 1]
        markers[i].set_data([x], [y])
    return markers

ani = animation.FuncAnimation(fig, update, frames=range(0, len(data), 100),
                              init_func=init, interval=20, blit=True)
plt.show()
