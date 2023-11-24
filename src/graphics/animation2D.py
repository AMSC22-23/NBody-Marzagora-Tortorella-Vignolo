# Import necessary libraries
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Number of bodies and time steps
N = 2
T = 100

# Generate a time array from 0 to 10 with a step of 0.01
t = np.arange(0, 10, 0.01)

# Generate x and y coordinates for each body
# Each body follows the same path in this example
x = [list(t) for _ in range(N)]
y = [list(t) for _ in range(N)]

# Create a new figure and axes
fig, ax = plt.subplots()

# Create a list of Line2D objects (one for each body)
# These will be used to plot the positions of the bodies
lines = [ax.plot([], [], marker='o')[0] for _ in range(len(x))]

# Initialization function
# This function is called once before the animation starts
# It clears the data of each line
def init():
    for line in lines:
        line.set_data([], [])
    return lines

# Animation function
# This function is called for each frame of the animation
# It updates the data of each line to plot only the current point
def animate(i):
    for j, line in enumerate(lines):
        line.set_data(x[j][i:i+1], y[j][i:i+1])  # Only plot the current point
    return lines

# Set the limits of x and y axes
ax.set_xlim([0, 3])
ax.set_ylim([0, 3])

# Create animation
# The FuncAnimation function creates an animation by repeatedly calling a function (in this case, animate)
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=2*T, interval=20, blit=True)

# Display the animation
plt.show()
