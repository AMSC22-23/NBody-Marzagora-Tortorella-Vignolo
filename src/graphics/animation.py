# Import necessary libraries
import matplotlib; matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
from matplotlib.widgets import Slider
import numpy as np
from matplotlib.widgets import Button
import random
import subprocess

#do a function that comnpile and execute a c++ program
def compile_and_execute():
    subprocess.call(["g++", "-std=c++11", "-fopenmp", "-o", "main", "../main/main.cpp"])
    subprocess.call(["./main"])

compile_and_execute()
duration = 5; # duration of the simulation
interval_millisec = 50;  # [milliseconds], duration of the frame
frames = int(duration/(interval_millisec*0.001)); # Calculate the number of frames

# Open the file for reading
with open('coordinates.txt', 'r') as f:
    # Read the number of the particles to simulate
    num_particles = int(next(f).strip())
    # Read the dimension of the simulation area
    dim = int(next(f).strip())
    # Read the radius values from the next num_particles lines of the file
    radius = [float(next(f).strip()) for _ in range(num_particles)]
    # Initialize the vectors of vectors to store the coordinates of the particles
    x = [[] for _ in range(num_particles)]
    y = [[] for _ in range(num_particles)]
    # Read the remaining lines of the file to get the particle coordinates
    for line in f:
        # Split the line into id, x, and y
        id, x_val, y_val = map(float, line.strip().split(','))
        # Append the x and y values to the appropriate vectors
        x[int(id)].append(x_val)
        y[int(id)].append(y_val)

# Create a new figure and axes
fig, ax = plt.subplots()

# Create a list of Circle objects (one for each particle)
circles = [Circle((0, 0), radius[i], color=(random.random(), random.random(), random.random())) for i in range(len(radius))]

# Add the circles to the axes
for circle in circles:
    ax.add_patch(circle)

# Create a list of Line2D objects (one for each particle)
lines = [ax.plot([], [], marker='o', markersize=radius[i])[0] for i in range(len(x))]

# Initialization function
# This function is called once before the animation starts
# It sets the initial position of each circle
def init():
    for i, circle in enumerate(circles):
        circle.center = (x[i][0], y[i][0])
    return circles

# Create a slider axes
# ax_slider = plt.axes([0.1, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
# slider = Slider(ax_slider, 'Time', 0, frames, valinit=0, valstep=1)

# Set the limits of x and y axes
ax.set_xlim([-dim, dim])
ax.set_ylim([-dim, dim])

# Animation function
# This function is called for each frame of the animation
# It updates the position of each circle
def animate(i):
    for j, circle in enumerate(circles):
        circle.center = (x[j][i], y[j][i])  # Update the position of the circle
    return circles

# Slider update function
# def update(i):
#     print("start from frame:")
#     print(i)
#     for j, line in enumerate(lines):
#         line.set_data(x[j][i:i+1], y[j][i:i+1])
#     return lines

def reset_animation(event):
    ani.frame_seq = ani.new_frame_seq()
    ani.event_source.start()

def stop_animation(event):
    ani.event_source.stop()

def resume_animation(event):
    ani.event_source.start

reser_button_ax = plt.axes([0.8, 0.9, 0.1, 0.05])
reser_button_ax = Button(reser_button_ax, 'Start', color='lightgoldenrodyellow', hovercolor='0.975')
reser_button_ax.on_clicked(reset_animation)

stop_button_ax = plt.axes([0.8, 0.025, 0.1, 0.04])
stop_button = Button(stop_button_ax, 'Stop', color='lightgoldenrodyellow', hovercolor='0.975')
stop_button.on_clicked(stop_animation)

# Create animation
# The FuncAnimation function creates an animation by repeatedly calling a function (in this case, animate)
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=1000, interval=10, blit=True)

# Connect the slider to the update function
#slider.on_changed(update)

# Display the animation
plt.show()
