# Import necessary libraries
import matplotlib; matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
import numpy as np
from matplotlib.widgets import Button

duration = 5; # 5 sec = duration of the simulation
interval_millisec = 20;  # 20 milliseconds = duration of the frame
frames = int(duration/(interval_millisec*0.001)); # calculate the number of frames given the duration of the simulation

# Open the file for reading
with open('coordinates.txt', 'r') as f:
    num_particles = int(next(f).strip())
    # Read the radius values from the first N lines of the file
    radius = [float(next(f).strip()) for _ in range(num_particles)]
    # Read the remaining lines of the file to get the particle coordinates
    print(num_particles)
    # Initialize the vectors of vectors
    x = [[] for _ in range(num_particles)]
    y = [[] for _ in range(num_particles)]
    print(x)
    # Read the rest of the file line by line
    for line in f:
        # Split the line into id, x, and y
        id, x_val, y_val = map(float, line.strip().split(','))
        # Append the x and y values to the appropriate vectors
        x[int(id)].append(x_val)
        y[int(id)].append(y_val)

# Create a new figure and axes
fig, ax = plt.subplots()

t = np.arange(0.0, duration, 0.01)

# Create a list of Line2D objects (one for each body)
lines = [ax.plot([], [], marker='o', markersize=radius[i])[0] for i in range(len(x))]

# Initialization function
# This function is called once before the animation starts
# It clears the data of each line
def init():
    for line in lines:
        line.set_data([], [])
    return lines

# Create a slider axes
# ax_slider = plt.axes([0.1, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
# slider = Slider(ax_slider, 'Time', 0, frames, valinit=0, valstep=1)

# Set the limits of x and y axes
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])


# Animation function
# This function is called for each frame of the animation
# It updates the data of each line to plot only the current point
def animate(i):
    for j, line in enumerate(lines):
        line.set_data(x[j][i:i+1], y[j][i:i+1])  # Only plot the current point
    return lines

# Slider update function
def update(i):
    print("start from frame:")
    print(i)
    for j, line in enumerate(lines):
        line.set_data(x[j][i:i+1], y[j][i:i+1])
    return lines

def reset_animation(event):
    ani.frame_seq = ani.new_frame_seq()
    ani.event_source.start()

def stop_animation(event):
    ani.event_source.stop()

def resume_animation(event):
    # your code here
    #ani.frame_seq = ani.new_frame_seq()
    ani.event_source.start

reser_button_ax = plt.axes([0.8, 0.9, 0.1, 0.05])
reser_button_ax = Button(reser_button_ax, 'Start', color='lightgoldenrodyellow', hovercolor='0.975')
reser_button_ax.on_clicked(reset_animation)

stop_button_ax = plt.axes([0.8, 0.025, 0.1, 0.04])
stop_button = Button(stop_button_ax, 'Stop', color='lightgoldenrodyellow', hovercolor='0.975')
stop_button.on_clicked(stop_animation)

# Create animation
# The FuncAnimation function creates an animation by repeatedly calling a function (in this case, animate)
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=interval_millisec, blit=True)

# Connect the slider to the update function
#slider.on_changed(update)

# Display the animation
plt.show()