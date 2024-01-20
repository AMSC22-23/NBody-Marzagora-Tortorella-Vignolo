# Import necessary libraries
import matplotlib;
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
import numpy as np
from matplotlib.widgets import Button
import random
import subprocess
import sys


#write a funcion that compiles and runs the simulation
def run_simulation():
    subprocess.call(["g++", "-fopenmp", "-I../utils", "../main/main.cpp", "-o", "main"])
    subprocess.call(["./main"])
#if len(sys.argv) != 2:
#    subprocess.call(["g++", "-fopenmp", "-I../utils", "../main/main.cpp", "-o", "main"])
#    subprocess.call(["./main"])
#else:
    #given a string bounded by "" in the command line, the string is passed as a parameter
if len(sys.argv) != 2:
    subprocess.call(["g++", "-fopenmp", "-I../utils", "../main/main.cpp", "-o", "main"])
    result = subprocess.call(["./main"])
else:
    option = sys.argv[1]
    subprocess.call(["g++", "-fopenmp", "-I../utils", "../main/main.cpp", "-o", "main"])
    command = ["./main"] + option.split(" ")
    result = subprocess.call(command)



#print(command)
#process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#subprocess.call(["./main", option])


#run_simulation()
interval_millisec = 20;  # [milliseconds], duration of the frame

# Open the file for reading
with open('Info.txt', 'r') as f:

    # Read the number of the particles to simulate
    num_particles = int(next(f).strip())

    # Read the dimension of the simulation area
    dim = int(next(f).strip())

    # Read the dimension of the simulation (2D or 3D)
    dimension = int(next(f).strip())

    # Read frames
    frames = int(next(f).strip())

    # Read the number of files
    num_files = int(next(f).strip())

    # Read the radius values from the next num_particles lines of the file
    radius = [float(next(f).strip()) for _ in range(num_particles)]

    f.close()

# Initialize the vectors of vectors to store the coordinates of the particles
x = [[] for _ in range(num_particles)]
y = [[] for _ in range(num_particles)]
z = [[] for _ in range(num_particles)]

# Create a list of file names to open in this 
file_names = ['Coordinates_' + str(i) + '.txt' for i in range(num_files)]
# Create an empty list to store file objects
file_objs = []

# Loop through the list of file names to open each file and append the file object to the file_objs list
for file_name in file_names:
    file_obj = open(file_name, 'r')
    file_objs.append(file_obj)

for i in range (num_files): 
    # Read the remaining lines of the file to get the particle coordinates
    if(dimension == 2):
        for line in file_objs[i]:
            # Split the line into id, x, and y
            id, x_val, y_val = map(float, line.strip().split(','))
            x[int(id)].append(x_val)
            y[int(id)].append(y_val)
    else:
        for line in file_objs[i]:
            # Split the line into id, x, y, anx z
            id, x_val, y_val, z_val = map(float, line.strip().split(','))
            x[int(id)].append(x_val)
            y[int(id)].append(y_val)
            z[int(id)].append(z_val)

#print(x)

# Close all the file objects
for file_obj in file_objs:
    file_obj.close()

if(dimension == 2):
    print("2D animation")

    #Create a new figure and axes
    fig, ax = plt.subplots()

    # Create a list of Circle objects (one for each particle)
    circles_2D = [Circle((0, 0), radius[i], color=(random.random(), random.random(), random.random())) for i in range(len(radius))]

    # Add the circles to the axes
    for circle in circles_2D:
        ax.add_patch(circle)

    # Create a list of Line2D objects (one for each particle)
    lines = [ax.plot([], [], marker='o', markersize=radius[i])[0] for i in range(len(x))]

    def init():
        for i, circle in enumerate(circles_2D):
            circle.center = (x[i][0], y[i][0])
        return circles_2D

    # Set the limits of x and y axes
    ax.set_xlim([-dim, dim])
    ax.set_ylim([-dim, dim])

    def animate(i):
        for j, circle in enumerate(circles_2D):
            circle.center = (x[j][i], y[j][i])  # Update the position of the circle
        return circles_2D

    def reset_animation(event):
        ani.frame_seq = ani.new_frame_seq()
        ani.event_source.start()

    def stop_animation(event):
        ani.event_source.stop()

    def resume_animation(event):
        ani.event_source.start

    # Draw restart buttons
    reser_button_ax = plt.axes([0.8, 0.9, 0.1, 0.05])
    reser_button_ax = Button(reser_button_ax, 'Restart', color='lightgoldenrodyellow', hovercolor='0.975')
    reser_button_ax.on_clicked(reset_animation)

    # Draw stop buttons
    stop_button_ax = plt.axes([0.8, 0.025, 0.1, 0.04])
    stop_button = Button(stop_button_ax, 'Stop', color='lightgoldenrodyellow', hovercolor='0.975')
    stop_button.on_clicked(stop_animation)

    # Create animation
    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=interval_millisec, blit=True)
    plt.show()

else:

    print("3D animation")
    fig_3D = plt.figure()
    ax_3D = plt.axes(projection="3d")

    # Set the limits of x, y, and z axes
    ax_3D.set_xlim([-dim, dim])
    ax_3D.set_ylim([-dim, dim])
    ax_3D.set_zlim([-dim, dim])
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Function to plot a sphere
    def plot_sphere(ax, center, radius, color):
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x1 = center[0] + radius * np.outer(np.cos(u), np.sin(v))
        y2 = center[1] + radius * np.outer(np.sin(u), np.sin(v))
        z3 = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
        return ax.plot_surface(x1,y2,z3,color=color)

    def init_3D():
        return fig_3D

    surface_color = "tab:blue"
    def animate_3D(i):
        ax_3D.clear()
        ax_3D.set_xlim([-dim, dim])
        ax_3D.set_ylim([-dim, dim])
        ax_3D.set_zlim([-dim, dim])
        for j in range(len(radius)):
            x1 = x[j][i] + radius[j] * np.outer(np.cos(u), np.sin(v))
            y1 = y[j][i] + radius[j] * np.outer(np.sin(u), np.sin(v))
            z1 = z[j][i] + radius[j] * np.outer(np.ones(np.size(u)), np.cos(v))
            ax_3D.plot_surface(x1, y1, z1, color=surface_color)
        return fig_3D

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
    ani = animation. FuncAnimation(fig_3D, animate_3D, init_func=init_3D, frames=frames, interval=interval_millisec)
    plt.show()
