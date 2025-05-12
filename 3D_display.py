"""" """
import numpy as np
import sys
import datetime
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation


def generate_sphere_points(n_points):
    """
    Generates approximately uniform points on a unit sphere using the Fibonacci lattice.
    """
    points = []
    offset = 2.0 / n_points
    increment = np.pi * (3.0 - np.sqrt(5.0))
    
    for i in range(n_points):
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - y * y)
        phi = (i * increment) % (2 * np.pi)
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        points.append([x, y, z])
        
    return np.array(points)

def plot_sphere_points(points, radius=1.0):
    """
    Plots the given points on a 3D sphere.
    """
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Scale the points
    points = radius * points
    x, y, z = points[:, 0], points[:, 1], points[:, 2]

    # Plot the points
    ax.scatter(x, y, z, color='blue', s=20)
    ax.set_box_aspect([1, 1, 1])
    ax.set_title(f"{len(points)} Surface Points on Sphere (Radius = {radius})")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

def animate_sphere(points, radius=1.0, frames=120, interval=100):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.set_title(f"{len(points)} Surface Points on Sphere (Radius = {radius})")

    # Scale points
    points = radius * points
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    scatter = ax.scatter(x, y, z, color='blue', s=20)

    def update(frame):
        ax.view_init(elev=20, azim=frame * 3)  # rotate azimuthally
        return scatter,

    ani = FuncAnimation(fig, update, frames=frames, interval=interval, blit=False)
    plt.show()

if __name__ == "__main__":

    # Display the current date of run
    today = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    start_process = time.time()

    NUMBER_OF_POINTS = 92

    # Generate
    print(f"{today}Generating points...")
    points = generate_sphere_points(NUMBER_OF_POINTS)
    end_generate_point = time.time()
    
    plot_opt = int(input('Select an option :\
    \n1 = plot \n2 = animated plot \nOption : '))
    start_display = time.time()
    if plot_opt == 1:
        # Plot 92 points on a sphere
        print("'Plot' option selected")
        plot_sphere_points(points)
    elif plot_opt == 2:
        # Animate 92 points
        print("'Animate' option selected")
        animate_sphere(points)
    else :
        print("Invalid selection. Exiting...")
        sys.exit()
    end_display = time.time()
    print(f"Time to generate {NUMBER_OF_POINTS} points : {end_generate_point - start_process} seconds\
    \nTime of display : {end_display - start_display} seconds")
