import numpy as np
import sys
import matplotlib.pyplot as plt
from PIL import Image

sys.path.insert(0, '/workspaces/408-Project/src')
import systems

D = 2
N = 300
q = 100 #unit charge
dt = 0.01
radius = 1
omega = np.pi

"""
# create system
dipole_system = systems.Spinning_Dipole(q, radius, omega, dt)

#setup
x = np.linspace(-5, 5, 100)  # Generate x values
y = np.linspace(-5, 5, 100)  # Generate y values

# Create a grid of (x, y) values
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)

#iterate over time using point charge objects
total_time = .1
current_time = 0
while current_time < total_time:
    dipole_system.step()
    potential_field_func = dipole_system.get_potential_field(X,Y,Z)  # Create a callable potential field function
    heatmap = plt.contourf(X,Y,potential_field_func, levels=50, cmap='viridis')
    plt.colorbar(heatmap, label = 'Heat Map Value')
    current_time += dt



#export figures as a gif

#export frames
iterations = (int) (total_time/dt)
for i in range(iterations):
    plt.figure(i + 1)
    plt.savefig(f"frames/frame_{i:03d}.png")
    plt.close()

#stitch frames
frames = []
for i in range(iterations):
    frame = Image.open(f"frames/frame_{i:03d}.png")
    frames.append(frame)

# Save as a GIF
frames[0].save("animation.gif", save_all=True, append_images=frames[1:], duration=200, loop=0)


"""


#iterate over time using continuous distribution
epsilon = 0.001

#https://en.wikipedia.org/wiki/Dirac_delta_function
weak_dirac_delta = f"exp(-((x-(x0))**2 + (y-(y0))**2)/(2*{epsilon}))/({epsilon}*sqrt(2*pi))"

# + exp(-((x-(x1))**2 + (y-(y1))**2)/(2*{epsilon}))/({epsilon}*sqrt(2*pi))

dynamics = [[f"{radius}*cos({omega}*t)",f"{radius}*sin({omega}*t)"], [f"-{radius}*cos({omega}*t)",f"-{radius}*sin({omega}*t)"]]

#using $ as place holders for positions of particles
continuous_dipole_system = systems.Continuous_Distribution_System(weak_dirac_delta, dynamics, "cartesian", dt)

continuous_dipole_system.get_potential_pypde()