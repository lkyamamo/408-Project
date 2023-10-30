import numpy as np
import point_charge_systems as pcs
import matplotlib.pyplot as plt
from PIL import Image

D = 2
N = 300
q = 1 #unit charge
dt = 0.01

# create system
dipole_system = pcs.Spinning_Dipole(q, 1, np.pi, dt)

#setup
x = np.linspace(-5, 5, 100)  # Generate x values
y = np.linspace(-5, 5, 100)  # Generate y values

# Create a grid of (x, y) values
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)

#iterate over time
total_time = .1
current_time = 0
while current_time < total_time:
    dipole_system.step()
    #potential_field_func = dipole_system.get_potential_field()  # Create a callable potential field function
    #potential_field_values = potential_field_func(X,Y,Z)
    plt.contour([X,Y],1, cmap='viridis')
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