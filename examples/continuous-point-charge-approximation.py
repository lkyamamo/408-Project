import numpy as np
import sys
import matplotlib.pyplot as plt
from PIL import Image
import pde

sys.path.insert(0, '/workspaces/408-Project/src')
import systems

plt.ioff()

D = 2
N = 300
q = 100 #unit charge
dt = 0.01
radius = 1
omega = 3* np.pi

#iterate over time using continuous distribution
epsilon = 0.000001

#https://en.wikipedia.org/wiki/Dirac_delta_function
weak_dirac_delta = f"{q}*exp(- ((x-(x0))**2 + (y-(y0))**2)/(2*{epsilon}))/({epsilon}*sqrt(2*pi)) + {-1*q}*exp(- ((x-(x1))**2 + (y-(y1))**2)/(2*{epsilon}))/({epsilon}*sqrt(2*pi))"

#second term to create spinning dipole however there is a bug where including this causes the evaluation to overflow
# + exp(- ((x-(x1))**2 + (y-(y1))**2)/(2*{epsilon}))/({epsilon}*sqrt(2*pi))

#circular motion
dynamics = [[f"{radius}*cos({omega}*t)",f"{radius}*sin({omega}*t)"], [f"-{radius}*cos({omega}*t)",f"-{radius}*sin({omega}*t)"]]

continuous_dipole_system = systems.Continuous_Distribution_System(weak_dirac_delta, dynamics, "cartesian", dt)

#iterate over time saving plots of 
total_time = .1
current_time = 0
collection = []
i = 0
while current_time < total_time:
    continuous_dipole_system.step()
    result = continuous_dipole_system.get_potential_pypde()
    result.plot(kind="image", action="close", filename=f"frames/frame_{i:03d}.png")
    current_time += dt
    i += 1

#stitch frames
iterations = (int) (total_time/dt)
frames = []
for i in range(iterations):
    frame = Image.open(f"frames/frame_{i:03d}.png")
    frames.append(frame)

# Save as a GIF
frames[0].save("continuous-point-charge-distribution-animation.gif", save_all=True, append_images=frames[1:], duration=200, loop=0)