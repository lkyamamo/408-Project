import numpy as np
import matplotlib.pyplot as plt

#setup
x = np.linspace(-5, 5, 100)  # Generate x values
y = np.linspace(-5, 5, 100)  # Generate y values

# Create a grid of (x, y) values
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)

print(X)
position = np.array([1,0,0])

V = np.zeros_like(X)

X0 = np.full_like(X, position[0])
Y0 = np.full_like(Y, position[1])
Z0 = np.full_like(Z, position[2])
r =  np.array([X-X0, Y-Y0, Z-Z0])



# Create a contour map
contour = plt.contour(X, Y, 1, levels=20, cmap='viridis')

# Add a color bar
plt.colorbar(contour)

# Add labels and a title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Contour Map Example')
plt.savefig(f"Contour_Map_Example.png")