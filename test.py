import numpy as np
import matplotlib.pyplot as plt

#setup
x = np.linspace(-5, 5, 100)  # Generate x values
y = np.linspace(-5, 5, 100)  # Generate y values

# Create a grid of (x, y) values
X, Y = np.meshgrid(x, y)

Z = np.sin(np.sqrt(X**2 + Y**2))

# Create a contour map
contour = plt.contour(X, Y, Z, levels=20, cmap='viridis')

# Add a color bar
plt.colorbar(contour)

# Add labels and a title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Contour Map Example')
plt.savefig(f"Contour_Map_Example.png")