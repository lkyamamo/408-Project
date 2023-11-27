import numpy as np

#setup
x = np.linspace(-5, 5, 100)  # Generate x values
y = np.linspace(-5, 5, 100)  # Generate y values

# Create a grid of (x, y) values
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)

print(X[0][0])