import pde
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import h5py

size = 32


grid = pde.CartesianGrid([[0, 4*np.pi], [0, 2*np.pi]], [128, 32])
field = pde.ScalarField.from_expression(grid, 'sin(x) * cos(y)')

laplace_dir = field.laplace({'value': 0})
laplace_dir.plot(title='Laplacian of field with Dirichlet boundary conditions', colorbar=True);