import numpy as np
import matplotlib.pyplot as plt
from fenics import *

class Constants:
    @property
    def Coulomb_Constant(self):
        return 8.99 * 10**9
    
    @property
    def Electric_Constant(self):
        return 8.85 * 10**-12
    
    @property
    def Magnetic_Constant(self):
        return 1.26 * 10**-6
    
dr = 0.1
dt = 0.01 
size = 1000
q = 1 #Coulombs

mesh = UnitSquareMesh(size,size)
V = FunctionSpace(mesh, "Lagrange", 1)

def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

# Create mesh and define function space
nx, ny = 100, 100  # adjust as needed
mesh = UnitSquareMesh(nx, ny)
V = VectorFunctionSpace(mesh, 'P', 1)

# Define Maxwell's equations
E = TrialFunction(V)
H = TestFunction(V)

mu_0 = Constant(4 * pi * 1e-7)
epsilon_0 = Constant(8.854e-12)

# Magnetic field H formulation (curl-curl)
a_H = dot(curl(H), curl(E)) * dx
L_H = Constant(0) * H * dx

# Electric field E formulation (curl-curl)
a_E = dot(curl(E), curl(H)) * dx
L_E = Constant(0) * H * dx

# Assemble matrices and vectors
A_H, rhs_H = assemble_system(a_H, L_H, bc)
A_E, rhs_E = assemble_system(a_E, L_E, bc)

# Solve the linear systems
H_solution = Function(V)
E_solution = Function(V)

solve(A_H, H_solution.vector(), rhs_H)
solve(A_E, E_solution.vector(), rhs_E)


plot(E_solution, title="Electric Field (E)")
plt.show()

# Plot the magnetic field solution
plot(H_solution, title="Magnetic Field (H)")
plt.show()

# Plot vector fields using quiver
def plot_vector_field(solution, title):
    mesh = solution.function_space().mesh()
    X = mesh.coordinates()
    U = solution.compute_vertex_values(mesh)

    # Extract x and y components of the vector field
    Ux = U[0::2]
    Uy = U[1::2]

    # Create a quiver plot
    fig, ax = plt.subplots()
    ax.quiver(X[:, 0], X[:, 1], Ux, Uy, scale=20, color='b', width=0.005)
    ax.set_aspect('equal')
    ax.set_title(title)
    plt.show()

# Plot electric field vector
plot_vector_field(E_solution, title="Electric Field (E)")

# Plot magnetic field vector
plot_vector_field(H_solution, title="Magnetic Field (H)")
