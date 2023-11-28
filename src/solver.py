import numpy as np
from dolfin import *

def get_potential_superposition(self, X, Y, Z):
    V_total = np.zeros_like(X)
    for charge in self.charges:
        V_total += self.constants.Coulomb_Constant * charge.charge / np.sqrt((X - charge.position[0])**2 + (Y - charge.position[1])**2 + (Z - charge.position[2])**2)
    return V_total

#finds solution to the poisson equation given distribution rho
def get_potential_fenics(self, rho):
    # Create mesh and define function space
    mesh = UnitSquareMesh(32, 32)
    V = FunctionSpace(mesh, "Lagrange", 1)

    # Define Dirichlet boundary (x = 0 or x = 1)
    def boundary(x):
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

    # Define boundary condition
    u0 = Constant(0.0)
    bc = DirichletBC(V, u0, boundary)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Expression(rho, degree=2)
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)

def get_potential_multipole(self):
    