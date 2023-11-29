import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
from pde import CartesianGrid, ScalarField, solve_poisson_equation
from sympy import *


class Constants:
    @property
    def Coulomb_Constant(self):
        return 8.99 * 10**9
    
size = 10

class Point_Charge:
    def __init__(self, charge, position):
        self.charge = charge
        self.position = position
        
class System:
    def __init__(self):
        self.constants = Constants()


class Point_Charge_System(System):
    def __init__(self, charges):
        self.charges = charges
        

    def get_electric_field(self, X, Y, Z):
        def electric_field(r):
            E = np.array([0.0,0.0,0.0])
            for q in self.charges:
                for i in range(3):
                    try:
                        result = self.constants.Coulomb_Constant * q.charge / (r[i] - q.position[i])**2
                    except:
                        result = 0
                    E[i] += result
            return E
        return electric_field
    
    #exports a numpy array with the same shape as X which spans the whole space
    #directly calculate 
    def get_potential_superposition(self, X, Y, Z):
        V_total = np.zeros_like(X)
        for charge in self.charges:
            V_total += self.constants.Coulomb_Constant * charge.charge / np.sqrt((X - charge.position[0])**2 + (Y - charge.position[1])**2 + (Z - charge.position[2])**2)
        return V_total

    #only getting up to dipole term
    def get_potential_multipole(self):
        dipole = [(q.charge * q.position) for q in self.charges]

    #

    #directly calculate
    def get_electric_field_superposition(self, X,Y,Z):
        #create array for components of E field
        E = [0,0,0]
        for i in range(len(E)):
            E[i] = np.zeros_like(X)
            for charge in self.charges:
                E[i] += self.constants.Coulomb_Constant * charge.charge / ((X - charge.position[0])**2 + (Y - charge.position[1])**2 + (Z - charge.position[2])**2)
        return E

    #calculate by taking 
    def get_electric_field_from_potential(self, X, Y, Z):
        i = 1

    
    """
    #returns magnetic field as a function of a numpy 
    def get_magnetic_field(self):
        def magnetic_field(r):
            B = np.array([0.0,0.0,0.0])
            for q in self.charges:
                for i in range(3):
                    try:
                        result = 
                    except:
                        result = 0
                    B
            return B
    """
            
class Spinning_Dipole(Point_Charge_System):
    def __init__(self, charge, radius, omega, dt):
        super().__init__([Point_Charge(charge, np.array([radius, 0,0])), Point_Charge(-charge, np.array([-radius, 0, 0]))])
        self.omega = omega
        self.dt = dt
        self.r = radius
    
    def step(self):
        rotation_matrix = np.array([[np.cos(self.omega * self.dt), -np.sin(self.omega * self.dt), 0],
                            [np.sin(self.omega * self.dt), np.cos(self.omega * self.dt), 0],
                            [0,0,1]])
        for charge in self.charges:
            charge.position = np.dot(rotation_matrix, charge.position)
        plt.figure()
        plt.xlim(-size/2,size/2)
        plt.ylim(-size/2,size/2)
        current_positions = self.get_positions_2D()
        plt.scatter(current_positions[0], current_positions[1])

    def get_positions_2D(self):
        positions = [[],[]]
        for i in range(len(self.charges)):
            positions[0].append(self.charges[i].position[0])
            positions[1].append(self.charges[i].position[1])
        return positions

class Dynamics:
    def __init__(self, expression, basis):
        self.expression = expression
        self.basis = basis

    def evaluate(self):
        
        

class Distribution:
    def __init__(self, expression, dynamics, basis):
        self.expression = expression
        self.basis = basis
        self.dynamics = Dynamics(dynamics, basis)

        instance = expression.find('$')
        while(instance != -1):
            self.positions.append(None)
    
    #returns expression 
    def get_expression(self):
        if self.positions[0] == None:
            raise "input necessary for expression"
        
        i = 0
        instance = self.expression.find('$')
        while(instance != -1):
            self.expression[instance] = self.positions[i]
            i += 1
            instance = self.expression.find('$')

class Continuous_Distribution_System(System):
    def __init__(self, distribution, dynamics, basis):
        self.distribution = Distribution(distribution, dynamics, basis)
        self.basis = basis

    #finds solution to the poisson equation given distribution rho
    def get_potential_fenics(self):
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
        f = Expression(self.distribution, degree=2)
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx

        # Compute solution
        u = Function(V)
        solve(a == L, u, bc)


    def get_potential_pypde(self):
        #https://en.wikipedia.org/wiki/Dirac_delta_function
        weak_delta = self.distribution

        grid = CartesianGrid([[0, 1],[0, 1]], 50, periodic=False)
        field = ScalarField.from_expression(grid, weak_delta)
        result = solve_poisson_equation(field, bc=[[{"value": 0}, {"value": 0}], [{"value": 0}, {"value": 0}]])

        result.plot()
