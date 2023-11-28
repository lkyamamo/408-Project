import numpy as np
import matplotlib.pyplot as plt
import solver as sol

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
    #each 
    def get_potential_field(self, X, Y, Z):
        result = sol.get_potential_superposition(self, X, Y, Z)
        return result
        
    
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

class Continuous_Distribution_System(System):
    def __init__(self, distribution, basis):
        self.distribution = distribution
        self.basis = basis