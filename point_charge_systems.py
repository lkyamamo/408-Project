import numpy as np
import matplotlib.pyplot as plt

class Constants:
    @property
    def Coulomb_Constant(self):
        return 8.99 * 10**9
    
size = 10

class Point_Charge:
    def __init__(self, charge, position):
        self.charge = charge
        self.position = position
        

class Point_Charge_System:
    def __init__(self, charges):
        self.charges = charges
        self.constants = Constants()

    def get_electric_field(self):
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
    
    def get_potential_field(self):
        def potential_field(x,y,z):
            r = np.array([x,y,z])
            V = 0
            for q in self.charges:
                distance = np.linalg.norm(r - q.position)
                result = 0
                if distance != 0:
                    result = self.constants.Coulomb_Constant * q.charge / distance
                V += result
            return V
        return potential_field
    
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
