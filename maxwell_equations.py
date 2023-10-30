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
size = 10
q = 1 #Coulombs

x = np.arange(-size/2,size/2,dr)
y = np.arange(-size/2,size/2,dr)
z = np.arange(-size/2,size/2,dr)

X,Y,Z = np.meshgrid(x,y,z)

#Gauss's Law of Electric Fields
rho = np.zeros_like(X)
rho[X==1 and Y == 0 and Z == 0] = q
rho[X==-1 and Y == 0 and Z == 0] = -q