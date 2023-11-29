from pde import CartesianGrid, ScalarField, solve_poisson_equation
from dolfin import *

epsilon = 0.00001
#https://en.wikipedia.org/wiki/Dirac_delta_function
weak_delta = f"exp(- ((x-0.25)**2 + (y-0.25)**2)/(2*{epsilon}))/({epsilon}*sqrt(2*pi)) + exp(- ((x-0.75)**2 + (y-0.75)**2)/(2*{epsilon}))/({epsilon}*sqrt(2*pi))"

grid = CartesianGrid([[0, 1],[0, 1]], 50, periodic=False)
field = ScalarField.from_expression(grid, weak_delta)
result = solve_poisson_equation(field, bc=[[{"value": 0}, {"value": 0}], [{"value": 0}, {"value": 0}]])

result.plot()