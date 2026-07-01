from orbitlab import solve_kepler
import math
E=solve_kepler(0.5,0.3)
assert abs((E-0.3*math.sin(E))-0.5)<1e-9
