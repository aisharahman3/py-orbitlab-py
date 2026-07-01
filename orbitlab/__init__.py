import math
from .constants import G

def gravitational_force(m1, m2, r):
    return G * m1 * m2 / r ** 2

def period(a, m_central):
    return 2 * math.pi * math.sqrt(a ** 3 / (G * m_central))

def circular_velocity(r, m_central):
    return math.sqrt(G * m_central / r)

def escape_velocity(r, m_central):
    return math.sqrt(2 * G * m_central / r)

def vis_viva(r, a, m_central):
    return math.sqrt(G * m_central * (2.0 / r - 1.0 / a))

def solve_kepler(M, e, tol=1e-12, max_iter=60):
    """Solve M = E - e*sin(E) for the eccentric anomaly via Newton-Raphson."""
    E = M if e < 0.8 else math.pi
    for _ in range(max_iter):
        f = E - e * math.sin(E) - M
        E -= f / (1 - e * math.cos(E))
        if abs(f) < tol:
            break
    return E

def true_anomaly(M, e):
    E = solve_kepler(M, e)
    return 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2),
                          math.sqrt(1 - e) * math.cos(E / 2))

def hohmann_dv(r1, r2, m_central):
    mu = G * m_central
    a_t = (r1 + r2) / 2
    dv1 = math.sqrt(mu / r1) * (math.sqrt(2 * r2 / (r1 + r2)) - 1)
    dv2 = math.sqrt(mu / r2) * (1 - math.sqrt(2 * r1 / (r1 + r2)))
    return dv1, dv2, dv1 + dv2
