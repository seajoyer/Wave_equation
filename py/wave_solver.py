import numpy as np
from typing import Callable, Tuple

class WaveEquationSolver:
    def __init__(
            self,
            x_points: int,
            t_points: int,
            order: int = 2,
            x_max: float = 1.0,
            t_max: float = 1.0
    ):
        self.x_points = x_points
        self.t_points = t_points
        self.x_max = x_max
        self.t_max = t_max
        self.order = order

        self.dx = x_max / (x_points - 1)
        self.dt = t_max / (t_points - 1)

        self.x = np.linspace(0, x_max, x_points)
        self.t = np.linspace(0, t_max, t_points)

        self.u = np.zeros((t_points, x_points))

    def set_initial_conditions(
            self,
            initial_displacement: Callable[[float], float],
            initial_velocity: Callable[[float], float]
    ):
        # initial displacement
        for j in range(self.x_points):
            self.u[0, j] = initial_displacement(self.x[j])

        # first time step
        for j in range(1, self.x_points-1):
            self.u[1, j] = self.u[0, j] + self.dt * initial_velocity(self.x[j])

            # correction based on order
            if self.order >= 2:
                # Second-order correction
                d2u_dx2 = (self.u[0, j+1] - 2*self.u[0, j] + self.u[0, j-1]) / (self.dx**2)
                source = self.source_term(self.x[j], 0)
                self.u[1, j] += 0.5 * self.dt**2 * (d2u_dx2 + source)

        # Apply boundary conditions
        self.apply_boundary_conditions(1)

    def apply_boundary_conditions(self, t_idx: int):
        t = self.t[t_idx]

        if self.order == 1:
            # First-order boundary conditions
            self.u[t_idx, 0] = self.u[t_idx, 1] - self.dx*(1 + np.tan(t/2))       # Left

            rhs = t + 1/(1 + np.cos(t-1))
            self.u[t_idx, -1] = (self.u[t_idx, -2] + self.dx*rhs) / (1 + self.dx) # Right
        else:
            # Second-order boundary conditions
            self.u[t_idx, 0] = (4*self.u[t_idx, 1] - self.u[t_idx, 2]) / 3 - \
                2*self.dx*(1 + np.tan(t/2)) / 3                                   # Left

            rhs = t + 1/(1 + np.cos(t-1))
            self.u[t_idx, -1] = (4*self.u[t_idx, -2] - self.u[t_idx, -3] + 2*self.dx*rhs) / (3 + 2*self.dx) # Right

    def source_term(self, x: float, t: float) -> float:
        denom = 1 + np.cos(t - x)
        return (2 + 2*np.cos(t - x) + x*np.sin(t - x)) / (denom * denom)

    def solve(self):
        r = (self.dt / self.dx)**2
        epsilon = 0.01

        for i in range(1, self.t_points-1):
            for j in range(1, self.x_points-1):
                # wave equation discretization
                u_next = 2*self.u[i, j] - self.u[i-1, j] + \
                    r*(self.u[i, j+1] - 2*self.u[i, j] + self.u[i, j-1])

                u_next += self.dt**2 * self.source_term(self.x[j], self.t[i])

                u_next += epsilon * (self.u[i, j+1] - 2*self.u[i, j] + self.u[i, j-1])

                self.u[i+1, j] = u_next

            self.apply_boundary_conditions(i+1)

    def get_exact_solution(self, x: float, t: float) -> float:
        return t + x + x * np.tan((t - x) / 2)

    def get_error_norms(self) -> Tuple[float, float, float]:
        exact_solution = np.array([[self.get_exact_solution(x, t)
                                  for x in self.x]
                                 for t in self.t])
        error = self.u - exact_solution

        l2_norm = np.sqrt(np.mean(error**2))
        l1_norm = np.mean(np.abs(error))
        linf_norm = np.max(np.abs(error))

        return l2_norm, l1_norm, linf_norm
