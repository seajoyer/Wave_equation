import numpy as np
import matplotlib.pyplot as plt
from wave_solver import WaveEquationSolver

def run_solver(order):
    h = 0.05
    tau = 0.05
    x_points = int(1/h) + 1
    t_points = int(1/tau) + 1

    solver = WaveEquationSolver(x_points, t_points, order=order)

    def initial_displacement(x):
        return x - x * np.tan(x/2)

    def initial_velocity(x):
        return 1 + x/(1 + np.cos(x))

    solver.set_initial_conditions(initial_displacement, initial_velocity)
    solver.solve()

    return solver

def main():
    # first and second orders
    solver1 = run_solver(1)
    solver2 = run_solver(2)

    fig = plt.figure(figsize=(15, 15))

    X, T = np.meshgrid(solver1.x, solver1.t)
    exact_solution = np.array([[solver1.get_exact_solution(x, t)
                              for x in solver1.x]
                              for t in solver1.t])

    # Plot solutions
    plt.subplot(321)
    plt.pcolormesh(X, T, solver1.u, shading='auto')
    plt.colorbar(label='u(x,t)')
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title('First-order Numerical Solution')

    plt.subplot(322)
    plt.pcolormesh(X, T, solver2.u, shading='auto')
    plt.colorbar(label='u(x,t)')
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title('Second-order Numerical Solution')

    # Plot errors
    error1 = np.abs(solver1.u - exact_solution)
    error2 = np.abs(solver2.u - exact_solution)

    vmax = max(np.max(error1), np.max(error2))

    plt.subplot(323)
    plt.pcolormesh(X, T, error1, shading='auto', vmin=0, vmax=vmax)
    plt.colorbar(label='|Error|')
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title('First-order Error')

    plt.subplot(324)
    plt.pcolormesh(X, T, error2, shading='auto', vmin=0, vmax=vmax)
    plt.colorbar(label='|Error|')
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title('Second-order Error')

    t_index = solver1.t_points // 2
    t_value = solver1.t[t_index]

    plt.subplot(313)
    plt.plot(solver1.x, solver1.u[t_index, :], 'b-', label='First-order')
    plt.plot(solver2.x, solver2.u[t_index, :], 'g-', label='Second-order')
    plt.plot(solver1.x, exact_solution[t_index, :], 'r--', label='Exact')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.title(f'Solution comparison at t = {t_value:.2f}')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
