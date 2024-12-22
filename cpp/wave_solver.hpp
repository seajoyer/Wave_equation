#pragma once

#include <vector>
#include <functional>
#include <cmath>
#include <fstream>
#include <string>
#include <tuple>
#include <algorithm>

class WaveEquationSolver {
public:
    WaveEquationSolver(int x_points, int t_points, int order = 2,
                      double x_max = 1.0, double t_max = 1.0)
        : x_points(x_points), t_points(t_points), order(order),
          x_max(x_max), t_max(t_max) {

        dx = x_max / (x_points - 1);
        dt = t_max / (t_points - 1);

        // Initialize grids
        x.resize(x_points);
        t.resize(t_points);
        u.resize(t_points, std::vector<double>(x_points, 0.0));

        // Set up coordinate vectors
        for (int i = 0; i < x_points; ++i) {
            x[i] = i * dx;
        }
        for (int i = 0; i < t_points; ++i) {
            t[i] = i * dt;
        }
    }

    void set_initial_conditions(
        const std::function<double(double)>& initial_displacement,
        const std::function<double(double)>& initial_velocity) {

        // Set initial displacement
        for (int j = 0; j < x_points; ++j) {
            u[0][j] = initial_displacement(x[j]);
        }

        // Set first time step
        for (int j = 1; j < x_points-1; ++j) {
            u[1][j] = u[0][j] + dt * initial_velocity(x[j]);

            // Add correction based on order
            if (order >= 2) {
                double d2u_dx2 = (u[0][j+1] - 2*u[0][j] + u[0][j-1]) / (dx*dx);
                double source = source_term(x[j], 0);
                u[1][j] += 0.5 * dt*dt * (d2u_dx2 + source);
            }
        }

        apply_boundary_conditions(1);
    }

    void solve() {
        double r = (dt/dx) * (dt/dx);
        double epsilon = 0.01;  // artificial viscosity coefficient

        for (int i = 1; i < t_points-1; ++i) {
            // Interior points
            for (int j = 1; j < x_points-1; ++j) {
                // Standard wave equation discretization
                double u_next = 2*u[i][j] - u[i-1][j] +
                              r*(u[i][j+1] - 2*u[i][j] + u[i][j-1]);

                // Add source term
                u_next += dt*dt * source_term(x[j], t[i]);

                // Add artificial viscosity
                u_next += epsilon * (u[i][j+1] - 2*u[i][j] + u[i][j-1]);

                u[i+1][j] = u_next;
            }

            apply_boundary_conditions(i+1);
        }
    }

    std::tuple<double, double, double> get_error_norms() const {
        double l2_sum = 0.0;
        double l1_sum = 0.0;
        double linf_norm = 0.0;
        int total_points = x_points * t_points;

        for (int i = 0; i < t_points; ++i) {
            for (int j = 0; j < x_points; ++j) {
                double exact = get_exact_solution(x[j], t[i]);
                double error = std::abs(u[i][j] - exact);

                l2_sum += error * error;
                l1_sum += error;
                linf_norm = std::max(linf_norm, error);
            }
        }

        double l2_norm = std::sqrt(l2_sum / total_points);
        double l1_norm = l1_sum / total_points;

        return {l2_norm, l1_norm, linf_norm};
    }

    double get_exact_solution(double x, double t) const {
        return t + x + x * tan((t - x) / 2);
    }

    void save_solution_data(const std::string& prefix) const {
        // Save numerical solution
        std::ofstream num_file(prefix + "_numerical.dat");
        for (int i = 0; i < t_points; ++i) {
            for (int j = 0; j < x_points; ++j) {
                num_file << x[j] << " " << t[i] << " " << u[i][j] << "\n";
            }
            num_file << "\n";  // Empty line between t slices for gnuplot
        }

        // Save exact solution
        std::ofstream exact_file(prefix + "_exact.dat");
        for (int i = 0; i < t_points; ++i) {
            for (int j = 0; j < x_points; ++j) {
                exact_file << x[j] << " " << t[i] << " "
                          << get_exact_solution(x[j], t[i]) << "\n";
            }
            exact_file << "\n";
        }

        // Save error data
        std::ofstream error_file(prefix + "_error.dat");
        for (int i = 0; i < t_points; ++i) {
            for (int j = 0; j < x_points; ++j) {
                double exact = get_exact_solution(x[j], t[i]);
                double error = std::abs(u[i][j] - exact);
                error_file << x[j] << " " << t[i] << " " << error << "\n";
            }
            error_file << "\n";
        }

        // Save solution comparison at middle time
        int mid_t = t_points / 2;
        std::ofstream comp_file(prefix + "_comparison.dat");
        for (int j = 0; j < x_points; ++j) {
            double exact = get_exact_solution(x[j], t[mid_t]);
            comp_file << x[j] << " " << u[mid_t][j] << " " << exact << "\n";
        }
    }

private:
    void apply_boundary_conditions(int t_idx) {
        double t_val = t[t_idx];

        if (order == 1) {
            // First-order boundary conditions
            // Left boundary (Neumann)
            u[t_idx][0] = u[t_idx][1] - dx*(1 + tan(t_val/2));

            // Right boundary (Robin)
            double rhs = t_val + 1/(1 + cos(t_val-1));
            u[t_idx][x_points-1] = (u[t_idx][x_points-2] + dx*rhs) / (1 + dx);
        } else {
            // Second-order boundary conditions
            // Left boundary (Neumann)
            u[t_idx][0] = (4*u[t_idx][1] - u[t_idx][2]) / 3 -
                          2*dx*(1 + tan(t_val/2)) / 3;

            // Right boundary (Robin)
            double rhs = t_val + 1/(1 + cos(t_val-1));
            u[t_idx][x_points-1] = (4*u[t_idx][x_points-2] - u[t_idx][x_points-3] + 2*dx*rhs) /
                                  (3 + 2*dx);
        }
    }

    double source_term(double x, double t) const {
        double denom = 1 + cos(t - x);
        return (2 + 2*cos(t - x) + x*sin(t - x)) / (denom * denom);
    }

    int x_points;
    int t_points;
    int order;
    double x_max;
    double t_max;
    double dx;
    double dt;

    std::vector<double> x;
    std::vector<double> t;
    std::vector<std::vector<double>> u;
};
