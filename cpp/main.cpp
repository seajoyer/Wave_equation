#include "wave_solver.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

// Function to execute system command and get output
std::string exec_command(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

// Function to check if gnuplot is available
bool is_gnuplot_available() {
    try {
        std::string result = exec_command("which gnuplot");
        return !result.empty();
    } catch (const std::runtime_error&) {
        return false;
    }
}

WaveEquationSolver run_solver(int order) {
    double h = 0.05;  // spatial step
    double tau = 0.05;  // temporal step
    int x_points = static_cast<int>(1.0/h) + 1;
    int t_points = static_cast<int>(1.0/tau) + 1;

    WaveEquationSolver solver(x_points, t_points, order);

    auto initial_displacement = [](double x) {
        return x - x * tan(x/2);
    };

    auto initial_velocity = [](double x) {
        return 1 + x/(1 + cos(x));
    };

    solver.set_initial_conditions(initial_displacement, initial_velocity);
    solver.solve();

    return solver;
}

void create_gnuplot_script() {
    std::ofstream gnuplot_script("plot_solutions.gnu");
    gnuplot_script
    << "set term x11 persist size 1200,800\n"
    << "set multiplot layout 3,2 title 'Wave Equation Solutions'\n"
    << "set grid xtics ytics ls 3\n"
    << "set grid lt 1 lc rgb '#dddddd'\n"
    << "set style line 1 lc rgb '#0000FF' lt 1 lw 2\n"
    << "set style line 2 lc rgb '#FF0000' lt 1 lw 2\n"
    << "set style line 3 lc rgb '#00FF00' lt 1 lw 2\n"
    << "set xtics 0.2\n"
    << "set ytics 0.2\n"
    << "set xrange [0:1]\n"
    << "set title 'First-order Numerical Solution'\n"
    << "set ylabel 't'\n"
    << "set yrange [0:1]\n"
    << "set view map\n"
    << "set palette defined (0 '#000080', 0.5 '#00FF00', 1 '#FFFF00')\n"
    << "set cbrange [0:1.6]\n"
    << "plot 'first_order_numerical.dat' using 1:2:3 with image notitle\n"
    << "set title 'Second-order Numerical Solution'\n"
    << "plot 'second_order_numerical.dat' using 1:2:3 with image notitle\n"
    << "set title 'First-order Error'\n"
    << "set palette defined (0 '#000080', 0.5 '#00FF00', 1 '#FFFF00')\n"
    << "set cbrange [0:0.35]\n"
    << "plot 'first_order_error.dat' using 1:2:3 with image notitle\n"
    << "set title 'Second-order Error'\n"
    << "plot 'second_order_error.dat' using 1:2:3 with image notitle\n"
    << "unset view\n"
    << "set origin 0.0,0.0\n"
    << "set size 1.0,0.34\n"
    << "set title 'Second-order solution at t = 0.50'\n"
    << "set ylabel 'u(x,t)'\n"
    << "set yrange [0.4:1.3]\n"
    << "set key bottom right\n"
    << "plot 'first_order_comparison.dat' using 1:2 title 'First-order' with lines ls 1, \\\n"
    << "     'second_order_comparison.dat' using 1:2 title 'Second-order' with lines ls 2, \\\n"
    << "     'first_order_comparison.dat' using 1:3 title 'Exact' with lines dashtype 2 lc rgb '#FF0000'\n"
    << "unset multiplot\n";
}

int main() {
    // Check if gnuplot is available
    if (!is_gnuplot_available()) {
        std::cerr << "Error: gnuplot is not installed or not in PATH\n";
        return 1;
    }

    // Run both solvers
    WaveEquationSolver solver1 = run_solver(1);
    WaveEquationSolver solver2 = run_solver(2);

    // Save data for plotting
    solver1.save_solution_data("first_order");
    solver2.save_solution_data("second_order");

    // Print error norms
    auto [l2_1, l1_1, linf_1] = solver1.get_error_norms();
    auto [l2_2, l1_2, linf_2] = solver2.get_error_norms();

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "First-order error norms:\n"
              << "L2 norm: " << l2_1 << "\n"
              << "L1 norm: " << l1_1 << "\n"
              << "L∞ norm: " << linf_1 << "\n\n";

    std::cout << "Second-order error norms:\n"
              << "L2 norm: " << l2_2 << "\n"
              << "L1 norm: " << l1_2 << "\n"
              << "L∞ norm: " << linf_2 << "\n";

    // Create and run gnuplot script
    create_gnuplot_script();

    if (std::system("gnuplot plot_solutions.gnu") != 0) {
        std::cerr << "Error: Failed to run gnuplot\n";
        return 1;
    }

    return 0;
}
