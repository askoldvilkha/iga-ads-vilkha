// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "co2_sequestration_2d.hpp"
#include <iostream>
#include <lyra/lyra.hpp>

// UI example for parsing arguments - look into 'iga-ads/examples/cg/main.cpp'

/*
int main() {
    ads::dim_config dim{2, 50}; // {P, N} 2 - number of B-splines (P); higher P - more precision; 2-3 is okay; 40 - size of the mesh (N)
    ads::timesteps_config steps{10000, 1e-4};
    int ders = 1; // order of derivative 

    ads::config_2d c{dim, dim, steps, ders};
    ads::problems::co2_sequestration_2d sim{c};
    sim.run();
}
*/

int main(int argc, char* argv[]) {

    double mesh_x;
    double mesh_y;
    int num_steps;
    double timestep_size;

    bool show_help = false;
    bool verbose = false;

    double mu_w = 1.25;
    double mu_g = 2;
    double K = 1;
    double phi = 1;
    double rho_w = 2;
    double rho_g = 1;
    double g = 1;

    auto cli = lyra::help(show_help)  //
               | lyra::arg(mesh_x, "mesh_size - x")("mesh resolution in the x direction").required()  //
               | lyra::arg(mesh_y, "mesh_size - y")("mesh resolution in the y direction").required()  //
               | lyra::arg(num_steps, "num_steps")("number of time steps").required()  //
               | lyra::arg(timestep_size, "timestep_size")("size of the timestep").required()  //
               | lyra::opt(mu_w, "mu_w")["--mu_w"]("mu_w (parameter) - brine viscosity")  //
               | lyra::opt(mu_g, "mu_g")["--mu_g"]("mu_g (parameter) - gas viscosity")  //
               | lyra::opt(K, "K")["--K"]("K (parameter) - permeability tensor ") //
               | lyra::opt(phi, "phi")["--phi"]("phi (parameter) - porosity") //
               | lyra::opt(rho_w, "rho_w")["--rho_w"]("rho_w (parameter) - density of the brine") //
               | lyra::opt(rho_g, "rho_g")["--rho_g"]("rho_g (parameter) - density of the gas") // 
               | lyra::opt(g, "g")["--g"]("g (parameter) - gravitational acceleration") //
               | lyra::opt(verbose)["--verbose"];

    auto const result = cli.parse({argc, argv});

    if (!result) {
        std::cerr << "Error: " << result.errorMessage() << std::endl;
        std::cerr << cli << std::endl;
        std::exit(1);
    }

    if (show_help) {
        std::cout << cli << std::endl;
        std::exit(0);
    }

    if (verbose) {
        std::cout << "Argument parsing is complete" << std::endl;
        // print the values:
        std::cout << "mesh_size x: " << mesh_x << std::endl;
        std::cout << "mesh_size y: " << mesh_y << std::endl;
        std::cout << "num_steps: " << num_steps << std::endl;
        std::cout << "timestep_size: " << timestep_size << std::endl;
        std::cout << "mu_w: " << mu_w << std::endl;
        std::cout << "mu_g: " << mu_g << std::endl;
        std::cout << "K: " << K << std::endl;
        std::cout << "phi: " << phi << std::endl;
        std::cout << "rho_w: " << rho_w << std::endl;
        std::cout << "rho_g: " << rho_g << std::endl;
        std::cout << "g: " << g << std::endl;

    }

    int n_elem_x;
    int n_elem_y;

    n_elem_x = static_cast<int>(mesh_x);
    n_elem_y = static_cast<int>(mesh_y);

    ads::dim_config dim_x{2, n_elem_x, 0, mesh_x}; // {P, N} 2 - number of B-splines (P); higher P - more precision; 2-3 is okay
    ads::dim_config dim_y{2, n_elem_y, 0, mesh_y}; 
    ads::timesteps_config steps{num_steps, timestep_size};
    int ders = 1; // order of derivative 

    ads::config_2d c{dim_x, dim_y, steps, ders};
    auto sim = ads::problems::co2_sequestration_2d(c, mu_w, mu_g, K, phi, rho_w, rho_g, g, verbose);
    sim.run();
}
