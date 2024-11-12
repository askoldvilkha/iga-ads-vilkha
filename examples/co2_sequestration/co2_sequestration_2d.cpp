// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "co2_sequestration_2d.hpp"

// UI example for parsing arguments - look into 'iga-ads/examples/cg/main.cpp'

int main() {
    ads::dim_config dim{2, 50}; // {P, N} 2 - number of B-splines (P); higher P - more precision; 2-3 is okay; 40 - size of the mesh (N)
    ads::timesteps_config steps{10000, 1e-4};
    int ders = 1; // order of derivative 

    ads::config_2d c{dim, dim, steps, ders};
    ads::problems::co2_sequestration_2d sim{c};
    sim.run();
}

/*
#include "co2_sequestration_2d.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <mesh_size> <num_steps> <timestep_size>\n";
        return 1;
    }

    int mesh_size = std::stoi(argv[1]);
    int num_steps = std::stoi(argv[2]);
    double timestep_size = std::stod(argv[3]);

    ads::dim_config dim{2, mesh_size}; // {P, N} 2 - number of B-splines (P); higher P - more precision; 2-3 is okay
    ads::timesteps_config steps{num_steps, timestep_size};
    int ders = 1; // order of derivative 

    ads::config_2d c{dim, dim, steps, ders};
    ads::problems::co2_sequestration_2d sim{c};
    sim.run();
}
*/