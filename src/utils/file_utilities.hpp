#include <vector>
#include <fstream>
#include <iostream>
#include <random>
#include "particle.hpp" 


/**
 * @brief Template function that writes on file the total number of particles
 *and the size of the area of the simulation and then the initial state of the
 *particles in a file.
 *
 * @tparam Dimension Number of dimensions of the simulation
 * @param particles Reference to the vector of particles that are to be printed
 *on the file
 * @param dim Dimension of the simulation area
 * @param fileName Name of the file in which the initial states of the particles
 *are going to be written
 * @param file Reference to the file in which the initial states of the
 *particles are going to be written
 * @param it Number of iterations
 **/

template <size_t Dimension>
void printInitialStateOnFile(std::vector<Particle<Dimension>>* particles, int dim, std::string fileName,
                             std::ofstream& file, int it, int speedUp) {
    if (file.is_open()) {
        file << (*particles).size() << std::endl;
        file << dim << std::endl;
        file << Dimension << std::endl;
        file << it / speedUp << std::endl;

        for (const Particle<Dimension>& p : (*particles)) file << p.getRadius() << std::endl;
        for (const Particle<Dimension>& p : (*particles)) {
            file << p.getId() << ",";
            const auto& pos = p.getPos();
            for (size_t i = 0; i < Dimension; ++i) {
                file << pos[i];
                if (i < Dimension - 1) file << ",";
            }
            file << std::endl;
        }
    } else {
        std::cout << "Unable to open file in printInitialStateOnFile";
    }
}

// file utilities

template <size_t Dimension>
void writeParticlePositionsToFile(const std::vector<Particle<Dimension>>& particles, std::ofstream& file) {
    for (const auto& q : particles) {
        if (file.is_open()) {
            file << q.getId() << ",";
            const auto& pos = q.getPos();
            for (size_t i = 0; i < Dimension; ++i) {
                file << pos[i];
                if (i < Dimension - 1) {
                    file << ",";
                }
            }
            file << std::endl;
        } else {
            std::cout << "Unable to open file";
            break; 
        }
    }
}