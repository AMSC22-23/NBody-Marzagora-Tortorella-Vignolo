#include <vector>
#include "../../utils/particle.hpp"
#include "../../utils/force.hpp"
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <omp.h>

// Function that randomly generates particles:

template<size_t Dimension>
std::vector<Particle<Dimension>> generateRandomParticlesSerial(int N, int posBoundary = 100, int minProperty = 1, int maxProperty = 99, int maxVx = 100, int maxVy = 100, int minRadius = 0, int maxRadius = 15, bool type = false) {
    std::vector<Particle<Dimension>> particles;

    // Initialize random seed
    srand(time(0));

    bool uniquePosition = false;
    double x, y;
    int r;

    for (int i = 0; i < N; i++) {

        bool uniquePosition = false;
        //double x, y;
        std::array<double, Dimension> pos;

        while (!uniquePosition) {
            uniquePosition = true;
            
            // Generate random radius between minRadius and maxRadius
            r = rand() % (maxRadius - minRadius + 1) + minRadius;

            // Generate random position between -posBoundary+r and +posBoundary-r
            //x = -posBoundary + r + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*posBoundary - 2*r)));
            //y = -posBoundary + r + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*posBoundary - 2*r)));
            for(size_t i=0; i<Dimension; ++i)
                pos[i] = -posBoundary + r + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*posBoundary - 2*r)));

            // Check if the position is unique
            for (const Particle<Dimension> &p : particles) {
                //TODO: replace method pow(...) with manual multiplication 
                double distance = sqrt(p.getPos()[0] - pos[0] *  p.getPos()[0] - pos[0] + p.getPos()[1] - pos[1] * p.getPos()[1] - pos[1]);
                if (distance < p.getRadius() + r) {
                    uniquePosition = false;
                    break;
                }
            }
        }

        // Generate random value of property between 1 and 100
        // should use the `random` library instead of `rand()`
        // see: https://en.cppreference.com/w/cpp/numeric/random
        double property = minProperty + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(maxProperty-1)));

        std::array<double, Dimension> vel;
        // Generate random velocity between -100 and 100
        for(size_t i=0; i<Dimension; ++i)
        //double vx = -maxVx + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVx)));
        //double vy = -maxVy + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVy)));
            vel[i] = -maxVy + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVy)));

        // Create a new particle with the random value of property, position, and velocity
        Particle<Dimension> p(i, property, pos, vel, r, type);

        // Add the particle to the vector
        particles.push_back(p);
    }

    return particles;
}


template<size_t Dimension>
void parallelSimulation(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f){
    
    // Start of the simulation
    //#pragma omp parallel
    {
        // Set the number of threads
        omp_set_num_threads(omp_get_max_threads());
        //create a vector of thread that contains a vector of forces of dimension Dimension
        std::vector<std::vector<std::array<double, Dimension>>> local_forces(omp_get_max_threads(), std::vector<std::array<double, Dimension>>(particles->size()));
        std::array<double, Dimension> force;

        for (int z = 0; z < it; ++z){

            //calculate forces and add it to a local array for each thread
            #pragma omp for
            for (int i = 0; i < (*particles).size(); ++i) {
                for(int y = 0; y < Dimension; ++y)
                    force[y] = 0;
                Particle<Dimension> &k = (*particles)[i];
                for(int j = i+1; j < (*particles).size(); ++j){
                    Particle<Dimension> &q = (*particles)[j];
                    force = f.calculateForce(k, q);
                    for (int y = 0; y < Dimension; ++y) {
                        local_forces[omp_get_thread_num()][i][y] += force[y];
                        local_forces[omp_get_thread_num()][j][y] -= force[y];
                    }
                }
            }
            #pragma omp barrier

            //sum all the forces calculated by each thread
            #pragma omp for
            for (int i = 0; i < (*particles).size(); ++i) {
                Particle<Dimension> &q = (*particles)[i];
                q.resetForce();
                for (int j = 0; j < omp_get_max_threads(); ++j) {
                    q.addForce(local_forces[j][i]);
                }
            }
            #pragma omp barrier

            #pragma omp for
            for(int j = 0; j < (*particles).size(); ++j){
                for (int y = 0; y < Dimension; ++y) {
                    local_forces[omp_get_thread_num()][j][y] = 0.0;
                }
            }
            #pragma omp barrier

            //update the position of the particles
            #pragma omp single
            for (int i = 0; i < (*particles).size(); ++i) {
                Particle<Dimension> &q = (*particles)[i];
                q.update(delta_t);
            }
            #pragma omp barrier

            //don't do this in parallel
            #pragma omp single
            {
                // Write on file the updates after delta_t
                for (int i = 0; i < (*particles).size(); i++) {
                    Particle<Dimension> &q = (*particles)[i];
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
                        std::cout << "Unable to open file in parallelSimulation ";
                    }
                }
            }
            #pragma omp barrier
        }
    }
} 

int main() {
    
    // Define of simulation variables
    const int d = 2; //2D or 3D
    const double delta_t = 0.01; // in seconds
    const double dim = 50; // Dimension of the simulation area
    int it = 1000; // number of iteration
    int n = 20; // number of particles
    double softening = 0.7; // Softening parameter
    time_t start, end;
    std::vector<Particle<d>> particles; // Create a vector of particles
    Force<d>* f = new CustomForce<d>(2000); // Create force
    std::string fileName = "../../graphics/coordinates.txt"; // File name


    // Generate n random particles
    //particles = generateRandomParticlesSerial<d>(n, dim, 1, 99, 50, 50, 1, 10, false);

    //generate two particles with velocity null
    std::array<double, d> pos1 = {0, 0};
    std::array<double, d> pos2 = {20, 0};
    std::array<double, d> vel1 = {0, 0};
    std::array<double, d> vel2 = {0, 100};

    //Particle(int id, double p, std::array<double, Dimension> pos, std::array<double, Dimension> v, double radius, bool type)
    Particle<d> p1(0, 100, pos1, vel1, 1, false);
    Particle<d> p2(1, 1, pos2, vel2, 1, false);
    particles.push_back(p1);
    particles.push_back(p2);


    // Print on file the initial state of the particles
    std::ofstream file(fileName); // Open file
    if (file.is_open()) {

        // Write on file the total number of particles and the size of the area of the simulation
        file << particles.size() << std::endl;
        file << dim << std::endl;

        // Write on file the radius of the particles and the initial state
        for (const Particle<d> & p : particles) 
            file << p.getRadius() << std::endl;
        for (const Particle<d> &p : particles) {
            file << p.getId() << ",";
            const auto& pos = p.getPos();
            for (size_t i = 0; i < d; ++i) {
                file << pos[i];
                if (i < d - 1) file << ",";
            }
            file << std::endl;
        }
    }else {
        std::cout << "Unable to open file";
        return 0;
    }


    // Parallel simulation
    time(&start);
    parallelSimulation<d>(it, &particles, dim, softening, delta_t, fileName, file, *f);
    time(&end);
    printf("Time taken by the parallel implementation: %ld seconds\n", end - start);


    // Print the final state of the particles
    std::cout << "--------------------------------------------\n";
    std::cout << "Final state:\n";
    for (const Particle<d> &p : particles) {
        p.printStates();
    }
    file.close();


    return 0;
    
}
