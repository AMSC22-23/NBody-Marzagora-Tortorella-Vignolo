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
std::vector<Particle<Dimension>> generateRandomParticlesParallel(int N = 10, int posBoundary = 100, int minMass = 1, int maxMass = 99, int maxVx = 100, int maxVy = 100, int minRadius = 0, int maxRadius = 15, bool type = false) {
    std::vector<Particle<Dimension>> particles;

    // Initialize random seed
    unsigned int seed;

    // TODO: make it a parameter that we can choose outside of this function
    int num_threads = 4, r;
    double mass, charge, x, y, vx, vy;
    bool uniquePosition = false;

    // Initialize random seeds for each thread beacuse if multiple threads call rand_r() with the same seed, it can lead to data races and unpredictable behavior.
    unsigned int seeds[num_threads];
    #pragma omp parallel
    for (int i = 0; i < num_threads; i++) {
        seeds[i] = 1234 + i;
    }

    // Parallelize the outer for loop
    #pragma omp parallel for num_threads(num_threads) private(x, y, r, mass, charge, vx, vy) shared(particles, seed)
    for (int i = 0; i < N; i++) {

        // get the thread seed
        seed = seeds[omp_get_thread_num()];
        uniquePosition = false;
        while (!uniquePosition) {
            uniquePosition = true;

            // Generate random radius between 0 and 
            r = rand() % (maxRadius - minRadius + 1) + minRadius;

            // Generate random position between -posBoundary and posBoundary
            // Use rand_r() with the seed to get thread-safe random numbers
            x = -posBoundary + r + static_cast<double>(rand_r(&seed)) / (static_cast<double>(RAND_MAX/(2*posBoundary - 2*r)));
            y = -posBoundary + r + static_cast<double>(rand_r(&seed)) / (static_cast<double>(RAND_MAX/(2*posBoundary - 2*r)));

            // Check if the position is unique
            // Use of critical section to access the shared vector of particles:
            // This is necessary because multiple threads are writing to the same memory location, which can lead to data races.
            #pragma omp critical
            {
                // Check if the position is unique and that the particles are not overlapping
                for (const Particle& p : particles) {
                    double distance = sqrt(pow(p.getPos()[0] - x, 2) + pow(p.getPos()[1] - y, 2));
                    if (distance < p.getRadius() + r) {
                        uniquePosition = false;
                        break;
                    }
                }
            }
        }

        // Generate random mass between minMass and maxMass
        mass = minMass + static_cast<double>(rand_r(&seed)) / (static_cast<double>(RAND_MAX/(maxMass-1)));
        charge = 0.0;

        // Generate random velocity between -maxVx and maxVx
        vx = -maxVx + static_cast<double>(rand_r(&seed)) / (static_cast<double>(RAND_MAX/(2*maxVx)));
        vy = -maxVy + static_cast<double>(rand_r(&seed)) / (static_cast<double>(RAND_MAX/(2*maxVy)));

        // Create a new particle with the random mass, position, and velocity
        Particle p(i, mass, charge, x, y, vx, vy, r);

        // Add the particle to the vector
        // Use an atomic operation to avoid data corruption
        #pragma omp critical
        particles.push_back(p);
    }

    return particles;
}



template<size_t Dimension>
std::vector<Particle<Dimension>> generateRandomParticlesSerial(int N, int posBoundary = 100, int minProperty = 1, int maxProperty = 99, int maxVx = 50, int maxVy = 50, int minRadius = 0, int maxRadius = 15, bool type = false) {
    std::vector<Particle<Dimension>> particles;

    // Initialize random seedgit pull
    srand(time(0));

    bool uniquePosition = false;
    double x, y;
    int r;

    for (int i = 0; i < N; i++) {
        uniquePosition = false;
        while (!uniquePosition) {
            uniquePosition = true;

            // Generate random radius between minRadius and maxRadius
            r = rand() % (maxRadius - minRadius + 1) + minRadius;

            // Generate random position between -posBoundary+r and +posBoundary-r
            x = -posBoundary + r + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*posBoundary - 2*r)));
            y = -posBoundary + r + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*posBoundary - 2*r)));

            // Check if the position is unique and that the particles are not overlapping
            for (const Particle& p : particles) {
                double distance = sqrt(pow(p.getPos()[0] - x, 2) + pow(p.getPos()[1] - y, 2));
                if (distance < p.getRadius() + r) {
                    uniquePosition = false;
                    break;
                }
            }
        }

        // Generate random value of property between 1 and 100
        double property = minProperty + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(maxProperty-1)));

        // Generate random velocity between -100 and 100
        double vx = -maxVx + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVx)));
        double vy = -maxVy + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVy)));

        // Create a new particle with the random value of property, position, and velocity
        Particle p(i, property, x, y, vx, vy, r, type);

        // Add the particle to the vector
        particles.push_back(p);
    }

    return particles;
}


int main() {
    // Define of variables
    const double delta_t = 0.01; // in seconds
    const double dim = 500; // Dimension of the simulation area

    // number of iteration
    int it = 1000;
    
    //2D or 3D
    const int d = 2;

    // Create a vector of particles
    std::vector<Particle<d>> particles;
    Force<d>* f = new CustomForce<d>();
    std::string fileName = "../graphics/coordinates.txt";
    std::ofstream file(fileName); // Open file
    
    //generate 100 particles
    int n = 50;
    particles = generateRandomParticlesParallel<d>(n, dim, 1, 99, 50, 50, 1, 10, false);

    double softening = 0.7; // Softening parameter

    // Print the initial state of the particles
    std::cout << "Initial state:\n";
    if (file.is_open()) {

        // Write on file the total number of particles
        file << n << std::endl;

        // Write on file the size of the area of the simulation
        file << dim << std::endl;

        // Write on file the radius of the particles and the initial state
        for (const Particle<d> & p : particles) 
            file << p.getRadius() << std::endl;
        for (const Particle<d> &p : particles) {
            //file << p.getId() << "," << p.getPos()[0] << "," <<  p.getPos()[1] << std::endl;
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
    }

    // Start of the simulation
    for (int z = 0; z < it; ++z){
        for (int i = 0; i < particles.size(); ++i) {
            Particle<d> &q = particles[i];

            // Check if the particle hits the bounday

            if(q.hitsBoundary(dim)){
                q.manageCollision(q, dim);
            }
            for (int j = i + 1; j < particles.size(); j++) {
                Particle<d> &k = particles[j];

                // check collisions between particles
                if(q.squareDistance(k) < (((q.getRadius() + k.getRadius())*(q.getRadius() + k.getRadius())))*softening){

                    // Call the collision method
                    q.manageCollision(k, 0.0);
                }

                std::array<double, d> force_qk;
                //force_qk = f.calculateForce(k,q);
                q.addForce(k, *f);
            }
            z==it-1? q.update(delta_t):q.updateAndReset(delta_t);
        }
        
        // Write on file the updates after delta_t
        for (int i = 0; i < particles.size(); i++) {
            Particle<d> &q = particles[i];
            // Write on file the updates after delta_t
            if (file.is_open()) {

                //file << q.getId() << "," << q.getPos()[0] << "," <<  q.getPos()[1] << std::endl;
                file << q.getId() << ",";

                const auto& pos = q.getPos();
                for (size_t i = 0; i < d; ++i) {
                    file << pos[i];
                    if (i < d - 1) {
                        file << ",";
                    }
                }

                file << std::endl;

            } else {
                std::cout << "Unable to open file";
            }
        }
    }

    // Print the final state of the particles
    std::cout << "--------------------------------------------\n";
    std::cout << "Final state:\n";
    for (const Particle<d> &p : particles) {
        p.printStates();

        double power = 0.0;
        for (size_t i = 0; i < d; ++i) power = power + (p.getPos()[i] * p.getPos()[i]);
        std::cout << "Distance from origin: " << sqrt(power) << "\n";
    }
    file.close();
    return 0;
}