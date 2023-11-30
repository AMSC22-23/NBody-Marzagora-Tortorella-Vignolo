#include <vector>
#include "../utils/particle.hpp" 
#include "../utils/force.hpp"
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

// Function that randomly generates particles:
std::vector<Particle> generateRandomParticles(int N = 10, int posBoundary = 100, int minMass = 1, int maxMass = 99, int maxVx = 100, int maxVy = 100, int minRadius = 0, int maxRadius = 15) {
    std::vector<Particle> particles;

    // Initialize random seed
    unsigned int seed = 1234;

    // Define num threads
    int num_threads = 4;

    bool uniquePosition = false;
    double x, y;
    int r;

    // Parallelize the outer for loop with OpenMP
    #pragma omp parallel for num_threads(num_threads) private(x, y, r, mass, charge, radius, vx, vy, p) shared(particles, seed)
    for (int i = 0; i < N; i++) {
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
            // Use of critical section to access the shared vector of particles
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
        double mass = minMass + static_cast<double>(rand_r(&seed)) / (static_cast<double>(RAND_MAX/(maxMass-1)));
        double charge = 0.0;

        // Generate random radius between minRadius and maxRadius
        double radius = minRadius + static_cast<double>(rand_r(&seed)) / (static_cast<double>(RAND_MAX/(maxRadius-minRadius)));

        // Generate random velocity between -maxVx and maxVx
        double vx = -maxVx + static_cast<double>(rand_r(&seed)) / (static_cast<double>(RAND_MAX/(2*maxVx)));
        double vy = -maxVy + static_cast<double>(rand_r(&seed)) / (static_cast<double>(RAND_MAX/(2*maxVy)));

        // Create a new particle with the random mass, position, and velocity
        Particle p(i, mass, charge, x, y, vx, vy, r);

        // Add the particle to the vector
        // Use an atomic operation to avoid data corruption
        #pragma omp atomic
        particles.push_back(p);
    }


    return particles;
}


int main() {

    // Definition of variables
    const double delta_t = 0.01; // [sec]
    int dim = 250; // Dimension of the simulation area
    int it = 1000; // Number of iteration
    int n = 50; // Number of particles
    std::vector<Particle> particles;  // Create a vector of particles
    CustomForce f;
    std::string fileName = "../graphics/coordinates.txt";
    std::ofstream file(fileName); // Open file
    
    // Generation of n random particles
    particles  = generateRandomParticles(n, dim);

    // Print the initial state of the particles
    std::cout << "Initial state:\n";
    if (file.is_open()) {

        // Write on file the total number of particles
        file << n << std::endl;

        // Write on file the size of the area of the simulation
        file << dim << std::endl;

        // Write on file the radius of the particles and the initial state
        for (const Particle& p : particles) 
            file << p.getRadius() << std::endl;
        for (const Particle& p : particles) {
            file << p.getId() << "," << p.getPos()[0] << "," <<  p.getPos()[1] << std::endl;
    }
    }else {
        std::cout << "Unable to open file";
    }

    // Start of the simulation
    for (int z = 0; z < it; ++z){
        for (int i = 0; i < particles.size(); ++i) {
            Particle &q = particles[i];

            // Check if the particle hits the bounday
            if(q.getPos()[0]+ q.getRadius() > dim || 
                q.getPos()[0] - q.getRadius() < -dim ||
                q.getPos()[1] + q.getRadius()> dim || 
                q.getPos()[1] - q.getRadius()< -dim){
                q.manage_collision(q, dim);
            }
            for (int j = i + 1; j < particles.size(); j++) {
                Particle &k = particles[j];

                // check collisions between particles
                if(q.square_distance(k) < ((q.getRadius() + k.getRadius())*(q.getRadius() + k.getRadius()))){

                    // Call the collision method
                    q.manage_collision(k, 0.0);
                }
                // Add force
                q.addForce(k, f);
            }
            z==it-1? q.update(delta_t):q.update_and_reset(delta_t);
        }
        
        // Write on file the updates after delta_t
        for (int i = 0; i < particles.size(); i++) {
            Particle &q = particles[i];
            if (file.is_open()) {
                file << q.getId() << "," << q.getPos()[0] << "," <<  q.getPos()[1] << std::endl;
            } else {
                std::cout << "Unable to open file";
            }
        }
    }

    // Print the final state of the particles
    std::cout << "--------------------------------------------\n";
    std::cout << "Final state:\n";
    for (const Particle& p : particles) {
        p.print_states();

        std::cout << "Distance from origin: " << sqrt(pow(p.getPos()[0],2) + pow(p.getPos()[1],2)) << "\n";
    }
    file.close();
    return 0;
}