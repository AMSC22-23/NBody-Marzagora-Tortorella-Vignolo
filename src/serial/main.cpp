#include <vector>
#include "../utils/particle.hpp" 
#include "../utils/force.hpp"
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

// TODO: sistemare generazione quando stalla
std::vector<Particle> generateRandomParticles(int N, int posBoundary = 100, int minProperty = 1, int maxProperty = 99, int maxVx = 50, int maxVy = 50, int minRadius = 0, int maxRadius = 15, bool type = false) {
    std::vector<Particle> particles;

    // Initialize random seed
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

    // Definition of variables
    const double delta_t = 0.01; // [sec]
    int dim = 600; // Dimension of the simulation area
    int it = 1000; // Number of iteration
    int n = 50; // Number of particles
    double softening = 0.7; // Softening parameter
    std::vector<Particle> particles;  // Create a vector of particles
    CustomForce f;
    std::string fileName = "../graphics/coordinates.txt";
    std::ofstream file(fileName); // Open file
    
    // Generation of n random particles
    //particles  = generateRandomParticles(n, dim);

    //call the generateRandomParticles function to generate the particles and put minRadius = 1
    particles  = generateRandomParticles(n, dim, 1, 99, 50, 50, 1, 10, false);

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
            if(q.getPos()[0] + q.getRadius() > dim || 
                q.getPos()[0] - q.getRadius() < -dim ||
                q.getPos()[1] + q.getRadius()> dim || 
                q.getPos()[1] - q.getRadius()< -dim){
                q.manage_collision(q, dim);
            }
            for (int j = i + 1; j < particles.size(); j++) {
                Particle &k = particles[j];

                // check collisions between particles
                if(q.square_distance(k) < (((q.getRadius() + k.getRadius())*(q.getRadius() + k.getRadius())))*softening){

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
