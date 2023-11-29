#include <vector>
#include "../utils/particle.hpp" 
#include "../utils/force.hpp"
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>
#include <fstream>

// TODO: sistemare generazione quando stalla
std::vector<Particle> generateRandomParticles(int N, int minMass = 1, int maxMass = 99, int posBoundary = 100, int maxVx = 100, int maxVy = 100, int minRadius = 0, int maxRadius = 15) {
    std::vector<Particle> particles;

    // Initialize random seed
    srand(time(0));

    for (int i = 0; i < N; i++) {

        bool uniquePosition = false;
        double x, y;

        while (!uniquePosition) {
            uniquePosition = true;

            // Generate random position between -100 and 100
            x = -posBoundary + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*posBoundary)));
            y = -posBoundary + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*posBoundary)));

            // Check if the position is unique
            for (const Particle& p : particles) {
                double distance = sqrt(pow(p.getPos()[0] - x, 2) + pow(p.getPos()[1] - y, 2));
                if (distance < p.getRadius() + maxRadius) {
                    uniquePosition = false;
                    break;
                }
            }
        }

        // Generate random radius between 0 and 
        int r = rand() % (maxRadius - minRadius + 1) + minRadius;

        // Generate random mass between 1 and 100
        double mass = minMass + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(maxMass-1)));
        double charge = 0.0;

        // Generate random radius between 1 and 10
        double radius = 1 + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(10-1)));

        // Generate random velocity between -100 and 100
        double vx = -maxVx + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVx)));
        double vy = -maxVy + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVy)));

        // Create a new particle with the random mass, position, and velocity
        Particle p(i, mass, charge, x, y, vx, vy, r);

        // Add the particle to the vector
        particles.push_back(p);
    }

    return particles;
}


int main() {
    // Define the gravitational constant and time step
    const double delta_t = 0.01; // in seconds
    const double dim = 250;

    // number of iteration
    int it = 1000; 

    // Create a vector of particles
    std::vector<Particle> particles;
    CustomForce f;

    std::ofstream file("../graphics/coordinates.txt");

    //generate 100 particles
    int n = 50;
    particles = generateRandomParticles(n);

    // Print the initial state of the particles
    std::cout << "Initial state:\n";
    if (file.is_open()) {
        file << n << std::endl;
        for (const Particle& p : particles) {
        p.print_states();
        // Write on file the radius of the particles
        if (file.is_open()) {
            file << p.getRadius() << std::endl;
        } else {
            std::cout << "Unable to open file";
        }
    }
    }else {
        std::cout << "Unable to open file";
    }

    for (const Particle& p : particles) {
        // Write on file the initial state
        if (file.is_open()) {
            file << p.getId() << "," << p.getPos()[0] << "," <<  p.getPos()[1] << std::endl;
        } else {
            std::cout << "Unable to open file";
        }
    }

    for (int z = 0; z < it; ++z){
        for (int i = 0; i < particles.size(); i++) {
            Particle &q = particles[i];
            // check if the particle hits the bounday
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
                    // call the collision method
                    q.manage_collision(k, 0.0);
                }

                std::array<double,2> force_qk;
                //force_qk = f.calculateForce(k,q);
                q.addForce(k, f);
            }
            z==it-1? q.update(delta_t):q.update_and_reset(delta_t);
        }
        
        for (int i = 0; i < particles.size(); i++) {
            Particle &q = particles[i];
            // Write on file the updates after delta_t
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
