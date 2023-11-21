#include <vector>
#include "particle.hpp" 
#include "force.hpp"
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>

// Function that randomly generates particles:
std::vector<Particle> generateRandomParticles(int N, int minMass = 1, int maxMass = 99, int posBoundary = 100, int maxVx = 100, int maxVy = 100) {
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
                if (p.getPos()[0] == x && p.getPos()[1] == y) {
                    uniquePosition = false;
                    break;
                }
            }
        }

        // Generate random mass between 1 and 100
        double mass = minMass + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(maxMass-1)));
        double charge = 0.0;

        // Generate random velocity between -100 and 100
        double vx = -maxVx + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVx)));
        double vy = -maxVy + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVy)));

        // Create a new particle with the random mass, position, and velocity
        Particle p(mass, charge, x, y, vx, vy);

        // Add the particle to the vector
        particles.push_back(p);
    }

    return particles;
}

int main() {
    // Define the gravitational constant and time step
    const double delta_t = 0.001; // in seconds

    // Create a vector of particles
    std::vector<Particle> particles;
    GravitationalForce f;
    
    Particle p1(20.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    particles.push_back(p1);
    Particle p2(70.0, 0.0, 20.0, 0.0, 0.0, 0.0);
    particles.push_back(p2);

    // Print the initial state of the particles
    std::cout << "Initial state:\n";
    for (const Particle& p : particles) {
        //std::cout << "Position: " << p.getPos()[0] << " " << p.getPos()[1] << "\n";
        //std::cout << "Mass: " << p.getMass() << std::endl;
        //std::cout << "Force: " << p.getForce()[0] << " " << p.getForce()[1] << "\n";
        //std::cout << "Velocity: " << p.getVel()[0] << " " << p.getVel()[1] << "\n";
        p.print_states();
    }

    // Perform the n-body simulation
    //for (Particle& q : particles) {
    //    q.resetForce();
    //}

    
    //For now it performs only two iterations, for loop with z needs to be changed 
    //for(int z=0; z<2; ++z){
        for (int i = 0; i < particles.size(); i++) {
            Particle &q = particles[i];
            
            q.resetForce();
            for (int j = i + 1; j < particles.size(); j++) {
                Particle &k = particles[j];

                std::array<double,2> force_qk;

                force_qk = f.calculateForce(k,q);

                // Newton's third law
                q.addForce(force_qk[0], force_qk[1]);
                k.addForce(-force_qk[0], -force_qk[1]);

                // Update position and velocity
                q.update(delta_t);
                k.update(delta_t);


            }
        }
        //std::cout << "--------------------------------------------\n";
        //    for (const Particle& p : particles) {
        //        std::cout << "Position: " << p.getPos()[0] << " " << p.getPos()[1] << "\n";
        //        std::cout << "Mass: " << p.getMass() << std::endl;
        //        std::cout << "Force: " << p.getForce()[0] << " " << p.getForce()[1] << "\n";
        //        std::cout << "Velocity: " << p.getVel()[0] << " " << p.getVel()[1] << "\n";
        //    }
        //std::cout << abs(particles[0].getPos()[0] - particles[1].getPos()[0])<< std::endl;
        //std::cout << particles[0].getVel()[0] << std::endl;
    //}

    // Print the final state of the particles
    std::cout << "--------------------------------------------\n";
    std::cout << "Final state:\n";
    for (const Particle& p : particles) {
        p.print_states();
    }

    return 0;
}
