#include <vector>
#include "particle.hpp" 
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>

std::vector<Particle> generateRandomParticles(int N) {
    std::vector<Particle> particles;

    // Initialize random seed
    srand(time(0));

    for (int i = 0; i < N; i++) {
        bool uniquePosition = false;
        double x, y;

        while (!uniquePosition) {
            uniquePosition = true;

            // Generate random position between -100 and 100
            x = -100 + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/200));
            y = -100 + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/200));

            // Check if the position is unique
            for (const Particle& p : particles) {
                if (p.getPos()[0] == x && p.getPos()[1] == y) {
                    uniquePosition = false;
                    break;
                }
            }
        }

        // Generate random mass between 1 and 100
        double mass = 1 + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/99));

        // Generate random velocity between -100 and 100
        double vx = -100 + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/200));
        double vy = -100 + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/200));

        // Create a new particle with the random mass, position, and velocity
        Particle p(mass, x, y, vx, vy);

        // Add the particle to the vector
        particles.push_back(p);
    }

    return particles;
}

int main() {
        std::cout << "Start" << "\n";
    // Define the gravitational constant and time step
    const double G = 6.67430e-11; // in m^3 kg^-1 s^-2
    const double delta_t = 0.01; // in seconds

    // Create a vector of particles
    std::vector<Particle> particles;

    particles = generateRandomParticles(2);

    // Print the initial state of the particles
    std::cout << "Initial state:\n";
    for (const Particle& p : particles) {
        std::cout << p.getPos() << "\n";
    }

    // Perform the n-body simulation
    for (Particle& q : particles) {
        q.resetForce();
    }

    for (int i = 0; i < particles.size(); i++) {
        Particle& q = particles[i];
        for (int j = i + 1; j < particles.size(); j++) {
            Particle& k = particles[j];

            double x_diff = q.getPos()[0] - k.getPos()[0];
            double y_diff = q.getPos()[1] - k.getPos()[1];
            double dist = sqrt(x_diff * x_diff + y_diff * y_diff);
            double dist_cubed = dist * dist * dist;

            double force_qk[2];
            force_qk[0] = G * q.getMass() * k.getMass() / dist_cubed * x_diff;
            force_qk[1] = G * q.getMass() * k.getMass() / dist_cubed * y_diff;

            // Newton's third law
            q.addForce(force_qk[0], force_qk[1]);
            k.addForce(-force_qk[0], -force_qk[1]);

            // Update position and velocity
            q.update(delta_t);
            k.update(delta_t);
        }
    }

    // Print the final state of the particles
    std::cout << "Final state:\n";
    for (const Particle& p : particles) {
        std::cout << p.getPos() << "\n";
    }

    return 0;
}
