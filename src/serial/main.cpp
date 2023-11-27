#include <vector>
#include "particle.hpp" 
#include "force.hpp"
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>
#include <fstream>

// Function that randomly generates particles:
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
                if (p.getPos()[0] == x && p.getPos()[1] == y) {
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
    const double delta_t = 0.1; // in seconds
    const double dim = 150;

    // Create a vector of particles
    std::vector<Particle> particles;
    CustomForce f;

    std::ofstream file("coordinates.txt");
    //
    //Particle p1(0, 2.0, 0.0, 0.0, 0.0, 0.0, -50.0, 2.0);
    //particles.push_back(p1);
    //Particle p2(1, 3.0, 0.0, 20.0, 0.0, 0.0, 0.0, 3.0);
    //particles.push_back(p2);

    //generate 100 particles
    int n = 100;
    particles = generateRandomParticles(n);

    // Print the initial state of the particles
    std::cout << "Initial state:\n";
    if (file.is_open()) {
        file << n << std::endl;
        for (const Particle& p : particles) {
    
        p.print_states();

        // Write on file the initial state
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
        
        p.print_states();

        // Write on file the initial state
        if (file.is_open()) {
            file << p.getId() << "," << p.getPos()[0] << "," <<  p.getPos()[1] << std::endl;
        } else {
            std::cout << "Unable to open file";
        }
    }

    // Perform the n-body simulation
    //for (Particle& q : particles) {
    //    q.resetForce();
    //}

    
    //For now it performs only two iterations, for loop with z needs to be changed 
    //for(int z=0; z<2; ++z){
//        for (int i = 0; i < particles.size(); i++) {
//            Particle &q = particles[i];
//            
//            for (int j = i + 1; j < particles.size(); j++) {
//                Particle &k = particles[j];
//
//                std::array<double,2> force_qk;
//
//                force_qk = f.calculateForce(k,q);
//
//                // Newton's third law
//                q.addForce(force_qk[0], force_qk[1]);
//                k.addForce(-force_qk[0], -force_qk[1]);
//
//            }
//            q.update(delta_t);
//            q.resetForce();
//        }
    //}
    int it = 1000;
    for (int z = 0; z < it; ++z){
        for (int i = 0; i < particles.size(); i++) {
            Particle &q = particles[i];
            if(q.getPos()[0]+ q.getRadius() > dim || 
                q.getPos()[0] - q.getRadius() < -dim ||
                q.getPos()[1] + q.getRadius()> dim || 
                q.getPos()[1] - q.getRadius()< -dim){
                std::cout << "Out of bounds for " << q.getId() << std::endl;
                q.manage_collision(q, dim);
            }
            for (int j = i + 1; j < particles.size(); j++) {
                Particle &k = particles[j];

                //check collisions
                if(q.square_distance(k) < q.getRadius() + k.getRadius()){
                    std::cout << "Collision between " << q.getId() << " and " << k.getId() << std::endl;
                    //collision method
                    q.manage_collision(k, 0.0);
                }

                std::array<double,2> force_qk;
                //force_qk = f.calculateForce(k,q);
                // Newton's third law
                q.addForce(k, f);
            }
            z==it-1? q.update(delta_t):q.update_and_reset(delta_t);
//          q.resetForce();
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

        //std::cout << "--------------------------------------------\n";
        //    for (const Particle& p : particles) {
        //        std::cout << "Position: " << p.getPos()[0] << " " << p.getPos()[1] << "\n";
        //        std::cout << "Mass: " << p.getMass() << std::endl;
        //        std::cout << "Force: " << p.getForce()[0] << " " << p.getForce()[1] << "\n";
        //        std::cout << "Velocity: " << p.getVel()[0] << " " << p.getVel()[1] << "\n";
        //    }
        //std::cout << abs(particles[0].getPos()[0] - particles[1].getPos()[0])<< std::endl;
        //std::cout << particles[0].getVel()[0] << std::endl;
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
