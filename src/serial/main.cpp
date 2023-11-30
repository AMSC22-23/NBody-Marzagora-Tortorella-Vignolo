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
template<size_t Dimension>
std::vector<Particle<Dimension>> generateRandomParticles(int N, int posBoundary = 100, int minMass = 1, int maxMass = 99, int maxVx = 100, int maxVy = 100, int minRadius = 0, int maxRadius = 15) {
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

            
            // Generate random radius between 0 and 
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

        // Generate random mass between 1 and 100
        double mass = minMass + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(maxMass-1)));
        double charge = 0.0;

        // Generate random radius between 1 and 10
        double radius = 1 + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(10-1)));


        std::array<double, Dimension> vel;
        // Generate random velocity between -100 and 100
        for(size_t i=0; i<Dimension; ++i)
        //double vx = -maxVx + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVx)));
        //double vy = -maxVy + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVy)));
            vel[i] = -maxVy + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVy)));

        // Create a new particle with the random mass, position, and velocity
        //Particle p(i, mass, charge, x, y, vx, vy, r);
        Particle p(i, mass, charge, pos, vel, r);

        // Add the particle to the vector
        particles.push_back(p);
    }

    return particles;
}


int main() {
    // Define of variables
    const double delta_t = 0.01; // in seconds
    const double dim = 250; // Dimension of the simulation area

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
    particles = generateRandomParticles<d>(n,dim);

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
            file << p.getId() << "," << p.getPos()[0] << "," <<  p.getPos()[1] << std::endl;
    }
    }else {
        std::cout << "Unable to open file";
    }

    // Start of the simulation
    for (int z = 0; z < it; ++z){
        for (int i = 0; i < particles.size(); ++i) {
            Particle<d> &q = particles[i];

            // Check if the particle hits the bounday
            if(q.getPos()[0]+ q.getRadius() > dim || 
                q.getPos()[0] - q.getRadius() < -dim ||
                q.getPos()[1] + q.getRadius()> dim || 
                q.getPos()[1] - q.getRadius()< -dim){
                q.manage_collision(q, dim);
            }
            for (int j = i + 1; j < particles.size(); j++) {
                Particle<d> &k = particles[j];

                // check collisions between particles
                if(q.square_distance(k) < ((q.getRadius() + k.getRadius())*(q.getRadius() + k.getRadius()))){

                    // Call the collision method
                    q.manage_collision(k, 0.0);
                }

                std::array<double, d> force_qk;
                //force_qk = f.calculateForce(k,q);
                q.addForce(k, *f);
            }
            z==it-1? q.update(delta_t):q.update_and_reset(delta_t);
        }
        
        // Write on file the updates after delta_t
        for (int i = 0; i < particles.size(); i++) {
            Particle<d> &q = particles[i];
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
    for (const Particle<d> &p : particles) {
        p.print_states();

        std::cout << "Distance from origin: " << sqrt(pow(p.getPos()[0],2) + pow(p.getPos()[1],2)) << "\n";
    }
    file.close();
    return 0;
}