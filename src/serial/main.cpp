#include <vector>
// include should have only the name
// add "-I../utils" to the compile flags 
#include "../utils/particle.hpp" 
#include "../utils/force.hpp" 
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

// TODO: sistemare generazione quando stalla
// `Particle` is a template class, it does not name a "concrete" type
// until you specify its dimension, so you have two options:
// 1. define the function as `std::vector<Particle<2>> generateRandomParticles(...)`
// or better
// 2. make the function a template, as follows
template<size_t Dimension>
std::vector<Particle<Dimension>> generateRandomParticles(int N, int posBoundary = 100, int minProperty = 1, int maxProperty = 99, int maxVx = 100, int maxVy = 100, int minRadius = 0, int maxRadius = 15, bool type = false) {
    std::vector<Particle<Dimension>> particles;

    // Initialize random seed
    srand(time(0));

    bool uniquePosition = false;
    double x, y;
    int r;
    double power = 0.0;

    for (int i = 0; i < N; i++) {

        bool uniquePosition = false;
        //double x, y;
        std::array<double, Dimension> pos;

        while (!uniquePosition) {
            uniquePosition = true;
            
            // Generate random radius between minRadius and maxRadius
            r = rand() % (maxRadius - minRadius + 1) + minRadius;

            // Generate random position between -posBoundary+r and +posBoundary-r
            for(size_t i=0; i<Dimension; ++i)
                pos[i] = -posBoundary + r + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*posBoundary - 2*r)));

            // Check if the position is unique
            for (const Particle<Dimension> &p : particles) {

                for(size_t i = 0; i < Dimension; ++i) 
                    power = power + ((p.getPos()[i] - pos[i])*(p.getPos()[i] - pos[i]));

                double distance = sqrt(power);
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
            vel[i] = -maxVy + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX/(2*maxVy)));

        // Create a new particle with the random value of property, position, and velocity
        Particle<Dimension> p(i, property, pos, vel, r, type);

        // Add the particle to the vector
        particles.push_back(p);
    }

    return particles;
}

template<size_t Dimension>
void serialSimulation(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f){
    std::array<double,Dimension> force_qk;
    
    // Start of the simulation
    for (int z = 0; z < it; ++z){
        for (int i = 0; i < (*particles).size(); ++i) {
            Particle<Dimension> &q = (*particles)[i];

            // Check if the particle hits the bounday
            if(q.hitsBoundary(dim)){
                q.manageCollision(q, dim);
            }
            for (int j = i + 1; j < (*particles).size(); j++) {
                Particle<Dimension> &k = (*particles)[j];

                // check collisions between particles
                if(q.squareDistance(k) < (((q.getRadius() + k.getRadius())*(q.getRadius() + k.getRadius())))*softening){

                    // Call the collision method
                    q.manageCollision(k, 0.0);
                }

                force_qk = f.calculateForce(q, k);
                q.addForce(force_qk);
                for(size_t i = 0; i < Dimension; ++i) force_qk[i] = -force_qk[i];
                k.addForce(force_qk);

            }
            z==it-1? q.update(delta_t):q.updateAndReset(delta_t);
        }
        
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
                std::cout << "Unable to open file";
            }
        }
    }
}

int main() {
    // Define of simulation variables
    const int d = 2; //2D or 3D
    const double delta_t = 0.01; // in seconds
    const double dim = 50; // Dimension of the simulation area
    int it = 1000; // number of iteration
    int n = 2; // number of particles
    double softening = 0.7; // Softening parameter
    time_t start, end; // Time variables
    std::vector<Particle<d>> particles; // Create a vector of particles
    Force<d>* f = new CustomForce<d>(2000); // Create force
    std::string fileName = "../graphics/coordinates.txt"; // File name
    std::ofstream file(fileName); // Open file

    // Generate n random particles
    //particles = generateRandomParticles<d>(n, dim, 1, 99, 50, 50, 1, 10, false);

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
    if (file.is_open()) {

        // Write on file the total number of particles and the size of the area of the simulation
        file << n << std::endl;
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

    // Start of simulation and compute the time taken
    start = time(NULL);
    serialSimulation<d>(it, &particles, dim, softening, delta_t, fileName, file, *f);
    end = time(NULL);
    std::cout << "Time taken by serial simulation: " << end - start << " seconds" << std::endl;

    // Print the final state of the particles
    //std::cout << "--------------------------------------------\n";
    //std::cout << "Final state:\n";
    //for (const Particle<d> &p : particles) {
    //    p.printStates();
    //    double power = 0.0;
    //    for (size_t i = 0; i < d; ++i) 
    //        power = power + (p.getPos()[i] * p.getPos()[i]);
    //    std::cout << "Distance from origin: " << sqrt(power) << "\n";
    //}
    file.close();
    return 0;
}