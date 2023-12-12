#include <vector>
#include "particle.hpp" 
#include "force.hpp" 
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>


template<size_t Dimension>
std::vector<Particle<Dimension>> generateRandomParticles(int N, int posBoundary = 100, int minProperty = 1, int maxProperty = 99, int maxVel = 100, int minRadius = 0, int maxRadius = 15, bool type = false) {

    // Initialize variables
    int maxRetry  = 15, counter = 0;
    bool overlapping = false;
    double x, y, r, property, squareDistance = 0.0;
    std::vector<Particle<Dimension>> particles;
    std::array<double, Dimension> pos, vel;

    // Create a random number generators using 'random' library
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> disProperty(minProperty, maxProperty);
    std::uniform_real_distribution<double> disVel(-maxVel, maxVel);
    std::uniform_real_distribution<double> disRadius(minRadius, maxRadius);
    std::uniform_real_distribution<double> disPos(-posBoundary+maxRadius, posBoundary-maxRadius);

    while (particles.size() < N) {

        // Check if the function struggles to generate to generate non-overlapping particles: 
        // if the function has to regenerate the particle's position more than maxRetry times in a row, then the program ends and print an error message
        if(counter == maxRetry)
            throw std::runtime_error("ERROR: the dimension of the simulation area is too little, please specify a bigger area or generate less particles.");

        overlapping = false;
                
        // Generate random radius between minRadius and maxRadius
        r = disRadius(gen);

        // Generate random position between -posBoundary+r and +posBoundary-r
        for(size_t i=0; i<Dimension; ++i)
            pos[i] = disPos(gen);

        // Check that it is not overlapping with any existing circle
        // Brute force approach
        for (const Particle<Dimension> &p : particles) {

            // Calculate the distance between the particle that we would like to create and the particle p in particles vector
            squareDistance = 0.0;
            for(size_t i = 0; i < Dimension; ++i) 
                squareDistance = squareDistance + ((p.getPos()[i] - pos[i])*(p.getPos()[i] - pos[i]));
            if (sqrt(squareDistance) < p.getRadius() + r) { // Particles are overlapped
                counter++;
                overlapping = true;
                break; // As soon as an overlap is detected, the break statement is executed. This immediately exits the loop, skipping the remaining particles in the vector.
            }
        }

        // Add valid circles to array
        if (!overlapping) {
            counter = 0;
            // Generate random value of property between minProperty and maxProperty
            property = disProperty(gen);

            // Generate random velocity between -maxVel and maxVel
            for(size_t i=0; i<Dimension; ++i)
                vel[i] = disVel(gen);

            // Create a new particle with the random value of property, position, and velocity
            Particle<Dimension> p(particles.size(), property, pos, vel, r, type);

            // Add the particle to the vector
            particles.push_back(p);
        }else{
            counter++;
        }
        
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

                // Check collisions between particles
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



template<size_t Dimension>
void parallelSimulation(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f, int u){

    particles->size() < omp_get_max_threads()? omp_set_num_threads(particles->size()):omp_set_num_threads(omp_get_max_threads());  

    // Start of the simulation
    std::vector<std::vector<std::array<double, Dimension>>> local_forces(omp_get_max_threads(), std::vector<std::array<double, Dimension>>(particles->size()));



    #pragma omp parallel shared(particles, f, it, delta_t, fileName, file, softening, dim)
    {
        //#pragma omp single
        //{
        //    std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
        //    std::cout << "Number of max threads: " << omp_get_max_threads() << std::endl;
        //}
 
        //create a vector of thread that contains a vector of forces of dimension Dimension
        std::array<double, Dimension> force;

        for (int z = 0; z < it; ++z){

            //calculate forces and add it to a local array for each thread
            //assign block thread in a cyclic way, not static
            #pragma omp for schedule(dynamic, u)
            for (int i = 0; i < (*particles).size(); ++i) {
                for(int y = 0; y < Dimension; ++y) force[y] = 0;
                
                Particle<Dimension> &k = (*particles)[i];

                // Check if the particle hits the boundary
                if(k.hitsBoundary(dim)){
                    k.manageCollision(k, dim);
                }
            
                for(int j = i+1; j < (*particles).size(); ++j){
                    Particle<Dimension> &q = (*particles)[j];

                    if(q.squareDistance(k) < (((q.getRadius() + k.getRadius())*(q.getRadius() + k.getRadius())))*softening){
                        q.manageCollision(k, 0.0);
                    }

                    force = f.calculateForce(k, q);
                    for (int y = 0; y < Dimension; ++y) {
                        local_forces[omp_get_thread_num()][i][y] += force[y];
                        local_forces[omp_get_thread_num()][j][y] -= force[y];
                    }
                }
            }
            //#pragma omp barrier

            //sum all the forces calculated by each thread
            #pragma omp for schedule(static, particles->size()/omp_get_num_threads())
            for (int i = 0; i < (*particles).size(); ++i) {
                Particle<Dimension> &q = (*particles)[i];
                q.resetForce();
                for (int j = 0; j < omp_get_num_threads(); ++j) {
                    q.addForce(local_forces[j][i]);
                }
            }
            //#pragma omp barrier

            #pragma omp for schedule(static, 1)
            for(int k = 0; k < omp_get_num_threads(); ++k){
                for (int j = 0; j < (*particles).size(); ++j) {
                    for (int y = 0; y < Dimension; ++y) {
                        local_forces[k][j][y] = 0.0;
                    }
                }
            }
            //#pragma omp barrier

            //update the position of the particles
            #pragma omp for schedule(static, particles->size()/omp_get_num_threads())
            for (int i = 0; i < (*particles).size(); ++i) {
                Particle<Dimension> &q = (*particles)[i];
                q.update(delta_t);
            }
            //#pragma omp barrier

            // Write on file the updates after delta_t
            #pragma omp single
            {
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
            //#pragma omp barrier
        }
    }
}



template<size_t Dimension>
void printAllParticlesStateAndDistance(std::vector<Particle<Dimension>>* particles){
    std::cout << "Print All Particles State And Distance:\n";
    double squareDistance = 0.0;
    for (int i = 0; i < (*particles).size(); ++i) {
        Particle<Dimension> &p = (*particles)[i];
        std::cout << "--------------------------------------\n";
        p.printStates();
        squareDistance = 0.0;
        for (size_t i = 0; i < Dimension; ++i) 
            squareDistance = squareDistance + (p.getPos()[i] * p.getPos()[i]);
        std::cout << "Distance from origin: " << sqrt(squareDistance) << "\n";
    }
}



template<size_t Dimension>
std::vector<Particle<Dimension>> generateOrbitTestParticles( double size, double costantForce){
    std::vector<Particle<Dimension>> particles; 
    double orbitRadius = size/2.0;
    double mass1 = 100;
    double mass2 = 1;
    double force = costantForce*(mass1*mass2)/(orbitRadius*orbitRadius);
    double velocity = sqrt(force*orbitRadius/mass2);
    std::array<double, Dimension> pos1;
    std::array<double, Dimension> pos2;
    std::array<double, Dimension> vel1;
    std::array<double, Dimension> vel2;

    for (int i=0; i<Dimension; ++i)
    {
        pos1[i] = 0;
        pos2[i] = i==0? orbitRadius:0;
        vel1[i] = 0;
        vel2[i] = i==(Dimension-1)? velocity:0;
    } 

    // Particle(int id, double p, std::array<double, Dimension> pos, std::array<double, Dimension> v, double radius, bool type)
    Particle<Dimension> p1(0, mass1, pos1, vel1, size/100, false);
    Particle<Dimension> p2(1, mass2, pos2, vel2, size/100, false);
    particles.push_back(p1);
    particles.push_back(p2);

    return particles;
}

template<size_t Dimension>
void printInitialStateOnFile(std::vector<Particle<Dimension>>* particles, int dim, std::string fileName, std::ofstream& file, int it){
    
    // Print on file the initial state of the particles
    if (file.is_open()) {

        // Write on file the total number of particles and the size of the area of the simulation
        file << (*particles).size() << std::endl;
        file << dim << std::endl;
        file << Dimension << std::endl;
        file << it << std::endl;

        // Write on file the radius of the particles and the initial state
        for (const Particle<Dimension> & p : (*particles)) 
            file << p.getRadius() << std::endl;
        for (const Particle<Dimension> &p : (*particles)) {
            file << p.getId() << ",";
            const auto& pos = p.getPos();
            for (size_t i = 0; i < Dimension; ++i) {
                file << pos[i];
                if (i < Dimension - 1) file << ",";
            }
            file << std::endl;
        }
    } else {
        std::cout << "Unable to open file";
    }
}


#ifndef SIMULATION_TYPE
    #define SIMULATION_TYPE 1 // Default value
#endif

#ifdef DIMENSION
    const size_t d = DIMENSION;
#else
    const size_t d = 2; // Default value
#endif

int main(int argc, char** argv) {

    // Simulation variables    
    const double delta_t = 0.01; // In seconds
    int dim = 200; // Dimension of the simulation area
    int it = 100; // Number of iteration
    int n = 20; // Number of particles
    int mass = 50; // Mass
    int maxVel = 50; // Maximum velocity
    int maxRadius = 50; // Maximum radius of the particles
    double softening = 0.7; // Softening parameter
    time_t start, end; // Time variables
    std::vector<Particle<d>> particles; // Create a vector of particles
    Force<d>* f = new CustomForce<d>(2000); // Create force
    std::string fileName = "../graphics/coordinates.txt"; // File name
    std::ofstream file(fileName); // Open file

    // Generate random particles
    start = time(NULL);
    particles = generateRandomParticles<d>(n, dim, 1, mass, maxVel, 1, maxRadius, false);
    end = time(NULL);
    std::cout << "Time taken by generateRandomParticles function: " << end - start << " seconds" << std::endl;
    

    // Print on file the initial state of the particles
    printInitialStateOnFile(&particles, dim, fileName, file, it);

    // Start simulation
    if(SIMULATION_TYPE == 0){
        start = time(NULL);
        serialSimulation<d>(it, &particles, dim, softening, delta_t, fileName, file, *f);
        end = time(NULL);
        std::cout << "Time taken by serial simulation: " << end - start << " seconds" << std::endl;

    }else if(SIMULATION_TYPE == 1){
        start = time(NULL);
        parallelSimulation<d>(it, &particles, dim, softening, delta_t, fileName, file, *f, 1);
        end = time(NULL);
        std::cout << "Time taken by parallel simulation: " << end - start << " seconds" << std::endl;
    } else {
        std::cout << "You did not specify correctly the type of the simulation [serial/parallel]" << std::endl;
    }

    // In order to print the final state of the particles use: printAllParticlesStateAndDistance(&particles);

    file.close();
    return 0;
    
}