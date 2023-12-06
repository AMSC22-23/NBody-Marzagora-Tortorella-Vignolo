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
std::vector<Particle<Dimension>> generateRandomParticlesSerial(int N, int posBoundary = 100, int minProperty = 1, int maxProperty = 99, int maxVx = 100, int maxVy = 100, int minRadius = 0, int maxRadius = 2.5, bool type = false) {
    std::vector<Particle<Dimension>> particles;

    // Initialize random seed
    srand(time(0));

    bool uniquePosition = false;
    int r;
    double power = 0.0;

    for (int i = 0; i < N; i++) {

        bool uniquePosition = false;
        std::array<double, Dimension> pos;

        while (!uniquePosition) {
            uniquePosition = true;
            
            // Generate random radius between minRadius and maxRadius
            //r = rand() % (maxRadius - minRadius + 1) + minRadius;
            r = 1;
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
void parallelSimulation(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f){
    
    // Start of the simulation
    std::vector<std::vector<std::array<double, Dimension>>> local_forces(omp_get_max_threads(), std::vector<std::array<double, Dimension>>(particles->size()));

    #pragma omp parallel shared(particles, f, it, delta_t, fileName, file, softening, dim)
    {
        // Set the number of threads
        omp_set_num_threads(omp_get_max_threads());
        //create a vector of thread that contains a vector of forces of dimension Dimension
        std::array<double, Dimension> force;

        for (int z = 0; z < it; ++z){
            //calculate forces and add it to a local array for each thread
            #pragma omp for
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
                        // Call the collision method
                        q.manageCollision(k, 0.0);
                    }

                    force = f.calculateForce(k, q);
                    for (int y = 0; y < Dimension; ++y) {
                        local_forces[omp_get_thread_num()][i][y] += force[y];
                        local_forces[omp_get_thread_num()][j][y] -= force[y];
                    }
                }
            }

            #pragma omp barrier

            //sum all the forces calculated by each thread
            #pragma omp for
            for (int i = 0; i < (*particles).size(); ++i) {
                Particle<Dimension> &q = (*particles)[i];
                q.resetForce();
                for (int j = 0; j < omp_get_max_threads(); ++j) {
                    q.addForce(local_forces[j][i]);
                }
            }
            #pragma omp barrier

            #pragma omp for
            for(int k = 0; k < omp_get_max_threads(); ++k){
                for (int j = 0; j < (*particles).size(); ++j) {
                    for (int y = 0; y < Dimension; ++y) {
                        local_forces[k][j][y] = 0.0;
                    }
                }
            }

            #pragma omp barrier

            //update the position of the particles
            #pragma omp for
            for (int i = 0; i < (*particles).size(); ++i) {
                Particle<Dimension> &q = (*particles)[i];
                q.update(delta_t);
            }
            #pragma omp barrier

            //print forces of particles
            //#pragma omp single
            //{
            //    for (int i = 0; i < (*particles).size(); ++i) {
            //        Particle<Dimension> &q = (*particles)[i];
            //        for(int j = 0; j < Dimension; ++j)
            //            std::cout << q.getForce()[j] << " ";
            //       std::cout << std::endl;
            //    }
            //}

            //don't do this in parallel
            #pragma omp single
            {
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
                        std::cout << "Unable to open file in parallelSimulation ";
                    }
                }
            }
            #pragma omp barrier
        }
    }
} 

int main() {
    
    // Define of simulation variables
    const int d = 2; //2D or 3D
    const double delta_t = 0.01; // in seconds
    const double dim = 100; // Dimension of the simulation area
    int it = 1000; // number of iteration
    int n = 100; // number of particles
    double softening = 0.7; // Softening parameter
    time_t start, end;
    std::vector<Particle<d>> particles; // Create a vector of particles
    Force<d>* f = new CustomForce<d>(20); // Create force
    std::string fileName = "../../graphics/coordinates.txt"; // File name


    // Generate n random particles
    particles = generateRandomParticlesSerial<d>(n, dim, 1, 99, 50, 50, 1, 10, false);

    //generate two particles with velocity null
    //std::array<double, d> pos1 = {0, 0};
    //std::array<double, d> pos2 = {20, 0};
    //std::array<double, d> vel1 = {0, 0};
    //std::array<double, d> vel2 = {0, 100};

    //Particle(int id, double p, std::array<double, Dimension> pos, std::array<double, Dimension> v, double radius, bool type)
//    Particle<d> p1(0, 100, pos1, vel1, 1, false);
//    Particle<d> p2(1, 1, pos2, vel2, 1, false);
//    particles.push_back(p1);
//    particles.push_back(p2);
//
//    ////generate onther particle
//    std::array<double, d> pos3 = {50, 50};
//    std::array<double, d> vel3 = {0, 0};
//    Particle<d> p3(2, 1, pos3, vel3, 1, false);
//    particles.push_back(p3);


    // Apri il file in modalità di scrittura
    std::ofstream outputFile("particles.txt");

    // Verifica se il file è stato aperto correttamente
    if (outputFile.is_open()) {
        // Scrivi le particelle nel file
        for (const auto& particle : particles) {
            outputFile << particle.getId() << "," << particle.getProperty() << ",";
            for (size_t i = 0; i < d; ++i ) outputFile << particle.getPos()[i] << ",";
            for (size_t i = 0; i < d; ++i ) outputFile << particle.getVel()[i] << ",";
            outputFile << particle.getRadius() << "," << particle.getType();
            outputFile << "\n";
        }

        // Chiudi il file dopo aver scritto tutte le particelle
        outputFile.close();

        std::cout << "Particelle scritte con successo nel file particles.txt" << std::endl;
    } else {
        // Se il file non può essere aperto, stampa un messaggio di errore
        std::cerr << "Impossibile aprire il file particles.txt" << std::endl;
        return 1; // Restituisci un codice di errore
    }



    // Print on file the initial state of the particles
    std::ofstream file(fileName); // Open file
    if (file.is_open()) {

        // Write on file the total number of particles and the size of the area of the simulation
        file << particles.size() << std::endl;
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


    // Parallel simulation
    time(&start);
    parallelSimulation<d>(it, &particles, dim, softening, delta_t, fileName, file, *f);
    time(&end);
    printf("Time taken by the parallel implementation: %ld seconds\n", end - start);


    // Print the final state of the particles
    std::cout << "--------------------------------------------\n";
    std::cout << "Final state:\n";
    for (const Particle<d> &p : particles) {
       p.printStates();
    }
    file.close();


    return 0;
    
}
