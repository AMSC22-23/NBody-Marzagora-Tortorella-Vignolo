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

/**
 * @brief Template function that generates particles in order to perform the simulation. 
 * First it creates randomly (using the random library) the distributions for all the parameters, then it checks if the function struggles to generate non-overlapping particles: 
 * if the function has to regenerate the particle's position more than maxRetry times in a row, then the program ends and print an error message. 
 * Otherwise, it generates randomly the radius and the position of the particles (inside the distributions generated before). After that it checks that the newly generated particle
 * is not overlapping with any existing circle: it calculates the distance between the particle that is being created and the particle p in particles vector. If the particles do not
 * overlap, the remaining properties of the particles are randomly generated and, lastly, the particle is added to the vector of particles; otherwise, variable counter is increased.
 * Finally, this method returns a vector of particles
 *   
 * @param N Number of particles to generate.
 * @param posBoundary Position boundary for particle positions (default: 100).
 * @param minProperty Minimum value of property for generated particles (default: 1).
 * @param maxProperty Maximum value of property for generated particles (default: 99).
 * @param maxVel Maximum velocity for particles (default: 100).
 * @param minRadius Minimum radius for particles (default: 0).
 * @param maxRadius Maximum radius for particles (default: 15).
 * @param type Type flag for particles to identify either they are particles of the gravitational or the Coulomb force (default: false, which refers to gravitational force).
 * @return A vector of particles with randomly generated properties.
 **/
template<size_t Dimension>
std::vector<Particle<Dimension>> generateRandomParticles(int N, int posBoundary = 100, int minProperty = 1, int maxProperty = 99, int maxVel = 100, int minRadius = 0, int maxRadius = 15, bool type = false) {

    int maxRetry  = 15, counter = 0;
    bool overlapping = false;
    double x, y, r, property, squareDistance = 0.0;
    std::vector<Particle<Dimension>> particles;
    std::array<double, Dimension> pos, vel;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> disProperty(minProperty, maxProperty);
    std::uniform_real_distribution<double> disVel(-maxVel, maxVel);
    std::uniform_real_distribution<double> disRadius(minRadius, maxRadius);
    std::uniform_real_distribution<double> disPos(-posBoundary+maxRadius, posBoundary-maxRadius);

    while (particles.size() < N) {

        if(counter == maxRetry)
            throw std::runtime_error("ERROR: the dimension of the simulation area is too little, please specify a bigger area or generate less particles.");

        overlapping = false;

        r = disRadius(gen);

        for(size_t i=0; i<Dimension; ++i) pos[i] = disPos(gen);

        for (const Particle<Dimension> &p : particles) {

            squareDistance = 0.0;
            for(size_t i = 0; i < Dimension; ++i) 
                squareDistance = squareDistance + ((p.getPos()[i] - pos[i])*(p.getPos()[i] - pos[i]));
            if (sqrt(squareDistance) < p.getRadius() + r) { 
                counter++;
                overlapping = true;
                break; 
            }
        }

        if (!overlapping) {
            counter = 0;

            property = disProperty(gen);

            for(size_t i=0; i<Dimension; ++i)
                vel[i] = disVel(gen);

            Particle<Dimension> p(particles.size(), property, pos, vel, r, type);

            particles.push_back(p);
        }else{
            counter++;
        }
        
    }

    return particles;
}


/**
 * @brief Template function that executes the NBody simulation serially for a given number of iterations. 
 * For each particle,  it checks first if the particle hits the boundary and manages that collision; then, it checks for collision between particles and call the manageCollision function to
 * take care of it. Finally it computed the force between the particles and adds it to a vector. 
 * After checking all the particles, it updates the values calculated before and resets the vector; finally, it writes the updated positions of the particles on the file after delta_t time.
 *   
 * @param it Number of iterations
 * @param particles Reference to the vector of particles
 * @param dim Number of dimensions of the simulation (2D,3D)
 * @param softening Overhead to avoid particles overlapping and fusing together
 * @param delta_t Time step after which the simulation is updated
 * @param fileName Name of tile in which the function writes the coordinates of the particles computed during the simulation
 * @param file Reference to the file in which the coordinates are written
 * @param f Reference to the Force object responsible for calculating particle interactions.
 **/
template<size_t Dimension>
void serialSimulation(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f){
    std::array<double,Dimension> force_qk;

    for (int z = 0; z < it; ++z){
        for (int i = 0; i < (*particles).size(); ++i) {
            Particle<Dimension> &q = (*particles)[i];

            if(q.hitsBoundary(dim)){
                q.manageCollision(q, dim);
            }

            for (int j = i + 1; j < (*particles).size(); j++) {
                Particle<Dimension> &k = (*particles)[j];

                if(q.squareDistance(k) < (((q.getRadius() + k.getRadius())*(q.getRadius() + k.getRadius())))*softening){

                    q.manageCollision(k, 0.0);
                }

                force_qk = f.calculateForce(q, k);
                q.addForce(force_qk);
                for(size_t i = 0; i < Dimension; ++i) force_qk[i] = -force_qk[i];
                k.addForce(force_qk);

            }
            z==it-1? q.update(delta_t):q.updateAndReset(delta_t);
        }
        
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

/**
 * @brief Template function that executes the NBody simulation in parallel for a given number of iterations, 
 * using OpenMP directives to parallelize the simulation loop, calculating forces and updating particle positions concurrently. 
 * Firstly, the function initializes the number of threads for the parallel section based on the minimum between the size of 
 * the particle vector or the maximum available threads. Then it starts the simulation, assigning blocks to threads dynamically:
 * collision among particles or between a particle and the boundary and computation of the forces are computed concurrently and then
 * summed up. The forces are then updated and then the local vector is reset; finally, the function periodically writes the 
 * updated positions inside a file.
 * 
 * @param it Number of iterations
 * @param particles Reference to the vector of particles
 * @param dim Number of dimensions of the simulation (2D,3D)
 * @param softening Overhead to avoid particles overlapping and fusing together
 * @param delta_t Time step after which the simulation is updated
 * @param fileName Name of tile in which the function writes the coordinates of the particles computed during the simulation
 * @param file Reference to the file in which the coordinates are written
 * @param f Reference to the Force object responsible for calculating particle interactions.
 * @param u 
 **/

template<size_t Dimension>
void parallelSimulation(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f, int u){

    particles->size() < omp_get_max_threads()? omp_set_num_threads(particles->size()):omp_set_num_threads(omp_get_max_threads());  
    int num_threads = omp_get_max_threads();

    std::vector<std::vector<std::array<double, Dimension>>> local_forces(omp_get_max_threads(), std::vector<std::array<double, Dimension>>(particles->size()));
    
    std::array<double,Dimension> force;

    for (int z = 0; z < it; ++z){
        #pragma omp parallel shared(particles, f, it, delta_t, fileName, file, softening, dim, num_threads) private(force)
        {
           #pragma omp for schedule(dynamic, u)
           for (int i = 0; i < (*particles).size(); ++i) {
               Particle<Dimension> &k = (*particles)[i];

               
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

           #pragma omp for schedule(static, particles->size()/omp_get_num_threads())
           for (int i = 0; i < (*particles).size(); ++i) {
               Particle<Dimension> &q = (*particles)[i];
               q.resetForce();
               for (int j = 0; j < num_threads; ++j) {
                   q.addForce(local_forces[j][i]);
               }
           }


           #pragma omp for schedule(static, 1)
           for(int k = 0; k < num_threads; ++k){
               for (int j = 0; j < (*particles).size(); ++j) {
                   for (int y = 0; y < Dimension; ++y) {
                       local_forces[k][j][y] = 0.0;
                   }
               }
           }

           #pragma omp for schedule(static, particles->size()/num_threads)
           for (int i = 0; i < (*particles).size(); ++i) {
               Particle<Dimension> &q = (*particles)[i];
               q.update(delta_t);
           }        
        }
        
        if(z%u==0){
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
    }
}


/**
* @brief Template function that prints the particles state
* @param particles Refence to the vector of particles whose states are to be printed
*/
template<size_t Dimension>
void printAllParticlesStateAndDistance(const std::vector<Particle<Dimension>>* particles){
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


/**
 * @brief Template function that generates two particles where one stays still and the other orbits around it.
 * For test purposes only.
 * 
 * @param size Size of the orbit
 * @param constantForce Constant of the force applied to the particles
 * @return Vector fo particles for a two-body orbit test
 **/
template<size_t Dimension>
std::vector<Particle<Dimension>> generateOrbitTestParticles(double size, double costantForce){
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

    Particle<Dimension> p1(0, mass1, pos1, vel1, size/100, false);
    Particle<Dimension> p2(1, mass2, pos2, vel2, size/100, false);
    particles.push_back(p1);
    particles.push_back(p2);

    return particles;
}

/**
 * @brief Template function that writes on file the total number of particles and the size of the area of the simulation and then the initial state of the particles in a file.
 * @param particles Reference to the vector of particles that are to be printed on the file
 * @param dim Dimension of the simulation area
 * @param fileName Name of the file in which the initial states of the particles are going to be written
 * @param file Reference to the file in which the initial states of the particles are going to be written
 * @param it Number of iterations
 **/

template<size_t Dimension>
void printInitialStateOnFile(std::vector<Particle<Dimension>>* particles, int dim, std::string fileName, std::ofstream& file, int it){
    
    if (file.is_open()) {

        
        file << (*particles).size() << std::endl;
        file << dim << std::endl;
        file << Dimension << std::endl;
        file << it << std::endl;

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


/**
 * @brief Template function that wraps main function for the 2D simulation 
 * 
 * 
 **/
template<size_t Dimension>
void main2DSimulation(int simType, char forceType){

    const double delta_t = 0.01; // In seconds
    int dim = 10; // Dimension of the simulation area
    int it = 1000; // Number of iteration
    int n = 50; // Number of particles
    int mass = 50; // Mass
    int maxVel = 50; // Maximum velocity
    int maxRadius = 50; // Maximum radius of the particles
    double softening = 0.7; // Softening parameter
    time_t start, end; // Time variables
    std::vector<Particle<Dimension>> particles; 
    std::vector<Particle<Dimension>> particlesSerial;// Create a vector of particles
    std::string fileName = "../graphics/coordinates.txt"; // File name
    std::ofstream file(fileName); // Open file
    Force<Dimension>* f;

    if(forceType == 'c') f = new CoulombForce<Dimension>();
    else f = new GravitationalForce<Dimension>(); // Create force
    

    // Generate random particles
    start = time(NULL);
    particles = generateRandomParticles<Dimension>(n, dim, 1, mass, maxVel, 1, maxRadius, false);
    particlesSerial = particles;
    end = time(NULL);
    std::cout << "Time taken by generateRandomParticles function: " << end - start << " seconds" << std::endl;
    

    // Print on file the initial state of the particles
    printInitialStateOnFile(&particles, dim, fileName, file, it);

    // Start simulation
    if(simType == 0){
        start = time(NULL);
        serialSimulation<Dimension>(it, &particles, dim, softening, delta_t, fileName, file, *f);
        end = time(NULL);
        std::cout << "Time taken by serial simulation: " << end - start << " seconds" << std::endl;

    }else if(simType == 1){
        start = time(NULL);
        parallelSimulation<Dimension>(it, &particles, dim, softening, delta_t, fileName, file, *f, 1);
        end = time(NULL);
        std::cout << "Time taken by parallel simulation: " << end - start << " seconds" << std::endl;
    } else {
        std::cout << "You did not specify correctly the type of the simulation [serial/parallel]" << std::endl;
    }

    // In order to print the final state of the particles use: printAllParticlesStateAndDistance(&particles);
    //printAllParticlesStateAndDistance(&particles);

    file.close();
}

template<size_t Dimension>
void main3DSimulation(int symType, char forceType){
    // Simulation variables    
    const double delta_t = 0.01; // In seconds
    int dim = 500; // Dimension of the simulation area
    int it = 1000; // Number of iteration
    int n = 50; // Number of particles
    int mass = 50; // Mass
    int maxVel = 50; // Maximum velocity
    int maxRadius = 50; // Maximum radius of the particles
    double softening = 0.7; // Softening parameter
    time_t start, end; // Time variables
    std::vector<Particle<Dimension>> particles; // Create a vector of particles
    std::string fileName = "../graphics/coordinates.txt"; // File name
    std::ofstream file(fileName); // Open file
    Force<Dimension>* f;

    if(forceType == 'c') f = new CoulombForce<Dimension>();
    else  f = new GravitationalForce<Dimension>(); // Create force

    // Generate random particles
    start = time(NULL);
    particles = generateRandomParticles<Dimension>(n, dim, 1, mass, maxVel, 1, maxRadius, false);
    end = time(NULL);
    std::cout << "Time taken by generateRandomParticles function: " << end - start << " seconds" << std::endl;
    

    // Print on file the initial state of the particles
    printInitialStateOnFile(&particles, dim, fileName, file, it);

    // Start simulation
    if(symType == 0){
        start = time(NULL);
        serialSimulation<Dimension>(it, &particles, dim, softening, delta_t, fileName, file, *f);
        end = time(NULL);
        std::cout << "Time taken by serial simulation: " << end - start << " seconds" << std::endl;

    }else if(symType == 1){
       start = time(NULL);
        parallelSimulation<Dimension>(it, &particles, dim, softening, delta_t, fileName, file, *f, 1);
        end = time(NULL);
        std::cout << "Time taken by parallel simulation v2: " << end - start << " seconds" << std::endl;
    } else {
        std::cout << "You did not specify correctly the type of the simulation [serial/parallel]" << std::endl;
    }

    // In order to print the final state of the particles use: printAllParticlesStateAndDistance(&particles);
    //printAllParticlesStateAndDistance(&particles);

    file.close();
}

int main(int argc, char** argv) {

    int dim = 0;
    int simType = 0;
    char forceType = 'c';
    
    std::cout << "Insert 0 for serial simulation and 1 for parallel simulation: " <<std::endl;
    std::cin >> simType;
    while(simType != 0 && simType != 1){
        std::cout << "Invalid simulation type! Insert 0 for serial simulation and 1 for parallel simulation: " << std::endl;
        std::cin >> simType;
    }
    
    std::cout << "Insert dimension for the simulation: " <<std::endl;
    std::cin >> dim;
    
    while(dim != 2 && dim != 3){
        std::cout << "No feasible dimension, insert another one!" << std::endl;
        std::cout << "Insert dimension for the simulation: " <<std::endl;
        std::cin >> dim;
    }

    std::cout << "Insert c for Coulomb force or g for gravitational force: " <<std::endl;
    std::cin >> forceType;
    while(simType != 'c' && simType != 'g'){
        std::cout << "Invalid force type! Insert c for Coulomb force or g for gravitational force: " << std::endl;
        std::cin >> forceType;
    }

    if(dim == 2) {
        main2DSimulation<2>(simType, forceType);
    } else if (dim == 3) {
        main3DSimulation<3>(simType, forceType);    
    }

    return 0;
    
}
