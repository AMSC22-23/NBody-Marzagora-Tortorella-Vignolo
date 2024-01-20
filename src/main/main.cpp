#include <vector>
#include "particle.hpp" 
#include "force.hpp" 
#include "treenode.hpp"
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
#include <cstring>

// Function to create the tree structure from a vector of particles
template<size_t Dimension>
TreeNode<Dimension>* createTree(std::vector<Particle<Dimension>>& particles, double dimSimulationArea) {
    if (particles.empty()) {
        return nullptr;
    }

    // Determine the simulation boundaries (this might need to be adjusted)
    double minX =  -1 * dimSimulationArea * 0.5;
    double maxX = dimSimulationArea * 0.5;
    double minY = minX;
    double maxY = maxX;

    std::cout<< dimSimulationArea<<"->"<< minX << ", " << maxX <<endl;

    // Create the root of the tree
    double width = std::max(maxX - minX, maxY - minY);
    TreeNode<Dimension>* root = new TreeNode<Dimension>(minX + width / 2, minY + width / 2, width);

    // Insert each particle into the tree
    for (auto& particle : particles) {
        root->insert(&particle);
    }

    return root;
}


/**
 * @brief Template function that generates particles in order to perform the simulation. 
 * First it creates randomly (using the random library) the distributions for all the parameters, then it checks if the function struggles to generate non-overlapping particles: 
 * if the function has to regenerate the particle's position more than maxRetry times in a row, then the program ends and print an error message. 
 * Otherwise, it generates randomly the radius and the position of the particles (inside the distributions generated before). After that it checks that the newly generated particle
 * is not overlapping with any existing circle: it calculates the distance between the particle that is being created and the particle p in particles vector. If the particles do not
 * overlap, the remaining properties of the particles are randomly generated and, lastly, the particle is added to the vector of particles; otherwise, variable counter is increased.
 * Finally, this method returns a vector of particles
 * 
 * @tparam Dimension Number of dimensions of the simulation  
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
 * @tparam Dimension Number of dimensions of the simulation   
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
void serialSimulation(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f, int speedup){
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
        
        if(z%speedup==0){
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
 * @tparam Dimension Number of dimensions of the simulation 
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
*
* @tparam Dimension Number of dimensions of the simulation 
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

// Funzione per generare un vettore di 7 particelle con posizioni specifiche
template<size_t Dimension>
std::vector<Particle<Dimension>> generateTreeTestParticles() {
    // Creazione di un vettore di 7 particelle con posizioni specifiche
    std::vector<Particle<Dimension>> particles;

    // Specifica delle posizioni per le 7 particelle
    std::array<double, Dimension> pos1 = {1.2237, 10.4216};
    std::array<double, Dimension> pos2 = {4.05444, -69.0231}; 
    std::array<double, Dimension> pos3 = {-55.1028, 64.0666};
    std::array<double, Dimension> pos4 = {57.4899, -71.0595}; 
    std::array<double, Dimension> pos5 = {24.387, 39.731};
    std::array<double, Dimension> pos6 = {32.3245, -20.8399}; 
    std::array<double, Dimension> pos7 = {-16.4472, 4.49095};

    // Aggiungi le particelle al vettore
    particles.push_back(Particle<Dimension>(0, 1.0, pos1, {1.0, 1.0}, 5.0, false));
    particles.push_back(Particle<Dimension>(1, 1.5, pos2, {2.0, 2.0}, 7.0, false));
    particles.push_back(Particle<Dimension>(2, 1.0, pos3, {1.0, 1.0}, 5.0, false));
    particles.push_back(Particle<Dimension>(3, 1.5, pos4, {2.0, 2.0}, 7.0, false));
    particles.push_back(Particle<Dimension>(4, 1.0, pos5, {1.0, 1.0}, 5.0, false));
    particles.push_back(Particle<Dimension>(5, 1.5, pos6, {2.0, 2.0}, 7.0, false));
    particles.push_back(Particle<Dimension>(6, 1.0, pos7, {1.0, 1.0}, 5.0, false));


    return particles;
}


/**
 * @brief Template function that generates two particles where one stays still and the other orbits around it.
 * For test purposes only.
 * 
 * @tparam Dimension Number of dimensions of the simulation 
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

//create a function that generates 7 particles with specific positions



/**
 * @brief Template function that writes on file the total number of particles and the size of the area of the simulation and then the initial state of the particles in a file.
 * 
 * @tparam Dimension Number of dimensions of the simulation 
 * @param particles Reference to the vector of particles that are to be printed on the file
 * @param dim Dimension of the simulation area
 * @param fileName Name of the file in which the initial states of the particles are going to be written
 * @param file Reference to the file in which the initial states of the particles are going to be written
 * @param it Number of iterations
 **/

template<size_t Dimension>
void printInitialStateOnFile(std::vector<Particle<Dimension>>* particles, int dim, std::string fileName, std::ofstream& file, int it, int speedUp){
    
    if (file.is_open()) {

        
        file << (*particles).size() << std::endl;
        file << dim << std::endl;
        file << Dimension << std::endl;
        file << it / speedUp << std::endl;

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
 * @brief Template function that wraps main function for the 2D simulation: calls the function that generates the particles, the one that prints the initial states on the file and then
 * calls the function which starts the simulation chosen by the user.
 * 
 * @tparam Dimension Number of dimensions of the simulation 
 * @param simType Simulation type: 0 for serial, 1 for parallel
 * @param forceType Type of the force: g for gravitational force, c for coulomb force
 * 
 **/
template<size_t Dimension>
void main2DSimulationBarnesHut(int forceType, int simType, double delta_t, int dimSimulationArea, int iterationNumber, int numParticles, int mass, int maxVel, int maxRadius, double softening, std::string fileName, int speedUp){
    
    time_t start, end;
    std::vector<Particle<Dimension>> particles; 
    numParticles = 7;
    dimSimulationArea = 80;
    
    Force<Dimension>* f;
    if(forceType == 1 ) f = new CoulombForce<Dimension>();
    else f = new GravitationalForce<Dimension>();
    
    std::ofstream file(fileName);

    //start = time(NULL);
    //particles = generateRandomParticles<Dimension>(numParticles, dimSimulationArea, (forceType)? -mass:1, mass, maxVel, 1, maxRadius, forceType);
    //end = time(NULL);
    //std::cout << "Time taken by generateRandomParticles function: " << end - start << " seconds" << std::endl;

    particles = generateTreeTestParticles<Dimension>();

    std::cout << "\nPosition of particles generated: "<< std::endl;
    for(auto p:particles){
        std::cout<< p.getPos()[0]<< ", "<< p.getPos()[1]<<endl;
    }

    double deltaTime = 0.01; // Δt
    double totalTime = 0.1;  // Durata totale della simulazione
    double theta = 0.5;

    TreeNode<Dimension>* treeRoot = createTree(particles, 2*dimSimulationArea);
    
    std::cout << "\n\nTree structure: "<< std::endl;
    //After building the tree
    if (treeRoot != nullptr) {
        treeRoot->printTree(); // Display the tree structure
    }

    // TEST calculateNetForce on a particle: 

    // Choose a particle for which you want to calculate the net force
    Particle<Dimension>* targetParticle = &particles[0]; // example: first particle in the list

    // Reset force on the target particle before calculation
    targetParticle->resetForce();

    // Optionally, you can print the state of the particle to see the net force
    std::cout << "\n\nInitial state of targetParticle: "<< std::endl;
    targetParticle->printStates();

    // Calculate the net force on the target particle
    std::cout << "\n\nStart of calculateNetForce: "<< std::endl;
    calculateNetForce(treeRoot, targetParticle, theta);

    // Optionally, you can print the state of the particle to see the net force
    std::cout << "\n\nFinal state of targetParticle after calculateNetForce: "<< std::endl;
    targetParticle->printStates();
    
    // Clean up the tree
    delete treeRoot;

    //for (double currentTime = 0; currentTime < totalTime; currentTime += deltaTime) {
    //    // Costruisci l'albero di Barnes-Hut con le particelle attuali
    //    
    //    // After building the tree
    //    if (treeRoot != nullptr) {
    //        treeRoot->printTree(); // Display the tree structure
    //    }
    //
    //    // Calcola le forze su ogni particella
    //    calculateForces(treeRoot, 0.5);
    //
    //    // Aggiorna la posizione e la velocità di ogni particella
    //    for (auto& particle : particles) {
    //        UpdateParticle(particle, deltaTime);
    //    }
    //
    //    // Eventualmente, salva lo stato attuale per l'analisi o la visualizzazione
    //}

}



/*
To calculate the net force acting on body b, use the following recursive procedure, starting with the root of the quad-tree:

1)If the current node is an external node (and it is not body b), calculate the force exerted by the current node on b, and add this amount to b’s net force.
2)Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal node as a single body, and calculate the force it exerts on body b, and add this amount to b’s net force.
3)Otherwise, run the procedure recursively on each of the current node’s children.

*/
template <size_t Dimension>
void calculateNetForce(TreeNode<Dimension>* node, Particle<Dimension>* b, double theta) {
    std::cout << "Calculating net force for particle " << b->getId() << std::endl;

    if (node == nullptr) {
        std::cout << "Node is null. Exiting function." << std::endl;
        return;
    }else{
        if(node->getParticle() != nullptr)
            std::cout << "node:  " << node->getParticle()->getId() << std::endl;

    }

    // Check if the node is a leaf and not the same as the particle we're calculating the force on
    if (node->isLeaf() && node->getParticle() != nullptr && node->getParticle() != b) {
        std::cout << "Node is a leaf and contains a different particle. Applying direct force." << std::endl;
        ApplyDirectForce(b, node->getParticle()); // Apply direct force
    } 
    else if (!node->isLeaf()) {
        std::cout << "Node is not a leaf. Checking if approximation can be used." << std::endl;

        // Compute s/d ratio (size of the region represented by the node divided by distance)
        double s = node->getWidth(); // Assuming getWidth() returns the size of the node region

        // TODO: check memory access here (segmentation fault problem occours here)
        double d = std::sqrt(b->squareDistance(*node->getParticle())); // Using squareDistance from Particle class


        // Check for division by zero or very small d
        if (d == 0) {
            std::cout << "Distance is zero, avoiding division by zero." << std::endl;
            return;
        }
        
        std::cout << "s/d ratio: " << s/d << " (Threshold: " << theta << ")" << std::endl;
        if (s / d < theta) {
            std::cout << "Using approximation to apply force." << std::endl;
            ApplyApproximateForce(b, node); // Apply approximate force
        } 
        else {
            std::cout << "s/d ratio too high. Recursing into children nodes." << std::endl;
            // Recursively calculate force on children
            for (auto& child : node->getChildren()) {
                calculateNetForce(child, b, theta);
            }
        }
    }
        std::cout << "WTF" << std::endl;
}

template <size_t Dimension>
double computeRatio(const std::array<double, Dimension>& pos, const std::array<double, Dimension>& centerMass, double width) {
    double sum = 0.0;

    for (size_t i = 0; i < Dimension; ++i) {
        sum += std::pow(pos[i] - centerMass[i], 2);
    }

    return width/std::sqrt(sum);
}

template <size_t Dimension>
void ApplyApproximateForce(Particle<Dimension>* particle, TreeNode<Dimension>* node) {
    double G = 6.67430e-11; // Gravitational constant

    // Calculate the distance between the particle and the node's center of mass
    std::array<double, Dimension> diff;
    for (size_t i = 0; i < Dimension; ++i) {
        diff[i] = node->getTotalCenter()[i] - particle->getPos()[i];
    }
    double distSquared = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);

    // Avoid division by zero or very small distances
    if (distSquared < 1e-10) {
        return;
    }

    // Calculate force magnitude
    double forceMagnitude = G * (particle->getProperty() * node->getTotalMass()) / distSquared;

    // Normalize the direction and scale by force magnitude
    double length = std::sqrt(distSquared);
    std::array<double, Dimension> force;
    for (size_t i = 0; i < Dimension; ++i) {
        force[i] = (diff[i] / length) * forceMagnitude;
    }

    // Add the calculated force to the particle
    particle->addForce(force);
}

template <size_t Dimension>
void ApplyDirectForce(Particle<Dimension>* particle1, Particle<Dimension>* particle2) {
    std::cout<<"ApplyDirectForce"<<endl;
    double G = 6.67430e-11; // Gravitational constant

    double distSquared = particle1->squareDistance(*particle2);

    // Avoid division by zero or very small distances
    if (distSquared < 1e-10) {
        return;
    }

    // Calculate force magnitude
    double forceMagnitude = G * (particle1->getProperty() * particle2->getProperty()) / distSquared;

    // Calculate force direction
    std::array<double, Dimension> forceDirection;
    for (size_t i = 0; i < Dimension; ++i) {
        forceDirection[i] = particle2->getPos()[i] - particle1->getPos()[i];
    }

    // Normalize the direction and scale by force magnitude
    double length = std::sqrt(std::inner_product(forceDirection.begin(), forceDirection.end(), forceDirection.begin(), 0.0));
    std::array<double, Dimension> force;
    for (size_t i = 0; i < Dimension; ++i) {
        force[i] = (forceDirection[i] / length) * forceMagnitude;
    }

    // Add the calculated force to the particles
    particle1->addForce(force);
    particle2->addForce(force);
}


/**
 * @brief Template function that wraps main function for the 2D simulation: calls the function that generates the particles, the one that prints the initial states on the file and then
 * calls the function which starts the simulation chosen by the user.
 * 
 * @tparam Dimension Number of dimensions of the simulation 
 * @param simType Simulation type: 0 for serial, 1 for parallel
 * @param forceType Type of the force: g for gravitational force, c for coulomb force
 * 
 **/
template<size_t Dimension>
void main2DSimulation(int forceType, int simType, double delta_t, int dimSimulationArea, int iterationNumber, int numParticles, int mass, int maxVel, int maxRadius, double softening, std::string fileName, int speedUp){
    
    time_t start, end;
    std::vector<Particle<Dimension>> particles; 
    
    Force<Dimension>* f;
    if(forceType == 1 ) f = new CoulombForce<Dimension>();
    else f = new GravitationalForce<Dimension>();
    
    std::ofstream file(fileName);

    start = time(NULL);
    particles = generateRandomParticles<Dimension>(numParticles, dimSimulationArea, (forceType)? -mass:1, mass, maxVel, 1, maxRadius, forceType);
    end = time(NULL);
    std::cout << "Time taken by generateRandomParticles function: " << end - start << " seconds" << std::endl;

    printInitialStateOnFile(&particles, dimSimulationArea, fileName, file, iterationNumber, speedUp);

    if(simType == 0){
        start = time(NULL);
        serialSimulation<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, fileName, file, *f, speedUp);
        end = time(NULL);
        std::cout << "Time taken by serial simulation: " << end - start << " seconds" << std::endl;

    }else if(simType == 1){
        start = time(NULL);
        parallelSimulation<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, fileName, file, *f, speedUp);
        end = time(NULL);
        std::cout << "Time taken by parallel simulation: " << end - start << " seconds" << std::endl;
    }

    file.close();
}

/**
 * @brief Template function that wraps main function for the 3D simulation: calls the function that generates the particles, the one that prints the initial states on the file and then
 * calls the function which starts the simulation chosen by the user.
 * 
 * @tparam Dimension Number of dimensions of the simulation 
 * @param simType Simulation type: 0 for serial, 1 for parallel
 * @param forceType Type of the force: g for gravitational force, c for coulomb force
 * 
 **/
template<size_t Dimension>
void main3DSimulation(int forceType, int symType, double delta_t, int dimSimulationArea, int iterationNumber, int numParticles, int mass, int maxVel, int maxRadius, int softening, std::string fileName, int speedUp){
    
    time_t start, end; 
    std::vector<Particle<Dimension>> particles;
    Force<Dimension>* f;
    if(forceType == 1 ) f = new CoulombForce<Dimension>();
    else  f = new GravitationalForce<Dimension>(); 

    std::ofstream file(fileName); 

    start = time(NULL);
    particles = generateRandomParticles<Dimension>(numParticles, dimSimulationArea, (forceType)? -mass:1, mass, maxVel, 1, maxRadius, forceType);
    end = time(NULL);
    std::cout << "Time taken by generateRandomParticles function: " << end - start << " seconds" << std::endl;
    
    printInitialStateOnFile(&particles, dimSimulationArea, fileName, file, iterationNumber, speedUp);

    if(symType == 0){
        start = time(NULL);
        serialSimulation<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, fileName, file, *f, speedUp);
        end = time(NULL);
        std::cout << "Time taken by serial simulation: " << end - start << " seconds" << std::endl;

    }else if(symType == 1){
       start = time(NULL);
        parallelSimulation<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, fileName, file, *f, speedUp);
        end = time(NULL);
        std::cout << "Time taken by parallel simulation: " << end - start << " seconds" << std::endl;
    } 

    file.close();
}

void showHelp() {
    std::cout << "Change the following parameters if you don't want to run the default simualtion: " <<std::endl;
    std::cout << "      -h : prints helper " << std::endl ;
    std::cout << "      -dim <int> : number of dimension of the simulation (2D,3D) " << std::endl;
    std::cout << "      -simT <int> : simulation type (0 for serial, 1 for parallel) " << std::endl;
    std::cout << "      -force <int> : type of the force of the simulation " << std::endl;
    std::cout << "      -delta <double>: time step of the simulation " << std::endl;
    std::cout << "      -simA <int> : dimension of the simulation area " << std::endl;
    std::cout << "      -it <int> : iteration number " << std::endl;
    std::cout << "      -numP <int> : number of particles " << std::endl;
    std::cout << "      -maxPr <int> : maximum property of the particle(mass, charge) " << std::endl;
    std::cout << "      -maxVel <int> : maximum velocity of the particles " << std::endl;
    std::cout << "      -maxR <int> : maximum radius of the particles " << std::endl;
    std::cout << "      -soft <double> : softener of the particles " << std::endl;
    std::cout << "      -spUp <int> : speedup of the simulation  " << std::endl;
    std::cout << "      -file <std::string> : file in which the output is written  " << std::endl;
}

/**
 * @brief Main function which asks the users for the simulation type, the dimensions and the force type and then runs the simulation accordingly
*/

int main(int argc, char** argv) {

    int dim = 2; 
    int simType = 0;
    int forceType = 0;
    double delta_t = 0.01;
    int dimSimulationArea = 500; 
    int iterationNumber = 1000; 
    int numParticles = 50;
    int mass = 50; 
    int maxVel = 50; 
    int maxRadius = 5; 
    double softening = 0.7;
    int speedUp = 1;
    std::string fileName = "../graphics/coordinates.txt";

     if (argc < 2) {
        char a;
        std::cout << "Enter 'd' to run the default simulation: " <<std::endl;
        std::cin >> a;
        if(a == 'd') dim = 2; 
        else {
            showHelp();
            return 1;
        }
    }

    for(int i = 1; i < argc; i++){
        if (strcmp(argv[i], "-dim") == 0) {
             if (++i < argc) {
                dim = atoi(argv[i]);
                if(dim != 2 && dim != 3) {
                    std::cout << "No feasible dimension." << std::endl;
                    return 1;
                }
            } else {
                std::cout <<" Error: flag -dim requires values 2 or 3 to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-simT") == 0) {
             if (++i < argc) {
                simType = atoi(argv[i]);
                if(simType != 0 && simType != 1) {
                    std::cout << "No feasible simulation." << std::endl;
                    return 1;
                }
            } else {
                std::cout <<" Error: flag -simT requires values 0 for serial or 1 for parallel to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-force") == 0) {
             if (++i < argc) {
                forceType = atoi(argv[i]);
                if(forceType != 0 && forceType != 1) {
                    std::cout << "No feasible force." << std::endl;
                    return 1;
                }
            } else {
                std::cout <<" Error: flag -force requires values 0 for gravitational force or 1 for Coulomb force to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-delta") == 0) {
             if (++i < argc) {
                double delta = std::__cxx11::stof(argv[i]);
                if(delta < 0){
                    std::cout << "No feasible delta t." << std::endl;
                    return 1;
                }
                delta_t = delta;
            } else {
                std::cout <<" Error: flag -delta requires a positive value to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-simA") == 0) {
             if (++i < argc) {
                int simArea= atoi(argv[i]);
                if(simArea < 0){
                    std::cout << "No feasible simulation area." << std::endl;
                    return 1;
                }
                dimSimulationArea = simArea;
            } else {
                std::cout <<" Error: flag -simA requires a positive value to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-it") == 0) {
             if (++i < argc) {
                int it = atoi(argv[i]);
                if(it < 0) {
                    std::cout << "No feasible number of iterations." << std::endl;
                    return 1;
                }
                iterationNumber = it;
            } else {
                std::cout <<" Error: flag -it requires a positive value to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-numP") == 0) {
             if (++i < argc) {
                int numP = atoi(argv[i]);
                if(numP < 0) {
                    std::cout << "No feasible number of particles." << std::endl;
                    return 1;
                }
                numParticles = numP;
            } else {
                std::cout <<" Error: flag -numP requires a positive value to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-maxPr") == 0) {
             if (++i < argc) {
                int maxPr = atoi(argv[i]);
                if(maxPr < 0){
                    std::cout << "No feasible value of maximum property." << std::endl;
                    return 1;
                }
                mass = maxPr;
            } else {
                std::cout <<" Error: flag -maxPr requires a positive value to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-maxVel") == 0) {
             if (++i < argc) {
                int maxV = atoi(argv[i]);
                if(maxV < 0){
                    std::cout << "No feasible value of radius of the particles." << std::endl;
                    return 1;
                }
                maxVel = maxV;
            } else {
                std::cout <<" Error: flag -maxVEl requires a positive value to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-maxR") == 0) {
             if (++i < argc) {
                int maxR = atoi(argv[i]);
                if(maxR < 0) {
                    std::cout << "No feasible value of radius of the particles." << std::endl;
                    return 1;
                }
                maxRadius = maxR;
            } else {
                std::cout <<" Error: flag -maxR requires a positive value to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-soft") == 0) {
             if (++i < argc) {
                double soft = std::__cxx11::stof(argv[i]);
                if(soft < 0){ 
                    std::cout << "No feasible value of softening." << std::endl;
                    return 1;
                }
                softening = soft;
            } else {
                std::cout <<" Error: flag -soft requires a positive value to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-spUp") == 0) {
             if (++i < argc) {
                int spUp = atoi(argv[i]);
                if(spUp < 0) {
                    std::cout << "No feasible value of speedup." << std::endl;
                    return 1;
                }
                speedUp = spUp;
            } else {
                std::cout <<" Error: flag -spUp requires a positive value to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-file") == 0) {
             if (++i < argc) {
                std::string file = argv[i];
                fileName = file;
            } else {
                std::cout <<" Error: flag -file requires a valid file name to work. " <<std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-h") == 0) {
            showHelp(); 
            return 0;
        }
        
    }
    

    if(dim == 2) {
        main2DSimulationBarnesHut<2>(forceType, simType, delta_t, dimSimulationArea, iterationNumber, numParticles, mass, maxVel, maxRadius, softening, fileName, speedUp); 
        //main2DSimulation<2>(forceType, simType, delta_t, dimSimulationArea, iterationNumber, numParticles, mass, maxVel, maxRadius, softening, fileName, speedUp);
    } else if (dim == 3) {
        main3DSimulation<3>(forceType, simType, delta_t, dimSimulationArea, iterationNumber, numParticles, mass, maxVel, maxRadius, softening, fileName, speedUp);    
    }

    return 0;
    
}
