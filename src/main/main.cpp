#include <vector>
#include "../utils/particle.hpp"
#include "../utils/force.hpp"
#include "../utils/quadtreeNode.hpp"
#include <cstdlib> 
#include <ctime>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
#include <cstring>
#include <memory>

template<size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> createQuadTree(std::vector<Particle<Dimension>>& particles, double dimSimulationArea) {
    if (particles.empty()) {
        std::cout<<"empty"<<std::endl;
        return nullptr;
    }

    // Determine the simulation boundaries (this might need to be adjusted)
    double minX =  -1 * dimSimulationArea * 0.5;
    double maxX = dimSimulationArea * 0.5;
    double minY = minX;
    double maxY = maxX;

    //std::cout<< dimSimulationArea<<"->"<< minX << ", " << maxX <<endl;
    //std::cout<<"num particles"<< particles.size() <<endl;

    // Create the root of the tree
    double width = std::max(maxX - minX, maxY - minY);
    std::unique_ptr<QuadtreeNode<Dimension>> root = std::make_unique<QuadtreeNode<Dimension>>(minX + width / 2, minY + width / 2, width, 20);

    // Insert each particle into the tree
    for (auto& particle : particles) {
        std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
        root->insert(p);
    }

    return root;
}

template<size_t Dimension>
void assignRegions(double regionWidth, double dim, double x, double y, std::vector<std::array<double, Dimension>>* centers){
    if(dim == regionWidth){
        centers->push_back({x-regionWidth/2, y+regionWidth/2});
        centers->push_back({x+regionWidth/2, y+regionWidth/2});
        centers->push_back({x-regionWidth/2, y-regionWidth/2});
        centers->push_back({x+regionWidth/2, y-regionWidth/2});

    }
    else{
        assignRegions(regionWidth, dim/2, x-dim/2, y+dim/2, centers);
        assignRegions(regionWidth, dim/2, x+dim/2, y+dim/2, centers);
        assignRegions(regionWidth, dim/2, x-dim/2, y-dim/2, centers);
        assignRegions(regionWidth, dim/2, x+dim/2, y-dim/2, centers);
    }

}



template<size_t Dimension>
std::vector<std::array<double, Dimension>> getSubRegionsCoordinates(int num_threads, double dim) {
    std::vector<std::array<double, Dimension>> regions;
    
    int rows = std::sqrt(num_threads);
    while (num_threads % rows != 0) {
        rows--;
    }
    int cols = num_threads / rows;

    double regionWidth = dim / cols;
    double regionHeight = dim / rows;

    /*for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double centerX = (j + 0.5) * regionWidth - dim / 2;
            double centerY = (i + 0.5) * regionHeight - dim / 2;
            regions.push_back({centerX, centerY});
        }
    }*/
    // funzione ricorsiva
    assignRegions(regionWidth, dim/2, 0.0, 0.0,&regions);

    return regions;
}

template<size_t Dimension>
std::array<double, Dimension> getSubRegionsDimension(int num_threads, double dim) {
    std::array<double, Dimension> regionsDim;
    
    int rows = std::sqrt(num_threads);
    while (num_threads % rows != 0) {
        rows--;
    }
    int cols = num_threads / rows;

    double regionWidth = dim / cols;
    double regionHeight = dim / rows;
    regionsDim[0] = regionWidth;
    regionsDim[1] = regionHeight;

    return regionsDim;
}



template<size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> merging(std::vector<std::unique_ptr<QuadtreeNode<Dimension>>>* roots, double subRegionsDimension){
    
    int size = (*roots).size(); // 4^n elements
    int k = 4;
    int y = 1;



    std::array<std::unique_ptr<QuadtreeNode<Dimension>>, 4> temp_roots;
    int num_threads = omp_get_max_threads();
    
    double x = 0.0;
    double z = 0.0;

    int l = 0;

    while (size >= 4) {

        if(size < num_threads){
            num_threads = size;
            omp_set_num_threads(num_threads);
        }

        #pragma omp parallel for private(temp_roots) shared(k,y) schedule(static, size/num_threads)
        for (int i = 0; i < size/4; i++) {

            for(int j = 0; j<4 ; j++){
                int index = k*i + j*y;
                swap(temp_roots[j],(*roots)[index]);

            }

            std::cout<<"subRegionsDimension : "<<subRegionsDimension*(l+1)<<std::endl;

            (*roots).at(k*i) = std::move(mergeRoots(&temp_roots, subRegionsDimension*(l+1)));

        }

        k = k*4;
        y = y*4;
        l++;
        size /= 4; // Reduce the size of the vector by a factor of 4
    }
    
    if((*roots)[0] != nullptr)
        (*roots)[0]->printTree();

    return std::move((*roots)[0]);
}


template<size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> mergeRoots(std::array<std::unique_ptr<QuadtreeNode<Dimension>>, 4>* roots, double subRegionsDimension){
    
    double x,y;

    std::unique_ptr<QuadtreeNode<Dimension>> intNode;
    double totx = 0.0;
    double toty = 0.0;

    std::cout<<"............................"<<std::endl;
    for(int i=0;i<4;i++){
        totx += (*roots)[i]->getCenter()[0];
        toty += (*roots)[i]->getCenter()[1];
         std::cout<<"totx "<<totx<<std::endl;
    }
    std::cout<<"............................"<<std::endl;

    x = totx/4.0;
    y = toty/4.0;

    intNode = std::make_unique<QuadtreeNode<Dimension>>(x, y, 2*subRegionsDimension, 20);

    std::cout<<"\n x = "<<x<<std::endl;
    std::cout<<" y = "<<y<<std::endl;

    for(int j = 0; j<4; j++){
        std::cout<<" \n --------------------------------- "<<std::endl;
        if((*roots)[j]->isLeaf()){
            if((*roots)[j]->getParticle() != nullptr){
                std::cout<<"Leaf node. Num particles: 1"<<std::endl;
                intNode->insert((*roots)[j]->getParticle());
            }else{
                std::cout<<"Leaf node. Num particles: 0"<<std::endl;
            }
        } 
        else{
            std::cout<<"Not a leaf node, Num particles: "<<(*roots)[j]->getCount()<<std::endl; 
            intNode->insertNode(std::move((*roots)[j]));
            std::cout<<"insertNode done"<<std::endl; 
        }
    }

    //if(intNode != nullptr)
    //    intNode->printTree();

    return intNode;
}

template<size_t Dimension>
int generateTreeParallel(std::vector<Particle<Dimension>>* particles, double dimSimulationArea){

    std::cout<<"start generateTreeParallel"<<std::endl;

    int num_threads = omp_get_max_threads(); // Example value

    std::cout<<"num_threads "<<num_threads<<std::endl;
    int num_quad = pow(4, ceil(log(num_threads)/log(4)));

    std::cout<<"num_quad "<<num_quad<<std::endl;

    std::vector<std::array<double, 2>> regions = getSubRegionsCoordinates<Dimension>(num_quad, dimSimulationArea);
    std::array<double, 2> subRegionsDimension = getSubRegionsDimension<Dimension>(num_quad, dimSimulationArea);

    for (const auto& region : regions) {
        std::cout << "Center: (" << region[0] << ", " << region[1] 
                << "); " << std::endl;
    }
    std::cout<<subRegionsDimension[0]<<" : "<<subRegionsDimension[1]<<std::endl;
    std::cout<<num_quad<<std::endl;

    std::vector<std::unique_ptr<QuadtreeNode<Dimension>>> roots;
    std::unique_ptr<QuadtreeNode<Dimension>> roots_;
    roots.reserve(num_quad);
    for(int i = 0; i<num_quad; i++) {
        roots.push_back(nullptr);
    }

    #pragma omp parallel for shared(particles, regions)
    for(int i = 0; i<num_quad; i++){
        //roots.push_back(nullptr);
        std::cout<<i<<std::endl;
        (roots).at(i) = std::move(createQuadTreeParallel(*particles, regions[i], subRegionsDimension[0], dimSimulationArea));
        //std::cout<<num_quad<std::endl;
    }
    
    for(int i = 0; i<num_quad; i++){
        std::cout<<"------------printTree------------"<<std::endl;
        roots[i]->printTree();
    }

    std::unique_ptr<QuadtreeNode<Dimension>> treeRoot;
    
    treeRoot = merging(&roots, subRegionsDimension[0]);

    // merge
    std::cout<<"end generateTreeParallel"<<std::endl;
    return 1;
}

template<size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> createQuadTreeParallel(std::vector<Particle<Dimension>>& particles, std::array<double, Dimension> center, double width, double dimSimulationArea) {
    
    std::cout<<"start createQuadTreeParallel"<<std::endl;

    if (particles.empty()) {
        std::cout<<"empty"<<std::endl;
        return nullptr;
    }

    // Determine the simulation boundaries (this might need to be adjusted)
    double minX = center[0] - width * 0.5;
    double maxX = center[0] + width * 0.5;
    double minY = center[1] - width * 0.5;
    double maxY = center[1] + width * 0.5;

    std::cout<< width <<" -> "<< minX << ", " << maxX <<endl;
    //std::cout<<"num particles"<< particles.size() <<endl;

    // Create the root of the tree
    std::unique_ptr<QuadtreeNode<Dimension>> root = std::make_unique<QuadtreeNode<Dimension>>(center[0], center[1], width, 20);

    // Insert each particle into the tree
    if(maxX == dimSimulationArea/2 && minY == -dimSimulationArea/2){
        std::cout<<"caso particolarissimo"<<std::endl;
        for (auto& particle : particles) {
            if(particle.getPos()[0] >= minX && particle.getPos()[0] <= maxX && particle.getPos()[1] <= maxY && particle.getPos()[1] >= minY){
                std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
                root->insert(p);
            }
        }

    } else if(minY == -dimSimulationArea/2){
        std::cout<<"minY == -dimSimulationArea/2"<<std::endl;
        for (auto& particle : particles) {
            if(particle.getPos()[0] >= minX && particle.getPos()[0] < maxX && particle.getPos()[1] <= maxY && particle.getPos()[1] >= minY){
                std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
                root->insert(p);
            }
        }
        
    } else if(maxX == dimSimulationArea/2){
        std::cout<<"maxX == dimSimulationArea/2"<<std::endl;
        for (auto& particle : particles) {
            if(particle.getPos()[0] >= minX && particle.getPos()[0] <= maxX && particle.getPos()[1] <= maxY && particle.getPos()[1] > minY){
                std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
                root->insert(p);
            }
        }
    } else {
        std::cout<<"caso normale"<<std::endl;
        for (auto& particle : particles) {
            if(particle.getPos()[0] >= minX && particle.getPos()[0] < maxX && particle.getPos()[1] <= maxY && particle.getPos()[1] > minY){
                std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
                root->insert(p);
            }
        }
    }

    std::cout<<"end  ------------"<<std::endl;
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
        //std::cout<<z<<std::endl;
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
                    std::cout << "e to open file";
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
    std::array<double, Dimension> pos1 = {2.0, 10.0};
    std::array<double, Dimension> pos2 = {4.0, -70.0}; 
    std::array<double, Dimension> pos3 = {-55.0, 65.0};
    std::array<double, Dimension> pos4 = {55.0, -70.0}; 
    std::array<double, Dimension> pos5 = {25.0, 40.0};
    std::array<double, Dimension> pos6 = {30.0, -20.0}; 
    std::array<double, Dimension> pos7 = {-15.0, 5.0};

    // Aggiungi le particelle al vettore
    particles.push_back(Particle<Dimension>(0, 10.0, pos1, {10.0, 10.0}, 5.0, false));
    particles.push_back(Particle<Dimension>(1, 10.5, pos2, {20.0, 20.0}, 7.0, false));
    particles.push_back(Particle<Dimension>(2, 10.0, pos3, {10.0, 10.0}, 5.0, false));
    particles.push_back(Particle<Dimension>(3, 10.5, pos4, {20.0, 20.0}, 7.0, false));
    particles.push_back(Particle<Dimension>(4, 10.0, pos5, {10.0, 10.0}, 5.0, false));
    particles.push_back(Particle<Dimension>(5, 10.5, pos6, {20.0, 20.0}, 7.0, false));
    particles.push_back(Particle<Dimension>(6, 10.0, pos7, {10.0, 10.0}, 5.0, false));


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
void serialBarnesHut(int iterationNumber, std::vector<Particle<Dimension>>* particles, int dimSimulationArea, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f, int speedUp){
    
    time_t start, end;
    //std::vector<Particle<Dimension>> particles; 

    double theta = 0.5;
    std::unique_ptr<QuadtreeNode<Dimension>> quadtree;
    int numParticles = (*particles).size();

    // Inizio dell'algoritmo di Barnes-Hut
    
    start = time(NULL);
    for (int iter = 0; iter < iterationNumber; ++iter) {
        
        quadtree = createQuadTree(*particles, 2*dimSimulationArea);

        // Calcolo delle forze per ogni particella
        for (int i = 0; i < numParticles; ++i) {

            //calculateNetForce(treeRoot, &particles[i], theta, *f);
            std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>((*particles)[i]);
            
            calculateNetForceQuadtree(quadtree, p, theta, f, dimSimulationArea, softening);
            if(p != nullptr)
                (*particles)[i].addForce(p->getForce());
                (*particles)[i].setVel(p->getVel());
        }

        // Aggiornamento delle posizioni delle particelle
        for (size_t i = 0; i < numParticles; ++i) { 
             (*particles)[i].updateAndReset(delta_t);
        }

        for (size_t i = 0; i < numParticles; ++i) {
            Particle<Dimension>& q = (*particles)[i];
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

        quadtree.reset();
    }
    
    end = time(NULL);
    std::cout << "Time taken by generateRandomParticles function: " << end - start << " seconds" << std::endl;
    file.close();

}


template<size_t Dimension>
void parallelBarnesHut(int iterationNumber, std::vector<Particle<Dimension>>* particles, int dimSimulationArea, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f, int speedUp){

    particles->size() < omp_get_max_threads()? omp_set_num_threads(particles->size()):omp_set_num_threads(omp_get_max_threads());  
    int num_threads = omp_get_max_threads();
    time_t start, end, start1, end1;

    double theta = 0.5;
    std::unique_ptr<QuadtreeNode<Dimension>> quadtree;

    double totalTreeTime = 0.0;

    if(num_threads % 4 != 0){
        omp_set_num_threads(num_threads - num_threads%4);
    }

    // Inizio dell'algoritmo di Barnes-Hut
    
    start1 = time(NULL);
    iterationNumber = 1;
    for (int iter = 0; iter < iterationNumber; ++iter) {
        //std::cout<<"fanculo"<<std::endl;
        start = time(NULL);;
        int a;
        a = generateTreeParallel<Dimension>(particles, 2*dimSimulationArea);
        //quadtree = createQuadTreeParallel(particles, 2*dimSimulationArea, num_threads);
        //end = time(NULL);
        //totalTreeTime = totalTreeTime + end - start;
        //#pragma omp parallel shared(particles, f, iter, delta_t, file, softening, dimSimulationArea, num_threads, quadtree)
        //{
//
        //    // Calcolo delle forze per ogni particella
        //    #pragma omp for schedule(static, particles.size()/omp_get_num_threads())
        //    for (size_t i = 0; i < numParticles; ++i) {
//
        //        //calculateNetForce(treeRoot, &particles[i], theta, *f);
        //        std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particles[i]);
        //        
        //        calculateNetForceQuadtree(quadtree, p, theta, *f, dimSimulationArea, softening);
        //        if(p != nullptr)
        //            particles[i].addForce(p->getForce());
        //            particles[i].setVel(p->getVel());
        //    }
//
        //    // Aggiornamento delle posizioni delle particelle
        //    #pragma omp for schedule(static, particles.size()/omp_get_num_threads())
        //    for (auto& particle : particles) {  
        //        particle.updateAndReset(delta_t);
        //    }
        //}
//
        //for (size_t i = 0; i < numParticles; ++i) {
        //    Particle<Dimension> q = particles[i];
        //    if (file.is_open()) {
        //        file << q.getId() << ",";
        //        const auto& pos = q.getPos();
        //        for (size_t i = 0; i < Dimension; ++i) {
        //            file << pos[i];
        //            if (i < Dimension - 1) {
        //                file << ",";
        //            }
        //        }
        //        file << std::endl;
        //    } else {
        //        std::cout << "Unable to open file";
        //    }
        //}
//
        //quadtree.reset();      

    }
    std::cout << "Time taken by createQuadTree function: " << totalTreeTime << " seconds" << std::endl;
    
    end1 = time(NULL);
    std::cout << "Time taken by Barnes parallel function: " << end1 - start1 << " seconds" << std::endl;
    file.close();
}

/*
To calculate the net force acting on body b, use the following recursive procedure, starting with the root of the quad-tree:

1)If the current node is an external node (and it is not body b), calculate the force exerted by the current node on b, and add this amount to b’s net force.
2)Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal node as a single body, and calculate the force it exerts on body b, and add this amount to b’s net force.
3)Otherwise, run the procedure recursively on each of the current node’s children.

*/
template <size_t Dimension>
void calculateNetForceQuadtree(const std::unique_ptr<QuadtreeNode<Dimension>>& node, std::shared_ptr<Particle<Dimension>> p, double theta, Force<Dimension>& f, double dimSimulationArea, double softening) {
   
    double s, d;
    
    std::array<double, Dimension> force_qk;
    std::unique_ptr<Particle<Dimension>> approxParticle;

    if(node == nullptr) {
        return;
    }


    if(p->hitsBoundary(dimSimulationArea)){
        p->manageCollision(*p, dimSimulationArea);
    }   

    //std::cout << "is leaf: " << node->isLeaf() << std::endl;
    if(node->isLeaf() && node->getParticle() != nullptr && node->getParticle()->getId() != p->getId()){

        if(node->getParticle()->squareDistance(*p) < (((node->getParticle()->getRadius() + p->getRadius())*(node->getParticle()->getRadius() + p->getRadius())))*softening){
            node->getParticle()->manageCollision(*p, 0.0);
        }

        //std::cout << "Node is a leaf and contains a different particle. Applying direct force." << std::endl;
        force_qk = f.calculateForce(*p, *node->getParticle());
        //std::cout << "force_qk: " << force_qk[0]<< " : " <<force_qk[1]<<std::endl;
        p->addForce(force_qk);
    }
    else if(!node->isLeaf()){
        s = node->getWidth();
        approxParticle = node->createApproximateParticle();
        d = std::sqrt(p->squareDistance(*approxParticle));
        if (d == 0) {
            return;
        }
        if(s/d < theta){
            force_qk = f.calculateForce(*p, *approxParticle);
            p->addForce(force_qk);
        }
        else {
            for(auto& child : node->getChildren()) {
                
                calculateNetForceQuadtree(child, p, theta, f, dimSimulationArea, softening);
            }
        }
    }
}

template <size_t Dimension>
double computeRatio(const std::array<double, Dimension>& pos, const std::array<double, Dimension>& centerMass, double width) {
    double sum = 0.0;

    for (size_t i = 0; i < Dimension; ++i) {
        sum += std::pow(pos[i] - centerMass[i], 2);
    }

    return width/std::sqrt(sum);
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
    //particles = generateOrbitTestParticles<Dimension>(dimSimulationArea, 10000);
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
    }else if(simType == 2){
        start = time(NULL);
        serialBarnesHut<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, fileName, file, *f, speedUp);
        end = time(NULL);
        std::cout << "Time taken by serial Barnes Hut: " << end - start << " seconds" << std::endl;
    }else if(simType == 3){
        start = time(NULL);
        parallelBarnesHut<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, fileName, file, *f, speedUp);
        end = time(NULL);
        std::cout << "Time taken by parallel Barnes Hut: " << end - start << " seconds" << std::endl;
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
/*template<size_t Dimension>
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
}*/

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
    int simType = 3;
    int forceType = 0;
    double delta_t = 0.01;
    int dimSimulationArea = 200; 
    int iterationNumber = 1000; 
    int numParticles = 10;
    int mass = 50; 
    int maxVel = 50; 
    int maxRadius = 15; 
    double softening = 0.7;
    int speedUp = 1;
    std::string fileName = "../graphics/coordinates.txt";
    std::string fileName2 = "../graphics/coordinates2.txt";

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
                if(simType != 0 && simType != 1 && simType!= 2 && simType != 3) {
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
        main2DSimulation<2>(forceType, simType, delta_t, dimSimulationArea, iterationNumber, numParticles, mass, maxVel, maxRadius, softening, fileName, speedUp);
    } else if (dim == 3) {
        //main3DSimulation<3>(forceType, simType, delta_t, dimSimulationArea, iterationNumber, numParticles, mass, maxVel, maxRadius, softening, fileName, speedUp);
    }

    return 0;
    
}
