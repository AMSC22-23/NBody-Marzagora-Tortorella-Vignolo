#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "force.hpp"
#include "particle.hpp"
#include "quadtreeNode.hpp"

#include "file_utilities.hpp"
#include "simulation_utilities.hpp"
#include "tree_utilities.hpp"


template <size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> merging(std::vector<std::unique_ptr<QuadtreeNode<Dimension>>>* roots,
                                                 double subRegionsDimension, int num_threads) {
    int size = (*roots).size();  // always 4^n elements
    int k = 4;
    int y = 1;
    double x = 0.0;
    double z = 0.0;
    int l = 0;

    std::array<std::unique_ptr<QuadtreeNode<Dimension>>, 4> temp_roots;

    while (size >= 4) {
        if (size < num_threads) {
            num_threads = size;
            omp_set_num_threads(num_threads);
        }

        #pragma omp parallel for private(temp_roots) shared(k, y) schedule(static, size / num_threads)
        for (int i = 0; i < size / 4; i++) {
            for (int j = 0; j < 4; j++) {
                int index = k * i + j * y;
                swap(temp_roots[j], (*roots)[index]);
            }

            (*roots).at(k * i) = std::move(mergeRoots(&temp_roots, subRegionsDimension * (l + 1)));
        }

        k = k * 4;
        y = y * 4;
        l++;
        size /= 4;  // Reduce the size of the vector by a factor of 4
    }


    return std::move((*roots)[0]);
}

template <size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> mergeRoots(std::array<std::unique_ptr<QuadtreeNode<Dimension>>, 4>* roots,
                                                    double subRegionsDimension) {
    double x, y;

    std::unique_ptr<QuadtreeNode<Dimension>> intNode;
    double totx = 0.0;
    double toty = 0.0;

    for (int i = 0; i < 4; i++) {
        totx += (*roots)[i]->getCenter()[0];
        toty += (*roots)[i]->getCenter()[1];
    }

    x = totx / 4.0;
    y = toty / 4.0;

    intNode = std::make_unique<QuadtreeNode<Dimension>>(x, y, 2 * subRegionsDimension, 20);

    for (int j = 0; j < 4; j++) {
        if ((*roots)[j]->isLeaf()) {
            if ((*roots)[j]->getParticle() != nullptr) {
                intNode->insert((*roots)[j]->getParticle());
            }
        } else {
            intNode->insertNode(std::move((*roots)[j]));
        }
    }

    // if(intNode != nullptr)
    //     intNode->printTree();

    return intNode;
}

template <size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> generateTreeParallel(std::vector<Particle<Dimension>>* particles,
                                                              double dimSimulationArea, int num_threads) {

    int num_quad = pow(4, ceil(log(num_threads) / log(4)));
    std::unique_ptr<QuadtreeNode<Dimension>> treeRoot;

    std::vector<std::array<double, 2>> regions = getSubRegionsCoordinates<Dimension>(num_quad, dimSimulationArea);
    std::array<double, 2> subRegionsDimension = getSubRegionsDimension<Dimension>(num_quad, dimSimulationArea);

    std::vector<std::unique_ptr<QuadtreeNode<Dimension>>> roots;
    roots.reserve(num_quad);

    for (int i = 0; i < num_quad; i++) {
        roots.push_back(nullptr);
    }

    #pragma omp parallel for shared(particles, regions)
    for (int i = 0; i < num_quad; i++) {
        (roots).at(i) =
            std::move(createQuadTreeParallel(*particles, regions[i], subRegionsDimension[0], dimSimulationArea));
    }

    treeRoot = merging(&roots, subRegionsDimension[0], num_threads);

    omp_set_num_threads(num_threads);

    return treeRoot;
}

template <size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> createQuadTreeParallel(std::vector<Particle<Dimension>>& particles,
                                                                std::array<double, Dimension> center, double width,
                                                                double dimSimulationArea) {
    if (particles.empty()) {
        return nullptr;
    }

    // Determine the simulation boundaries (this might need to be adjusted)
    double minX = center[0] - width * 0.5;
    double maxX = center[0] + width * 0.5;
    double minY = center[1] - width * 0.5;
    double maxY = center[1] + width * 0.5;

    // Create the root of the tree
    std::unique_ptr<QuadtreeNode<Dimension>> root =
        std::make_unique<QuadtreeNode<Dimension>>(center[0], center[1], width, 20);

    // Insert each particle into the tree
    if (maxX == dimSimulationArea / 2 && minY == -dimSimulationArea / 2) {
        for (auto& particle : particles) {
            if (particle.getPos()[0] >= minX && particle.getPos()[0] <= maxX && particle.getPos()[1] <= maxY &&
                particle.getPos()[1] >= minY) {
                std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
                root->insert(p);
            }
        }

    } else if (minY == -dimSimulationArea / 2) {
        for (auto& particle : particles) {
            if (particle.getPos()[0] >= minX && particle.getPos()[0] < maxX && particle.getPos()[1] <= maxY &&
                particle.getPos()[1] >= minY) {
                std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
                root->insert(p);
            }
        }

    } else if (maxX == dimSimulationArea / 2) {
        for (auto& particle : particles) {
            if (particle.getPos()[0] >= minX && particle.getPos()[0] <= maxX && particle.getPos()[1] <= maxY &&
                particle.getPos()[1] > minY) {
                std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
                root->insert(p);
            }
        }
    } else {  // standard case
        for (auto& particle : particles) {
            if (particle.getPos()[0] >= minX && particle.getPos()[0] < maxX && particle.getPos()[1] <= maxY &&
                particle.getPos()[1] > minY) {
                std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
                root->insert(p);
            }
        }
    }

    return root;
}

/*
To calculate the net force acting on body b, use the following recursive
procedure, starting with the root of the quad-tree:

1)If the current node is an external node (and it is not body b), calculate the
force exerted by the current node on b, and add this amount to b’s net force.
2)Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal node as a
single body, and calculate the force it exerts on body b, and add this amount to
b’s net force. 3)Otherwise, run the procedure recursively on each of the current
node’s children.
*/
template <size_t Dimension>
void calculateNetForceQuadtree(const std::unique_ptr<QuadtreeNode<Dimension>>& node,
                               std::shared_ptr<Particle<Dimension>> p, double theta, Force<Dimension>& f,
                               double dimSimulationArea, double softening) {
    double s, d;

    std::array<double, Dimension> force_qk;
    std::unique_ptr<Particle<Dimension>> approxParticle;

    if (node == nullptr) {
        return;
    }

    if (p->hitsBoundary(dimSimulationArea)) {
        p->manageCollision(*p, dimSimulationArea);
    }

    if (node->isLeaf() && node->getParticle() != nullptr && node->getParticle()->getId() != p->getId()) {
        if (node->getParticle()->squareDistance(*p) < (((node->getParticle()->getRadius() + p->getRadius()) *
                                                        (node->getParticle()->getRadius() + p->getRadius()))) *
                                                          softening) {
            node->getParticle()->manageCollision(*p, 0.0);
        }

        force_qk = f.calculateForce(*p, *node->getParticle());
        p->addForce(force_qk);
    } else if (!node->isLeaf()) {
        s = node->getWidth();
        approxParticle = node->createApproximateParticle();
        d = std::sqrt(p->squareDistance(*approxParticle));
        if (d == 0) {
            return;
        }
        if (s / d < theta) {
            force_qk = f.calculateForce(*p, *approxParticle);
            p->addForce(force_qk);
        } else {
            for (auto& child : node->getChildren()) {
                calculateNetForceQuadtree(child, p, theta, f, dimSimulationArea, softening);
            }
        }
    }
}

template <size_t Dimension>
double computeRatio(const std::array<double, Dimension>& pos, const std::array<double, Dimension>& centerMass,
                    double width) {
    double sum = 0.0;

    for (size_t i = 0; i < Dimension; ++i) {
        sum += std::pow(pos[i] - centerMass[i], 2);
    }

    return width / std::sqrt(sum);
}


/**
 serialBarnesHut
 **/
template <size_t Dimension>
void serialBarnesHut(int iterationNumber, std::vector<Particle<Dimension>>* particles, int dimSimulationArea,
                     double softening, double delta_t, Force<Dimension>& f,int speedUp) {

    double theta = 0.5;
    std::unique_ptr<QuadtreeNode<Dimension>> quadtree;
    int numParticles = (*particles).size();
    std::FILE* file;
    std::string fileName = "../graphics/Coordinates_0.txt";
    const char* re = fileName.c_str();
    file = fopen(re, "w");

    if (file == NULL) {
        std::cout << "Error opening file!" << std::endl;
        exit(1);
    }

    // Inizio dell'algoritmo di Barnes-Hut

    for (int iter = 0; iter < iterationNumber; ++iter) {
        quadtree = createQuadTree(*particles, 2 * dimSimulationArea);

        // Calcolo delle forze per ogni particella
        for (int i = 0; i < numParticles; ++i) {
            // calculateNetForce(treeRoot, &particles[i], theta, *f);
            std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>((*particles)[i]);

            calculateNetForceQuadtree(quadtree, p, theta, f, dimSimulationArea, softening);
            if (p != nullptr) (*particles)[i].addForce(p->getForce());
            (*particles)[i].setVel(p->getVel());
        }

        // Aggiornamento delle posizioni delle particelle
        for (size_t i = 0; i < numParticles; ++i) {
            //(*particles)[i].updateAndReset(delta_t);
            (*particles)[i].velocityVerletUpdate(delta_t);
            (*particles)[i].resetForce();
        }

        writeParticlePositionsToFile(*particles, file);

        quadtree.reset();
    }
    fclose(file);
}


/**
 parallelBarnesHut
 **/
template <size_t Dimension>
void parallelBarnesHut(int iterationNumber, std::vector<Particle<Dimension>>* particles, int dimSimulationArea,
                       double softening, double delta_t, Force<Dimension>& f,int speedUp, size_t numFilesAndThreads) {
    int num_threads = numFilesAndThreads;
    time_t start, end;
    double theta = 0.5;
    double totalTreeTime = 0.0;
    std::vector<std::vector<std::array<double, Dimension>>> local_forces(omp_get_max_threads(), std::vector<std::array<double, Dimension>>(particles->size()));
    std::vector<std::ofstream> coordinateFiles(numFilesAndThreads);
    std::array<double,Dimension> force;
    std::ofstream file;
    std::vector<std::FILE*> files(numFilesAndThreads);
    //initialize all files
    for (int i = 0; i < numFilesAndThreads; ++i) {
        std::string fileName = "../graphics/Coordinates_" + std::to_string(i) + ".txt";
        const char* re = fileName.c_str();
        files[i] = fopen(re, "w");
    }
    int numParticles = (*particles).size();
    int idThread;

    //check if all files are open
    for (int i = 0; i < numFilesAndThreads; ++i) {
        if (files[i] == NULL) {
            std::cout << "Error opening file!" << std::endl;
            exit(1);
        }
    }

    //std::shared_ptr<Particle<Dimension>> p;
    std::unique_ptr<QuadtreeNode<Dimension>> quadtree;

    
    
    for (int iter = 0; iter < iterationNumber; ++iter) {
        start = time(NULL);
        quadtree = generateTreeParallel<Dimension>(particles, 2 * dimSimulationArea, num_threads);
        end = time(NULL);
        totalTreeTime = totalTreeTime + end - start;
        #pragma omp parallel shared(numParticles, particles, f, delta_t, file, softening, dimSimulationArea, num_threads, quadtree), private(idThread)
        {
            idThread = omp_get_thread_num();
            // Calcolo delle forze per ogni particella
            #pragma omp for schedule(static, numParticles / num_threads)
            for (int i = 0; i < numParticles; ++i) {
                std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>((*particles)[i]);
                calculateNetForceQuadtree(quadtree, p, theta, f, dimSimulationArea, softening);
                if (p != nullptr) {
                    (*particles)[i].addForce(p->getForce());
                    (*particles)[i].setVel(p->getVel());
                }
            }
            // Aggiornamento delle posizioni delle particelle
            #pragma omp for schedule(static, numParticles / num_threads)
            for (auto& particle : (*particles)) {
                //particle.updateAndReset(delta_t);
                particle.velocityVerletUpdate(delta_t);
                particle.resetForce();
            }
    
            if(iter%speedUp==0){
                #pragma omp for schedule(static, numParticles/num_threads)
                for (int i = 0; i < numParticles; i++) {
                    fprintf(files[idThread], "%d,", (*particles)[i].getId());
                    const auto& pos = (*particles)[i].getPos();
                    for (size_t i = 0; i < Dimension; ++i) {
                        fprintf(files[idThread], "%f", pos[i]);
                        if (i < Dimension - 1) {
                            fprintf(files[idThread], ",");
                        }   
                    }
                    fprintf(files[idThread], "\n");
                }
            }
        }
        quadtree.reset();
    }
    #pragma omp for schedule(static, 1)
    for(int i = 0; i < numFilesAndThreads; ++i){
        fclose(files[i]);
    }
    std::cout << "Time taken by createQuadTree function: " << totalTreeTime << " seconds" << std::endl;
}



/**
 * @brief Template function that executes the NBody simulation serially for a
 *given number of iterations. For each particle,  it checks first if the
 *particle hits the boundary and manages that collision; then, it checks for
 *collision between particles and call the manageCollision function to take care
 *of it. Finally it computed the force between the particles and adds it to a
 *vector. After checking all the particles, it updates the values calculated
 *before and resets the vector; finally, it writes the updated positions of the
 *particles on the file after delta_t time.
 *
 * @tparam Dimension Number of dimensions of the simulation
 * @param it Number of iterations
 * @param particles Reference to the vector of particles
 * @param dim Number of dimensions of the simulation (2D,3D)
 * @param softening Overhead to avoid particles overlapping and fusing together
 * @param delta_t Time step after which the simulation is updated
 * @param fileName Name of tile in which the function writes the coordinates of
 *the particles computed during the simulation
 * @param file Reference to the file in which the coordinates are written
 * @param f Reference to the Force object responsible for calculating particle
 *interactions.
 **/
template<size_t Dimension>
void serialSimulation(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, Force<Dimension>& f, int speedup){
    std::array<double,Dimension> force_qk;
    //std::ofstream file("../graphics/Coordinates_0.txt");
    std::FILE* file;
    std::string fileName = "../graphics/Coordinates_0.txt";
    const char* re = fileName.c_str();
    file = fopen(re, "w");

    if (file == NULL) {
        std::cout << "Error opening file!" << std::endl;
        exit(1);
    }
    
    for (int z = 0; z < it; ++z) {
        // std::cout<<z<<std::endl;
        for (int i = 0; i < (*particles).size(); ++i) {
            Particle<Dimension>& q = (*particles)[i];

            if (q.hitsBoundary(dim)) {
                q.manageCollision(q, dim);
            }

            for (int j = i + 1; j < (*particles).size(); j++) {
                Particle<Dimension>& k = (*particles)[j];

                if (q.squareDistance(k) <
                    (((q.getRadius() + k.getRadius()) * (q.getRadius() + k.getRadius()))) * softening) {
                    q.manageCollision(k, 0.0);
                }

                force_qk = f.calculateForce(q, k);
                q.addForce(force_qk);
                for (size_t i = 0; i < Dimension; ++i) force_qk[i] = -force_qk[i];
                k.addForce(force_qk);
            }
            z == it - 1 ? q.update(delta_t) : q.updateAndReset(delta_t);
        }

        if (z % speedup == 0) {
            writeParticlePositionsToFile(*particles, file);
        }
    }
    fclose(file);
}

/**
 * @brief Template function that executes the NBody simulation in parallel for a
 *given number of iterations, using OpenMP directives to parallelize the
 *simulation loop, calculating forces and updating particle positions
 *concurrently. Firstly, the function initializes the number of threads for the
 *parallel section based on the minimum between the size of the particle vector
 *or the maximum available threads. Then it starts the simulation, assigning
 *blocks to threads dynamically: collision among particles or between a particle
 *and the boundary and computation of the forces are computed concurrently and
 *then summed up. The forces are then updated and then the local vector is
 *reset; finally, the function periodically writes the updated positions inside
 *a file.
 *
 * @tparam Dimension Number of dimensions of the simulation
 * @param it Number of iterations
 * @param particles Reference to the vector of particles
 * @param dim Number of dimensions of the simulation (2D,3D)
 * @param softening Overhead to avoid particles overlapping and fusing together
 * @param delta_t Time step after which the simulation is updated
 * @param fileName Name of tile in which the function writes the coordinates of
 *the particles computed during the simulation
 * @param file Reference to the file in which the coordinates are written
 * @param f Reference to the Force object responsible for calculating particle
 *interactions.
 * @param u
 **/

template<size_t Dimension>
void parallelSimulation(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, Force<Dimension>& f, int u, size_t numFilesAndThreads){

    int id_thread;
    int num_particles = (*particles).size();
    std::array<double,Dimension> force_qk;
    bool filesOpen = true;
    int num_threads = numFilesAndThreads;
    std::vector<std::vector<std::array<double, Dimension>>> local_forces(omp_get_max_threads(), std::vector<std::array<double, Dimension>>(particles->size()));
    std::vector<std::ofstream> coordinateFiles(numFilesAndThreads);
    std::array<double,Dimension> force;
    std::ofstream file;
    std::vector<std::FILE*> files(numFilesAndThreads);
    //initialize all files
    for (int i = 0; i < numFilesAndThreads; ++i) {
        std::string fileName = "../graphics/Coordinates_" + std::to_string(i) + ".txt";
        const char* re = fileName.c_str();
        files[i] = fopen(re, "w");
    }

    //check if all files are open
    for (int i = 0; i < numFilesAndThreads; ++i) {
        if (files[i] == NULL) {
            std::cout << "Error opening file!" << std::endl;
            exit(1);
        }
    }

    Particle<Dimension> *q;
    Particle<Dimension> *k;

        #pragma omp parallel shared(particles, f, it, delta_t, coordinateFiles, softening, dim, num_threads, files, num_particles) private(force, id_thread, file, k, q)
    {
        id_thread = omp_get_thread_num();
    {
        for (int z = 0; z < it; ++z){
               #pragma omp for schedule(dynamic, u)
               for (int i = 0; i < num_particles; ++i) {
                    k = &(*particles)[i];


                   if(k->hitsBoundary(dim)){
                       k->manageCollision(*k, dim);
                   }

                   for(int j = i+1; j < num_particles; ++j){
                       q = &(*particles)[j];

                       if(q->squareDistance(*k) < (((q->getRadius() + k->getRadius())*(q->getRadius() + k->getRadius())))*softening){
                           q->manageCollision(*k, 0.0);
                       }

                       force = f.calculateForce(*k, *q);
                       for (int y = 0; y < Dimension; ++y) {
                           local_forces[id_thread][i][y] += force[y];
                           local_forces[id_thread][j][y] -= force[y];
                       }
                   }
               }

               #pragma omp for schedule(static, particles->size()/omp_get_num_threads())
               for (int i = 0; i < num_particles; ++i) {
                   q = &(*particles)[i];
                   q->resetForce();
                   for (int j = 0; j < num_threads; ++j) {
                       q->addForce(local_forces[j][i]);
                   }
               }

               #pragma omp for schedule(static, 1)
               for(int k = 0; k < num_threads; ++k){
                   for (int j = 0; j < num_particles; ++j) {
                       for (int y = 0; y < Dimension; ++y) {
                           local_forces[k][j][y] = 0.0;
                       }
                   }
               }

               #pragma omp for schedule(static, particles->size()/num_threads)
               for (int i = 0; i < num_particles; ++i) {
                   q = &(*particles)[i];
                   q->update(delta_t);
                }

                if(z%u==0){
                    #pragma omp for schedule(static, particles->size()/num_threads)
                    for (int i = 0; i < num_particles; i++) {
                        q = &(*particles)[i];
                        fprintf(files[id_thread], "%d,", q->getId());
                        const auto& pos = q->getPos();
                        for (size_t i = 0; i < Dimension; ++i) {
                            fprintf(files[id_thread], "%f", pos[i]);
                            if (i < Dimension - 1) {
                                fprintf(files[id_thread], ",");
                            }   
                        }
                        fprintf(files[id_thread], "\n");
                    }
                }    
        }
        #pragma omp for schedule(static, 1)
        for(int i = 0; i < numFilesAndThreads; ++i){
            fclose(files[i]);
        }
    }
    }
}

/**
 * @brief Template function that wraps main function for the 2D simulation:
 *calls the function that generates the particles, the one that prints the
 *initial states on the file and then calls the function which starts the
 *simulation chosen by the user.
 *
 * @tparam Dimension Number of dimensions of the simulation
 * @param simType Simulation type: 0 for serial, 1 for parallel
 * @param forceType Type of the force: g for gravitational force, c for coulomb
 *force
 *
 **/
template <size_t Dimension>
void main2DSimulation(int forceType, int simType, double delta_t, int dimSimulationArea, int iterationNumber,
                      int numParticles, int mass, int maxVel, int maxRadius, double softening, std::string fileName,
                      int speedUp) {
    time_t start, end;
    std::vector<Particle<Dimension>> particles;
    size_t numFilesAndThreads;
    Force<Dimension>* f;
    if (forceType == 1)
        f = new CoulombForce<Dimension>();
    else
        f = new GravitationalForce<Dimension>();

    start = time(NULL);
    particles = generateRandomParticles<Dimension>(numParticles, dimSimulationArea, (forceType) ? -mass : 1, mass,
                                                   maxVel, 1, maxRadius, forceType);

    end = time(NULL);
    std::cout << "Time taken by generateRandomParticles function: " << end - start << " seconds" << std::endl;

    if(!simType) numFilesAndThreads = 1;
    else{
        particles.size() < omp_get_max_threads()? omp_set_num_threads(particles.size()):omp_set_num_threads(omp_get_max_threads());  
        numFilesAndThreads = omp_get_max_threads();
    }
    printInitialStateOnFile(&particles, dimSimulationArea, fileName, iterationNumber, speedUp, numFilesAndThreads);

    if (simType == 0) {
        start = time(NULL);
        serialSimulation<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, *f, speedUp);
        end = time(NULL);
        std::cout << "Time taken by serial simulation: " << end - start << " seconds" << std::endl;
    } else if (simType == 1) {
        start = time(NULL);
        parallelSimulation<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, *f, speedUp, numFilesAndThreads);
        end = time(NULL);
        std::cout << "Time taken by parallel simulation: " << end - start << " seconds" << std::endl;
    } else if (simType == 2) {
        start = time(NULL);
        serialBarnesHut<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, *f, speedUp);
        end = time(NULL);
        std::cout << "Time taken by serial Barnes Hut: " << end - start << " seconds" << std::endl;
    } else if (simType == 3) {
        start = time(NULL);
        parallelBarnesHut<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, *f, speedUp, numFilesAndThreads);
        end = time(NULL);
        std::cout << "Time taken by parallel Barnes Hut: " << end - start << " seconds" << std::endl;
    }

    //file.close();
}

/**
 * @brief Template function that wraps main function for the 3D simulation:
 *calls the function that generates the particles, the one that prints the
 *initial states on the file and then calls the function which starts the
 *simulation chosen by the user.
 *
 * @tparam Dimension Number of dimensions of the simulation
 * @param simType Simulation type: 0 for serial, 1 for parallel
 * @param forceType Type of the force: g for gravitational force, c for coulomb
 *force
 *
 **/
template <size_t Dimension>
void main3DSimulation(int forceType, int symType, double delta_t, int dimSimulationArea, int iterationNumber,
                      int numParticles, int mass, int maxVel, int maxRadius, int softening, std::string fileName,
                      int speedUp) {
    time_t start, end;
    std::vector<Particle<Dimension>> particles;
    Force<Dimension>* f;
    size_t numFilesAndThreads;
    if (forceType == 1)
        f = new CoulombForce<Dimension>();
    else
        f = new GravitationalForce<Dimension>();

    start = time(NULL);
    particles = generateRandomParticles<Dimension>(numParticles, dimSimulationArea, (forceType) ? -mass : 1, mass, maxVel, 1, maxRadius, forceType);
    end = time(NULL);
    std::cout << "Time taken by generateRandomParticles function: " << end - start << " seconds" << std::endl;

    if(!symType) numFilesAndThreads = 1;
    else{
        particles.size() < omp_get_max_threads()? omp_set_num_threads(particles.size()):omp_set_num_threads(omp_get_max_threads());  
        numFilesAndThreads = omp_get_max_threads();
    }
    printInitialStateOnFile(&particles, dimSimulationArea, fileName, iterationNumber, speedUp, numFilesAndThreads);

    if (symType == 0) {
        start = time(NULL);
        serialSimulation<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, *f, speedUp);
        end = time(NULL);
        std::cout << "Time taken by serial simulation: " << end - start << "seconds" << std::endl;

    } else if (symType == 1) {
        start = time(NULL);
        parallelSimulation<Dimension>(iterationNumber, &particles, dimSimulationArea, softening, delta_t, *f, speedUp, numFilesAndThreads);
        end = time(NULL);
        std::cout << "Time taken by parallel simulation: " << end - start << " seconds" << std::endl;
    }

    //file.close();
}

void showHelp() {
    std::cout << "Change the following parameters if you don't want to run the default simualtion: "<< std::endl;
    std::cout << "      -h : prints helper " << std::endl;
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
 * @brief Main function which asks the users for the simulation type, the
 * dimensions and the force type and then runs the simulation accordingly
 */

int main(int argc, char** argv) {
    int dim = 2;
    int simType = 3;
    int forceType = 0;
    double delta_t = 0.1;
    int dimSimulationArea = 10000;
    int iterationNumber = 1000;
    int numParticles = 1000;
    int mass = 50;
    int maxVel = 10;
    int maxRadius = 15;
    double softening = 0.7;
    int speedUp = 10;
    std::string fileName = "../graphics/Info.txt";

    if (argc < 2) {
        char a;
        std::cout << "Enter 'd' to run the default simulation: " << std::endl;
        std::cin >> a;
        if (a == 'd')
            dim = 2;
        else {
            showHelp();
            return 1;
        }
    }

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-dim") == 0) {
            if (++i < argc) {
                dim = atoi(argv[i]);
                if (dim != 2 && dim != 3) {
                    std::cout << "No feasible dimension." << std::endl;
                    return 1;
                }
            } else {
                std::cout << " Error: flag -dim requires values 2 or 3 to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-simT") == 0) {
            if (++i < argc) {
                simType = atoi(argv[i]);
                if (simType != 0 && simType != 1 && simType != 2 && simType != 3) {
                    std::cout << "No feasible simulation." << std::endl;
                    return 1;
                }
            } else {
                std::cout << " Error: flag -simT requires values 0 for serial "
                             "or 1 for parallel to work. "
                          << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-force") == 0) {
            if (++i < argc) {
                forceType = atoi(argv[i]);
                if (forceType != 0 && forceType != 1) {
                    std::cout << "No feasible force." << std::endl;
                    return 1;
                }
            } else {
                std::cout << " Error: flag -force requires values 0 for "
                             "gravitational force or 1 for Coulomb force to work. "
                          << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-delta") == 0) {
            if (++i < argc) {
                double delta = std::__cxx11::stof(argv[i]);
                if (delta < 0) {
                    std::cout << "No feasible delta t." << std::endl;
                    return 1;
                }
                delta_t = delta;
            } else {
                std::cout << " Error: flag -delta requires a positive value to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-simA") == 0) {
            if (++i < argc) {
                int simArea = atoi(argv[i]);
                if (simArea < 0) {
                    std::cout << "No feasible simulation area." << std::endl;
                    return 1;
                }
                dimSimulationArea = simArea;
            } else {
                std::cout << " Error: flag -simA requires a positive value to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-it") == 0) {
            if (++i < argc) {
                int it = atoi(argv[i]);
                if (it < 0) {
                    std::cout << "No feasible number of iterations." << std::endl;
                    return 1;
                }
                iterationNumber = it;
            } else {
                std::cout << " Error: flag -it requires a positive value to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-numP") == 0) {
            if (++i < argc) {
                int numP = atoi(argv[i]);
                if (numP < 0) {
                    std::cout << "No feasible number of particles." << std::endl;
                    return 1;
                }
                numParticles = numP;
            } else {
                std::cout << " Error: flag -numP requires a positive value to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-maxPr") == 0) {
            if (++i < argc) {
                int maxPr = atoi(argv[i]);
                if (maxPr < 0) {
                    std::cout << "No feasible value of maximum property." << std::endl;
                    return 1;
                }
                mass = maxPr;
            } else {
                std::cout << " Error: flag -maxPr requires a positive value to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-maxVel") == 0) {
            if (++i < argc) {
                int maxV = atoi(argv[i]);
                if (maxV < 0) {
                    std::cout << "No feasible value of radius of the particles." << std::endl;
                    return 1;
                }
                maxVel = maxV;
            } else {
                std::cout << " Error: flag -maxVEl requires a positive value "
                             "to work. "
                          << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-maxR") == 0) {
            if (++i < argc) {
                int maxR = atoi(argv[i]);
                if (maxR < 0) {
                    std::cout << "No feasible value of radius of the particles." << std::endl;
                    return 1;
                }
                maxRadius = maxR;
            } else {
                std::cout << " Error: flag -maxR requires a positive value to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-soft") == 0) {
            if (++i < argc) {
                double soft = std::__cxx11::stof(argv[i]);
                if (soft < 0) {
                    std::cout << "No feasible value of softening." << std::endl;
                    return 1;
                }
                softening = soft;
            } else {
                std::cout << " Error: flag -soft requires a positive value to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-spUp") == 0) {
            if (++i < argc) {
                int spUp = atoi(argv[i]);
                if (spUp < 0) {
                    std::cout << "No feasible value of speedup." << std::endl;
                    return 1;
                }
                speedUp = spUp;
            } else {
                std::cout << " Error: flag -spUp requires a positive value to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-file") == 0) {
            if (++i < argc) {
                std::string file = argv[i];
                fileName = file;
            } else {
                std::cout << " Error: flag -file requires a valid file name to work. " << std::endl;
                return 1;
            }
        }

        if (strcmp(argv[i], "-h") == 0) {
            showHelp();
            return 0;
        }
    }

    if (dim == 2) {
        main2DSimulation<2>(forceType, simType, delta_t, dimSimulationArea, iterationNumber, numParticles, mass, maxVel,
                            maxRadius, softening, fileName, speedUp);
    } else if (dim == 3) {
        // main3DSimulation<3>(forceType, simType, delta_t, dimSimulationArea,
        // iterationNumber, numParticles, mass, maxVel, maxRadius, softening,
        // fileName, speedUp);
    }

    return 0;
}
