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
    int count = 1;
    double x, y, squareDistance = 0.0;
    double r;
    std::vector<Particle<Dimension>> particles;
    std::array<double, Dimension> pos;
    std::array<double, Dimension> vel;

    // Create a random number generators using 'random' library
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> disProperty(minProperty, maxProperty);
    std::uniform_real_distribution<double> disVel(-maxVel, maxVel);
    std::uniform_real_distribution<double> disRadius(minRadius, maxRadius);

    for (int i = 0; i < N; i++) {
        count = 0;
        do{ 
            // Check if the function struggles to generate to generate non-overlapping particles: 
            // if the function has to regenerate the particle's position more than 10 times in a row, then the program ends and print an error message
            if (count >= 10){
                std::cout << "ERROR: the dimension of the simulation area is too little, please specify a bigger area or generate less particles." << count << "\n";
                exit(0);
            }
            else{

                // Generate random radius between minRadius and maxRadius
                r = disRadius(gen);
                std::uniform_real_distribution<double> disPos(-posBoundary+r, posBoundary-r);

                // Generate random position between -posBoundary+r and +posBoundary-r
                for(size_t i=0; i<Dimension; ++i)
                    pos[i] = disPos(gen);
                
                // Loop all over the particles previously generated: 
                // in the loop the program check if the newly generated particle overlaps with any of the previously generated particles
                for (const Particle<Dimension> &p : particles) {

                    // Calculate the distance between the particle just created and the particle p
                    squareDistance = 0.0;
                    for(size_t i = 0; i < Dimension; ++i) 
                        squareDistance = squareDistance + ((p.getPos()[i] - pos[i])*(p.getPos()[i] - pos[i]));
                    if (sqrt(squareDistance) < p.getRadius() + r) { 
                        // Particles are overlapped
                        // Update 'count' variable, used to check if the funciton struggles to generate non-overlapping particles
                        count = count + 1;
                        break;
                    }else{ 
                        // Particles are NOT overlapped
                        count = -1;
                    }
                }
            }
        }while (count >= 0 && !particles.empty());

        // Generate random value of property between 1 and 100
        double property = disProperty(gen);

        // Generate random velocity between -maxVel and maxVel
        for(size_t i=0; i<Dimension; ++i)
            vel[i] = disVel(gen);

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
    
    std::array<double,Dimension> force;
    bool go;

    //set in prgama a local variable force
    #pragma omp parallel shared(particles, f, it, delta_t, fileName, file, softening, dim, go) private(force)
    
    {
        //#pragma omp single
        //{
        //    std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
        //    std::cout << "Number of max threads: " << omp_get_max_threads() << std::endl;
        //}
 
        //create a vector of thread that contains a vector of forces of dimension Dimension
        ///std::array<double, Dimension> force;

        for (int z = 0; z < it; ++z){

            //calculate forces and add it to a local array for each thread
            //assign block thread in a cyclic way, not static and i want
            #pragma omp for schedule(dynamic, u) 
                for (int i = 0; i < (*particles).size(); ++i) {
                //if(i==0)
                //    go = false;
                //i==0? std::cout << "INIZIO CALCOLO" << std::endl:std::cout << "" ;
                //for(int y = 0; y < Dimension; ++y) force[y] = 0;
                
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
                //i==(*particles).size()-1? std::cout << "FINE CALCOLO" << std::endl:std::cout << "" ;
                //if(i==(*particles).size()-1)
                //    go = true;
            }

            //while(!go){}



            //sum all the forces calculated by each thread
            #pragma omp for schedule(static, particles->size()/omp_get_num_threads())
            for (int i = 0; i < (*particles).size(); ++i) {
                Particle<Dimension> &q = (*particles)[i];
                q.resetForce();
                for (int j = 0; j < omp_get_num_threads(); ++j) {
                    q.addForce(local_forces[j][i]);
                }
                //i==(*particles).size()/omp_get_num_threads()*2-1? std::cout << "FINE CALCOLO" << std::endl:std::cout << "" ;
            }
            //#pragma omp barrier

            #pragma omp for schedule(static, 1) nowait
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
            #pragma omp single nowait
            {
                //std::cout << "INIZIO SCRITTURA" << std::endl;
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
                //std::cout << "FINITO SCRITTURA " << omp_get_thread_num << std::endl;
            }
            //#pragma omp barrier
        }
    }
}

template<size_t Dimension>
void parallelSimulation_v2(int it, std::vector<Particle<Dimension>>* particles, int dim, double softening, double delta_t, std::string fileName, std::ofstream& file, Force<Dimension>& f, int u){

    particles->size() < omp_get_max_threads()? omp_set_num_threads(particles->size()):omp_set_num_threads(omp_get_max_threads());  
    int num_threads = omp_get_max_threads();

    // Start of the simulation
    std::vector<std::vector<std::array<double, Dimension>>> local_forces(omp_get_max_threads(), std::vector<std::array<double, Dimension>>(particles->size()));
    
    std::array<double,Dimension> force;

    for (int z = 0; z < it; ++z){
    //set in prgama a local variable force
    #pragma omp parallel shared(particles, f, it, delta_t, fileName, file, softening, dim, num_threads) private(force)
    
    {
        //#pragma omp single
        //{
        //    std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
        //    std::cout << "Number of max threads: " << omp_get_max_threads() << std::endl;
        //}
 
        //create a vector of thread that contains a vector of forces of dimension Dimension
        ///std::array<double, Dimension> force;

        ///for (int z = 0; z < it; ++z){

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
                for (int j = 0; j < num_threads; ++j) {
                    q.addForce(local_forces[j][i]);
                }
            }
            //#pragma omp barrier

            #pragma omp for schedule(static, 1)
            for(int k = 0; k < num_threads; ++k){
                for (int j = 0; j < (*particles).size(); ++j) {
                    for (int y = 0; y < Dimension; ++y) {
                        local_forces[k][j][y] = 0.0;
                    }
                }
            }
            //#pragma omp barrier

            //update the position of the particles
            #pragma omp for schedule(static, particles->size()/num_threads)
            for (int i = 0; i < (*particles).size(); ++i) {
                Particle<Dimension> &q = (*particles)[i];
                q.update(delta_t);
            }
            //#pragma omp barrier

            // Write on file the updates after delta_t
            ///#pragma omp single
            ///{
            ///    for (int i = 0; i < (*particles).size(); i++) {
            ///        Particle<Dimension> &q = (*particles)[i];
            ///        if (file.is_open()) {
            ///            file << q.getId() << ",";
            ///            const auto& pos = q.getPos();
            ///            for (size_t i = 0; i < Dimension; ++i) {
            ///                file << pos[i];
            ///                if (i < Dimension - 1) {
            ///                    file << ",";
            ///                }
            ///            }
            ///            file << std::endl;
            ///        } else {
            ///            std::cout << "Unable to open file in parallelSimulation ";
            ///        }
            ///    }
            ///}
            //#pragma omp barrier
        ///}
    }
    if(z%1==0)
    {for (int i = 0; i < (*particles).size(); i++) {
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
    }}
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

int main() {

    //prepocessor instruction that set dimension of the simulation
    #ifdef DIMENSION
        const int d = DIMENSION;
    #else
        const int d = 2;
    #endif

    // Simulation variables
    //const int d = 3; //2D or 3D
    const double delta_t = 1; // In seconds
    const double dim = 10000; // Dimension of the simulation area
    int it = 1000; // Number of iteration
    int n = 1000; // Number of particles
    int mass = 50; // Mass
    int maxVel = 10; // Maximum velocity
    int maxRadius = 10; // Maximum radius of the particles
    double softening = 0.7; // Softening parameter
    time_t start, end; // Time variables
    std::vector<Particle<d>> particles; // Create a vector of particles
    Force<d>* f = new CustomForce<d>(2000); // Create force
    std::string fileName = "../graphics/coordinates.txt"; // File name
    std::ofstream file(fileName); // Open file

    // Generate random particles
    // In order to launch an easy test use: generateTestParticles<d>();
    particles = generateRandomParticles<d>(n, dim, 1, mass, maxVel, 1, maxRadius, false);
    //particles = generateOrbitTestParticles<d>(dim, 2000);

    // Print on file the initial state of the particles
    if (file.is_open()) {

        // Write on file the total number of particles and the size of the area of the simulation
        file << particles.size() << std::endl;
        file << dim << std::endl;
        file << d << std::endl;
        file << it-1 << std::endl;

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
    } else {
        std::cout << "Unable to open file";
        return 0;
    }

    // Start of simulation
    for(int i=1; i<=1; ++i)
    {
        start = time(NULL);
        parallelSimulation<d>(it, &particles, dim, softening, delta_t, fileName, file, *f, i);
        end = time(NULL);
        std::cout << "Time taken by parallel simulation: " << end - start << " seconds" << std::endl;
    }
        for(int i=1; i<=1; ++i)
    {
        start = time(NULL);
        parallelSimulation_v2<d>(it, &particles, dim, softening, delta_t, fileName, file, *f, i);
        end = time(NULL);
        std::cout << "Time taken by parallel simulation v2: " << end - start << " seconds" << std::endl;
    }
    
    //start = time(NULL);
    //serialSimulation<d>(it, &particles, dim, softening, delta_t, fileName, file, *f);
    //end = time(NULL);
    //std::cout << "Time taken by serial simulation: " << end - start << " seconds" << std::endl;

    // In order to print the final state of the particles use: printAllParticlesStateAndDistance(&particles);
    //printAllParticlesStateAndDistance(&particles);

    file.close();
    return 0;
}