#include <vector>
#include <fstream>
#include <iostream>
#include <random>
#include "particle.hpp" 
#include "quadtreeNode.hpp"
#include "force.hpp"

/**
 * @brief Template function that generates particles in order to perform the
 *simulation. First it creates randomly (using the random library) the
 *distributions for all the parameters, then it checks if the function struggles
 *to generate non-overlapping particles: if the function has to regenerate the
 *particle's position more than maxRetry times in a row, then the program ends
 *and print an error message. Otherwise, it generates randomly the radius and
 *the position of the particles (inside the distributions generated before).
 *After that it checks that the newly generated particle is not overlapping with
 *any existing circle: it calculates the distance between the particle that is
 *being created and the particle p in particles vector. If the particles do not
 * overlap, the remaining properties of the particles are randomly generated
 *and, lastly, the particle is added to the vector of particles; otherwise,
 *variable counter is increased. Finally, this method returns a vector of
 *particles
 *
 * @tparam Dimension Number of dimensions of the simulation
 * @param N Number of particles to generate.
 * @param posBoundary Position boundary for particle positions (default: 100).
 * @param minProperty Minimum value of property for generated particles
 *(default: 1).
 * @param maxProperty Maximum value of property for generated particles
 *(default: 99).
 * @param maxVel Maximum velocity for particles (default: 100).
 * @param minRadius Minimum radius for particles (default: 0).
 * @param maxRadius Maximum radius for particles (default: 15).
 * @param type Type flag for particles to identify either they are particles of
 *the gravitational or the Coulomb force (default: false, which refers to
 *gravitational force).
 * @return A vector of particles with randomly generated properties.
 **/

template <size_t Dimension>
std::vector<Particle<Dimension>> generateRandomParticles(int N, int posBoundary = 100, int minProperty = 1,
                                                         int maxProperty = 99, int maxVel = 100, int minRadius = 0,
                                                         int maxRadius = 15, bool type = false) {
    int maxRetry = 15, counter = 0;
    bool overlapping = false;
    double x, y, r, property, squareDistance = 0.0;
    std::vector<Particle<Dimension>> particles;
    std::array<double, Dimension> pos, vel;

    std::random_device rd;
    std::mt19937 gen(rd());//@note nice use of random_device
    std::uniform_real_distribution<double> disProperty(minProperty, maxProperty);
    std::uniform_real_distribution<double> disVel(-maxVel, maxVel);
    std::uniform_real_distribution<double> disRadius(minRadius, maxRadius);
    std::uniform_real_distribution<double> disPos(-posBoundary + maxRadius, posBoundary - maxRadius);

    while (particles.size() < N) {
        if (counter == maxRetry)
            throw std::runtime_error(
                "ERROR: the dimension of the simulation area is too little, "
                "please specify a bigger area or generate less particles.");

        overlapping = false;

        r = disRadius(gen);

        for (size_t i = 0; i < Dimension; ++i) pos[i] = disPos(gen);

        for (const Particle<Dimension>& p : particles) {
            squareDistance = 0.0;
            for (size_t i = 0; i < Dimension; ++i)
                squareDistance = squareDistance + ((p.getPos()[i] - pos[i]) * (p.getPos()[i] - pos[i]));
            if (sqrt(squareDistance) < p.getRadius() + r) {
                counter++;
                overlapping = true;
                break;
            }
        }

        if (!overlapping) {
            counter = 0;

            property = disProperty(gen);

            for (size_t i = 0; i < Dimension; ++i) vel[i] = disVel(gen);

            Particle<Dimension> p(particles.size(), property, pos, vel, r, type);

            particles.push_back(p);
        } else {
            counter++;
        }
    }

    return particles;
}



/**
 * @brief Template function that prints the particles state
 *
 * @tparam Dimension Number of dimensions of the simulation
 * @param particles Refence to the vector of particles whose states are to be
 * printed
 */
template <size_t Dimension>
void printAllParticlesStateAndDistance(const std::vector<Particle<Dimension>>* particles) {
    std::cout << "Print All Particles State And Distance:\n";
    double squareDistance = 0.0;
    for (int i = 0; i < (*particles).size(); ++i) {
        Particle<Dimension>& p = (*particles)[i];
        std::cout << "--------------------------------------\n";
        p.printStates();
        squareDistance = 0.0;
        for (size_t i = 0; i < Dimension; ++i) squareDistance = squareDistance + (p.getPos()[i] * p.getPos()[i]);
        std::cout << "Distance from origin: " << sqrt(squareDistance) << "\n";
    }
}


/**
 * @brief Template function that generates two particles where one stays still
 *and the other orbits around it. For test purposes only.
 *
 * @tparam Dimension Number of dimensions of the simulation
 * @param size Size of the orbit
 * @param constantForce Constant of the force applied to the particles
 * @return Vector fo particles for a two-body orbit test
 **/
template <size_t Dimension>
std::vector<Particle<Dimension>> generateOrbitTestParticles(double size, double costantForce) {
    std::vector<Particle<Dimension>> particles;
    double orbitRadius = size / 2.0;
    double mass1 = 100;
    double mass2 = 1;
    double force = costantForce * (mass1 * mass2) / (orbitRadius * orbitRadius);
    double velocity = sqrt(force * orbitRadius / mass2);
    std::array<double, Dimension> pos1;
    std::array<double, Dimension> pos2;
    std::array<double, Dimension> vel1;
    std::array<double, Dimension> vel2;

    for (int i = 0; i < Dimension; ++i) {
        pos1[i] = 0;
        pos2[i] = i == 0 ? orbitRadius : 0;
        vel1[i] = 0;
        vel2[i] = i == (Dimension - 1) ? velocity : 0;
    }

    Particle<Dimension> p1(0, mass1, pos1, vel1, size / 100, false);
    Particle<Dimension> p2(1, mass2, pos2, vel2, size / 100, false);
    particles.push_back(p1);
    particles.push_back(p2);

    return particles;
}