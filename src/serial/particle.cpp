// A class for a particle in an n-body simulation
#include "particle.hpp"

double Particle::getMass() const {
    return mass;
}

std::array<double, 2> Particle::getPos() const{
    return pos;
}


std::array<double, 2> Particle::getVel() const{
    return vel;
}

std::array<double, 2> Particle::getForce() const{
    return force;
}

void Particle::resetForce() {
    //force.resize(2);
    force[0] = 0.0;
    force[1] = 0.0;
}

void Particle::addForce(double fx, double fy) {
    force[0] += fx;
    force[1] += fy;
}

void Particle::update(double delta_t) {
    // Euler integration for position and velocity update
    pos[0] += vel[0] * delta_t;
    pos[1] += vel[1] * delta_t;
    
    vel[0] += (force[0] / mass) * delta_t;
    vel[1] += (force[1] / mass) * delta_t;
}

void Particle::print_states() const{
    std::cout << "Position: " << pos[0] << " " << pos[1] << std::endl;
    std::cout << "Mass: " << mass << std::endl;
    std::cout << "Force: " << force[0] << " " << force[1] << std::endl;
    std::cout << "Velocity: " << vel[0] << " " << vel[1] << std::endl;
}
