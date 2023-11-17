// A class for a particle in an n-body simulation
#include "particle.hpp"

double Particle::getMass() {
    return mass;
}

double* Particle::getPos() {
    return pos;
}

double* Particle::getVel() {
    return vel;
}

double* Particle::getForce() {
    return force;
}

void Particle::resetForce() {
    force[0] = 0.0;
    force[1] = 0.0;
}

void Particle::addForce(double fx, double fy) {
    force[0] += fx;
    force[1] += fy;
}

void Particle::update(double delta_t) {
    // Euler integration for position and velocity update
    vel[0] += (force[0] / mass) * delta_t;
    vel[1] += (force[1] / mass) * delta_t;

    pos[0] += vel[0] * delta_t;
    pos[1] += vel[1] * delta_t;
}
