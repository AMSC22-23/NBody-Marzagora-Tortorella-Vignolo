// A class for a particle in an n-body simulation
#include "force.hpp"
#include "particle.hpp"
#include <cmath>

double Particle::getMass() const {
    return mass;
}

double Particle::getCharge() const {
    return charge;
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

int Particle::getId() const{
    return id;
}

void Particle::resetForce() {
    //force.resize(2);
    force[0] = 0.0;
    force[1] = 0.0;
}

//void Particle::addForce(double fx, double fy) {
//    force[0] += fx;
//    force[1] += fy;
//}

void Particle::addForce(Particle &k, const Force& f) {
    std::array<double,2> force_qk;
    force_qk = f.calculateForce(k, *this);
    force[0] += force_qk[0];
    force[1] += force_qk[1];
    k.force[0] -= force_qk[0];
    k.force[1] -= force_qk[1];
}

void Particle::update(double delta_t) {
    // Euler integration for position and velocity update
    pos[0] += vel[0] * delta_t;
    pos[1] += vel[1] * delta_t;
    
    vel[0] += (force[0] / mass) * delta_t;
    vel[1] += (force[1] / mass) * delta_t;
}

void Particle::update_and_reset(const double delta_t) {
    // Euler integration for position and velocity update
    pos[0] += vel[0] * delta_t;
    pos[1] += vel[1] * delta_t;
    
    vel[0] += (force[0] / mass) * delta_t;
    vel[1] += (force[1] / mass) * delta_t;
    resetForce();
}

void Particle::print_states() const{
    std::cout << "Id: " << id << std::endl;
    std::cout << "Position: " << pos[0] << " " << pos[1] << std::endl;
    std::cout << "Mass: " << mass << std::endl;
    std::cout << "Force: " << force[0] << " " << force[1] << std::endl;
    std::cout << "Velocity: " << vel[0] << " " << vel[1] << std::endl;
}

//not to sure about this, maybe it's better to take it as input
//depends on the interactions (ask!!!)
bool Particle::returnType() const{
    if (charge != 0.0) return true; 
    else return false;
}
