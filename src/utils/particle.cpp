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

double Particle::square_distance(const Particle &p) const{
    const auto& k_pos = getPos();
    const auto& p_pos = p.getPos();

    double x_diff = k_pos[0] - p_pos[0];
    double y_diff = k_pos[1] - p_pos[1];

    double square_dist = x_diff * x_diff + y_diff * y_diff;
    return square_dist;
}

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

//method that returns the radius of the particle
double Particle::getRadius() const{
    return radius;
}

//not to sure about this, maybe it's better to take it as input
//depends on the interactions (ask!!!)
bool Particle::returnType() const{
    if (charge != 0.0) return true; 
    else return false;
}

void Particle::setVel(double vx, double vy){
    vel[0] = vx;
    vel[1] = vy;
}

void Particle::setMass(double m){
    mass = m;
}

void Particle::manage_collision(Particle &p, double dim){
    if(dim){
        //manages collision with the borders of the simulation (visual purposes only)
        if(pos[0] + radius > dim || pos[0] - radius < -dim){
            vel[0] = -vel[0];
        }
        if(pos[1] + radius > dim || pos[1] - radius < -dim){
            vel[1] = -vel[1];
        }
    }
    else{ //manages collision between particles
        //elastic collision
        
        //if(elastic){
            double xvel_prev, yvel_prev, xvel_new, yvel_new;
            xvel_prev = vel[0];
            yvel_prev = vel[1];
            vel[0] = ((mass - p.getMass())*vel[0]+2*p.getMass()*p.getVel()[0]) / (mass + p.getMass());
            vel[1] = ((mass - p.getMass())*vel[1]+2*p.getMass()*p.getVel()[1]) / (mass + p.getMass());
            xvel_new = ((p.getMass() - mass)*p.getVel()[0]+2*mass*xvel_prev) / (mass + p.getMass());
            yvel_new = ((p.getMass() - mass)*p.getVel()[1]+2*mass*yvel_prev) / (mass + p.getMass());
            p.setVel(xvel_new, yvel_new);
        /*} else {
            //inelastic collision -- da rivedere 
            double xvel_prev, yvel_prev, xvel_new, yvel_new;
            vel[0] = (mass * vel[0] + p.getMass() * p.getVel()[0]) / (mass + p.getMass());
            vel[1] = (mass * vel[1] + p.getMass() * p.getVel()[1]) / (mass + p.getMass());
            mass = mass + p.getMass();
            radius = radius + p.getRadius(); 
            p.setMass(0.0);

           */ 
        //}
    }
}