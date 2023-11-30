// A class for a particle in an n-body simulation
#include "force.hpp"
#include "particle.hpp"
#include <cmath>

double Particle::getProperty() const {
    return property;
}

template<size_t Dimension>
std::array<double, Dimension> Particle<Dimension>::getPos() const{
    return pos;
}

template<size_t Dimension>
std::array<double, Dimension> Particle<Dimension>::getVel() const{
    return vel;
}

template<size_t Dimension>
std::array<double, Dimension> Particle<Dimension>::getForce() const{
    return force;
}

template<size_t Dimension>
int Particle<Dimension>::getId() const{
    return id;
}

template<size_t Dimension>
void Particle<Dimension>::resetForce() {
    for (size_t i = 0; i < Dimension; ++i ) {
        force[i] = 0.0; 
    }
    //force[0] = 0.0;
    //force[1] = 0.0;
}

template<size_t Dimension>
double Particle<Dimension>::square_distance(const Particle<Dimension> &p) const{
    const auto& k_pos = getPos();
    const auto& p_pos = p.getPos();

    //double x_diff = k_pos[0] - p_pos[0];
    //double y_diff = k_pos[1] - p_pos[1];
    std::array<double, Dimension> diff;
    for(size_t i=0; i < Dimension; ++i) diff[i] = k_pos[i] - p_pos[i];
    //double square_dist = x_diff * x_diff + y_diff * y_diff;
    
    double square_dist = 0.0;
    for(const auto& d : diff) square_dist += d * d;
    
    return square_dist;
}

template<size_t Dimension>
void Particle<Dimension>::addForce(Particle<Dimension> &k, const Force<Dimension>& f) {
    //std::array<double,2> force_qk;
    std::array<double,Dimension> force_qk;
    force_qk = f.calculateForce(k, *this);

    //force[0] += force_qk[0];
    //force[1] += force_qk[1];
    // k.force[0] -= force_qk[0];
    // k.force[1] -= force_qk[1];

    for(size_t i = 0; i < Dimension; ++i){
        force[i] += force_qk[i];
        k.force[i] -= force_qk[i];
    }
}

template<size_t Dimension>
void Particle<Dimension>::update(double delta_t) {
    // Euler integration for position and velocity update
    
    //pos[0] += vel[0] * delta_t;
    //pos[1] += vel[1] * delta_t;
    
    vel[0] += (force[0] / property) * delta_t;
    vel[1] += (force[1] / property) * delta_t;
}

void Particle::update_and_reset(const double delta_t) {
    // Euler integration for position and velocity update
    //pos[0] += vel[0] * delta_t;
    //pos[1] += vel[1] * delta_t;
    
    vel[0] += (force[0] / property) * delta_t;
    vel[1] += (force[1] / property) * delta_t;
    resetForce();
}

template<size_t Dimension>
void Particle<Dimension>::print_states() const{
    std::cout << "Id: " << id << std::endl;
    std::cout << "Position: " << pos[0] << " " << pos[1] << std::endl;
    (!type? std::cout << "Mass: " : std::cout << "Charge: ");
    std::cout << property << std::endl;
    std::cout << "Force: " << force[0] << " " << force[1] << std::endl;
    std::cout << "Velocity: " << vel[0] << " " << vel[1] << std::endl;
}

//method that returns the radius of the particle
template<size_t Dimension>
double Particle<Dimension>::getRadius() const{
    return radius;
}

//not to sure about this, maybe it's better to take it as input
//depends on the interactions (ask!!!)
bool Particle::getType() const{
    return type;
}

template<size_t Dimension> 
void Particle<Dimension>::setVel(std::array<double, Dimension> v){ //original method took as arguments (double vx, double vy)
    //vel[0] = vx;
    //vel[1] = vy;
    for(size_t i=0; i < Dimension; ++i) vel[i] = v[i];
}

void Particle::setProperty(double p){
    property = p;
}

template<size_t Dimension> 
void Particle<Dimension>::manage_collision(Particle<Dimension> &p, double dim){
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
            /*double xvel_prev, yvel_prev, xvel_new, yvel_new;
            xvel_prev = vel[0];
            yvel_prev = vel[1];
            vel[0] = ((property - p.getProperty())*vel[0]+2*p.getProperty()*p.getVel()[0]) / (property + p.getProperty());
            vel[1] = ((property - p.getProperty())*vel[1]+2*p.getProperty()*p.getVel()[1]) / (property + p.getProperty());
            xvel_new = ((p.getProperty() - property)*p.getVel()[0]+2*property*xvel_prev) / (property + p.getProperty());
            yvel_new = ((p.getProperty() - property)*p.getVel()[1]+2*property*yvel_prev) / (property + p.getProperty());
            p.setVel(xvel_new, yvel_new);
        /*} else {
            //inelastic collision -- da rivedere 
            double xvel_prev, yvel_prev, xvel_new, yvel_new;
            vel[0] = (property * vel[0] + p.getProperty() * p.getVel()[0]) / (property + p.getProperty());
            vel[1] = (property * vel[1] + p.getProperty() * p.getVel()[1]) / (property + p.getProperty());
            property = property + p.getProperty();
            radius = radius + p.getRadius(); 
            p.setproperty(0.0);

           */ 
        //}
    }
}