// A class for a particle in an n-body simulation
#include "force.hpp"
#include "particle.hpp"
#include <cmath>

template<size_t Dimension>
double Particle<Dimension>::getMass() const {
    return mass;
}

template<size_t Dimension>
double Particle<Dimension>::getCharge() const {
    return charge;
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
    
    // vel[0] += (force[0] / mass) * delta_t;
    // vel[1] += (force[1] / mass) * delta_t;
    for(size_t i = 0; i < Dimension; ++i){
        pos[i] += vel[i] * delta_t;
        vel[i] += (force[i] / mass) * delta_t;
    }
}

template<size_t Dimension>
void Particle<Dimension>::update_and_reset(const double delta_t) {
    // Euler integration for position and velocity update
    //pos[0] += vel[0] * delta_t;
    //pos[1] += vel[1] * delta_t;
    
    //vel[0] += (force[0] / mass) * delta_t;
    //vel[1] += (force[1] / mass) * delta_t;

    for(size_t i = 0; i < Dimension; ++i){
        pos[i] += vel[i] * delta_t;
        vel[i] += (force[i] / mass) * delta_t;
    }
    resetForce();
}

template<size_t Dimension>
void Particle<Dimension>::print_states() const{
    std::cout << "Id: " << id << std::endl;
    std::cout << "Position: " << pos[0] << " " << pos[1] << std::endl;
    std::cout << "Mass: " << mass << std::endl;
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
template<size_t Dimension>
bool Particle<Dimension>::returnType() const{
    if (charge != 0.0) return true; 
    else return false;
}

template<size_t Dimension> 
void Particle<Dimension>::setVel(std::array<double, Dimension> v){ //original method took as arguments (double vx, double vy)
    //vel[0] = vx;
    //vel[1] = vy;
    for(size_t i=0; i < Dimension; ++i) vel[i] = v[i];
}

template<size_t Dimension> 
void Particle<Dimension>::setMass(double m){
    mass = m;
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
            vel[0] = ((mass - p.getMass())*vel[0]+2*p.getMass()*p.getVel()[0]) / (mass + p.getMass());
            vel[1] = ((mass - p.getMass())*vel[1]+2*p.getMass()*p.getVel()[1]) / (mass + p.getMass());
            xvel_new = ((p.getMass() - mass)*p.getVel()[0]+2*mass*xvel_prev) / (mass + p.getMass());
            yvel_new = ((p.getMass() - mass)*p.getVel()[1]+2*mass*yvel_prev) / (mass + p.getMass());
            p.setVel(xvel_new, yvel_new);*/

            std::array<double, Dimension> prev_vel;
            std::array<double, Dimension> new_vel;

            for(size_t i = 0; i < Dimension; ++i){
                prev_vel[i] = vel[i];
            }

            for(size_t i = 0; i < Dimension; ++i){
                vel[i] = ((mass - p.getMass())*vel[i]+2*p.getMass()*p.getVel()[i]) / (mass + p.getMass());
            }

            for(size_t i = 0; i < Dimension; ++i){
                new_vel[i] = ((p.getMass() - mass)*p.getVel()[i]+2*mass*prev_vel[i]) / (mass + p.getMass());
            }

            p.setVel(new_vel);

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