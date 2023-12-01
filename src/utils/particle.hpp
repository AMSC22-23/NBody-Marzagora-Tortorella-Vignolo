#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <iostream>
#include <array>
//#include "force.hpp"

template<size_t Dimension>
class Force;

template<size_t Dimension>
class Particle {
    public:
        // A constructor that initializes the mass, position, and velocity of the particle
        Particle(int id, double p, std::array<double, Dimension> pos, std::array<double, Dimension> v, double radius, bool type)
            : id(id), property(p), pos{pos}, vel{v}, force{}, type(type), radius(radius) {}

        //setter methods for the attributes of the class    

        void setVel(const std::array<double, Dimension> &v){ //original method took as arguments (double vx, double vy)
            //vel[0] = vx;
            //vel[1] = vy;
            for(size_t i=0; i < Dimension; ++i) vel[i] = v[i];
        }

        void setProperty(double p){
            property = p;
        }

        //getter methods for the attributes of the class
        double getProperty() const {
            return property;
        }

        std::array<double, Dimension> getPos() const{
            return pos;
        }

        std::array<double, Dimension> getVel() const{
            return vel; 
        }


        std::array<double, Dimension> getForce() const{
        return force;
        }


        int getId() const{
            return id;
        }

        
        double getRadius() const{
            return radius;
        }

        bool getType() const{
            return type;
        }

        //method that resets total force for next implementation
        void resetForce() {
            for (size_t i = 0; i < Dimension; ++i ) {
                force[i] = 0.0; 
            }
        }

        //method that calculates the square of the distance
        double squareDistance(const Particle<Dimension> &p) const{
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

        //method that adds fo 
        void addForce(Particle<Dimension> &k, const Force<Dimension>& f) {
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

        //method that updates positions and velocities using Euler integration
        void update(double delta_t) {
            
            //pos[0] += vel[0] * delta_t;
            //pos[1] += vel[1] * delta_t;

            for(size_t i = 0; i < Dimension; ++i) pos[i] += vel[i] * delta_t;
            
            //vel[0] += (force[0] / property) * delta_t;
            //vel[1] += (force[1] / property) * delta_t;

            for(size_t i = 0; i < Dimension; ++i) vel[i] += (force[i] / property) * delta_t;
        }

        //method that updates positions and velocities using Euler integration + resets force
        void updateAndReset(const double delta_t) {
            
            //pos[0] += vel[0] * delta_t;
            //pos[1] += vel[1] * delta_t;

            for(size_t i = 0; i < Dimension; ++i) pos[i] += vel[i] * delta_t;
            
            // vel[0] += (force[0] / property) * delta_t;
            // vel[1] += (force[1] / property) * delta_t;

            for(size_t i = 0; i < Dimension; ++i) vel[i] += (force[i] / property) * delta_t;

            resetForce();
        }

        //method that prints info about the particles
        void printStates() const{
            std::cout << "Id: " << id << std::endl;
            
            //std::cout << "Position: " << pos[0] << " " << pos[1] << std::endl;
            std::cout << "Position: ";
            for (size_t i = 0; i < Dimension; ++i) {
                std::cout << pos[i];
                if (i < Dimension - 1) {
                    std::cout << " ";
                }
            }
            std::cout << std::endl;

            (!type? std::cout << "Mass: " : std::cout << "Charge: ");
            std::cout << property << std::endl;

            //std::cout << "Force: " << force[0] << " " << force[1] << std::endl;
            std::cout << "Force: ";
            for (size_t i = 0; i < Dimension; ++i) {
                std::cout << force[i];
                if (i < Dimension - 1) {
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
            
            //std::cout << "Velocity: " << vel[0] << " " << vel[1] << std::endl;
            std::cout << "Velocity: ";
            for (size_t i = 0; i < Dimension; ++i) {
                std::cout << vel[i];
                if (i < Dimension - 1) {
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
        }

        //method that manages collision
        void manageCollision(Particle<Dimension> &p, double dim){
            if(dim){
                //manages collision with the borders of the simulation (visual purposes only)
                for(size_t i = 0; i < Dimension; ++i){
                    if (pos[i] + radius > dim || pos[i] - radius < -dim) vel[i] = -vel[i];
                }
                
                /*if(pos[0] + radius > dim || pos[0] - radius < -dim){
                    vel[0] = -vel[0];
                }
                if(pos[1] + radius > dim || pos[1] - radius < -dim){
                    vel[1] = -vel[1];
                }*/
            }
            else{ //manages collision between particles
                //elastic collision
                //if(elastic){
                    //double xvel_prev, yvel_prev, xvel_new, yvel_new;
                    std::array<double, Dimension> prev_vel, new_vel;
                    //xvel_prev = vel[0];
                    //yvel_prev = vel[1];
                    for(size_t i = 0; i < Dimension; ++i) prev_vel[i] = vel[i];
                    for(size_t i = 0; i < Dimension; ++i) vel[i] = ((property - p.getProperty())*vel[i]+2*p.getProperty()*p.getVel()[i]) / (property + p.getProperty());
                    //vel[0] = ((property - p.getProperty())*vel[0]+2*p.getProperty()*p.getVel()[0]) / (property + p.getProperty());
                    //vel[1] = ((property - p.getProperty())*vel[1]+2*p.getProperty()*p.getVel()[1]) / (property + p.getProperty());
                    for(size_t i = 0; i < Dimension; ++i) new_vel[i] = ((p.getProperty() - property)*p.getVel()[i]+2*property*prev_vel[i]) / (property + p.getProperty());
                    //xvel_new = ((p.getProperty() - property)*p.getVel()[0]+2*property*xvel_prev) / (property + p.getProperty());
                    //yvel_new = ((p.getProperty() - property)*p.getVel()[1]+2*property*yvel_prev) / (property + p.getProperty());
                    p.setVel(new_vel);
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

        bool hitsBoundary(double dim){
            bool bound_touched = false;
            for(size_t i = 0; i < Dimension; ++i){
                if(getPos()[i] + getRadius() > dim || getPos()[i] - getRadius() < -dim ){ 
                    bound_touched = true;
                    return bound_touched;
                }
            } 
            
            return bound_touched;
        }

        ~Particle(){}
    private:
        // If of the particle
        int id;
        // The property (mass, charge...) of the particle
        double property;
        // The position of the particle as a two-dimensional vector
        //std::array<double, 2> pos;
        std::array<double, Dimension> pos;
        // The velocity of the particle as a two-dimensional vector
        //std::array<double, 2> vel;
        std::array<double, Dimension> vel;
        // The force acting on the particle as a two-dimensional vector
        std::array<double, Dimension> force;
        //type of the property (mass, charge...)
        bool type;
        // radius of the particle
        double radius;
};

#endif