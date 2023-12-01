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

        void setVel(const std::array<double, Dimension> &v){
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

            std::array<double, Dimension> diff;
            for(size_t i=0; i < Dimension; ++i) diff[i] = k_pos[i] - p_pos[i];
            
            double square_dist = 0.0;
            for(const auto& d : diff) square_dist += d * d;
            
            return square_dist;
        }

        //method that adds force 
        void addForce(Particle<Dimension> &k, const Force<Dimension>& f) {
            std::array<double,Dimension> force_qk;
            force_qk = f.calculateForce(k, *this);

            for(size_t i = 0; i < Dimension; ++i){
                force[i] += force_qk[i];
                k.force[i] -= force_qk[i];
            }
        }

        //method that updates positions and velocities using Euler integration
        void update(double delta_t) {;

            for(size_t i = 0; i < Dimension; ++i) pos[i] += vel[i] * delta_t;
            for(size_t i = 0; i < Dimension; ++i) vel[i] += (force[i] / property) * delta_t;
        }

        //method that updates positions and velocities using Euler integration + resets force
        void updateAndReset(const double delta_t) {

            for(size_t i = 0; i < Dimension; ++i) pos[i] += vel[i] * delta_t;

            for(size_t i = 0; i < Dimension; ++i) vel[i] += (force[i] / property) * delta_t;

            resetForce();
        }

        //method that prints info about the particles
        void printStates() const{
            std::cout << "Id: " << id << std::endl;
            
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

            std::cout << "Force: ";
            for (size_t i = 0; i < Dimension; ++i) {
                std::cout << force[i];
                if (i < Dimension - 1) {
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
            
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
            }
            else{ //manages collision between particles
                std::array<double, Dimension> prev_vel, new_vel;
            
                for(size_t i = 0; i < Dimension; ++i) prev_vel[i] = vel[i];
                for(size_t i = 0; i < Dimension; ++i) vel[i] = ((property - p.getProperty())*vel[i]+2*p.getProperty()*p.getVel()[i]) / (property + p.getProperty());
                for(size_t i = 0; i < Dimension; ++i) new_vel[i] = ((p.getProperty() - property)*p.getVel()[i]+2*property*prev_vel[i]) / (property + p.getProperty());
                p.setVel(new_vel);
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