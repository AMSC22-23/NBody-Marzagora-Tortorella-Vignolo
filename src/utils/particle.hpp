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

        /*//getter methods for class attributes
        double getProperty() const;
        int getId() const;
        
        //std::array<double, 2> getPos() const;
        //std::array<double, 2> getVel() const;
        //std::array<double, 2> getForce() const;

        //getter methods for class template
        std::array<double, Dimension> getPos() const;
        std::array<double, Dimension> getVel() const;
        std::array<double, Dimension> getForce() const;

        
        void resetForce();
        //method that adds new values of a force
        void addForce(Particle<Dimension> &k, const Force<Dimension>& f);
        // method that updates current values of a force
        void update(double delta_t);
        //method that updates current values of a force and resets it for next implementation
        void update_and_reset(const double delta_t);
        //method that prints current states of the particles 
        void print_states() const;
        //method that returns type of the particles to instatiate the right kind of force
        bool getType() const;
        //method that returns the distance between two particles
        double square_distance(const Particle<Dimension> &p) const;
        //method that returns the radius of the particle
        double getRadius() const;
        //method that manage collisions
        void manage_collision(Particle<Dimension> &p, double dim);
        //method that set the velocity of the particle
        // pass by const reference, not copy!
        void setVel(const std::array<double, Dimension> &v);
        //method that sets the property of the particle (only to manage the inelastic collisions)
        void setProperty(const double m);*/

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
        double square_distance(const Particle<Dimension> &p) const{
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
        void update_and_reset(const double delta_t) {
            
            //pos[0] += vel[0] * delta_t;
            //pos[1] += vel[1] * delta_t;

            for(size_t i = 0; i < Dimension; ++i) pos[i] += vel[i] * delta_t;
            
            // vel[0] += (force[0] / property) * delta_t;
            // vel[1] += (force[1] / property) * delta_t;

            for(size_t i = 0; i < Dimension; ++i) vel[i] += (force[i] / property) * delta_t;

            resetForce();
        }

        //method that prints info about the particles
        void print_states() const{
            std::cout << "Id: " << id << std::endl;
            std::cout << "Position: " << pos[0] << " " << pos[1] << std::endl;
            (!type? std::cout << "Mass: " : std::cout << "Charge: ");
            std::cout << property << std::endl;
            std::cout << "Force: " << force[0] << " " << force[1] << std::endl;
            std::cout << "Velocity: " << vel[0] << " " << vel[1] << std::endl;
        }

        //method that manages collision
        void manage_collision(Particle<Dimension> &p, double dim){
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