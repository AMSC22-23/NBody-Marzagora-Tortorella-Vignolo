#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <iostream>
#include <array>
//#include "force.hpp"

class Force;

class Particle {
    public:
        /* Old constructor, badly implemented
        Particle(double m, double c, double x, double y, double vx, double vy) {
            mass = m;
            charge = c;
            pos[0] = x;
            pos[1] = y;
            vel[0] = vx;
            vel[1] = vy;
            force[0] = 0.0;
            force[1] = 0.0;
            type = 0;
        }*/

        // A constructor that initializes the mass, position, and velocity of the particle
        Particle(int id, double p, double x, double y, double vx, double vy, double radius, bool type)
            : id(id), property(p), pos{x, y}, vel{vx, vy}, force{0.0, 0.0}, type(type), radius(radius) {}

        //getter methods for class attributes
        double getProperty() const;
        int getId() const;
        std::array<double, 2> getPos() const;
        std::array<double, 2> getVel() const;
        std::array<double, 2> getForce() const;
        //method that resets total force for next implementation
        void resetForce();
        //method that adds new values of a force
        //void addForce(double fx, double fy);
        void addForce(Particle &k, const Force& f);
        // method that updates current values of a force
        void update(double delta_t);
        //method that updates current values of a force and resets it for next implementation
        void update_and_reset(const double delta_t);
        //method that prints current states of the particles 
        void print_states() const;
        //method that returns type of the particles to instatiate the right kind of force
        bool getType() const;
        //method that returns the distance between two particles
        double square_distance(const Particle &p) const;
        //method that returns the radius of the particle
        double getRadius() const;
        //method that manage collisions
        void manage_collision(Particle &p, double dim);
        //method that set the velocity of the particle
        void setVel(double vx, double vy);
        //method that sets the property of the particle (only to manage the inelastic collisions)
        void setProperty(double m);

        ~Particle(){}
    private:
        // If of the particle
        int id;
        // The property (mass, charge...) of the particle
        double property;
        // The position of the particle as a two-dimensional vector
        std::array<double, 2> pos;
        // The velocity of the particle as a two-dimensional vector
        std::array<double, 2> vel;
        // The force acting on the particle as a two-dimensional vector
        std::array<double, 2> force;
        //type of the property (mass, charge...)
        bool type;
        // radius of the particle
        double radius;
};

#endif 