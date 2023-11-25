#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <iostream>
#include <array>

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
        Particle(double m, double c, double x, double y, double vx, double vy)
            : mass(m), charge(c), pos{x, y}, vel{vx, vy}, force{0.0, 0.0}, type(0) {}

        //getter methods for class attributes
        double getMass() const;
        double getCharge() const;
        std::array<double, 2> getPos() const;
        std::array<double, 2> getVel() const;
        std::array<double, 2> getForce() const;
        //method that resets total force for next implementation
        void resetForce();
        //method that adds new values of a force
        void addForce(double fx, double fy);
        // method that updates current values of a force
        void update(double delta_t);
        //method that prints current states of the particles 
        void print_states() const;
        //method that returns type of the particles to instatiate the right kind of force
        bool returnType() const;
        ~Particle(){}
    private:
        // The mass of the particle
        double mass;
        // The position of the particle as a two-dimensional vector
        std::array<double, 2> pos;
        // The velocity of the particle as a two-dimensional vector
        std::array<double, 2> vel;
        // The force acting on the particle as a two-dimensional vector
        std::array<double, 2> force;
        // Charge of the particle(for the Coulomb force)
        double charge;
        // type of particle: 1 if Coulomb, 0 otherwise 
        bool type;
        // variables to compute distance
        double x_diff, y_diff;
};

#endif 