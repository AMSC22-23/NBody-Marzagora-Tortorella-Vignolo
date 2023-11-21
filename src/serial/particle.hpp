#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <iostream>
#include <array>

class Particle {
    public:
        // A constructor that initializes the mass, position, and velocity of the particle
        //to do: make it better
        Particle(double m, double c, double x, double y, double vx, double vy) {
            mass = m;
            charge = c;
            //pos.resize(2);
            pos[0] = x;
            pos[1] = y;
            //vel.resize(2);
            vel[0] = vx;
            vel[1] = vy;
            //force.resize(2);
            force[0] = 0.0;
            force[1] = 0.0;
        }

        double getMass() const;
        double getCharge() const;
        std::array<double, 2> getPos() const;
        std::array<double, 2> getVel() const;
        std::array<double, 2> getForce() const;
        void resetForce();
        void addForce(double fx, double fy);
        void update(double delta_t);
        void print_states() const;
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
};

#endif 