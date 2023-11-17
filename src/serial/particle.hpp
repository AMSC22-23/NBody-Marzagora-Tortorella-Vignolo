#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <iostream>

class Particle {
    public:
        // A constructor that initializes the mass, position, and velocity of the particle
        Particle(double m, double x, double y, double vx, double vy) {
            mass = m;
            pos.resize(2);
            pos[0] = x;
            pos[1] = y;
            vel.resize(2);
            vel[0] = vx;
            vel[1] = vy;
            force.resize(2);
            force[0] = 0.0;
            force[1] = 0.0;
        }

        double getMass() const;
        //double* getPos();
        //double* getVel();
        //double* getForce();
        std::vector<double> getPos() const;
        std::vector<double> getVel() const;
        std::vector<double> getForce() const;
        void resetForce();
        void addForce(double fx, double fy);
        void update(double delta_t);
        void print_states() const;
    private:
        // The mass of the particle
        double mass;
        // The position of the particle as a two-dimensional vector
        std::vector<double> pos;
        //double pos[2];
        // The velocity of the particle as a two-dimensional vector
        std::vector<double> vel;
        //double vel[2];
        // The force acting on the particle as a two-dimensional vector
        std::vector<double> force;
        //double force[2];
};

#endif 