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
        //Particle(int id, double m, double c, double x, double y, double vx, double vy, double radius)
        //    : id(id), mass(m), charge(c), pos{x, y}, vel{vx, vy}, force{0.0, 0.0}, type(0), radius(radius) {}

        // Constructor for template class
         Particle(int id, double m, double c, std::array<double, Dimension> p, std::array<double, Dimension> v, double r)
            : id(id), mass(m), charge(c), pos(p), vel(v), force({}), type(0), radius(r) {}

        //getter methods for class attributes
        double getMass() const;
        double getCharge() const;
        int getId() const;
        
        //std::array<double, 2> getPos() const;
        //std::array<double, 2> getVel() const;
        //std::array<double, 2> getForce() const;

        //getter methods for class template
        std::array<double, Dimension> getPos() const;
        std::array<double, Dimension> getVel() const;
        std::array<double, Dimension> getForce() const;

        //method that resets total force for next implementation
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
        bool returnType() const;
        //method that returns the distance between two particles
        double square_distance(const Particle<Dimension> &p) const;
        //method that returns the radius of the particle
        double getRadius() const;
        //method that manage collisions
        void manage_collision(Particle<Dimension> &p, double dim);
        //method that set the velocity of the particle
        //void setVel(double vx, double vy);
        void setVel(std::array<double, Dimension> v);
        //method that sets the mass of the particle (only to manage the inelastic collisions)
        void setMass(double m);

        ~Particle(){}
    private:
        // If of the particle
        int id;
        // The mass of the particle
        double mass;
        // The position of the particle as a two-dimensional vector
        //std::array<double, 2> pos;
        std::array<double, Dimension> pos;
        // The velocity of the particle as a two-dimensional vector
        //std::array<double, 2> vel;
        std::array<double, Dimension> vel;
        // The force acting on the particle as a two-dimensional vector
        //std::array<double, 2> force;
        std::array<double, Dimension> force;
        // Charge of the particle(for the Coulomb force)
        double charge;
        // type of particle: 1 if Coulomb, 0 otherwise 
        bool type;
        // radius of the particle
        double radius;
};

#endif 