#ifndef FORCE_H
#define FORCE_H

#include <array>
#include <cmath>
#include "particle.hpp"

class Force{
    public:
        //Constructor
        Force(){}

        // Pure virtual function for calculating force
        virtual std::array<double,2> calculateForce(const Particle &p1, const Particle &p2) const = 0;

        //Destructor 
        virtual ~Force() {}
};

class GravitationalForce : public Force{
    public:
        std::array<double,2> calculateForce(const Particle &k, const Particle &q) const override;
    private:
        double const G = 1; //true value = 6.67430e-11
};

class CoulombForce : public Force{
    public:
        std::array<double,2> calculateForce(const Particle &k, const Particle &q) const override;
    private:
        double const K = 8.987e-09;
};

//Other forces to be added

#endif