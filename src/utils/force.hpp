#ifndef FORCE_H
#define FORCE_H

#include <array>
#include <cmath>
//#include "particle.hpp"

class Particle;

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
        double const G = 6.67430e-11;
};

class CoulombForce : public Force{
    public:
        std::array<double,2> calculateForce(const Particle &k, const Particle &q) const override;
    private:
        double const K = 8.987e-09;
};

class NuclearForce : public Force{
    public: 
        std::array<double,2> calculateForce(const Particle &k, const Particle &q) const override;
    private:
        double const r0 = 1.0;
};

class RepulsiveForce : public Force {
    public: 
        std::array<double,2> calculateForce(const Particle &k, const Particle &q) const override;
    private:
        double const kp = 1.0;
};

class CustomForce : public Force{
    public:
        std::array<double,2> calculateForce(const Particle &k, const Particle &q) const override;
    private:
        double const G = 2000; 
};

#endif
