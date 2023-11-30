#ifndef FORCE_H
#define FORCE_H

#include <array>
#include <cmath>
//#include "particle.hpp"

template<size_t Dimension>
class Particle;

template<size_t Dimension>
class Force{
    public:
        //Constructor
        Force(){}

        // Pure virtual function for calculating force
        //virtual std::array<double,2> calculateForce(const Particle &p1, const Particle &p2) const = 0;
        virtual std::array<double,Dimension> calculateForce(const Particle<Dimension> &p1, const Particle<Dimension> &p2) const = 0;

        //Destructor 
        virtual ~Force() {}
};


template<size_t Dimension>
class GravitationalForce : public Force<Dimension>{
    public:
        std::array<double,Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override;
    private:
        double const G = 6.67430e-11;
};

template<size_t Dimension>
class CoulombForce : public Force<Dimension>{
    public:
        std::array<double,Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override;
    private:
        double const K = 8.987e-09;
};

template<size_t Dimension>
class NuclearForce : public Force<Dimension>{
    public: 
        std::array<double,Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override;
    private:
        double const r0 = 1.0;
};

template<size_t Dimension>
class RepulsiveForce : public Force<Dimension>{
    public: 
        std::array<double, Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override;
    private:
        double const kp = 1.0;
};

template<size_t Dimension>
class CustomForce : public Force<Dimension>{
    public:
        std::array<double,Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override;
    private:
        double const G = 2000; 
};

#endif
