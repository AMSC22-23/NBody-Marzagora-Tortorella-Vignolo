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
        virtual std::array<double,Dimension> calculateForce(const Particle<Dimension> &p1, const Particle<Dimension> &p2) const = 0;

        //Destructor 
        virtual ~Force() {}
};


template<size_t Dimension>
class GravitationalForce : public Force<Dimension>{
    public:
        std::array<double,Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override{
            const auto& k_pos = k.getPos();
            const auto& q_pos = q.getPos();

            std::array<double, Dimension> pos_diff;
            for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = k_pos[i] - q_pos[i];

            double dist_squared = 0.0;
            for(size_t i = 0; i < Dimension; ++i) dist_squared = pos_diff[i] * pos_diff[i] + dist_squared;
            double dist = std::sqrt(dist_squared);
            double dist_cubed = dist * dist * dist;

            double gravF = G * q.getProperty() * k.getProperty() / dist_cubed;

            std::array<double, Dimension> force_qk;
            for(size_t i = 0; i < Dimension; ++i) force_qk[i] = gravF * pos_diff[i];

            return force_qk; 
        }
    private:
        double const G = 6.67430e-11;
};

template<size_t Dimension>
class CoulombForce : public Force<Dimension>{
    public:
        std::array<double,Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override{
            const auto& k_pos = k.getPos();
            const auto& q_pos = q.getPos();

            std::array<double, Dimension> pos_diff;
            for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = k_pos[i] - q_pos[i];

            double dist_squared = 0.0;
            for(size_t i = 0; i < Dimension; ++i) dist_squared = pos_diff[i] * pos_diff[i] + dist_squared;
            double dist = std::sqrt(dist_squared);
            double dist_cubed = dist * dist * dist;

            double coulF = K * q.getProperty() * k.getProperty() / dist_cubed;

            std::array<double, Dimension> force_qk;
            for(size_t i = 0; i < Dimension; ++i) force_qk[i] = coulF * pos_diff[i];

            return force_qk; 
        }
    private:
        double const K = 8.987e-09;
};

template<size_t Dimension>
class NuclearForce : public Force<Dimension>{
    public: 
        std::array<double,Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override{
            const auto& k_pos = k.getPos();
            const auto& q_pos = q.getPos();

            std::array<double, Dimension> pos_diff;
            for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = k_pos[i] - q_pos[i];

            double dist = 0.0;
            for(size_t i = 0; i < Dimension; ++i) dist = pos_diff[i] * pos_diff[i] + dist;

            double nukeF = -std::exp(-dist/r0)/dist;

            std::array<double, Dimension> force_qk;
        
            for(size_t i = 0; i < Dimension; ++i) force_qk[i] = nukeF * pos_diff[i];

            return force_qk;
        }
    private:
        double const r0 = 1.0;
};

template<size_t Dimension>
class RepulsiveForce : public Force<Dimension>{
    public: 
        std::array<double, Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override{
            const auto& k_pos = k.getPos();
            const auto& q_pos = q.getPos();

            std::array<double, Dimension> pos_diff;
            for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = k_pos[i] - q_pos[i];

            double dist = 0.0;
            for(size_t i = 0; i < Dimension; ++i) dist = pos_diff[i] * pos_diff[i] + dist;

            double repF = kp / dist;
            
            std::array<double, Dimension> force_qk;
            
            for(size_t i = 0; i < Dimension; ++i) force_qk[i] = repF * pos_diff[i] / dist;

            return force_qk;
        }
    private:
        double const kp = 1.0;
};

template<size_t Dimension>
class CustomForce : public Force<Dimension>{
    public:
        CustomForce(double G) : G(G) {}
        std::array<double,Dimension> calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const override{
            const auto& k_pos = k.getPos();
            const auto& q_pos = q.getPos();

            std::array<double, Dimension> pos_diff;
            for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = q_pos[i] - k_pos[i];

            double dist_squared = 0.0;
            for(size_t i = 0; i < Dimension; ++i) dist_squared = pos_diff[i] * pos_diff[i] + dist_squared;
            double dist = std::sqrt(dist_squared);
            double dist_cubed = dist * dist * dist;

            double customF = G * q.getProperty() * k.getProperty() / dist_cubed;

            std::array<double, Dimension> force_qk;
            for(size_t i = 0; i < Dimension; ++i) force_qk[i] = customF*pos_diff[i];

            return force_qk; 
        }
    private:
        double const G; 
};

#endif
