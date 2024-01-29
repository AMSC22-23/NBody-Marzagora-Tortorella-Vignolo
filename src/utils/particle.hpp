#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <iostream>
#include <array>


template<size_t Dimension>
class Force;

/**
 * @brief Particle class that models particles
 * @tparam Dimension Number of dimensions of the simulation 
 */
template<size_t Dimension>
class Particle {
    public:
        /**
        * @brief Constructor that initializes the mass, position, and velocity of the particle
        * @param id Id of the particle
        * @param p Property of particle: mass for gravitational force, charge for Coulomb force
        * @param pos Array of positions of the particle
        * @param v Array of velocities of the particle
        * @param radius Radius of the particle
        * @param type Type of the particle: 0 for gravitational, 1 for Coulomb 
        * */
        Particle(int id, double p, std::array<double, Dimension> pos, std::array<double, Dimension> v, double radius, bool type)
            : id(id), property(p), pos{pos}, vel{v}, force{}, type(type), radius(radius) {
                for(size_t j = 0; j < Dimension; ++j) accel[j] = 0.0;
            }

        /**
         * @brief Setter method for velocity of the particle
         * @param v Array of velocities of the particle
         * 
        */
        void setVel(const std::array<double, Dimension> &v){
            for(size_t i=0; i < Dimension; ++i) vel[i] = v[i];
        }

        /**
         * @brief Setter method for property of the particle
         * @param p Property of particle: mass for gravitational force, charge for Coulomb force
        */
        void setProperty(double p){
            property = p;
        }

        /**
         * @brief Getter method for property of the particle
         * @return property Property of particle: mass for gravitational force, charge for Coulomb force
        */
        double getProperty() const {
            return property;
        }

        std::array<double, Dimension> getAccel() const {
            return accel;
        }

        /**
         * @brief Getter method for positions of the particle
         * @return pos Array of positions of the particle
        */
        std::array<double, Dimension> getPos() const{
            return pos;
        }

        /**
         * @brief Getter method for velocity of the particle
         * @return vel Array of velocity of the particle
        */
        std::array<double, Dimension> getVel() const{
            return vel; 
        }

        /**
         * @brief Getter method for force of the particle
         * @return force Array of force of the particle
        */
        std::array<double, Dimension> getForce() const{
        return force;
        }


        /**
         * @brief Getter method for the ID of the particle
         * @return id ID of the particle
        */
        int getId() const{
            return id;
        }

        /**
         * @brief Getter method for the radius of the particle
         * @return radius Radius of the particle
        */
        double getRadius() const{
            return radius;
        }

        /**
         * @brief Getter method for the ID of the particle
         * @return type Type of the particle: 0 for gravitational, 1 for Coulomb 
        */
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
        double squareDistance(const Particle<Dimension> &p) const{
            const auto& k_pos = getPos();
            const auto& p_pos = p.getPos();

            std::array<double, Dimension> diff;
            for(size_t i=0; i < Dimension; ++i) diff[i] = k_pos[i] - p_pos[i];
            
            double square_dist = 0.0;
            for(const auto& d : diff) square_dist += d * d;
            
            return square_dist;
        }

        /**
         * @brief Method that adds force
         * @param force_qk Array of components of the force between particles q and k
         * 
        */
        void addForce(const std::array<double, Dimension> &force_qk) {

            for(size_t i = 0; i < Dimension; ++i){
                force[i] += force_qk[i];
            }

        }

        /**
         * @brief Method that updates positions and velocities using Euler integration
         * @param delta_t Time step after which the state of the particle is being updated
         * 
        */
        void update(const double delta_t) {
            for(size_t i = 0; i < Dimension; ++i) 
            {
                pos[i] += vel[i] * delta_t;
                vel[i] += (force[i] / ((property<0)? -property:property)) * delta_t;
            }
        }

        void velocityVerletUpdate(const double delta_t) {
            std::array<double, Dimension> prevAccel;
            for(size_t i=0; i<Dimension; ++i)
            {
                prevAccel[i] = accel[i];
                accel[i] = (force[i] / ((property<0)? -property:property));
                pos[i] += vel[i] * delta_t + 0.5 * prevAccel[i] * delta_t * delta_t;
                vel[i] += 0.5 * (prevAccel[i] + accel[i]) * delta_t;
            } 
        }

        /**
         * @brief Method that updates positions and velocities using Euler integration and resets force
         * @param delta_t Time step after which the state of the particle is being updated
         */
        void updateAndReset(const double delta_t) {

            for(size_t i = 0; i < Dimension; ++i) pos[i] += vel[i] * delta_t;
            for(size_t i = 0; i < Dimension; ++i) vel[i] += (force[i] / property) * delta_t;

            resetForce();
        }

        /**
         * @brief Method that prints info about the particles
         */
        void printStates() const{
            std::cout << "Id: " << id << std::endl;
            
            std::cout << "Position: ";
            for (size_t i = 0; i < Dimension; ++i) {
                std::cout << pos[i];
                if (i < Dimension - 1) {
                    std::cout << " ";
                }
            }
            std::cout << std::endl;

            (!type? std::cout << "Mass: " : std::cout << "Charge: ");
            std::cout << property << std::endl;

            std::cout << "Force: ";
            for (size_t i = 0; i < Dimension; ++i) {
                std::cout << force[i];
                if (i < Dimension - 1) {
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
            
            std::cout << "Velocity: ";
            for (size_t i = 0; i < Dimension; ++i) {
                std::cout << vel[i];
                if (i < Dimension - 1) {
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
        }

        /**
         * @brief Method that manages collisions between two particles and between a particle and the borders of the simulation (visual purposes only)
         * @param p Particle which has collided with this particle
         * @param dim Number of dimensions of the simulation
         */
        void manageCollision(Particle<Dimension> &p, double dim){
            if(dim){
                for(size_t i = 0; i < Dimension; ++i){
                    if (pos[i] + radius > dim){
                        vel[i] = -vel[i];
                        if(pos[i] > dim)
                            pos[i] = dim;
                    }
                    else if ( pos[i] - radius < -dim)
                    {
                        vel[i] = -vel[i];
                        if(pos[i] < -dim)
                            pos[i] = -dim;
                    }
                }
            }
            else{ 
                std::array<double, Dimension> prev_vel, new_vel;
            
                for(size_t i = 0; i < Dimension; ++i) prev_vel[i] = vel[i];
                for(size_t i = 0; i < Dimension; ++i) vel[i] = ((property - p.getProperty())*vel[i]+2*p.getProperty()*p.getVel()[i]) / (property + p.getProperty());
                for(size_t i = 0; i < Dimension; ++i) new_vel[i] = ((p.getProperty() - property)*p.getVel()[i]+2*property*prev_vel[i]) / (property + p.getProperty());
                p.setVel(new_vel);
            }
        }

        /**
         * @brief Method that checks if the boundary has been touched
         * @param dim Number of dimensions of the simulation
         * @return bound_touched true if the border has been touched, false otherwise
         */
        bool hitsBoundary(double dim){
            bool bound_touched = false;
            for(size_t i = 0; i < Dimension; ++i){
                if(getPos()[i] + getRadius() > dim || getPos()[i] - getRadius() < -dim ){ 
                    bound_touched = true;
                    return bound_touched;
                }
            } 
            
            return bound_touched;
        }

        /**
         * @brief Default destructor of Particle class
        */
        ~Particle(){}
    private:
        int id;
        double property;
        std::array<double, Dimension> pos;
        std::array<double, Dimension> vel;
        std::array<double, Dimension> force;
        std::array<double, Dimension> accel;
        bool type;
        double radius;
};

#endif