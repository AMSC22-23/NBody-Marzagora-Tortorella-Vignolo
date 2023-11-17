#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {
    public:
        // A constructor that initializes the mass, position, and velocity of the particle
        Particle(double m, double x, double y, double vx, double vy) {
            mass = m;
            pos[0] = x;
            pos[1] = y;
            vel[0] = vx;
            vel[1] = vy;
            force[0] = 0.0;
            force[1] = 0.0;
        }

        double getMass();
        double* getPos();
        double* getVel();
        double* getForce();
        void resetForce();
        void addForce(double fx, double fy);
        void update(double delta_t);
    private:
        // The mass of the particle
        double mass;
        // The position of the particle as a two-dimensional vector
        double pos[2];
        // The velocity of the particle as a two-dimensional vector
        double vel[2];
        // The force acting on the particle as a two-dimensional vector
        double force[2];
};

#endif 