// A class for a particle in an n-body simulation
class Particle {
    private:
        // The mass of the particle
        double mass;
        // The position of the particle as a two-dimensional vector
        double pos[2];
        // The velocity of the particle as a two-dimensional vector
        double vel[2];
        // The force acting on the particle as a two-dimensional vector
        double force[2];
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

        // A method that returns the mass of the particle
        double getMass() {
            return mass;
        }

        // A method that returns the position of the particle as a pointer to a two-dimensional array
        double* getPos() {
            return pos;
        }

        // A method that returns the velocity of the particle as a pointer to a two-dimensional array
        double* getVel() {
            return vel;
        }

        // A method that returns the force acting on the particle as a pointer to a two-dimensional array
        double* getForce() {
            return force;
        }

        // A method that sets the force acting on the particle to zero
        void resetForce() {
            force[0] = 0.0;
            force[1] = 0.0;
        }

        // A method that adds a given force to the force acting on the particle
        void addForce(double fx, double fy) {
            force[0] += fx;
            force[1] += fy;
        }

        // A method that updates the position and velocity of the particle using Euler's method
        void update(double delta_t) {
            // Update position
            pos[0] += delta_t * vel[0];
            pos[1] += delta_t * vel[1];
            // Update velocity
            vel[0] += delta_t / mass * force[0];
            vel[1] += delta_t / mass * force[1];
        }
};