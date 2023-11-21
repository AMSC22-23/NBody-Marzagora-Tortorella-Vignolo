#include "force.hpp"

std::array<double,2> GravitationalForce::calculateForce(const Particle &k, const Particle &q) const{
    double x_diff = k.getPos()[0] - q.getPos()[0];
    double y_diff = k.getPos()[1] - q.getPos()[1];
    double dist = sqrt(x_diff * x_diff + y_diff * y_diff);
    double dist_cubed = dist * dist * dist;

    std::array<double, 2> force_qk;
    force_qk[0] = G * q.getMass() * k.getMass() / dist_cubed * x_diff;
    force_qk[1] = G * q.getMass() * k.getMass() / dist_cubed * y_diff;

    return force_qk; 
}

std::array<double,2> CoulombForce::calculateForce(const Particle &k, const Particle &q) const{
    double x_diff = k.getPos()[0] - q.getPos()[0];
    double y_diff = k.getPos()[1] - q.getPos()[1];
    double dist = sqrt(x_diff * x_diff + y_diff * y_diff);
    double dist_cubed = dist * dist * dist;

    std::array<double, 2> force_qk;
    force_qk[0] = K * q.getCharge() * k.getCharge() / dist_cubed * x_diff;
    force_qk[1] = K * q.getCharge() * k.getCharge() / dist_cubed * y_diff; 

    return force_qk;
}

//Other forces to be added

/*
// Example of usage
int main() {
    Force* force = new GravitationalForce();
    std::vector<double> gravitationalForce = force->calculateForce();
    delete force;

    force = new AnotherForce();
    std::vector<double> anotherForce = force->calculateForce();
    delete force;

    return 0;
}
*/
