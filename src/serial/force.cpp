#include "force.hpp"
#include "particle.hpp"

std::array<double,2> GravitationalForce::calculateForce(const Particle &k, const Particle &q) const{
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    double x_diff = k_pos[0] - q_pos[0];
    double y_diff = k_pos[1] - q_pos[1];

   /*double x_diff = k.getPos()[0] - q.getPos()[0];
    double y_diff = k.getPos()[1] - q.getPos()[1];
    double dist = sqrt(x_diff * x_diff + y_diff * y_diff);*/
    double dist = std::hypot(x_diff, y_diff);
    double dist_cubed = dist * dist * dist;

    double gravF = G * q.getMass() * k.getMass() / dist_cubed;

    std::array<double, 2> force_qk;
    force_qk[0] = gravF * x_diff;
    force_qk[1] = gravF * y_diff;
    return force_qk; 
}

std::array<double,2> CustomForce::calculateForce(const Particle &k, const Particle &q) const{
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    double x_diff = k_pos[0] - q_pos[0];
    double y_diff = k_pos[1] - q_pos[1];

   /*double x_diff = k.getPos()[0] - q.getPos()[0];
    double y_diff = k.getPos()[1] - q.getPos()[1];
    double dist = sqrt(x_diff * x_diff + y_diff * y_diff);*/
    double dist = std::hypot(x_diff, y_diff);
    double dist_cubed = dist * dist * dist;

    double gravF = G * q.getMass() * k.getMass() / dist_cubed;

    std::array<double, 2> force_qk;
    force_qk[0] = gravF * x_diff;
    force_qk[1] = gravF * y_diff;
    return force_qk; 
}

std::array<double,2> CoulombForce::calculateForce(const Particle &k, const Particle &q) const{
    /*double x_diff = k.getPos()[0] - q.getPos()[0];
    double y_diff = k.getPos()[1] - q.getPos()[1];
    double dist = sqrt(x_diff * x_diff + y_diff * y_diff);*/
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    double x_diff = k_pos[0] - q_pos[0];
    double y_diff = k_pos[1] - q_pos[1];

    double dist = std::hypot(x_diff, y_diff);
    double dist_cubed = dist * dist * dist;

    double coulF = K * q.getCharge() * k.getCharge() / dist_cubed;

    std::array<double, 2> force_qk;
    force_qk[0] =  coulF * x_diff;
    force_qk[1] = coulF * y_diff;

    return force_qk;
}

std::array<double,2> NuclearForce::calculateForce(const Particle &k, const Particle &q) const {
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    double x_diff = k_pos[0] - q_pos[0];
    double y_diff = k_pos[1] - q_pos[1];

    double dist = std::hypot(x_diff, y_diff);

    double nukeF = -std::exp(-dist/r0)/dist;

    std::array<double, 2> force_qk;
    force_qk[0] = nukeF * x_diff / dist;
    force_qk[1] = nukeF * y_diff / dist;

    return force_qk;
}

std::array<double,2> RepulsiveForce::calculateForce(const Particle &k, const Particle &q) const{
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    double x_diff = k_pos[0] - q_pos[0];
    double y_diff = k_pos[1] - q_pos[1];

    double dist = std::hypot(x_diff, y_diff);

    double repF = kp / dist;
    
    std::array<double, 2> force_qk;
    force_qk[0] = repF * x_diff / dist;
    force_qk[1] = repF * x_diff / dist;

    return force_qk;

}
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
