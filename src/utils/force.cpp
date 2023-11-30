#include "force.hpp"
#include "particle.hpp"

template<size_t Dimension>
std::array<double,Dimension> GravitationalForce<Dimension>::calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const{
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    //double x_diff = k_pos[0] - q_pos[0];
    //double y_diff = k_pos[1] - q_pos[1];

    std::array<double, Dimension> pos_diff;
    for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = k_pos[i] - q_pos[i];

    //double dist = std::hypot(x_diff, y_diff);
    double dist = 0.0;
    for(size_t i = 0; i < Dimension; ++i) dist = pos_diff[i] * pos_diff[i] + dist;
    double dist_cubed = dist * dist * dist;

    double gravF = G * q.getProperty() * k.getProperty() / dist_cubed;

    std::array<double, Dimension> force_qk;
    //force_qk[0] = gravF * x_diff;
    //force_qk[1] = gravF * y_diff;

    for(size_t i = 0; i < Dimension; ++i) force_qk[i] = pos_diff[i];

    return force_qk; 
}

template<size_t Dimension>
std::array<double,Dimension> CustomForce<Dimension>::calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const{
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    //double x_diff = k_pos[0] - q_pos[0];
    //double y_diff = k_pos[1] - q_pos[1];

    std::array<double, Dimension> pos_diff;
    for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = k_pos[i] - q_pos[i];


    //double dist = std::hypot(x_diff, y_diff);
    double dist = 0.0;
    for(size_t i = 0; i < Dimension; ++i) dist = pos_diff[i] * pos_diff[i] + dist;
    double dist_cubed = dist * dist * dist;

    double gravF = G * q.getProperty() * k.getProperty() / dist_cubed;

    std::array<double, Dimension> force_qk;
    //force_qk[0] = gravF * x_diff;
    //force_qk[1] = gravF * y_diff;
    for(size_t i = 0; i < Dimension; ++i) force_qk[i] = pos_diff[i];

    return force_qk; 
}

template<size_t Dimension>
std::array<double,Dimension> CoulombForce<Dimension>::calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const{
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    //double x_diff = k_pos[0] - q_pos[0];
   // double y_diff = k_pos[1] - q_pos[1];

    std::array<double, Dimension> pos_diff;
    for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = k_pos[i] - q_pos[i];

    //double dist = std::hypot(x_diff, y_diff);
    double dist = 0.0;
    for(size_t i = 0; i < Dimension; ++i) dist = pos_diff[i] * pos_diff[i] + dist;
    double dist_cubed = dist * dist * dist;

    double coulF = K * q.getProperty() * k.getProperty() / dist_cubed;

    std::array<double, Dimension> force_qk;
    //force_qk[0] =  coulF * x_diff;
    //force_qk[1] = coulF * y_diff;
    for(size_t i = 0; i < Dimension; ++i) force_qk[i] = coulF * pos_diff[i];

    return force_qk;
}

template<size_t Dimension>
std::array<double,Dimension> NuclearForce<Dimension>::calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const {
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    //double x_diff = k_pos[0] - q_pos[0];
    //double y_diff = k_pos[1] - q_pos[1];

    std::array<double, Dimension> pos_diff;
    for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = k_pos[i] - q_pos[i];

    //double dist = std::hypot(x_diff, y_diff);
    double dist = 0.0;
    for(size_t i = 0; i < Dimension; ++i) dist = pos_diff[i] * pos_diff[i] + dist;

    double nukeF = -std::exp(-dist/r0)/dist;

    std::array<double, Dimension> force_qk;
    //force_qk[0] = nukeF * x_diff / dist;
    //force_qk[1] = nukeF * y_diff / dist;
    for(size_t i = 0; i < Dimension; ++i) force_qk[i] = nukeF * pos_diff[i];

    return force_qk;
}

template<size_t Dimension>
std::array<double,Dimension> RepulsiveForce<Dimension>::calculateForce(const Particle<Dimension> &k, const Particle<Dimension> &q) const{
    const auto& k_pos = k.getPos();
    const auto& q_pos = q.getPos();

    //double x_diff = k_pos[0] - q_pos[0];
    //double y_diff = k_pos[1] - q_pos[1];

    std::array<double, Dimension> pos_diff;
    for(size_t i = 0; i < Dimension; ++i) pos_diff[i] = k_pos[i] - q_pos[i];

    //double dist = std::hypot(x_diff, y_diff);
    double dist = 0.0;
    for(size_t i = 0; i < Dimension; ++i) dist = pos_diff[i] * pos_diff[i] + dist;

    double repF = kp / dist;
    
    std::array<double, Dimension> force_qk;
    //force_qk[0] = repF * x_diff / dist;
    //force_qk[1] = repF * x_diff / dist;
    for(size_t i = 0; i < Dimension; ++i) force_qk[i] = repF * pos_diff[i] / dist;

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
