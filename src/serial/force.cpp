#include <vector>

// Abstract base class 
class Force {
public:
    virtual ~Force() {}
    // Pure virtual function for calculating force
    virtual std::vector<double> calculateForce() = 0;
};

// Derived class for gravitational force
class GravitationalForce : public Force {
public:
    std::vector<double> calculateForce() override {
        std::vector<double> force;
        // Calculate gravitational force here
        return force;
    }
};

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
