#ifndef TREENODE_H
#define TREENODE_H

#include <array>
#include "particle.hpp"  

template<size_t Dimension>
class Particle;

template <size_t Dimension>
class TreeNode {
public:
    // Constructor for the TreeNode class
    // Initializes a TreeNode with specified x, y coordinates and width (w)
    TreeNode(float x, float y, float w) 
        : x(x), y(y), w(w), leaf(true), particle(nullptr), 
          totalCenter(Vector<Dimension>()), center(nullptr), 
          totalMass(0.0f), count(0) {
        children.fill(nullptr);  // Initialize all children pointers to nullptr
    }

    // Destructor for the TreeNode class
    // Responsible for cleaning up dynamically allocated child TreeNodes
    ~TreeNode() {
        for (auto& child : children) {
            delete child;  // Delete each dynamically allocated child
        }
    }

    // Method to split the TreeNode into four child nodes
    // This is typically used in quadtree implementations
    void split() {
        // Implement the logic to split the node into four children
    }

    // Method to determine which child node a given vector belongs to
    // This is used to decide in which quadrant a point lies
    int which(const Vector<Dimension>& v) {
        // Implement the logic to determine the appropriate quadrant
        return 0;  // Placeholder return value
    }

    // Method to insert a new Particle into the TreeNode
    // Handles the logic of placing the particle in the correct node
    void insert(Particle<Dimension>* newP) {
        // Implement the logic to insert a new Particle
    }

    // Method to display the TreeNode
    // This could be a graphical representation or textual output
    void display() {
        // Implement the logic for displaying the TreeNode
    }

private:
    float x, y, w;  // x, y coordinates and width of the TreeNode region
    std::array<TreeNode<Dimension>*, 4> children;  // Array of pointers to child TreeNodes
    bool leaf;  // Flag to indicate if the TreeNode is a leaf node
    Particle<Dimension>* particle;  // Pointer to a Particle stored in the TreeNode
    Vector<Dimension> totalCenter;  // Total center of mass for the TreeNode
    Vector<Dimension>* center;  // Center of mass for the TreeNode
    float totalMass;  // Total mass of all particles in the TreeNode
    int count;  // Count of particles in the TreeNode
};

#endif // TREENODE_H
