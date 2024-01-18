#ifndef TREENODE_H
#define TREENODE_H

#include <array>
#include <vector>
#include "particle.hpp"  

using namespace std;

template<size_t Dimension>
class Particle;

template <size_t Dimension>
class TreeNode {
public:
    // Constructor for the TreeNode class
    // Initializes a TreeNode with specified x, y coordinates and width (w)
    TreeNode(double x, double y, double w) 
        : x(x), y(y), w(w), leaf(true), particle(nullptr), 
          totalCenter(std::array<double, Dimension>()), center(nullptr), 
          totalMass(0.0), count(0) {
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
        float halfWidth = w / 2.0;
        float quarterWidth = w / 4.0;

        // Create four children
        children[0] = new TreeNode(x - quarterWidth, y - quarterWidth, halfWidth);
        children[1] = new TreeNode(x + quarterWidth, y - quarterWidth, halfWidth);
        children[2] = new TreeNode(x - quarterWidth, y + quarterWidth, halfWidth);
        children[3] = new TreeNode(x + quarterWidth, y + quarterWidth, halfWidth);

        leaf = false;

        // If there's already a particle in this node, reinsert it into the correct child
        if (particle != nullptr) {
            insert(particle);
            particle = nullptr; // Clear the particle pointer after reinsertion
        }
        // Implement the logic to split the node into four children
    }

    // Helper function to determine the quadrant index for a given particle
    int getQuadrantIndex(float x, float y, float centerX, float centerY) {
        bool isLeft = x < centerX;
        bool isTop = y < centerY;

        if (isLeft && isTop) return 0; // Top-left
        if (!isLeft && isTop) return 1; // Top-right
        if (isLeft && !isTop) return 2; // Bottom-left
        if (!isLeft && !isTop) return 3; // Bottom-right

        return -1; // Should not happen :D
    }

    // Method to insert a new Particle into the TreeNode
    // Handles the logic of placing the particle in the correct node
    void insert( Particle<Dimension>* newParticle) {
        updateAttributes(newParticle);
        if (leaf) {
            if (particle == nullptr) {
                // Node is empty, place the particle here
                particle = newParticle;
            } else {
                // Node is a leaf but already contains a particle
                split(); // Split the node

                // Reinsert the existing particle
                int existingParticleIndex = getQuadrantIndex(particle->getPos()[0], 
                                                            particle->getPos()[1], x, y);
                children[existingParticleIndex]->insert(particle);

                // Insert the new particle
                int newParticleIndex = getQuadrantIndex(newParticle->getPos()[0], 
                                                        newParticle->getPos()[1], x, y);
                children[newParticleIndex]->insert(newParticle);

                particle = nullptr; // Clear the particle pointer after reinsertion
            }
        } else {
            // Node is not a leaf, find the correct child for the new particle
            int index = getQuadrantIndex(newParticle->getPos()[0], 
                                        newParticle->getPos()[1], x, y);
            children[index]->insert(newParticle);
        }
    }

    // Method to update the node's mass and count attributes
    void updateAttributes(Particle<Dimension>* newParticle) {
        totalMass += newParticle->getProperty(); // Property of particle: mass for gravitational force, charge for Coulomb force.
        count++;
        // update: totalCenter & center
    }

    
    // Method to display the tree
    void display(int depth = 0) const {
        // Print the current node's details
        std::cout << std::string(depth * 4, ' ') << "Node at depth " << depth 
                  << ": [" << x << ", " << y << ", " << w << "], "
                  << "Mass: " << totalMass << ", Particles: " << count << "\n";

        if (leaf) {
            // If leaf, display the particle if it exists
            if (particle != nullptr) {
                std::cout << std::string((depth + 1) * 4, ' ') << "Particle: "
                          << particle->getId() << ") " << particle->getPos()[0] << ", " << particle->getPos()[1] << "\n"; // Assuming Particle class has toString method
            }
        } else {
            // If not a leaf, recursively display child nodes
            for (const auto& child : children) {
                if (child != nullptr) {
                    child->display(depth + 1);
                }
            }
        }
    }

private:
    double x, y, w;  // x, y coordinates and width of the TreeNode region
    std::array<TreeNode<Dimension>*, 4> children;  // Array of pointers to child TreeNodes
    bool leaf;  // Flag to indicate if the TreeNode is a leaf node
    Particle<Dimension>* particle;  // Pointer to a Particle stored in the TreeNode
    std::array<double, Dimension> totalCenter;  // Total center of mass for the TreeNode
    std::array<double, Dimension>* center;  // Center of mass for the TreeNode
    double totalMass;  // Total mass of all particles in the TreeNode
    int count;  // Count of particles in the TreeNode
};

#endif // TREENODE_H
