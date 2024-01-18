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

    // Method to determine which child node a given vector belongs to
    // This is used to decide in which quadrant a point lies
    int which(const Vector<Dimension>& v) {
        // Implement the logic to determine the appropriate quadrant
        return 0;  // Placeholder return value
    }

    // Method to insert a new Particle into the TreeNode
    // Handles the logic of placing the particle in the correct node
    void insert(Particle<Dimension>* newP) {
        updateAttributes(newParticle);
        if (leaf) {
            if (particle == nullptr) {
                // Node is empty, place the particle here
                particle = newParticle;
            } else {
                // Node is a leaf but already contains a particle
                split(); // Split the node

                // Reinsert the existing particle
                int existingParticleIndex = getQuadrantIndex(particle->getPosition()[0], 
                                                            particle->getPosition()[1], x, y);
                children[existingParticleIndex]->insert(particle);

                // Insert the new particle
                int newParticleIndex = getQuadrantIndex(newParticle->getPosition()[0], 
                                                        newParticle->getPosition()[1], x, y);
                children[newParticleIndex]->insert(newParticle);

                particle = nullptr; // Clear the particle pointer after reinsertion
            }
        } else {
            // Node is not a leaf, find the correct child for the new particle
            int index = getQuadrantIndex(newParticle->getPosition()[0], 
                                        newParticle->getPosition()[1], x, y);
            children[index]->insert(newParticle);
        }
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

    // Method to update the node's mass and count attributes
    void updateAttributes(Particle<Dimension>* newParticle) {
        totalMass += newParticle->getProperty(); // Property of particle: mass for gravitational force, charge for Coulomb force.
        count++;
    }
};

#endif // TREENODE_H
