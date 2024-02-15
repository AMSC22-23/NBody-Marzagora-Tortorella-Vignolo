#ifndef QUADTREE_NODE_HPP
#define QUADTREE_NODE_HPP

#include <array>
#include <memory>
#include <vector>

#include "particle.hpp"

using namespace std;

template <size_t Dimension>
class Particle;

/**
 * @brief QuadtreeNode class that models a node of the quadtree needed for the algorithm of Barnes-Hut.
 * @tparam Dimension Number of dimensions of the simulation 
 */
template <size_t Dimension>
class QuadtreeNode {
   public:
    /**
    * @brief Constructor that initializes the center of the node and its width
    * @param width Width of the region covered by the node
    * @param leaf Boolean value that indicates if the node is a leaf or not
    * @param totalCenter Coordinates of the center of mass of the particles of the node
    * @param center Coordinates of the center of the node
    * @param totalMass Total mass of the particles' of the node and its children
    * @param count Number of particles inside the node 
    * */
    QuadtreeNode(double x, double y, double w, int maxDepthValue = 10)
        : width(w),
          leaf(true), 
          totalCenter(std::array<double, Dimension>()),
          center({x, y}),
          totalMass(0.0),
          count(0),
          maxDepth(maxDepthValue) {
    }
    /**
     * @brief Default destructor for the quadtreeNode class
     * @note Not needed, the synthetised one is ok
    */
    ~QuadtreeNode() {}

    /**
     * @brief Getter function that returns the coordinates of the center of the node
     * 
     * @return Coordinates of the center 
    */
    const std::array<double, Dimension>& getCenter() const { return center; }

    /**
     * @brief Getter function that returns the width of the region 
     * 
     * @return Width of the region covered by the center
    */
    double getWidth() const { return width; }

    /**
     * @brief Boolean function that returns true if a node is leaf, false otherwise
     * 
     * @return True if a node is leaf, false otherwise
    */
    bool isLeaf() const { return leaf; }

    /**
     * @brief Getter function that returns the particle inside the node
     * 
     * @return Particle contained inside the node
    */
    const std::shared_ptr<Particle<Dimension>>& getParticle() const { return particle; }

    /**
     * @brief Getter function that returns the coordinates of the center of mass of the node
     * 
     * @return Coordinates of the center of mass of the node
    */
    const std::array<double, Dimension>& getTotalCenter() const { return totalCenter; }

    /**
     * @brief Getter function that returns the sum of the masses of the paticles inside the node
     * 
     * @return sum of the masses of the paticles inside the node
    */
    double getTotalMass() const { return totalMass; }

    /**
     * @brief Getter function that returns the number of particles inside a node and its children 
     * 
     * @return Number of particles
    */
    int getCount() const { return count; }

    /**
     * @brief Getter function that returns the children of a node
     * 
     * @return Array that contains the unique pointers to the children of a node
    */
    const std::array<std::unique_ptr<QuadtreeNode<Dimension>>, 4>& getChildren() const { return children; }

    /**
     * @brief Function that inserts a subtree as a child of a father node: if the current node is a leaf, it becomes an internal node and the 
     * subtree is inserted as a child of the father node. Otherwise, it finds the quadrant in which the subtree has to be inserted and calls:
     * if the quadrant node exists and the child node does not exist, it is initialized with the subtree, otherwise the subtree is inserted in 
     * the child node.
     * 
     * @param subTreeRoot Unique pointer to the subtree to be inserted
    */
    void insertNode(std::unique_ptr<QuadtreeNode<Dimension>> subTreeRoot, int depth = 0) {
        updateAttributes(subTreeRoot->createApproximateParticle(), subTreeRoot->getCount());

        if (!subTreeRoot) {
            cerr << "Error: the node  being inserted is null" << endl;
            return;
        }

        if (leaf) {
        
            split();
            leaf = false;
        }

        
        int index = getQuadrantIndex(subTreeRoot->getCenter());
        if (index == -1) {
            cerr << "Error: invalid index for the node to be inserted" << endl;
            return;
        }

        if (!children[index]->getParticle()) {
            children[index] = std::move(subTreeRoot);//@note nice use of move
        } else {
            children[index]->insertNode(std::move(subTreeRoot), depth + 1);
        }
    }

    /**
     * @brief Function that creates the children of a node: it creates the four children of a node and inserts them in the children array
    */
    void split() {
        leaf = false;

        double halfWidth = width / 2.0;
        double quarterWidth = width / 4.0;

        double x = center[0];
        double y = center[1];

        if (width == 0) return;

        children[0] = std::make_unique<QuadtreeNode<Dimension>>(x - quarterWidth, y - quarterWidth, halfWidth);
        children[1] = std::make_unique<QuadtreeNode<Dimension>>(x + quarterWidth, y - quarterWidth, halfWidth);
        children[2] = std::make_unique<QuadtreeNode<Dimension>>(x - quarterWidth, y + quarterWidth, halfWidth);
        children[3] = std::make_unique<QuadtreeNode<Dimension>>(x + quarterWidth, y + quarterWidth, halfWidth);
    }

    /**
     * @brief Function that returns the index of the quadrant in which the particle has to be inserted
     * 
     * @param pos Coordinates of the particle to be inserted
     * @return Index of the quadrant in which the particle has to be inserted
    */
    int getQuadrantIndex(const std::array<double, Dimension>& pos) {
        bool isLeft = pos[0] < center[0];
        bool isTop = pos[1] < center[1];

        if (isLeft && isTop) return 0;    // Top-left
        if (!isLeft && isTop) return 1;   // Top-right
        if (isLeft && !isTop) return 2;   // Bottom-left
        if (!isLeft && !isTop) return 3;  // Bottom-right

        return -1;
    }

    /**
     * @brief Function that updates the attributes of the node and then inserts a particle in the quadtree: if the node is a leaf and it is empty, 
     * the particle is inserted in the node. Otherwise, if a paricle is already inside the node, the particle is repositioned inside one of the 
     * children, based on its coordinates, and the method is recursively called for the new particle, in order to handle even the case where two 
     * particles should go in the same child.
     * 
     * @param p Particle to be inserted inside the node
    */
    void insert(std::shared_ptr<Particle<Dimension>> p, int depth = 0) {
        updateAttributes(p);

        if (leaf && particle == nullptr) {
            particle = p;
            return;
        }

        if (leaf && particle != nullptr) {
            split();
            int index = getQuadrantIndex(particle->getPos());
            if (index != -1)
                children[index]->insert(particle);
            else {
                cerr << "Error: quadrant index invalid while repositioning" << endl;
            }
            particle = nullptr;
        }

        int index = getQuadrantIndex(p->getPos());
        if (index != -1) {
            if (children[index] == nullptr) {
                cerr << "Error: child not initialised" << endl;
                return;
            }
            children[index]->insert(p, depth + 1);
        } else {
            cerr << "Error: invalid index for the new particle" << endl;
        }
    }

    /**
     * @brief Function that updates the attributes of the node: it updates the total mass of the node and its children and the coordinates of the
     * of the center of mass, in case it is needed to create the approximated particle for the evaluation of the force.
     * 
     * @param newParticle Particle to be inserted inside the node
     * @param numParticles Number of particles to be inserted the node (default is 1, needed to update the count of the particles inside the node)
    */
    void updateAttributes(std::shared_ptr<Particle<Dimension>> newParticle, int numParticles = 1) {
        double oldTotalMass = totalMass;
        totalMass += newParticle->getProperty();  // Property of particle: mass for gravitational force, charge for Coulomb force.
        count += numParticles;
        for (size_t i = 0; i < Dimension; ++i) {
            double sumWeights = totalCenter[i] * oldTotalMass;
            sumWeights += newParticle->getPos()[i] * newParticle->getProperty();
            totalCenter[i] = sumWeights / totalMass;
        }
    }

    /**
     * @brief Function that prints the current node's details: the level where the node is inside the tree, the coordinates of the center, the mass,
     * the number of particles inside the node and if the node is a leaf or an internal node. If the node is a leaf, it displays the particle and exits,
     * otherwise it recursively displays the child nodes.
    */
    void printTree(int depth = 0) const {
        // Print the current node's details
        std::cout << std::string(depth * 4, ' ') << "Node at depth " << depth << ": [" << center[0] << ", " << center[1]
                  << ", " << width << "], "
                  << "Mass: " << totalMass << ", Particles: " << count << ",";
        // Check if the node is a leaf or an internal node
        std::cout << (leaf ? " Leaf node\n" : " Internal node\n");

        if (leaf) {
            
            if (particle != nullptr) {
                std::cout << std::string((depth + 1) * 4, ' ') << "Particle: " << particle->getId() << ") "
                          << particle->getPos()[0] << ", " << particle->getPos()[1]
                          << ". mass: " << particle->getProperty()
                          << "\n"; 
            }
        } else {
            for (const auto& child : children) {
                if (child != nullptr) {
                    child->printTree(depth + 1);
                }
            }
        }
    }

    /**
     * @brief Function that creates the approximated particle which represent a group of particles identified by the center of mass of the particles inside 
     * the quadtree. It is used to evaluate the force between a particle a group of particles that are sufficiently far from each other, so
     * that the approximation is accurate enough. The new particle has id -1 in order not to confuse it for a real particle, the mass is the sum of the masses
     * of the particles inside the node and as center the center of mass of the particles (now initialized as 0.0).
     * 
     * @return Unique pointer to the approximated particle
    */
    std::unique_ptr<Particle<Dimension>> createApproximateParticle() {
        std::array<double, Dimension> pos;
        std::array<double, Dimension> vel;
        for (size_t i = 0; i < Dimension; ++i) {
            pos[i] = totalCenter[i];
            vel[i] = 0.0;
        }
        return std::make_unique<Particle<Dimension>>(-1, totalMass, pos, vel, 1.0, false);
    }

   private:
    double width;
    bool leaf;
    std::shared_ptr<Particle<Dimension>> particle;
    std::array<std::unique_ptr<QuadtreeNode<Dimension>>, 4> children;
    std::array<double, Dimension> totalCenter;
    std::array<double, Dimension> center;
    double totalMass;
    int count;
    int maxDepth;
};

#endif  // QUADTREE_NODE_H