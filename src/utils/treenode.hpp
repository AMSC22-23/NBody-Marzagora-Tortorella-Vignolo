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
          totalCenter(std::array<double, Dimension>()), center({x,y}), 
          totalMass(0.0), count(0) {
        children.fill(nullptr);  // Initialize all children pointers to nullptr
    }

    // Destructor for the TreeNode class
    // Responsible for cleaning up dynamically allocated child TreeNodes
    ~TreeNode() {
        for (auto& child : children) {
            if(child->approximatedParticle != nullptr) delete child->approximatedParticle;
            delete child;  // Delete each dynamically allocated child
        }
    }

    // Getters for private attributes
    double getX() const {
        return x;
    }

    double getY() const {
        return y;
    }

    double getWidth() const {
        return w;
    }

    bool isLeaf() const {
        return leaf;
    }

    Particle<Dimension>* getParticle() const {
        return particle;
    }

    const std::array<double, Dimension>& getTotalCenter() const {
        return totalCenter;
    }

    const std::array<double, Dimension> getCenter() const {
        return center;
    }

    double getTotalMass() const {
        return totalMass;
    }

    int getCount() const {
        return count;
    }

    const std::array<TreeNode<Dimension>*, 4>& getChildren() const {
        return children;
    }

    void split() {
        // Calcola la metà della larghezza e un quarto della larghezza per determinare i centri dei nuovi nodi figli
        double halfWidth = w / 2.0;
        double quarterWidth = w / 4.0;

        // Assicurati che i nodi figli siano nulli prima di crearne di nuovi
        for (auto& child : children) {
            if (child != nullptr) {
                delete child;
                child = nullptr;
            }
        }

        // Crea quattro nuovi nodi figli
        children[0] = new TreeNode(x - quarterWidth, y - quarterWidth, halfWidth); // Quadrante in alto a sinistra
        children[1] = new TreeNode(x + quarterWidth, y - quarterWidth, halfWidth); // Quadrante in alto a destra
        children[2] = new TreeNode(x - quarterWidth, y + quarterWidth, halfWidth); // Quadrante in basso a sinistra
        children[3] = new TreeNode(x + quarterWidth, y + quarterWidth, halfWidth); // Quadrante in basso a destra

        // Imposta il nodo corrente come non foglia
        leaf = false;
    }

    // Helper function to determine the quadrant index for a given particle
    int getQuadrantIndex(double x, double y, double centerX, double centerY) {
        bool isLeft = x < centerX;
        bool isTop = y < centerY;

        if (isLeft && isTop) return 0; // Top-left
        if (!isLeft && isTop) return 1; // Top-right
        if (isLeft && !isTop) return 2; // Bottom-left
        if (!isLeft && !isTop) return 3; // Bottom-right

        std::cout<<"NON VA BENE";
        return -1; // Should not happen :D
    }

    void insert(Particle<Dimension>* newParticle) {
        updateAttributes(newParticle);

        // Se il nodo è una foglia, ma non contiene ancora una particella
        if (leaf && particle == nullptr) {
            particle = newParticle; // Inserisci la particella qui
            approximatedParticle = createApproximatedParticle();
            return;
        }

        // Se il nodo è una foglia e contiene già una particella, bisogna dividerlo
        if (leaf && particle != nullptr) {
            split(); // Divide il nodo in 4 nodi figli

            // Reinserisci la particella esistente nei nuovi nodi figli
            int existingParticleIndex = getQuadrantIndex(particle->getPos()[0], particle->getPos()[1], x, y);
            if (existingParticleIndex != -1) {
                children[existingParticleIndex]->insert(particle);
                approximatedParticle = createApproximatedParticle();
            } else {
                cerr << "Errore: indice del quadrante non valido durante il reinserimento." << endl;
            }
            particle = nullptr; // Pulisci il puntatore della particella dopo il reinserimento
        }

        // Inserisci la nuova particella nel nodo figlio appropriato
        int newParticleIndex = getQuadrantIndex(newParticle->getPos()[0], newParticle->getPos()[1], x, y);
        if (newParticleIndex != -1) {
            if (children[newParticleIndex] == nullptr) {
                cerr << "Errore: nodo figlio non inizializzato." << endl;
                return;
            }
            children[newParticleIndex]->insert(newParticle);
        } else {
            cerr << "Errore: indice del quadrante non valido per la nuova particella." << endl;
        }
    }

    // Method to update the node's mass and count attributes
    void updateAttributes(Particle<Dimension>* newParticle) {
        double oldTotalMass = totalMass;
        totalMass += newParticle->getProperty(); // Property of particle: mass for gravitational force, charge for Coulomb force.
        count++;
        for(size_t i = 0; i < Dimension; ++i) {
            double sumWeights = totalCenter[i] * oldTotalMass;
            sumWeights += newParticle->getPos()[i] * newParticle->getProperty();
            totalCenter[i] = sumWeights / totalMass;
        }
    }

    
    // Method to display the tree
    void printTree(int depth = 0) const {
        // Print the current node's details
        std::cout << std::string(depth * 4, ' ') << "Node at depth " << depth 
                  << ": [" << x << ", " << y << ", " << w << "], " 
                  << "Mass: " << totalMass << ", Particles: " << count << ",";
        // Check if the node is a leaf or an internal node
        std::cout << (leaf ? " Leaf node\n" : " Internal node\n");
        
        if (leaf) {
            // If leaf, display the particle if it exists
            if (particle != nullptr) {
                std::cout << std::string((depth + 1) * 4, ' ') << "Particle: "
                          << particle->getId() << ") " << particle->getPos()[0] << ", " << particle->getPos()[1] <<". mass: "<< particle->getProperty()<<"\n"; // Assuming Particle class has toString method
            }
        } else {
            // If not a leaf, recursively display child nodes
            for (const auto& child : children) {
                if (child != nullptr) {
                    child->printTree(depth + 1);
                }
            }
        }
    }

    // Particle(int id, double p, std::array<double, Dimension> pos, std::array<double, Dimension> v, double radius, bool type
    Particle<Dimension>* createApproximatedParticle() {
        // Create a particle with the center of mass and total mass of the node
        // Return the particle
        return new Particle<Dimension>(count,  totalMass, totalCenter, {0.0, 0.0}, 1.0, false);
    }

    // Method that returns the approximated particle for the node
    Particle<Dimension>* getApproximatedParticle() const {
        return approximatedParticle;
    }

private:
    double x, y, w;  // x, y coordinates and width of the TreeNode region
    std::array<TreeNode<Dimension>*, 4> children;  // Array of pointers to child TreeNodes
    bool leaf;  // Flag to indicate if the TreeNode is a leaf node
    Particle<Dimension>* particle;  // Pointer to a Particle stored in the TreeNode
    std::array<double, Dimension> totalCenter;  // Total center of mass for the TreeNode
    std::array<double, Dimension> center;  // Center of mass for the TreeNode
    Particle<Dimension>* approximatedParticle;  // Approximated Particle for the TreeNode
    double totalMass;  // Total mass of all particles in the TreeNode
    int count;  // Count of particles in the TreeNode
};

#endif // TREENODE_H
