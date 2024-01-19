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

        return -1; // Should not happen :D
    }

    void insert(Particle<Dimension>* newParticle) {
        // Se il nodo è una foglia, ma non contiene ancora una particella
        if (leaf && particle == nullptr) {
            particle = newParticle; // Inserisci la particella qui
            return;
        }

        // Se il nodo è una foglia e contiene già una particella, bisogna dividerlo
        if (leaf && particle != nullptr) {
            split(); // Divide il nodo in 4 nodi figli

            // Reinserisci la particella esistente nei nuovi nodi figli
            int existingParticleIndex = getQuadrantIndex(particle->getPos()[0], particle->getPos()[1], x, y);
            if (existingParticleIndex != -1) {
                children[existingParticleIndex]->insert(particle);
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
