#ifndef QUADTREE_NODE_HPP
#define QUADTREE_NODE_HPP

#include <array>
#include <vector>
#include <memory>
#include "particle.hpp"  

using namespace std;

template<size_t Dimension>
class Particle;

template <size_t Dimension>
class QuadtreeNode{
    public:
    QuadtreeNode(double x, double y, double w, int maxDepthValue = 10) 
        : width(w), leaf(true), //particle(nullptr), 
          totalCenter(std::array<double, Dimension>()), center({x,y}), 
          totalMass(0.0), count(0), maxDepth(maxDepthValue)  {
          //children.fill(nullptr);
    }

    ~QuadtreeNode() {}

    const std::array<double, Dimension> & getCenter() const{
        return center;
    }

    double getWidth() const{
        return width;
    }

    bool isLeaf() const{
        return leaf;
    }

    const std::shared_ptr<Particle<Dimension>>& getParticle() const{
        return particle;
    }

    const std::array<double, Dimension>& getTotalCenter() const{
        return totalCenter;
    }   

    double getTotalMass() const{
        return totalMass;
    }

    int getCount() const{
        return count;
    }

    const std::array<std::unique_ptr<QuadtreeNode<Dimension>>, 4>& getChildren() const{
        return children;
    }

    void insertNode(std::unique_ptr<QuadtreeNode<Dimension>> subTreeRoot, int depth = 0) {
        std::cout<<"insertNode start"<<std::endl;
        if (!subTreeRoot) {
            cerr << "Errore: il nodo da inserire è nullo." << endl;
            return;
        }

        if (leaf) {
            std::cout<<"il nodo di intNode è una foglia"<<std::endl;
            // Se il nodo corrente è una foglia, lo trasformiamo in un nodo interno e inseriamo il sottoalbero
            split();
            std::cout<<"end split"<<std::endl;
            leaf = false;
        }

        // Trova il quadrante in cui inserire il sottoalbero
        int index = getQuadrantIndex(subTreeRoot->getCenter());
        std::cout<<"index: "<<index<<std::endl;
        if (index == -1) {
            cerr << "Errore: indice del quadrante non valido per il nodo da inserire." << endl;
            return;
        }

        if (!children[index]) {
            // Se il nodo figlio non esiste, lo inizializziamo con il sottoalbero
            children[index] = std::move(subTreeRoot);
        } else {
            // Altrimenti, inseriamo il sottoalbero nel nodo figlio esistente
            children[index]->insertNode(std::move(subTreeRoot), depth + 1);
            if(depth >5)
                return;
        }

        // Aggiorna gli attributi del nodo corrente in base al sottoalbero inserito
        //updateAttributesAfterNodeInsert(subTreeRoot);
    }

    void updateAttributesForSubtree(const QuadtreeNode<Dimension>* subtreeRoot) {
        if (!subtreeRoot) return;

        // Calcola la massa totale e il centro di massa del sottoalbero
        double subtreeMass = subtreeRoot->getTotalMass();
        const auto& subtreeCenter = subtreeRoot->getTotalCenter();

        // Aggiorna la massa totale del nodo corrente
        double oldTotalMass = totalMass;
        totalMass += subtreeMass;

        // Aggiorna il centro di massa del nodo corrente
        for (size_t i = 0; i < Dimension; ++i) {
            double sumWeights = totalCenter[i] * oldTotalMass;
            sumWeights += subtreeCenter[i] * subtreeMass;
            totalCenter[i] = sumWeights / totalMass;
        }

        // Incrementa il conteggio delle particelle
        count += subtreeRoot->getCount();
    }



    void split(){
        
        leaf = false;
        
        double halfWidth = width / 2.0;
        double quarterWidth = width / 4.0;

        double x = center[0];
        double y = center[1];

        if(width == 0) return;

        //std::cout<<"---------split"<<std::endl;
        //std::cout<<width<<std::endl;
        //std::cout<<x<<std::endl;
        //std::cout<<y<<std::endl;

        children[0] = std::make_unique<QuadtreeNode<Dimension>>(x - quarterWidth, y - quarterWidth, halfWidth);
        children[1] = std::make_unique<QuadtreeNode<Dimension>>(x + quarterWidth, y - quarterWidth, halfWidth);
        children[2] = std::make_unique<QuadtreeNode<Dimension>>(x - quarterWidth, y + quarterWidth, halfWidth);
        children[3] = std::make_unique<QuadtreeNode<Dimension>>(x + quarterWidth, y + quarterWidth, halfWidth);
    }

    int getQuadrantIndex(const std::array<double,Dimension>& pos){
        bool isLeft = pos[0] < center[0];
        bool isTop = pos[1] < center[1]; 

        if (isLeft && isTop) return 0; // Top-left
        if (!isLeft && isTop) return 1; // Top-right
        if (isLeft && !isTop) return 2; // Bottom-left
        if (!isLeft && !isTop) return 3; // Bottom-right

        return -1;
    }

    void insert(std::shared_ptr<Particle<Dimension>> p, int depth = 0){
        updateAttributes(p);

        if(leaf && particle == nullptr){
            particle = p;
            return;
        }

        if(leaf && particle != nullptr){
            //if (depth >= maxDepth) {
            //    //std::cout<<"MAXDEPTH"<<std::endl;
            //    // Gestisci la situazione, ad esempio ignorando l'inserimento o aggregando le particelle
            //    return;
            //}
            split();
            int index = getQuadrantIndex(particle->getPos());
            if(index != -1) children[index]->insert(particle);
            else{
                cerr << "Errore: indice del quadrante non valido durante il reinserimento." << endl;
            }
            particle = nullptr;
        }

        int index = getQuadrantIndex(p->getPos());
        if(index != -1){
            if(children[index] == nullptr){
                cerr << "Errore: nodo figlio non inizializzato." << endl;
                return;
            }
            children[index]->insert(p, depth + 1);
        }
        else{
            cerr << "Errore: indice del quadrante non valido per la nuova particella." << endl;
        }
    }

    void updateAttributes(std::shared_ptr<Particle<Dimension>> newParticle){
        double oldTotalMass = totalMass;
        totalMass += newParticle->getProperty(); // Property of particle: mass for gravitational force, charge for Coulomb force.
        count++;
        for(size_t i = 0; i < Dimension; ++i) {
            double sumWeights = totalCenter[i] * oldTotalMass;
            sumWeights += newParticle->getPos()[i] * newParticle->getProperty();
            totalCenter[i] = sumWeights / totalMass;
        }
    }

    void printTree(int depth = 0) const {
        // Print the current node's details
        std::cout << std::string(depth * 4, ' ') << "Node at depth " << depth 
                  << ": [" << center[0] << ", " << center[1] << ", " << width << "], " 
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

    std::unique_ptr<Particle<Dimension>> createApproximateParticle(){
        std::array<double, Dimension> pos;
        std::array<double, Dimension> vel;
        for(size_t i = 0; i < Dimension; ++i) {
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

#endif // QUADTREE_NODE_H