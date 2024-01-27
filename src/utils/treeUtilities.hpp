#include <vector>
#include <fstream>
#include <iostream>
#include <random>
#include "particle.hpp" 
#include "quadtreeNode.hpp"
#include "force.hpp"

/**
 * @brief Function that determines the simulation boundaries, creates the root of the tree and then inserts each particle into the tree.
 * It checks if the particle is inside the boundaries of the simulation area (standard case) and, if so, inserts it into the tree. 
 * It also handles the case in which the particle is on the boundary of the simulation area.
 * 
 * @tparam Dimension Number of dimensions of the simulation
 * @param particles Vector of particles to be inserted into the tree
 * @param dimSimulationArea Dimension of the simulation area
*/
template <size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> createQuadTree(std::vector<Particle<Dimension>>& particles,
                                                        double dimSimulationArea) {
    if (particles.empty()) {
        return nullptr;
    }

    double minX = -1 * dimSimulationArea * 0.5;
    double maxX = dimSimulationArea * 0.5;
    double minY = minX;
    double maxY = maxX;

    
    double width = std::max(maxX - minX, maxY - minY);
    std::unique_ptr<QuadtreeNode<Dimension>> root =
        std::make_unique<QuadtreeNode<Dimension>>(minX + width / 2, minY + width / 2, width, 20);

    for (auto& particle : particles) {
        std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
        root->insert(p);
    }

    return root;
}

/**
 * @brief Function that determines the centers of the subtrees based on the region width and the dimension of the simulation area.
 * 
 * @tparam Dimension Number of dimensions of the simulation
 * @param regionWidth Width of the subregion
 * @param dim Dimension of the simulation area
 * @param x X coordinate of the center of the subregion
 * @param y Y coordinate of the center of the subregion
 * @param centers Vector of arrays of doubles representing the centers of the subregions
*/
template <size_t Dimension>
void assignRegions(double regionWidth, double dim, double x, double y,
                   std::vector<std::array<double, Dimension>>* centers) {
    if (dim == regionWidth) {
        centers->push_back({x - regionWidth / 2, y + regionWidth / 2});
        centers->push_back({x + regionWidth / 2, y + regionWidth / 2});
        centers->push_back({x - regionWidth / 2, y - regionWidth / 2});
        centers->push_back({x + regionWidth / 2, y - regionWidth / 2});

    } else {
        assignRegions(regionWidth, dim / 2, x - dim / 2, y + dim / 2, centers);
        assignRegions(regionWidth, dim / 2, x + dim / 2, y + dim / 2, centers);
        assignRegions(regionWidth, dim / 2, x - dim / 2, y - dim / 2, centers);
        assignRegions(regionWidth, dim / 2, x + dim / 2, y - dim / 2, centers);
    }
}

/**
 * @brief Getter function that returns the centers of the subtrees based on the number of threads and the dimension of the simulation area.
 * 
 * @tparam Dimension Number of dimensions of the simulation
 * @param num_threads Number of threads
 * @param dim Dimension of the simulation area
 * @return std::vector<std::array<double, Dimension>> Vector of arrays of doubles representing the centers of the subregions
*/
template <size_t Dimension>
std::vector<std::array<double, Dimension>> getSubRegionsCoordinates(int num_threads, double dim) {
    std::vector<std::array<double, Dimension>> regions;

    int rows = std::sqrt(num_threads);
    while (num_threads % rows != 0) {
        rows--;
    }
    int cols = num_threads / rows;

    double regionWidth = dim / cols;
    double regionHeight = dim / rows;

    assignRegions(regionWidth, dim / 2, 0.0, 0.0, &regions);

    return regions;
}

/**
 * @brief Getter function that returns the dimension of the subtrees based on the number of threads and the dimension of the simulation area.
 * 
 * @tparam Dimension Number of dimensions of the simulation
 * @param num_threads Number of threads
 * @param dim Dimension of the simulation area
 * @return std::array<double, Dimension> Array of doubles representing the dimension of the subregions
*/
template <size_t Dimension>
std::array<double, Dimension> getSubRegionsDimension(int num_threads, double dim) {
    std::array<double, Dimension> regionsDim;

    int rows = std::sqrt(num_threads);

    while (num_threads % rows != 0) {
        rows--;
    }

    int cols = num_threads / rows;

    double regionWidth = dim / cols;
    double regionHeight = dim / rows;
    regionsDim[0] = regionWidth;
    regionsDim[1] = regionHeight;

    return regionsDim;
}



