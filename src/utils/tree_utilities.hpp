#include <vector>
#include <fstream>
#include <iostream>
#include <random>
#include "particle.hpp" 
#include "quadtreeNode.hpp"
#include "force.hpp"

template <size_t Dimension>
std::unique_ptr<QuadtreeNode<Dimension>> createQuadTree(std::vector<Particle<Dimension>>& particles,
                                                        double dimSimulationArea) {
    if (particles.empty()) {
        return nullptr;
    }

    // Determine the simulation boundaries (this might need to be adjusted)
    double minX = -1 * dimSimulationArea * 0.5;
    double maxX = dimSimulationArea * 0.5;
    double minY = minX;
    double maxY = maxX;

    // Create the root of the tree
    double width = std::max(maxX - minX, maxY - minY);
    std::unique_ptr<QuadtreeNode<Dimension>> root =
        std::make_unique<QuadtreeNode<Dimension>>(minX + width / 2, minY + width / 2, width, 20);

    // Insert each particle into the tree
    for (auto& particle : particles) {
        std::shared_ptr<Particle<Dimension>> p = std::make_shared<Particle<Dimension>>(particle);
        root->insert(p);
    }

    return root;
}

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



