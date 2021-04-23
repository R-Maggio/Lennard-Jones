#pragma once

#include<vector>

#include "Vector2D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"

class Cell; // forward declaration

class Grid
{
private:
    unsigned int gridSize; // gridSize x gridSize = number of cells
    real_t domainSize; // size of the domain
    real_t cellSize; // cellSize = domainSize / gridSize;
    std::vector<Cell> cells;

    /**
     * @brief private function that initialize the grid
     * 
     */
    void initGrid();

public:
    // constructors:
    Grid(real_t domainSize, unsigned int gridSize);
    //TODO: remove this constructor. gridSize will be determine in the LJSimulation class.
    Grid(real_t domainSize, unsigned int gridSize, real_t dc); // gridSize is determined by the cut-off distance => gridSize should be left empty.

    // getters and setters:

    unsigned int getGridSize();

    real_t getDomainSize();

    real_t getCellSize();

    std::vector<Cell>& getCells();
    const std::vector<Cell>& getCells() const;

    // functions:

    /**
     * @brief place a partice in the good cell
     * 
     * @param particle 
     */
    void placeParticle(Particle* particle);

    /**
     * @brief return the list of cells where we need the search (dc is the cut-off distance)
     * 
     * @param position 
     * @param dc dc is the cutt-off distance. Dc MUST be smaller or equal to the size of the cell.
     * @return std::vector<Cell*>
     */
    std::vector<Cell*> findCells(const Vector2D& position, real_t dc);
};
