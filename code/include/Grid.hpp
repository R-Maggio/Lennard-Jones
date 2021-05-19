#pragma once

#include<vector>

#include "Vector2D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"

class Cell; // forward declaration

class Grid
{
private:
    unsigned int gridSizeX; // gridSizeX x gridSizeY = number of cells
    unsigned int gridSizeY;
    Vector2D domainSize; // size of the domain
    Vector2D cellSize; // cellSize = domainSize / gridSize;
    std::vector<Cell> cells;
    LJBoundary boundary; // type of boundary

    /**
     * @brief private function that initialize the grid
     * 
     */
    void initGrid();

    /**
     * @brief init the Cells Neighbors
     * 
     */
    void initCellsNeighbors();

    /**
     * @brief find the neighbors of one cell at position (i, j)
     * 
     * @param i
     * @param j 
     * @return std::vector<Cell*> 
     */
    std::vector<Cell*> findNeighbors(size_t i, size_t j);

public:
    // constructors:
    Grid(real_t domainSize, unsigned int gridSize, LJBoundary boundary = LJBoundary::PERIODIC, real_t dc = 1.); // gridSize is determined by the cut-off distance => gridSize should be left empty.

        // constructors:
    Grid(const Vector2D& domainSize, const Vector2D& gridSize, LJBoundary boundary = LJBoundary::PERIODIC, real_t dc = 1.); // gridSize is determined by the cut-off distance => gridSize should be left empty.

    // getters and setters:

    Vector2D getGridSize();

    real_t getGridSizeX();

    real_t getGridSizeY();
    
    Vector2D getDomainSize();

    Vector2D getCellSize();

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
