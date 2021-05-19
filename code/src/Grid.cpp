
#include "Grid.hpp"
#include "Vector2D.hpp"
#include "Cell.hpp"


Grid::Grid(real_t domainSize, unsigned int gridSize, LJBoundary boundary, real_t dc) {
    this->boundary = boundary;
    this->domainSize = domainSize;
    if (gridSize == 0) {
        // since we know the cut-off distance, we can compute the size of the domain:
        this->gridSize = domainSize / dc;
    } else {
        this->gridSize = gridSize;
    }
    this->cellSize = domainSize / this->gridSize;
    initGrid();
}

// getters and setters:
unsigned int Grid::getGridSize() {return gridSize;}

real_t Grid::getDomainSize() {return domainSize;}

real_t Grid::getCellSize() {return cellSize;} 

std::vector<Cell>& Grid::getCells() {return cells;}
const std::vector<Cell>& Grid::getCells() const {return cells;}

void Grid::initGrid() {
    cells = std::vector<Cell>(gridSize*gridSize); // gridSize^2 cells
    // fill the cells:
    for (size_t j = 0; j < gridSize; j++) // y direction
    {
        for (size_t i = 0; i < gridSize; i++) // x direction
        {
            cells[j*gridSize + i].set(Vector2D(i*cellSize, j*cellSize), Vector2D((i+1)*cellSize, (j+1)*cellSize), this);
        }
        
    }
    initCellsNeighbors(); // init the neighbors of each cell
}

void Grid::placeParticle(Particle* particle) {
    auto position = particle->getPosition() % (gridSize * cellSize);
    // we need to find the cell according to the position of the particle
    // i is the index of the x-axis and j the index of the y-axis
    // we have that: i = int(px/cellSize), j = int(px/cellSize)
    unsigned int i = ((int) (position.x / cellSize) % gridSize); // % gridSize to avoid the case where position.x == gridsize*cellSize.... => index out of range
    unsigned int j = ((int) (position.y / cellSize) % gridSize);
    // add the particle tho the corresponding cell:
    cells[j*gridSize + i].addParticle(particle);
}

std::vector<Cell*> Grid::findCells(const Vector2D& position, real_t dc) {
    //! WARNING: we suppose that dc <= cellSize.
    //! dc is unused...
    // empty vector:
    std::vector<Cell*> nearCells;
    // index of the cell:
    unsigned int i = (((int) (position.x / cellSize)) + gridSize)%gridSize; // %gridSize because of the peridical boundaries. + gridSize to behave like a "true" modulus
    unsigned int j = (((int) (position.y / cellSize)) + gridSize)%gridSize;
    Cell* currentCell = &cells[j*gridSize + i];

    nearCells = currentCell->getNeighbors();
    nearCells.push_back(currentCell);
    return nearCells;
}

std::vector<Cell*> Grid::findNeighbors(size_t i, size_t j) {
    // empty vector:
    std::vector<Cell*> neighbors;
    // depend on the boundaries:
    if (boundary == LJBoundary::PERIODIC) // Periodic boundaries
    {
        // 8 neighboring cells in total:
        for (int y_shift = -1; y_shift < 2; y_shift++)
        {
            for (int x_shift = -1; x_shift < 2; x_shift++)
            {
                // %gridSize because of the peridical boundaries
                if (x_shift != 0 || y_shift != 0) // to avoid including itself
                    neighbors.push_back(&cells[((j + y_shift + gridSize)%gridSize)*gridSize + ((i + x_shift + gridSize)%gridSize)]); // + gridSize to behave like a "true" modulus
            }
        }
    } else if (boundary == LJBoundary::POISSEUILLE) { // boundaries for the Poisseuille flow experiment
        // the code is repeated for performance reasons (test boundary only once)
        for (int y_shift = -1; y_shift < 2; y_shift++)
        {
            for (int x_shift = -1; x_shift < 2; x_shift++)
            {
                // %gridSize because of the peridical boundaries
                if ((x_shift != 0 || y_shift != 0) && (j + y_shift < gridSize && j + y_shift >= 0 && i + x_shift < gridSize && i + x_shift >= 0)) // to avoid including itself + Poisseuille boundaries
                    neighbors.push_back(&cells[((j + y_shift + gridSize)%gridSize)*gridSize + ((i + x_shift + gridSize)%gridSize)]); // + gridSize to behave like a "true" modulus
            }
        } 
    }
    return neighbors;
}

void Grid::initCellsNeighbors() {
    for (size_t j = 0; j < gridSize; j++) // y direction
    {
        for (size_t i = 0; i < gridSize; i++) // x direction
        {
            cells[j*gridSize + i].setNeighbors(findNeighbors(i, j)); // init the neighbors
        }
        
    }
    
    
}
