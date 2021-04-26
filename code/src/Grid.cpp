
#include "Grid.hpp"
#include "Vector2D.hpp"
#include "Cell.hpp"

Grid::Grid(real_t domainSize, unsigned int gridSize) {
    this->domainSize = domainSize;
    this->gridSize = gridSize;
    this->cellSize = domainSize / gridSize;
    initGrid();
}

//TODO: remove this fct. gridSize will be determine in the LJSimulation class.
Grid::Grid(real_t domainSize, unsigned int gridSize, real_t dc) {
    this->domainSize = domainSize;
    // since we know the cut-off distance, we can compute the size of the domain:
    this->gridSize = domainSize / dc;
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
    //! bottleneck
    // TODO: each cell must know its neighbors
    // index of the cell:
    unsigned int i = (int) (position.x / cellSize);
    unsigned int j = (int) (position.y / cellSize);
    // empty vector:
    std::vector<Cell*> nearCells(9);
    // 9 neighboring cells in total:
    for (int y_shift = 0; y_shift < 3; y_shift++)
    {
        for (int x_shift = 0; x_shift < 3; x_shift++)
        {
            // %gridSize because of the peridical boundaries
            nearCells[x_shift*3 + y_shift] = &cells[((j + y_shift - 1 + gridSize)%gridSize)*gridSize + ((i + x_shift - 1 + gridSize)%gridSize)]; // + gridSize to behave like a "true" modulus
        }
        
    }
    return nearCells;
}
