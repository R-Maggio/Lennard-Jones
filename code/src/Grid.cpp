
#include "Grid.hpp"
#include "Vector2D.hpp"
#include "Cell.hpp"


Grid::Grid(real_t domainSize, unsigned int gridSize, LJBoundary boundary, real_t dc) {
    this->boundary = boundary;
    this->domainSize = {domainSize, domainSize};
    if (gridSize == 0) {
        // since we know the cut-off distance, we can compute the size of the domain:
        this->gridSizeX = domainSize / dc;
    } else {
        this->gridSizeX = gridSize;
    }
    this->gridSizeY = this->gridSizeX;
    this->cellSize = this->domainSize / gridSize;
    initGrid();
}

Grid::Grid(const Vector2D& domainSize, const Vector2D& gridSize, LJBoundary boundary, real_t dc) {
    this->boundary = boundary;
    this->domainSize = domainSize;
    this->gridSizeX = gridSize.x;
    this->gridSizeY = gridSize.y;
    this->cellSize = {this->domainSize.x / this->gridSizeX, this->domainSize.y / this->gridSizeY,};
    initGrid();
}

// getters and setters:
Vector2D Grid::getGridSize() {return {this->gridSizeX, this->gridSizeY};}

real_t Grid::getGridSizeX() {return this->gridSizeX;}

real_t Grid::getGridSizeY() {return this->gridSizeY;}

Vector2D Grid::getDomainSize() {return domainSize;}

Vector2D Grid::getCellSize() {return cellSize;} 

std::vector<Cell>& Grid::getCells() {return cells;}
const std::vector<Cell>& Grid::getCells() const {return cells;}

void Grid::initGrid() {
    cells = std::vector<Cell>(gridSizeX*gridSizeY); // gridSize^2 cells
    // fill the cells:
    for (size_t j = 0; j < gridSizeY; j++) // y direction
    {
        for (size_t i = 0; i < gridSizeX; i++) // x direction
        {
            cells[j*gridSizeX + i].set(Vector2D(i*cellSize.x, j*cellSize.y), Vector2D((i+1)*cellSize.x, (j+1)*cellSize.y), this);
        }
        
    }
    initCellsNeighbors(); // init the neighbors of each cell
}

void Grid::placeParticle(Particle* particle) {
    auto position = particle->getPosition() % (Vector2D(gridSizeX * cellSize.x, gridSizeY* cellSize.y));
    // we need to find the cell according to the position of the particle
    // i is the index of the x-axis and j the index of the y-axis
    // we have that: i = int(px/cellSize), j = int(px/cellSize)
    unsigned int i = ((int) (position.x / cellSize.x) % gridSizeX); // % gridSize to avoid the case where position.x == gridsize*cellSize.... => index out of range
    unsigned int j = ((int) (position.y / cellSize.y) % gridSizeY);
    // add the particle tho the corresponding cell:
    cells[j*gridSizeX + i].addParticle(particle);
}

std::vector<Cell*> Grid::findCells(const Vector2D& position, real_t dc) {
    //! WARNING: we suppose that dc <= cellSize.
    //! dc is unused...
    // empty vector:
    std::vector<Cell*> nearCells;
    // index of the cell:
    unsigned int i = (((int) (position.x / cellSize.x)) + gridSizeX)%gridSizeX; // %gridSize because of the peridical boundaries. + gridSize to behave like a "true" modulus
    unsigned int j = (((int) (position.y / cellSize.y)) + gridSizeY)%gridSizeY;
    Cell* currentCell = &cells[j*gridSizeX + i];

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
                    neighbors.push_back(&cells[((j + y_shift + gridSizeY)%gridSizeY)*gridSizeX + ((i + x_shift + gridSizeX)%gridSizeX)]); // + gridSize to behave like a "true" modulus
            }
        }
    } else if (boundary == LJBoundary::POISSEUILLE) { // boundaries for the Poisseuille flow experiment
        // the code is repeated for performance reasons (test boundary only once)
        for (int y_shift = -1; y_shift < 2; y_shift++)
        {
            for (int x_shift = -1; x_shift < 2; x_shift++)
            {
                // %gridSize because of the peridical boundaries
                if ((x_shift != 0 || y_shift != 0) && (j + y_shift < gridSizeY && j + y_shift >= 0 && i + x_shift < gridSizeX && i + x_shift >= 0)) // to avoid including itself + Poisseuille boundaries
                    neighbors.push_back(&cells[((j + y_shift + gridSizeY)%gridSizeY)*gridSizeX + ((i + x_shift + gridSizeX)%gridSizeX)]); // + gridSize to behave like a "true" modulus
            }
        } 
    }
    return neighbors;
}

void Grid::initCellsNeighbors() {
    for (size_t j = 0; j < gridSizeY; j++) // y direction
    {
        for (size_t i = 0; i < gridSizeX; i++) // x direction
        {
            cells[j*gridSizeX + i].setNeighbors(findNeighbors(i, j)); // init the neighbors
        }
        
    }
}
