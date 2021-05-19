#pragma once

#include <vector> 
#include <iterator>

#include "Vector2D.hpp"
#include "Particle.hpp"
#include "Grid.hpp"

// #include "typeDeclaration.hpp"

class Grid; // forward declaration

/**
 * @brief Class Cell
 * 
 */
class Cell
{

// bool isPositionInCell(const Vector2D& position) {
//     if (ll.x <= position.x && ll.y <= position.y && position.x <= ur.x && position.y <= ur.y)
//         return true;
//     return false;
// }
private:
    // lower left and upper right coordinates of the cell:
    Vector2D<> ll{};
    Vector2D<> ur{};
    std::vector<Cell*> neighbors; // list of neighbors cells
    // pointer on the grid:
    Grid* grid = nullptr;
    // list of particles inside the cell. we only store pointers:
    std::vector<Particle*> localParticles;

public:
    // constructors:
    Cell() {};
    Cell(const Vector2D<>& ll, const Vector2D<>& ur, Grid* grid): ll{ll}, ur{ur}, grid{grid} {}
    Cell(const Vector2D<>& ll, const Vector2D<>& ur, Grid* grid, std::vector<Particle*> localParticles): ll{ll}, ur{ur}, grid{grid}, localParticles{localParticles} {}
    Cell(const Vector2D<>& ll, const Vector2D<>& ur, Grid* grid, std::vector<Cell*> neighbors): ll{ll}, ur{ur}, grid{grid}, neighbors{neighbors} {}
    Cell(const Vector2D<>& ll, const Vector2D<>& ur, Grid* grid, std::vector<Particle*> localParticles, std::vector<Cell*> neighbors): ll{ll}, ur{ur}, grid{grid}, localParticles{localParticles}, neighbors{neighbors} {}

    // setters and getters:

    void set(const Vector2D<>& ll, const Vector2D<>& ur, Grid* grid);

    Vector2D<> getLl();
    void setLl(const Vector2D<>& ll);

    Vector2D<> getUr();
    void setUr(const Vector2D<>& ur);

    Grid* getGrid();
    void setGrid(Grid* grid);

    std::vector<Particle*>& getLocalParticles();
    const std::vector<Particle*>& getLocalParticles() const;
    void setLocalParticles(std::vector<Particle*> localParticles);

    std::vector<Cell*> getNeighbors();
    void setNeighbors(std::vector<Cell*> neighbors);
    // functions:

    /**
     * @brief add a pointer on a Particle to the localParticles list
     * 
     * @param particle 
     */
    void addParticle(Particle* particle);
    //! add removeParticle ? maybe inefficient...

    /**
     * @brief check if a point is located inside the cell or not
     * 
     * @param position 
     * @return true if the point is in the cell
     * @return false if the point is not in the cell
     */
    bool isPositionInCell(const Vector2D<>& position);
};

