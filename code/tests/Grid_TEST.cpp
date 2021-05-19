#include <iostream> // std::cout
#include <string> // std::string

#include "../include/Vector2D.hpp"
#include "../include/Particle.hpp"
#include "../include/Grid.hpp"
#include "../include/typeDefinitions.hpp"


int main(int argc, char const *argv[])
{
    //! TEST of the Grid/Cell class:
    
    Grid grid(100., 0, LJBoundary::PERIODIC, 2.6);
    std::cout << "domainSize:" << grid.getDomainSize() << ", gridSize:" << grid.getGridSize() << ", cellSize:" << grid.getCellSize() << '\n';

    std::cout << " grid.getCells()[39].getLl():   " << grid.getCells()[39].getLl() << '\n';

    Vector2D r = Vector2D(3, 3);
    Vector2D v = Vector2D(1, 1);
    Vector2D a = Vector2D(-1, 2);
    real_t m = 1.;
    Particle p1(r, v, a, m);

    grid.placeParticle(&p1);
    std::cout << " grid.getCells()[39].getLocalParticles()[0]->getPosition():   " << grid.getCells()[39].getLocalParticles()[0]->getPosition() << '\n';

    auto nearCells = grid.findCells(Vector2D(1., 1), 2.5);
    for (auto &e : nearCells)
    {
        std::cout << "grid.findCells(Vector2D(1, 1), 2.5)[i] :  " << e->getLl() << '\n';
    }
    std::cout << " grid.getCells()[0].isPositionInCell(Vector2D(2, 1)):   " << grid.getCells()[0].isPositionInCell(Vector2D(2, 1)) << '\n';
    std::cout << " grid.getCells()[0].isPositionInCell(Vector2D(1, -1)):   " << grid.getCells()[0].isPositionInCell(Vector2D(1, -1)) << '\n';

    return 0;
}