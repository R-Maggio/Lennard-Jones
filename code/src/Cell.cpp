#include <vector>

#include "Cell.hpp"
#include "Vector2D.hpp"


void Cell::set(const Vector2D& ll, const Vector2D& ur, Grid* grid) {
    this->ll = ll;
    this->ur = ur;
    this->grid = grid;
}

Vector2D Cell::getLl() {
    return this->ll;
}

void Cell::setLl(const Vector2D& ll) {
    this->ll = ll;
}

Vector2D Cell::getUr() {
    return this->ur;
}

void Cell::setUr(const Vector2D& ur) {
    this->ur = ur;
}

Grid* Cell::getGrid() {
    return this->grid;
}

void Cell::setGrid(Grid* grid) {
    this->grid = grid;
}

std::vector<Particle*>& Cell::getLocalParticles() {
    return this->localParticles;
}

const std::vector<Particle*>& Cell::getLocalParticles() const {
    return this->localParticles;
}

void Cell::setLocalParticles(std::vector<Particle*> localParticles) {
    this->localParticles = localParticles;
}

void Cell::addParticle(Particle* particle) {
    this->localParticles.push_back(particle);
}

bool Cell::isPositionInCell(const Vector2D& position) {
    if (ll.x <= position.x && ll.y <= position.y && position.x <= ur.x && position.y <= ur.y)
        return true;
    return false;
}