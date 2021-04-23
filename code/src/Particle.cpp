#include "Particle.hpp"

// getters and setters:

Vector2D Particle::getPosition() {
    return this->position;
}

void Particle::setPosition(const Vector2D& position) {
    this->position = position;
}

Vector2D Particle::getVelocity() {
    return this->velocity;
}

void Particle::setVelocity(const Vector2D& velocity) {
    this->velocity = velocity;
}

Vector2D Particle::getAcceleration() {
    return this->acceleration;
}

void Particle::setAcceleration(const Vector2D& acceleration) {
    this->acceleration = acceleration;
}

real_t Particle::getMass() {
    return mass;
}

void Particle::setMass(real_t mass) {
    this->mass = mass;
}

// return the eucledean distance between 2 particles a and b

real_t Particle::dist(const Particle& a, const Particle& b) {
    return Vector2D::dist(a.position, b.position);
}

// return the squared distance between 2 particles a and b

real_t Particle::squaredDist(const Particle& a, const Particle& b) {
    return Vector2D::squaredDist(a.position, b.position);
}

// compute the kinetic energy of a particle:

real_t Particle::computeKineticEnergy() {
    return 0.5 * mass * this->getVelocity().squaredNorm(); // 1/2 m v^2
}

// overload cout to print a particle

std::ostream& operator<<(std::ostream& os, const Particle& particle) {
    os << "[r:" << particle.position << ", v:" << particle.velocity
       << ", a:" << particle.acceleration << ", m:" << particle.mass << "]";
    return os;
}