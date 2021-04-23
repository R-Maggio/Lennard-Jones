#include <iostream> // std::cout
#include <string> // std::string

#include "Vector2D.hpp"
#include "Particle.hpp"

int main(int argc, char const *argv[])
{
    //! TEST of the Particle class:
    
    Vector2D r = Vector2D(1, -2);
    Vector2D v = Vector2D(1, 1);
    Vector2D a = Vector2D(-1, 2);
    real_t m = 1.;

    Particle p1(r, v, a, m);
    std::cout << p1 << '\n';

    std::cout << "p1 = Particle(r*2, v*2, a*2, m*2);" << '\n';
    p1 = Particle(r*2, v*2, a*2, m*2);
    std::cout << p1 << '\n';

    auto p2 = p1;
    std::cout << "auto p2 = p1.copy();" << '\n';
    std::cout << p2 << '\n';

    std::cout << "p1.getPosition(), p1.getVelocity(), p1.getAcceleration(), p1.getMass()" << '\n';
    std::cout << p1.getPosition() << p1.getVelocity() << p1.getAcceleration() << p1.getMass() << '\n';

    p2.setAcceleration(a);
    p2.setPosition(r);
    p2.setVelocity(v);
    p2.setMass(m);
    std::cout << "p2.setAcceleration(a); p2.setPosition(r); p2.setVelocity(v); p2.setMass(m);" << '\n';
    std::cout << "p2:" << p2 << '\n';

    std::cout << "p1, p2:" << p1 << ", " << p2 << '\n';
    std::cout << "Particle::dist(p1, p2)" << '\n';
    std::cout << Particle::dist(p1, p2) << '\n';

    std::cout << "Particle::squaredDist(p1, p2)" << '\n';
    std::cout << Particle::squaredDist(p1, p2) << '\n';

    std::cout << "p1.computeKineticEnergy()" << '\n';
    std::cout << p1.computeKineticEnergy() << '\n';

    return 0;
}