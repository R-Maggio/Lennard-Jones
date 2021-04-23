// class Particle
#pragma once

#include "typeDefinitions.hpp"
#include "Vector2D.hpp"

class Particle
{
private:
    Vector2D position {0., 0.};
    Vector2D velocity {0., 0.};
    Vector2D acceleration {0., 0.};
    //! Warning: mass is currently unused because we use dimensionless units
    real_t mass {1.};
public:
    Particle(const Vector2D& position, const Vector2D& velocity, const Vector2D& acceleration, real_t mass):
        position{position}, velocity{velocity}, mass{mass}, acceleration{acceleration} {}
    Particle(const Vector2D& position, const Vector2D& velocity, const Vector2D& acceleration):
        Particle(position, velocity, acceleration, 1.) {}

    // getters and setters:

    Vector2D getPosition();
    void setPosition(const Vector2D& position);

    Vector2D getVelocity();
    void setVelocity(const Vector2D& velocity);

    Vector2D getAcceleration();
    void setAcceleration(const Vector2D& acceleration);
    
    real_t getMass();
    void setMass(real_t mass);

    // other functions
    /**
     * @brief return the eucledean distance between 2 particles a and b (static function)
     * .This function is just a wrapper.
     * 
     * @param a 
     * @param b 
     * @return real_t 
     */
    static real_t dist(const Particle& a, const Particle& b);

    /**
     * @brief return the squared distance between 2 particles a and b (static function)
     * .This function is just a wrapper.
     * 
     * @param a 
     * @param b 
     * @return real_t 
     */
    static real_t squaredDist(const Particle& a, const Particle& b);

    /**
     * @brief compute the kinetic energy of a particle
     * 
     * @return real_t 
     */
    real_t computeKineticEnergy();

    /**
     * @brief overload cout to print a particle
     * 
     * @param os 
     * @param particle 
     * @return std::ostream& 
     */
    friend std::ostream& operator<<(std::ostream& os, const Particle& particle);

};
