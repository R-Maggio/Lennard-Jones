#pragma once

#include <vector>
#include <memory> // std::unique_ptr

#include "typeDefinitions.hpp"
#include "Vector2D.hpp"
#include "Particle.hpp"
#include "Cell.hpp"
#include "Grid.hpp"

/**
 * @brief class for the lennard-Jones potential simulation
 * 
 */
class LJSimulation
{
private:
    // list of all the particles in the simulation:
    std::vector<Particle> particles;
    // pointer on the grid:
    std::unique_ptr<Grid> grid;

    //* ----------------  Simulation parameters ("true" units):
    // duration of the simulations:
    //! simulation duration and nbSteps to set in the start function.
    real_t simulationDuration;
    // number of steps:
    unsigned int nbSteps;
    // sigma:
    real_t sigma;
    // mass of a particle:
    real_t mass;
    // eps:
    real_t eps;
    // cut-off distance:
    real_t dc;
    // domain size:
    real_t domainSize;

    // grid parameters:
    // size of the grid:
    unsigned int gridSize;
    //* ----------------
    //* Simulation parameters in dimensionless units (reduced units):
    /**
     * With the reduced units:
     * mass_r = 1
     * sigma_r = 1
     * eps_r = 1
     * domainSize_r = domainSize/sigma
     * time_r = time*sqrt(eps/(m*sigma^2))
     * 
    */
    real_t simulationDuration_r;
    real_t dc_r;
    real_t domainSize_r;
    //*------------------
    
    //* private functions:
    /**
     * @brief places all the particles in the grid.
     * 
     */
    void placeAllParticles();
    //*------------------

public:

    LJSimulation(real_t sigma, real_t mass, real_t eps, real_t dc, real_t domainSize, unsigned int gridSize);

    //TODO: add getters and setters

    // functions:

    /**
     * @brief compute the temperature inside the domaine
     * 
     * @return real_t 
     */
    real_t ComputeTemperature();

    /**
     * @brief compute the pressure inside the domain
     * 
     * @return real_t 
     */
    real_t ComputePressure();

    /**
     * @brief compute the new acceleration of all the particles.
     * 
     */
    void computeAccelerations();

    /**
     * @brief start the simulation
     * 
     */
    // void start();
};
