#pragma once

#include <vector>
#include <memory> // std::unique_ptr
#include <string>
#include <map>

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
    // time step = nbSteps/simulationDuration:
    real_t dt;
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

    // we can also apply a constant force on the particles:
    Vector2D constantForce;
    // or/and a constant acceleration:
    Vector2D constantAcceleration;
    
    // type of boundary:
    LJBoundary boundary;


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

    LJSimulation(real_t sigma, real_t mass, real_t eps, real_t dc, real_t domainSize, unsigned int gridSize, LJBoundary boundary = LJBoundary::PERIODIC, const Vector2D& constantForce = Vector2D(0., 0.), const Vector2D& constantAcceleration = Vector2D(0., 0.));
    //TODO: domainSize and gridSize should be vectors
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
     * @brief compute the next step of the simulation and update the velocities using the leap frog algorithm
     * 
     */
    void computeStep();

    /**
     * @brief compute the next step of the simulation and update the velocities using the leap frog algorithm
     * 
     * @param dt 
     */
    void computeStep(real_t dt);

    /**
     * @brief start the simulation
     * 
     * @param simulationDuration 
     * @param nbSteps
     */
    void start(real_t simulationDuration, unsigned int nbSteps);

    /**
     * @brief export the particles to an csv file
     * 
     * @param path 
     */
    void exportParticlesCSV(std::string path);

    /**
     * @brief export the configuration
     * 
     * @param path
     * @param customParam
     */
    void exportConfigJSON(std::string path, const std::map<std::string, std::string>& customParam = {});

    /**
     * @brief place randomly n particles inside the grid
     * 
     * @param nbParticles 
     */
    void placeRandomParticles(unsigned int nbParticles);

    void placeRandomParticles(unsigned int nbParticles, real_t muX, real_t muY, real_t sigmaX, real_t sigmaY, real_t corr);
};
