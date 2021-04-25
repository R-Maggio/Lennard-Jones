
#include <cmath> //std::sqrt, std::pow
#include <memory> // std::unique_ptr
#include <fstream> // write to file
#include <iostream>
#include <string>
#include <ctime> // srand, rand
#include "LJSimulation.hpp"
#include "Cell.hpp"
#include "constExpr.hpp"

LJSimulation::LJSimulation(real_t sigma, real_t mass, real_t eps, real_t dc, real_t domainSize, unsigned int gridSize) {
    this->sigma = sigma;
    this->mass = mass;
    this->eps = eps;
    this->dc = dc;
    this->domainSize = domainSize;
    this->gridSize = gridSize;
    // reduced units:
    /**
     * With the reduced units:
     * mass_r = 1
     * sigma_r = 1
     * eps_r = 1
     * domainSize_r = domainSize/sigma
     * time_r = time*sqrt(eps/(m*sigma^2))
     * 
    */
    this->dc_r = dc/sigma;
    this->domainSize_r = domainSize/sigma;
    // allocate memory for the grid:
    this->grid = std::make_unique<Grid>(this->domainSize_r, gridSize); // WARNING: we are using reduced units for the size of the domain.

    //TODO: call an init function for the particles (random for example)
    // places all the particles in the grid:
    placeAllParticles(); //! might be called in another place
}

void LJSimulation::placeAllParticles() {
    for (auto &p : particles) // puts all the particles in the grid
    {
        grid->placeParticle(&p);
    }
}

real_t LJSimulation::ComputeTemperature() {
    real_t E = 0; // kinetic energy
    for (auto &p : particles)
    {
        E += p.computeKineticEnergy();
    }
    E /= particles.size(); // Mean of the kineticEnergy
    /**
     * E* = E/eps = 3/2 Kb T/eps => T = 2 eps E* / 3Kb
    */
    return 2*eps*E/(3*Kb);
}

real_t LJSimulation::ComputePressure() {
    // PV = nRT => P = nRT/V
    return particles.size()*R_*ComputeTemperature()/(domainSize*domainSize);
}

void LJSimulation::computeAccelerations() {
    // TODO: make 2 different fct: 1 if all particles are the same, another if not.
    // TODO: rewrite the function with the non-reduced form
    // we need to iterate over the cells:
    for (auto &cell : grid->getCells())
    {
        for (auto &particle : cell.getLocalParticles())
        {
            // get the cells where we need to search:
            auto nearCells = grid->findCells(particle->getPosition(), dc_r);
            for (auto &nearCell : nearCells)
            {
                for (auto &otherParticle : nearCell->getLocalParticles())
                {
                    // otherParticle and particle are pointers.
                    // Fij = Fji => we only compute the force once.
                    // We don't want the interaction with itself (< and not <=).
                    if (otherParticle < particle)
                    {
                        // Vector_r is the difference between the positions of the 2 particles
                        const auto vector_r = particle->getPosition() - otherParticle->getPosition();
                        const real_t r_2 = vector_r.squaredNorm(vector_r); // r^2
                        const real_t r_inv2 = 1./r_2; // r^-2
                        const real_t r_inv6 = std::pow(r_inv2, 3); // r^-6
                        /*
                         * VLJ = 4(r^-12 - r^-6) : using reduced units
                         * Fij = dV/dr = 4(-12r^-13 + 6r^-7) = 48(-r_inv^12 + 0.5r_inv^6)*r_inv
                         * 
                         */
                        // the force will be: f = 4 * (0.5 * r_inv6 - r_inv6^2) * r_inv;
                        // in our case, f = a since m = 1
                        // to compute the vector of the acceleration, we can use the normalize vector distance:
                        // a = vector_r/r * f = vector_r * 4 * (0.5 * r_inv6 - r_inv6^2) * r_inv^2;
                        // this is why we only need to compute r_inv^2:
                        // f * r_inv is given by:
                        real_t force_r_inv = 4 * (0.5 * r_inv6 - r_inv6*r_inv6) * r_inv2;
                        // set the new accelerations:
                        particle->setAcceleration(particle->getAcceleration() + vector_r * force_r_inv);
                        otherParticle->setAcceleration(otherParticle->getAcceleration() + vector_r * -force_r_inv);
                    }
                    
                }
                
            }
            
        }
        
    }
    
}

void LJSimulation::computeStep(real_t dt) {
    
    /* the leap frog algorithm:
    * x(t+dt) = x(t) + dt * v(t + 1/2 dt)
    * v(t + 3/2 dt) = v(t + 1/2 dt) + dt * a(t + dt)
    */

    for (auto &p : particles)
    {
        // update the position
        p.setPosition((p.getPosition() + dt * p.getVelocity()) % this->domainSize_r);
        p.setAcceleration(Vector2D(0., 0.)); // reset the acceleration
    }

    // compute the new acceleration:
    computeAccelerations();

    // update the velocities
    for (auto &p : particles)
    {
        p.setVelocity(p.getVelocity() + dt * p.getAcceleration());
    }
    // replace the particles. We need to iterate over the cells:
    //! very unefficient !!!
    for (auto &cell : grid->getCells())
    {
        cell.setLocalParticles(std::vector<Particle *>()); // reset the cells
    }
    this->placeAllParticles(); // replace the particles
}

void LJSimulation::computeStep() {
    computeStep(this->dt);
}

void LJSimulation::start(real_t simulationDuration, unsigned int nbSteps) {
    // set the parameters for the simulation
    this->simulationDuration = simulationDuration;
    this->simulationDuration_r = simulationDuration / sigma;
    this->nbSteps = nbSteps;
    this->dt = this->simulationDuration_r / this->nbSteps;
    // compute the simulation
    for (int i = 0; i < nbSteps; i++)
    {
        computeStep();
    }
}

void LJSimulation::exportParticlesCSV(std::string path) {
    // open the file or create it
    std::ofstream file;
    file.open(path);
    if (file.is_open()) {
        // csv header:
        file << "px,py,vx,vy,ax,ay\n";
        for (auto &p : particles)
        {
            file << p.getPosition().x << ',' << p.getPosition().y << ','
                 << p.getVelocity().x << ',' << p.getVelocity().y << ','
                 << p.getAcceleration().x << ',' << p.getAcceleration().y << '\n';
        }
        
    } else {
        std::cerr << "Error: cannot open the file: " << path << '\n';
    }
    
    // close the file:
    file.close();
}

void LJSimulation::placeRandomParticles(unsigned int nbParticles) {
    std::srand(12345);

    for (int i = 0; i < nbParticles; i++)
    {
        real_t random_px = ((real_t) std::rand()) / (real_t) RAND_MAX * this->domainSize_r;
        real_t random_py = ((real_t) std::rand()) / (real_t) RAND_MAX * this->domainSize_r;
        real_t random_vx = ((real_t) std::rand()) / (real_t) RAND_MAX * 1;
        real_t random_vy = ((real_t) std::rand()) / (real_t) RAND_MAX * 1;
        particles.push_back(Particle(Vector2D(random_px, random_py), Vector2D(random_vx, random_vy), Vector2D()));
    }
    // put the particles inside the Grid:
    placeAllParticles();
}