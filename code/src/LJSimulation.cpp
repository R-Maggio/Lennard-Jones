
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
    // allocate memory for the grid:
    this->grid = std::make_unique<Grid>(this->domainSize, gridSize);

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
     * E = 3/2 Kb T => T = 2 E /(3*Kb)
    */
    return 2*E/(3*Kb);
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
            auto nearCells = grid->findCells(particle->getPosition(), dc);
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
                        // with periodic boudaries:
                        // if dx > X/2 => dx = X - dx
                        const auto vector_r = otherParticle->getPosition().PeriodicDiff(particle->getPosition(), Vector2D(this->domainSize, this->domainSize));
                        const real_t r_2 = vector_r.squaredNorm(vector_r); // r^2
                        const real_t r_inv2 = 1./r_2; // r^-2
                        const real_t r_inv6 = std::pow(r_inv2, 3); // r^-6
                        const real_t simga_6 = std::pow(sigma, 6);
                        /*
                         * VLJ = 4(r^-12 - r^-6) : using reduced units
                         * Fij = dV/dr = 4(-12r^-13 + 6r^-7) = 48(-r_inv^12 + 0.5r_inv^6)*r_inv
                         * 
                         * without reduced units:
                         * VLJ = 4 e0 (r^-12 - r^-6) : using reduced units
                         * Fij = dV/dr = 48(-12r^-13 + 6r^-7)
                         */
                        // the force will be: f = 4 * (0.5 * r_inv6 - r_inv6^2) * r_inv;
                        // in our case, f = a since m = 1
                        // to compute the vector of the acceleration, we can use the normalize vector distance:
                        // a = vector_r/r * f = vector_r * 4 * (0.5 * r_inv6 - r_inv6^2) * r_inv^2;
                        // this is why we only need to compute r_inv^2:
                        // f * r_inv is given by:
                        real_t force_r_inv = 48 * eps * (0.5 * r_inv6*simga_6 - r_inv6*r_inv6*simga_6*simga_6) * r_inv2;
                        // set the new accelerations:
                        particle->setAcceleration(particle->getAcceleration() + vector_r * (force_r_inv / particle->getMass()));
                        otherParticle->setAcceleration(otherParticle->getAcceleration() + vector_r * (-force_r_inv / particle->getMass()));
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
        p.setPosition((p.getPosition() + dt * p.getVelocity()) % this->domainSize);
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
    this->dt = this->simulationDuration / this->nbSteps;
    // compute the simulation
    for (size_t i = 0; i < nbSteps; i++)
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

void LJSimulation::exportConfigJSON(std::string path) {
    // open the file or create it
    std::ofstream file;
    file.open(path);
    if (file.is_open()) {
        // csv header:
        file << "{\n"
             << "\"domainSize\": " << this->domainSize
             << "\n}";
    } else {
        std::cerr << "Error: cannot open the file: " << path << '\n';
    }
    
    // close the file:
    file.close();
}

void LJSimulation::placeRandomParticles(unsigned int nbParticles) {
    std::srand(12345);

    for (size_t i = 0; i < nbParticles; i++)
    {
        // real_t px = ((real_t) std::rand()) / (real_t) RAND_MAX * this->domainSize;
        // real_t py = ((real_t) std::rand()) / (real_t) RAND_MAX * this->domainSize;
        real_t vx = (((real_t) std::rand()) / (real_t) RAND_MAX - 0.5) * 5;
        real_t vy = (((real_t) std::rand()) / (real_t) RAND_MAX -0.5) * 5;
        auto sqrt_nb = std::ceil(std::sqrt(nbParticles));
        real_t px = (i % (int) sqrt_nb + 0.5) * this->domainSize/sqrt_nb;
        real_t py = (i / (int) sqrt_nb + 0.5) * this->domainSize/sqrt_nb;
        particles.push_back(Particle(Vector2D(px, py), Vector2D(vx, vy), Vector2D(), this->mass));
    }
    // put the particles inside the Grid:
    placeAllParticles();
}