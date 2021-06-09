
#include <cmath> //std::sqrt, std::pow, std::fmod
#include <memory> // std::unique_ptr
#include <fstream> // write to file
#include <iostream>
#include <string>
#include <ctime> // srand, rand
#include <random> // std::normal_distribution
#include "LJSimulation.hpp"
#include "Cell.hpp"
#include "constExpr.hpp"

LJSimulation::LJSimulation(real_t sigma, real_t mass, real_t eps, real_t dc, const Vector2D& domainSize, const Vector2D& gridSize, LJBoundary boundary, const Vector2D& constantForce, const Vector2D& constantAcceleration, const Vector2D& maxVel, real_t simulationDuration) {
    this->sigma = sigma;
    this->mass = mass;
    this->eps = eps;
    this->dc = dc;
    this->domainSize = domainSize;
    this->gridSizeX = gridSize.x;
    this->gridSizeY = gridSize.y;
    this->boundary = boundary;
    this->constantForce = constantForce;
    this->constantAcceleration = constantAcceleration;
    this->maxVel = maxVel;
    this->simulationDuration = simulationDuration; // optional, can be 0
    // TODO: pass particle list
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
    this->grid = std::make_unique<Grid>(this->domainSize, gridSize, boundary, dc);

    //TODO: call an init function for the particles (random for example)
    // places all the particles in the grid:
    placeAllParticles(); //! might be called in another place
}

LJSimulation::LJSimulation(real_t sigma, real_t mass, real_t eps, real_t dc, real_t domainSize, unsigned int gridSize, LJBoundary boundary, const Vector2D& constantForce, const Vector2D& constantAcceleration, const Vector2D& maxVel, real_t simulationDuration):
    LJSimulation(sigma, mass, eps, dc, Vector2D(domainSize, domainSize), Vector2D(gridSize, gridSize), boundary, constantForce, constantAcceleration, maxVel, simulationDuration) { // call the other constructor
}

void LJSimulation::placeAllParticles() {
    for (auto &p : particles) // puts all the particles in the grid
    {
        grid->placeParticle(&p);
    }
    this->numberOfParticles = particles.size();
}

Vector2D LJSimulation::computeHelfandMoment() {
    Vector2D Gt = -viscosity_correction;
    for (auto &p : particles)
    {
        Gt += p.getMass() * p.getVelocity() * p.getPosition().flip();
    }
    return Gt;
    
}

real_t LJSimulation::computeMeanKineticEnergy() {
    real_t E = 0.;
    for (auto &p : particles) {
        E += p.getMass() * p.getVelocity().squaredNorm();
    }
    return 0.5 * E/particles.size();
}

real_t LJSimulation::ComputeTemperature() {
    auto E = computeMeanKineticEnergy(); // Mean of the kineticEnergy
    /**
     * E = 3/2 Kb T => T = 2 E /(3*Kb)
    */
    return 2*E/(3*Kb);
}

real_t LJSimulation::ComputePressure() {
    // PV = nRT => P = nRT/V
    return numberOfParticles * R_ * ComputeTemperature()/(domainSize.x*domainSize.y);
}

void LJSimulation::computeAccelerations() {
    // TODO: make 2 different fct: 1 if all particles are the same, another if not.
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
                        const auto vector_r = otherParticle->getPosition().PeriodicDiff(particle->getPosition(), domainSize);
                        const real_t r_2 = vector_r.squaredNorm(vector_r); // r^2
                        
                        if (r_2 <= dc) { // test if the distance is smaller than the cut-off distance
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
                            // to compute the force vector, we can use the normalize vector distance:
                            // vector_f = vector_r/r * f = vector_r * 4 * (0.5 * r_inv6 - r_inv6^2) * r_inv^2;
                            // example: rx/r = fx/f => fx = rx * f/r
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
    
}

void LJSimulation::computeStep(real_t dt) {
    
    /* the leap frog algorithm:
    * x(t+dt) = x(t) + dt * v(t + 1/2 dt)
    * v(t + 3/2 dt) = v(t + 1/2 dt) + dt * a(t + dt)
    */

    for (auto &p : particles)
    {
        // update the position
        if (this->boundary == LJBoundary::PERIODIC) { // Periodic boundaries
            const auto newPosition = p.getPosition() + dt * p.getVelocity();
            p.setPosition(newPosition % this->domainSize);
            // compute the size of the jump in the x direction
            // add the jump to the correction term:
            viscosity_correction += p.getMass() * p.getVelocity() * (p.getPosition() - newPosition).flip();
        } else if (this->boundary == LJBoundary::POISSEUILLE) { // Periodic boundaries
            const auto pxNoMod = p.getPosition().x  + dt * p.getVelocity().x; // new position without the modulus
            const auto px = std::fmod(std::fmod(pxNoMod, domainSize.x) + domainSize.x, domainSize.x);
            auto py = p.getPosition().y + dt * p.getVelocity().y;
            if (py > domainSize.y) {
                py = domainSize.y;
                // reverse velocity:
                p.setVelocity(-p.getVelocity());
            } else if (py < 0.) {
                py = 0.;
                // reverse velocity:
                p.setVelocity(-p.getVelocity()); 
            }
            p.setPosition({px, py});
            viscosity_correction += p.getMass() * p.getVelocity() * Vector2D(0., px - pxNoMod);
        }
        // reset the acceleration
        p.setAcceleration(this->constantAcceleration + (this->constantForce / p.getMass()));
    }

    // compute the new acceleration:
    computeAccelerations();

    // update the velocities
    for (auto &p : particles) {
        auto newVelocity = p.getVelocity() + dt * p.getAcceleration();
        newVelocity.x = newVelocity.x > maxVel.x ? maxVel.x : newVelocity.x;
        newVelocity.y = newVelocity.y > maxVel.y ? maxVel.y : newVelocity.y;

        newVelocity.x = newVelocity.x < -maxVel.x ? -maxVel.x : newVelocity.x;
        newVelocity.y = newVelocity.y < -maxVel.y ? -maxVel.y : newVelocity.y;

        p.setVelocity(newVelocity);
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
        file << "m,px,py,vx,vy,ax,ay\n";
        for (auto &p : particles)
        {
            file << p.getMass() << ','
                 << p.getPosition().x << ',' << p.getPosition().y << ','
                 << p.getVelocity().x << ',' << p.getVelocity().y << ','
                 << p.getAcceleration().x << ',' << p.getAcceleration().y << '\n';
        }
        
    } else {
        std::cerr << "Error: cannot open the file: " << path << '\n';
    }
    
    // close the file:
    file.close();
}

void LJSimulation::exportConfigJSON(std::string path, const std::map<std::string, std::string>& customParam) {
    // open the file or create it
    std::ofstream file;
    file.open(path);
    if (file.is_open()) {
        // csv header:
        file << "{\n"
             << "\"domainSizeX\": " << this->domainSize.x << ",\n"
             << "\"domainSizeY\": " << this->domainSize.y << ",\n"
             << "\"sigma\": " << this->sigma << ",\n";
        // import custom parameters:
        for (auto const& param : customParam)
        {
            file << '"' << param.first << '"'
                << ':'
                << param.second // param value
                << ',';
        }
        file << "\n}";
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
        real_t px = (i % (int) sqrt_nb + 0.5) * this->domainSize.x/sqrt_nb;
        real_t py = (i / (int) sqrt_nb + 0.5) * this->domainSize.y/sqrt_nb;
        particles.push_back(Particle(Vector2D(px, py), Vector2D(vx, vy), this->constantAcceleration + (this->constantForce / this->mass), this->mass));
    }
    // put the particles inside the Grid:
    placeAllParticles();

    // compute G0:
    G0 = computeHelfandMoment();
}

void LJSimulation::placeRandomParticles(unsigned int nbParticles, real_t muX, real_t muY, real_t sigmaX, real_t sigmaY, real_t corr) {
    // Print warning if velocity will be to big:
    if (simulationDuration > 0.) {
        // warning if 2 sigma +- mu > 0.01 * Dx/Dt
        auto AbsVelX = std::abs(muX) + 2*sigmaX;
        auto AbsVelY = std::abs(muY) + 2*sigmaY;
        if (0.01 * domainSize.x/simulationDuration < AbsVelX) {
            std::cerr << "\033[31m" << "Warning: particles velocity along X axis is too big" << "\033[0m" << "\n";
        }
        if (0.01 * domainSize.y/simulationDuration < AbsVelY) {
            std::cerr << "\033[31m" << "Warning: particles velocity along Y axis is too big" << "\033[0m" << "\n";
        }
    }

    std::default_random_engine randomGenerator;
    std::normal_distribution normalDistribution(0.0,1.0);
    for (size_t i = 0; i < nbParticles; i++) {
        real_t vx = muX + sigmaX * normalDistribution(randomGenerator);
        real_t vy;
        if (sigmaX == 0) {
            // TODO: check if correct
            vy = muY + sigmaY * normalDistribution(randomGenerator);
        } else {
            vy = (sigmaY/sigmaX) * corr * (vx - muX) + std::sqrt((1-corr*corr)*sigmaY*sigmaY) * normalDistribution(randomGenerator);
        }
        auto sqrt_nb = std::ceil(std::sqrt(nbParticles));
        real_t px = ((i % (int) sqrt_nb + 0.5 + 0.5 * ((i/(int) sqrt_nb)%2)) * this->domainSize.x/sqrt_nb);
        real_t py = (i / (int) sqrt_nb + 0.5) * this->domainSize.y/sqrt_nb;
        particles.push_back(Particle(Vector2D(px, py)%this->domainSize.x, Vector2D(vx, vy), this->constantAcceleration + (this->constantForce / this->mass), this->mass));
    }
    // put the particles inside the Grid:
    placeAllParticles();

    // compute G0:
    G0 = computeHelfandMoment();
}


real_t LJSimulation::computeDensity() {
    real_t rho = 0.;
    for (auto &p : particles) {
        rho += p.getMass();
    }
    return rho/(domainSize.x * domainSize.y);
}

real_t LJSimulation::computeDynamicViscosity(real_t t, bool positionDirectionX) {
    if (t == 0) t = simulationDuration;
    real_t coef = 3/(4 * computeDensity() * domainSize.x * domainSize.y * computeMeanKineticEnergy() * (particles.size() - 1) * t);
    auto Gt = computeHelfandMoment();
    if (positionDirectionX) {
        return coef * (Gt.y - G0.y)*(Gt.y - G0.y); // Gt.y because Gt.y contains Gyx(t) => x position
    }
    return coef * (Gt.x - G0.x)*(Gt.x - G0.x);
}