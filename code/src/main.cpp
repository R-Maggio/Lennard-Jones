#include <iostream>
#include <string>
#include <map>
#include "../include/LJSimulation.hpp"
#include "../include/typeDefinitions.hpp"

int main(int, char**) {

    std::string pathFolder = "../../printParticles/data/";
    std::string fileName = "_particles.csv";

    real_t sigma = 0.5;
    real_t mass = 1.;
    real_t eps = 1e-6;
    real_t dc = 2.5 * sigma;
    Vector2D domainSize = {15., 15.};
    Vector2D gridSize = {10, 10};
    auto constantAcceleration = Vector2D(0., 0.);
    auto constantForce = Vector2D(0., 0.);
    LJBoundary boundary = LJBoundary::PERIODIC;

    LJSimulation simulation = LJSimulation(sigma, mass, eps, dc, domainSize, gridSize, boundary, constantForce, constantAcceleration);
    unsigned int nbParticles = 300;
    simulation.placeRandomParticles(nbParticles, 0., 0., 0.5, 0.5, 0.5);
    // simulation.placeRandomParticles(500);
    constexpr int steps = 9000;
    const real_t dt = 0.001;
    real_t Dt = dt * steps;
    unsigned int nbFrames = 1000;
    unsigned int stepFrame = steps  / nbFrames;
    for (int i = 0; i < steps; i++)
    {
        simulation.computeStep(dt);
        if (i%stepFrame == 0)
            simulation.exportParticlesCSV(pathFolder + std::to_string(i/stepFrame) + fileName);
    }
    real_t dynamicViscosity = simulation.computeDynamicViscosity(Dt);
    // print viscosity:
    std::cout << "viscosity: " << dynamicViscosity << '\n';
    std::cout << "mean kinetic: " << simulation.computeMeanKineticEnergy() << '\n';
    std::cout << "rho: " << simulation.computeDensity() << '\n';
    std::cout << "number of particles: " << simulation.getNumberOfParticles() << '\n';
    simulation.exportConfigJSON(pathFolder + "config.json", {{"Dt", std::to_string(Dt)}});
}
