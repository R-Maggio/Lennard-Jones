#include <iostream>
#include <string>
#include <map>
#include "../include/LJSimulation.hpp"
#include "../include/typeDefinitions.hpp"

int main(int argc, char **argv) {

    if (argc > 1) {
        std::string input(argv[1]);
        if (input == "pineq")
        {
            std::string pathFolder = "../../pythonDataVisualization/dataPiNeq/";
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

            real_t maxCorr = 0.5; // maximum value for the correlation
            int nbOfSamples = 101;

            for (int i = 0; i < nbOfSamples; i++)
            {
                real_t corr = maxCorr * i / (nbOfSamples - 1);

                LJSimulation simulation = LJSimulation(sigma, mass, eps, dc, domainSize, gridSize, boundary, constantForce, constantAcceleration);
                unsigned int nbParticles = 225;
                simulation.placeRandomParticles(nbParticles, 0., 0., 0.25, 0.25, corr);
                simulation.exportParticlesCSV(pathFolder + std::to_string(i) + "_" + "in" + fileName);
                constexpr int steps = 20000;
                const real_t dt = 0.0005;
                real_t Dt = dt * steps;
                // start the simulation:
                simulation.start(Dt, steps);
                simulation.exportParticlesCSV(pathFolder + std::to_string(i) + "_" + "out" + fileName);

                real_t dynamicViscosity = simulation.computeDynamicViscosity(Dt);
                simulation.exportConfigJSON(pathFolder + std::to_string(i) + "_" + "config.json", {{"Dt", std::to_string(Dt)}, {"correlation", std::to_string(corr)}, {"visc", std::to_string(dynamicViscosity)}});
            }
            return 0;
        }          
    }

    std::string pathFolder = "../../pythonDataVisualization/data/";
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
    unsigned int nbParticles = 225;
    simulation.placeRandomParticles(nbParticles, 0., 0., 0.25, 0.25, 0.35);
    // simulation.placeRandomParticles(500);
    constexpr int steps = 20000;
    const real_t dt = 0.0005;
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

    return 0;
}
