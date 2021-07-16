#include <iostream>
#include <string>
#include <map>
#include <chrono>  // std::chrono
#include <iomanip>      // std::setprecision
#include "../include/LJSimulation.hpp"
#include "../include/typeDefinitions.hpp"

int main(int argc, char **argv) {
    // time measurement
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    if (argc = 2) {
        std::string input(argv[1]);
        if (input == "pineq")
        {
            std::string pathFolder = "../../pythonDataVisualization/dataPiNeq/";
            std::string fileName = "_particles.csv";

            real_t sigma = 0.5;
            real_t mass = 1.;
            real_t eps = 1e-10;
            real_t dc = 2.5 * sigma;
            Vector2D domainSize = {15., 15.};
            Vector2D gridSize = {10, 10};
            auto constantAcceleration = Vector2D(0., 0.);
            auto constantForce = Vector2D(0., 0.);
            LJBoundary boundary = LJBoundary::PERIODIC;
            unsigned int nbParticles = 3600;
            constexpr int steps = 2000;
            const real_t dt = 0.0005;
            real_t Dt = dt * steps;
            Vector2D maxVel = 0.5 * domainSize / dt;
            // real_t sigmaX = 0.005 * domainSize.x/Dt;
            // real_t sigmaY = 0.005 * domainSize.y/Dt;
            real_t sigmaX = 0.005 * domainSize.x;
            real_t sigmaY = 0.005 * domainSize.y;
            real_t maxCorr = 0.5; // maximum value for the correlation
            int nbOfSamples = 101;

            // #pragma omp parallel for schedule(static,1)
            #pragma omp parallel for
            for (int i = 0; i < nbOfSamples; i++)
            {
                real_t corr = maxCorr * i / (nbOfSamples - 1);

                LJSimulation simulation = LJSimulation(sigma, mass, eps, dc, domainSize, gridSize, boundary, constantForce, constantAcceleration, maxVel, Dt);
                simulation.placeRandomParticles(nbParticles, 0., 0., sigmaX, sigmaY, corr);
                simulation.exportParticlesCSV(pathFolder + std::to_string(i) + "_" + "in" + fileName);
                // start the simulation:
                simulation.start(Dt, steps);
                simulation.exportParticlesCSV(pathFolder + std::to_string(i) + "_" + "out" + fileName);

                real_t dynamicViscosity = simulation.computeDynamicViscosity(Dt);
                simulation.exportConfigJSON(pathFolder + std::to_string(i) + "_" + "config.json", {{"Dt", std::to_string(Dt)}, {"correlation", std::to_string(corr)}, {"visc", std::to_string(dynamicViscosity)}, {"steps", std::to_string(steps)}});
                // Print percentage:
                // std::cout << "\rprogress: " << std::setprecision(1) << 100*(i+1)/nbOfSamples << "% "; 
            }
            std::cout << '\n';
            return 0;
        }
            // error
            return 1;        
    } else {
        std::string pathFolder = "../../pythonDataVisualization/data/";
        std::string fileName = "_particles.csv";

        real_t sigma = 0.5;
        real_t mass = 1.;
        // real_t sigma = 0.005;
        // real_t mass = 0.0002784256472222222;
        real_t eps = 1e-10;
        real_t dc = 2.5 * sigma;
        Vector2D domainSize = {15, 15};
        Vector2D gridSize = {10, 10};
        auto constantAcceleration = Vector2D(0., 0.);
        auto constantForce = Vector2D(0., 0.);
        LJBoundary boundary = LJBoundary::PERIODIC;
        unsigned int nbParticles = 3600;
        // simulation.placeRandomParticles(500);
        constexpr int steps = 2300;
        const real_t dt = 0.005;
        real_t Dt = dt * steps;
        Vector2D maxVel = {10., 10.};
        unsigned int nbFrames = 230;
        unsigned int stepFrame = steps  / nbFrames;
        // real_t sigmaX = 0.5773866671701394;
        // real_t sigmaY = 0.5772862749939868;
        // real_t corr = -1.5574669947325287e-01;
        // real_t muX = 0.04846295;
        // real_t muY = -0.00560023;
        real_t sigmaX = 0.005;
        real_t sigmaY = 0.005;
        real_t corr = 0.5;
        real_t muX = 0;
        real_t muY = 0;
    
        // initialize the simulation
        LJSimulation simulation = LJSimulation(sigma, mass, eps, dc, domainSize, gridSize, boundary, constantForce, constantAcceleration, maxVel, Dt);
        simulation.placeRandomParticles(nbParticles, muX, muY, sigmaX, sigmaY, corr);

        simulation.exportParticlesCSV(pathFolder + std::to_string(0) + fileName);
        for (int i = 1; i < steps; i++)
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
        simulation.exportConfigJSON(pathFolder + "config.json", {{"Dt", std::to_string(Dt)}, {"correlation", std::to_string(corr)}, {"visc", std::to_string(dynamicViscosity)}, {"steps", std::to_string(steps)}});

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000. << " [s]" << std::endl;

        return 0;
    }
}
