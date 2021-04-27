#include <iostream>
#include <string>
#include "../include/LJSimulation.hpp"

int main(int, char**) {

    std::string pathFolder = "../../printParticles/data/";
    std::string fileName = "_particles.csv";

    LJSimulation simulation = LJSimulation(0.3, 10, 0.1, 10., 15., 10);
    simulation.placeRandomParticles(500);
    constexpr int steps = 1000;
    for (int i = 0; i < steps; i++)
    {
        simulation.computeStep(0.001);
        simulation.exportParticlesCSV(pathFolder + std::to_string(i) + fileName);
    }
    simulation.exportConfigJSON(pathFolder + "config.json");
}
