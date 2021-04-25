#include <iostream>
#include <string>
#include "../include/LJSimulation.hpp"

int main(int, char**) {

    std::string pathFolder = "../../printParticles/data/";
    std::string fileName = "_particles.csv";

    LJSimulation simulation = LJSimulation(1., 1., 1., 10., 10., 100);
    simulation.placeRandomParticles(1000);
    for (int i = 0; i < 100; i++)
    {
        simulation.computeStep(0.01);
        simulation.exportParticlesCSV(pathFolder + std::to_string(i) + fileName);
    }
}
