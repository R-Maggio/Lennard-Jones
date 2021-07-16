#include <iostream> // std::cout
#include <string> // std::string
#include <random> // std::normal_distribution
#include <vector>

#include "../include/Vector2D.hpp"
#include "../include/Particle.hpp"
#include "../include/Grid.hpp"
#include "../include/LJSimulation.hpp"


int main(int argc, char const *argv[])
{
    //! TEST of the LJSimulation class:
    // test 2D random distribution:
    int nbParticles = 3600;
    auto muX = 0.1;
    auto muY = -0.2;
    auto sigmaX = 0.3;
    auto sigmaY = 0.6;
    auto corr = 0.9;
    std::vector<real_t> Vx;
    std::vector<real_t> Vy;

    std::default_random_engine randomGenerator;
    std::normal_distribution<real_t> normalDistribution(0.0,1.0);
    for (size_t i = 0; i < nbParticles; i++) {
        real_t vx = muX + sigmaX * normalDistribution(randomGenerator);
        real_t vy;
        if (sigmaX == 0) {
            // TODO: check if correct
            vy = muY + sigmaY * normalDistribution(randomGenerator);
        } else {
            vy = (sigmaY/sigmaX) * corr * (vx - muX) + std::sqrt((1-corr*corr)*sigmaY*sigmaY) * normalDistribution(randomGenerator);
        }
        Vx.push_back(vx);
        Vy.push_back(vy);
    }

    std::cout << '[';
    for (const auto &v : Vx) {
        std::cout << v << ',';
    }
    std::cout << "]\n";
    
    std::cout << '[';
    for (const auto &v : Vy) {
        std::cout << v << ',';
    }
    std::cout << "]\n";
}