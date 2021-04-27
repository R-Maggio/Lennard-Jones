#include <iostream> // std::cout
#include <string> // std::string

#include "Vector2D.hpp"

int main(int argc, char const *argv[])
{
    //! TEST of Vector2D:
    
    std::cout << "Vector2D vector = Vector2D(1, 2);" << '\n';
    Vector2D vector = Vector2D(1, 2);
    std::cout << vector << '\n';

    std::cout << "Vector2D(vector);" << '\n';
    std::cout << Vector2D(vector) << '\n';

    std::cout << "vector;" << '\n';
    std::cout << vector << '\n';

    std::cout << "vector.set(4, 10);" << '\n';
    vector.set(4, 10);
    std::cout << vector << '\n';

    std::cout << "vector = Vector2D(2, 3)" << '\n';
    vector = Vector2D(2, 3);
    std::cout << vector << '\n';

    std::cout << "vector + Vector2D(1, 1)" << '\n';
    std::cout << vector + Vector2D(1, 1) << '\n';

    std::cout << "vector - Vector2D(1, 1)" << '\n';
    std::cout << vector - Vector2D(1, 1) << '\n';

    std::cout << "vector += Vector2D(1, 1);" << '\n';
    vector += Vector2D(1, 1);
    std::cout << vector << '\n';

    std::cout << "vector -= Vector2D(1, 1);" << '\n';
    vector -= Vector2D(1, 1);
    std::cout << vector << '\n';

    std::cout << "vector.dot(Vector2D(2, 2))" << '\n';
    std::cout << vector.dot(Vector2D(2, 2)) << '\n';

    std::cout << "0.5*vector" << '\n';
    std::cout << 0.5*vector << '\n';

    std::cout << "vector*0.5" << '\n';
    std::cout << vector*0.5 << '\n';

    vector *= 0.5;
    std::cout << "vector *= 0.5" << '\n';
    std::cout << vector << '\n';

    std::cout << "Vector2D(-12.7, 1500.9) % 10.5" << '\n';
    std::cout << Vector2D(-12.7, 1500.9) % 10.5 << '\n';

    std::cout << "vector.norm()" << '\n';
    std::cout << vector.norm() << '\n';

    std::cout << "Vector2D::norm(vector)" << '\n';
    std::cout << Vector2D::norm(vector) << '\n';

    std::cout << "Vector2D::squaredNorm(vector)" << '\n';
    std::cout << Vector2D::squaredNorm(vector) << '\n';

    std::cout << "vector.squaredNorm()" << '\n';
    std::cout << vector.squaredNorm() << '\n';

    std::cout <<  "Vector2D(2, 3).PeriodicDiff(Vector2D(8, 9), Vector2D(10., 10.))" << '\n';
    std::cout <<  Vector2D(2, 3).PeriodicDiff(Vector2D(8, 9), Vector2D(10., 10.)) << '\n';

    std::cout <<  "Vector2D(8, 9).PeriodicDiff(Vector2D(2, 3), Vector2D(10., 10.))" << '\n';
    std::cout <<  Vector2D(8, 9).PeriodicDiff(Vector2D(2, 3), Vector2D(10., 10.)) << '\n';

    std::cout << "Vector2D::dist(Vector2D(1, 2), Vector2D(5, 4))" << '\n';
    std::cout << Vector2D::dist(Vector2D(1, 2), Vector2D(5, 4)) << '\n';

    std::cout << "Vector2D::squaredDist(Vector2D(1, 2), Vector2D(5, 4))" << '\n';
    std::cout << Vector2D::squaredDist(Vector2D(1, 2), Vector2D(5, 4)) << '\n';

    return 0;
}