
#include <cmath> // sqrt
#include <iomanip> // std::setprecision

#include "Vector2D.hpp"

// setters:

void Vector2D::set(real_t x, real_t y) {
    this->x = x;
    this->y = y;
}

//sum of 2 vectors and returns a new vector

Vector2D Vector2D::operator+(const Vector2D& otherVector) const {
    return Vector2D(this->x + otherVector.x, this->y + otherVector.y);
}

// substraction:

Vector2D Vector2D::operator-(const Vector2D& otherVector) const {
    return Vector2D(this->x - otherVector.x, this->y - otherVector.y);
}
// add a vector to the current vector

Vector2D& Vector2D::operator+=(const Vector2D& otherVector) {
    this->x += otherVector.x;
    this->y += otherVector.y;
    return *this;
}

// -= operator overloading

Vector2D& Vector2D::operator-=(const Vector2D& otherVector) {
    this->x -= otherVector.x;
    this->y -= otherVector.y;
    return *this; 
}

// return the dot product

real_t Vector2D::dot(const Vector2D& otherVector) const {
    return this->x * otherVector.x + this->y * otherVector.y;
}

// return result of multiplying a vector by a scalar

Vector2D operator*(real_t scalar, const Vector2D& vector) {
    return Vector2D(vector.x * scalar, vector.y * scalar);
}
Vector2D operator*(const Vector2D& vector, real_t scalar) {
    return Vector2D(vector.x * scalar, vector.y * scalar);
}

//. *= operator overlading. multiply by scalar

Vector2D& Vector2D::operator*=(real_t scalar) {
    this->x *= scalar;
    this->y *= scalar;
    return *this; 
}

// return the true modulo operation

Vector2D operator%(const Vector2D& vector, real_t scalar) {
    const real_t x = std::fmod(vector.x, scalar);
    const real_t y = std::fmod(vector.y, scalar);
    x >= 0 ? x : x + scalar;
    y >= 0 ? y : y + scalar;
    return Vector2D(x, y);
}

// return the euclidean norm of a vector (static function)

real_t Vector2D::norm() const {
    return std::sqrt(this->x*this->x + this->y*this->y);
}

// return the euclidean norm of a vector

real_t Vector2D::norm(const Vector2D& vector) {
    return std::sqrt(vector.x*vector.x + vector.y*vector.y);
}

// return the absolute squared norm of a vector

real_t Vector2D::squaredNorm() const {
    return this->x*this->x + this->y*this->y;
}

// return the absolute squared norm of a vector (static function)

real_t Vector2D::squaredNorm(const Vector2D& vector) {
    return vector.x*vector.x + vector.y*vector.y;
}

// return the euclidean distance between 2 vectors a and b (static function)

real_t Vector2D::dist(const Vector2D& a, const Vector2D& b) {
    return Vector2D::norm(a - b);
}

// return the squared distance between 2 vectors a and b (static function)

real_t Vector2D::squaredDist(const Vector2D& a, const Vector2D& b) {
    return Vector2D::squaredNorm(a - b);
}

// overload cout to print a vector

std::ostream& operator<<(std::ostream& os, const Vector2D& vector) {
    os << "{x:" << std::setprecision(5) << vector.x << ", y:" << vector.y << "}";
    return os;
}