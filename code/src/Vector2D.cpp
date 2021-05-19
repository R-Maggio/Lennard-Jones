
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

Vector2D Vector2D::operator-() const {
    return Vector2D(-this->x, -this->y);
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

// division
Vector2D operator/(const Vector2D& vector, real_t scalar) {
    return Vector2D(vector.x / scalar, vector.y / scalar);
}

Vector2D operator/(real_t scalar, const Vector2D& vector) {
    return Vector2D(scalar / vector.x, scalar / vector.y);
}

Vector2D& Vector2D::operator/=(real_t scalar) {
    this->x /= scalar;
    this->y /= scalar;
    return *this; 
}

// return the true modulus operation

Vector2D operator%(const Vector2D& vector, real_t scalar) {
    const real_t x = std::fmod(std::fmod(vector.x, scalar) + scalar, scalar);
    const real_t y = std::fmod(std::fmod(vector.y, scalar) + scalar, scalar);
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

Vector2D Vector2D::PeriodicDiff(const Vector2D& otherVector, const Vector2D& domainSize) const {
    // we compute the "true" difference with periodic boundaries
    auto diff = *this - otherVector;
    if (diff.x > domainSize.x/2)
        diff.x = domainSize.x - diff.x;
    else if (diff.x < -domainSize.x/2)
        diff.x = domainSize.x + diff.x;
    if (diff.y > domainSize.y/2)
        diff.y = domainSize.y - diff.y;
    else if (diff.y < -domainSize.y/2)
        diff.y = domainSize.y + diff.y;
    return diff;
}

// overload cout to print a vector

std::ostream& operator<<(std::ostream& os, const Vector2D& vector) {
    os << "{x:" << std::setprecision(5) << vector.x << ", y:" << vector.y << "}";
    return os;
}