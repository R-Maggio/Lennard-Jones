
#include <cmath> // sqrt
#include <iomanip> // std::setprecision

#include "Vector2D.hpp"

// setters:

template <typename T>
void Vector2D<T>::set(T x, T y) {
    this->x = x;
    this->y = y;
}

//sum of 2 vectors and returns a new vector

template <typename T>
Vector2D<T> Vector2D<T>::operator+(const Vector2D& otherVector) const {
    return Vector2D(this->x + otherVector.x, this->y + otherVector.y);
}

// substraction:

template <typename T>
Vector2D<T> Vector2D<T>::operator-(const Vector2D& otherVector) const {
    return Vector2D(this->x - otherVector.x, this->y - otherVector.y);
}

template <typename T>
Vector2D<T> Vector2D<T>::operator-() const {
    return Vector2D(-this->x, -this->y);
}

// add a vector to the current vector

template <typename T>
Vector2D<T>& Vector2D<T>::operator+=(const Vector2D& otherVector) {
    this->x += otherVector.x;
    this->y += otherVector.y;
    return *this;
}

// -= operator overloading

template <typename T>
Vector2D<T>& Vector2D<T>::operator-=(const Vector2D& otherVector) {
    this->x -= otherVector.x;
    this->y -= otherVector.y;
    return *this; 
}

// return the dot product

template <typename T>
real_t Vector2D<T>::dot(const Vector2D& otherVector) const {
    return this->x * otherVector.x + this->y * otherVector.y;
}

// return result of multiplying a vector by a scalar

// TODO: if scalar is float => return float
template <typename T, typename T2>
Vector2D<real_t> operator*(T scalar, const Vector2D<T2>& vector) {
    return Vector2D(vector.x * scalar, vector.y * scalar);
}

template <typename T, typename T2>
Vector2D<real_t> operator*(const Vector2D<T>& vector, T2 scalar) {
    return Vector2D(vector.x * scalar, vector.y * scalar);
}

//. *= operator overlading. multiply by scalar

template <typename T>
Vector2D<T>& Vector2D<T>::operator*=(T scalar) {
    this->x *= scalar;
    this->y *= scalar;
    return *this;
}

// division
template <typename T, typename T2>
Vector2D<real_t> operator/(const Vector2D<T>& vector, T2 scalar) {
    return Vector2D(vector.x / scalar, vector.y / scalar);
}

template <typename T, typename T2>
Vector2D<real_t> operator/(T scalar, const Vector2D<T2>& vector) {
    return Vector2D(scalar / vector.x, scalar / vector.y);
}

template <typename T>
Vector2D<T>& Vector2D<T>::operator/=(T scalar) {
    this->x /= scalar;
    this->y /= scalar;
    return *this;
}

// return the true modulus operation

template <typename T, typename T2>
Vector2D<real_t> operator%(const Vector2D<T>& vector, T2 scalar) {
    const real_t x = std::fmod(std::fmod(vector.x, scalar) + scalar, scalar);
    const real_t y = std::fmod(std::fmod(vector.y, scalar) + scalar, scalar);
    return Vector2D(x, y);
}

Vector2D operator%(const Vector2D& vector, const Vector2D& otherVector) {
    const real_t x = std::fmod(std::fmod(vector.x, otherVector.x) + otherVector.x, otherVector.x);
    const real_t y = std::fmod(std::fmod(vector.y, otherVector.y) + otherVector.y, otherVector.y);
    return Vector2D(x, y);
}

// return the euclidean norm of a vector (static function)

template <typename T>
real_t Vector2D<T>::norm() const {
    return std::sqrt(this->x*this->x + this->y*this->y);
}

// return the euclidean norm of a vector

template <typename T>
real_t Vector2D<T>::norm(const Vector2D& vector) {
    return std::sqrt(vector.x*vector.x + vector.y*vector.y);
}

// return the absolute squared norm of a vector

template <typename T>
real_t Vector2D<T>::squaredNorm() const {
    return this->x*this->x + this->y*this->y;
}

// return the absolute squared norm of a vector (static function)

template <typename T>
real_t Vector2D<T>::squaredNorm(const Vector2D& vector) {
    return vector.x*vector.x + vector.y*vector.y;
}

// return the euclidean distance between 2 vectors a and b (static function)

template <typename T>
real_t Vector2D<T>::dist(const Vector2D& a, const Vector2D& b) {
    return Vector2D::norm(a - b);
}

// return the squared distance between 2 vectors a and b (static function)

template <typename T>
real_t Vector2D<T>::squaredDist(const Vector2D& a, const Vector2D& b) {
    return Vector2D::squaredNorm(a - b);
}

template <typename T>
Vector2D<T> Vector2D<T>::PeriodicDiff(const Vector2D& otherVector, const Vector2D& domainSize) const {
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

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector2D<T>& vector) {
    os << "{x:" << std::setprecision(5) << vector.x << ", y:" << vector.y << "}";
    return os;
}