#pragma once

#include <iostream>
#include "typeDefinitions.hpp"

/**
 * @brief 2D vector class
 * 
 */
class Vector2D
{
private:

public:
    // public x and y for better performances:
    real_t x {0.};
    real_t y {0.};
    //**----------
    Vector2D(): x{0.}, y{0.} {}
    Vector2D(real_t x, real_t y): x{x}, y{y} {}
    //TODO: remove copy constructor
    Vector2D(const Vector2D& v): x{v.x}, y{v.y} {}
    //**----------

    /**
     * @brief set the vector with values x and y
     * 
     * @param x 
     * @param y 
     */
    void set(real_t x, real_t y);

    /**
     * @brief sum of 2 vectors and returns a new vector
     * 
     * @param otherVector 
     * @return Vector2D
     */
    Vector2D operator+(const Vector2D& otherVector) const;

    /**
     * @brief substraction
     * 
     * @param otherVector 
     * @return Vector2D 
     */
    Vector2D operator-(const Vector2D& otherVector) const;

    /**
     * @brief add a vector to the current vector
     * 
     * @param otherVector 
     * @return Vector2D& 
     */ 
    Vector2D& operator+=(const Vector2D& otherVector);

    /**
     * @brief -= operator overloading
     * 
     * @param otherVector 
     * @return Vector2D& 
     */
    Vector2D& operator-=(const Vector2D& otherVector);

    /**
     * @brief return dot product
     * 
     * @param otherVector 
     * @return Vector2D 
     */
    real_t dot(const Vector2D& otherVector) const;

    /**
     * @brief return result of multiplying a vector by a scalar
     * 
     * @param vector 
     * @param scalar 
     * @return Vector2D 
     */
    friend Vector2D operator*(const Vector2D& vector, real_t scalar);

    /**
     * @brief return result of multiplying a scalar by a vecor
     * 
     * @param vector 
     * @param scalar 
     * @return Vector2D 
     */
    friend Vector2D operator*(real_t scalar, const Vector2D& vector);

    /**
     * @brief *= operator overlading. multiply by scalar
     * 
     * @param scalar
     * @return Vector2D& 
     */
    Vector2D& operator*=(real_t scalar);

    /**
     * @brief apply the true modulo operation on each element
     * 
     * @param vector 
     * @param scalar 
     * @return Vector2D 
     */
    friend Vector2D operator%(const Vector2D& vector, real_t scalar);

    /**
     * @brief return the euclidean norm of a vector (static function)
     * 
     * @return real_t 
     */
    static real_t norm(const Vector2D& vector);

    /**
     * @brief return the euclidean norm of a vector
     * 
     * @param vector 
     * @return real_t 
     */
    real_t norm() const;

    /**
     * @brief return the absolute squared norm of a vector (static function)
     * 
     * @param vector 
     * @return real_t 
     */
    static real_t squaredNorm(const Vector2D& vector);

    /**
     * @brief return the absolute squared norm of a vector
     * 
     * @return real_t 
     */
    real_t squaredNorm() const;

    /**
     * @brief return the euclidean distance between 2 vectors a and b (static function)
     * 
     * @param a 
     * @param b 
     * @return real_t 
     */
    static real_t dist(const Vector2D& a, const Vector2D& b);

    /**
     * @brief return the squared distance between 2 vectors a and b (static function)
     * 
     * @param a 
     * @param b 
     * @return real_t 
     */
    static real_t squaredDist(const Vector2D& a, const Vector2D& b);

    /**
     * @brief overload cout to print a vector
     * 
     * @param os 
     * @param vector 
     * @return std::ostream& 
     */
    friend std::ostream& operator<<(std::ostream& os, const Vector2D& vector);
};
