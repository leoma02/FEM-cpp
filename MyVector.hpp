#ifndef HEADER_MYVECTOR_HPP
#define HEADER_MYVECTOR_HPP

// This class implements a double vector type which is an extension of the standard array: the data structure is a standard vector, it has some comfortable members and operators like the
// norm method, the product between a vector and a scalar number and the dot product between two vectors. 
// N.B. The aim of this class is to extend the standard array, not the standard vector. As you can see, there is no member which allows the user to expand the data structure. For this
//      reason, the size is computed only when an instance of the class is created.

#include <vector>
#include <cmath>
#include <iostream>
#include <numeric>
#include <execution>

class MyVector {
private:
    // Data structure
    std::vector<double> coordinates;
    std::size_t size;

public:
    // Constructors
    MyVector() = default;
    MyVector(const std::vector<double> &c);
    MyVector(int d);

    // Subscript operators
    const double& operator()(std::size_t i) const;
    double&       operator()(std::size_t i);

    // Getters
    std::vector<double>::const_iterator cbegin() const;
    std::vector<double>::const_iterator cend()   const;
    const std::size_t& getSize() const;
    const std::vector<double>& get_vector() const;

    // Useful methods
    const double norm() const;
    void print() const;

};

// Operators
MyVector operator+(const MyVector& p1, const MyVector& p2);
MyVector operator*(const double& scalar, const MyVector& p);
double operator*(const MyVector& p1, const MyVector& p2);
MyVector operator-(const MyVector& p1, const MyVector& p2);

#endif