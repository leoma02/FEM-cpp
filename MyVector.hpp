#ifndef HEADER_MYVECTOR_HPP
#define HEADER_MYVECTOR_HPP

#include <vector>

class MyVector {
private:
    std::vector<double> coordinates;
    std::size_t size;

public:
    MyVector() = default;
    MyVector(const std::vector<double> &c);
    MyVector(int d);

    const double& operator()(std::size_t i) const;
    double& operator()(std::size_t i);
    const std::size_t &getSize() const;
    const std::vector<double>& get_vector() const;
    const double norm() const;
    void print() const;

};

MyVector operator+(const MyVector& p1, const MyVector& p2);

MyVector operator*(const double& scalar, const MyVector& p);
double operator*(const MyVector& p1, const MyVector& p2);

MyVector operator-(const MyVector& p1, const MyVector& p2);

#endif