#include "MyVector.hpp"
#include <cmath>
#include <iostream>

MyVector::MyVector(const std::vector<double> &c) : coordinates(c), size(c.size()) {};
MyVector::MyVector(int d) : coordinates(std::vector<double>(d,0.)), size(d) {};

const double& MyVector::operator()(std::size_t i) const {
    return coordinates[i];
}

double& MyVector::operator()(std::size_t i){
    return coordinates[i];
}

const std::size_t& MyVector::getSize() const{
    return size;
}

const double MyVector::norm() const{
    double sum = 0;
    for(auto i : coordinates) {
        sum += i*i;
    }
    return std::sqrt(sum);
}

MyVector operator+(const MyVector& p1, const MyVector& p2) {
    std::vector<double> c(p1.getSize());
    for (std::size_t i=0; i<p1.getSize(); i++) {
        c[i] = p1(i) + p2(i);
    }
    return MyVector(c);
}

MyVector operator*(const double& scalar, const MyVector& p) {
    std::vector<double> c(p.getSize());
    for (std::size_t i=0; i<p.getSize(); i++) {
        c[i] = scalar * p(i);
    }
    return MyVector(c);
}

double operator*(const MyVector& p1, const MyVector& p2) {
    double sum = 0;
    for (std::size_t i=0; i<p1.getSize(); i++) {
        sum += p1(i) * p2(i);
    }
    return sum;
}

MyVector operator-(const MyVector& p1, const MyVector& p2) {
    std::vector<double> c(p1.getSize());
    for (std::size_t i=0; i<p1.getSize(); i++) {
        c[i] = p1(i) - p2(i);
    }
    return MyVector(c);
}

void MyVector::print() const{
    for(auto i : coordinates){
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

const std::vector<double>& MyVector::get_vector() const{
    return coordinates;
}