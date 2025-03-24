#include "MyVector.hpp"

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
    auto sum    = [](double a, double b){return a+b;};
    auto square = [](double a){return a*a;};

    return std::sqrt(std::transform_reduce(coordinates.begin(), coordinates.end(), 0., sum, square));
}

MyVector operator+(const MyVector& p1, const MyVector& p2) {
    std::vector<double> res(p1.getSize());
    auto sum = [](double a, double b){return a+b;};
    std::transform(std::execution::par_unseq,p1.cbegin(), p1.cend(), p2.cbegin(), res.begin(),sum);
    return MyVector(res);
}

MyVector operator*(const double& scalar, const MyVector& p) {
    std::vector<double> res(p.getSize());
    auto prod = [scalar](double a){return a*scalar;};
    std::transform(std::execution::par_unseq, p.cbegin(), p.cend(), res.begin(), prod);
    return MyVector(res);
}

double operator*(const MyVector& p1, const MyVector& p2) {
    auto sum  = [](double a, double b){return a+b;};
    auto prod = [](double a, double b){return a*b;};
    return std::transform_reduce(p1.cbegin(), p1.cend(), p2.cbegin(), 0., sum, prod);
}

MyVector operator-(const MyVector& p1, const MyVector& p2) {
    std::vector<double> res(p1.getSize());
    auto sub = [](double a, double b){return a-b;};
    std::transform(p1.cbegin(), p1.cend(), p2.cbegin(), res.begin(), sub);
    return MyVector(res);
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

std::vector<double>::const_iterator MyVector::cbegin() const{
    return coordinates.cbegin();
}
std::vector<double>::const_iterator MyVector::cend() const{
    return coordinates.cend();
}