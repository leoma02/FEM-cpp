#ifndef HEADER_DENSE_MATRIX_HPP
#define HEADER_DENSE_MATRIX_HPP

#include <vector>
#include "MyVector.hpp"

class DenseMatrix {
private:
    std::vector<double> matrix;
    std::size_t rows;
    std::size_t cols;

public:
    DenseMatrix(std::size_t m, std::size_t n);
    DenseMatrix(const std::vector<double>& v);

    double determinant() const;
    DenseMatrix inverse() const;
    DenseMatrix transpose() const;
    void print() const;

    const double& operator()(std::size_t i, std::size_t j) const;
    double& operator()(std::size_t i, std::size_t j);
    DenseMatrix operator+(const DenseMatrix& other) const;
    DenseMatrix operator*(double scalar) const;
    MyVector operator*(const MyVector &v) const;
};

#endif