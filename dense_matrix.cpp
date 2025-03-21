#include "dense_matrix.hpp"
#include <iostream>

DenseMatrix::DenseMatrix(std::size_t m, std::size_t n) : rows(m), cols(n), matrix(n*m,0.) {}
DenseMatrix::DenseMatrix(const std::vector<double>& v) : matrix(v), rows(v.size()), cols(1) {};

double DenseMatrix::determinant() const{
    if(rows==2 && cols==2){
        return matrix[0]*matrix[3] - matrix[1]*matrix[2];
    }
    if(rows==3 && cols==3){
        return matrix[0] * matrix[4] * matrix[8] +  // Diagonale principale
               matrix[1] * matrix[5] * matrix[6] +
               matrix[2] * matrix[3] * matrix[7] -
               matrix[2] * matrix[4] * matrix[6] -  // Diagonale inversa
               matrix[0] * matrix[5] * matrix[7] -
               matrix[1] * matrix[3] * matrix[8];
    }

    return 0;
}

DenseMatrix DenseMatrix::inverse() const{
    double det = determinant();
    if(det==0){
        std::cerr << "Determinant is zero, matrix is not invertible" << std::endl;
    }
    
    if (rows == 2 && cols == 2) {
        DenseMatrix result(2, 2);
        result(0, 0) = matrix[3] / det;
        result(0, 1) = -matrix[1] / det;
        result(1, 0) = -matrix[2] / det;
        result(1, 1) = matrix[0] / det;
        return result;
    }

    if (rows == 3 && cols == 3) {
        DenseMatrix inv(rows, cols);;

        // Calcolo dei cofattori
        double a = (*this)(0, 0), b = (*this)(0, 1), c = (*this)(0, 2);
        double d = (*this)(1, 0), e = (*this)(1, 1), f = (*this)(1, 2);
        double g = (*this)(2, 0), h = (*this)(2, 1), i = (*this)(2, 2);
        
        // Calcolo dei cofattori con trasposizione (matrice aggiunta)
        inv(0, 0) = (e * i - f * h) / det;
        inv(0, 1) = (c * h - b * i) / det;
        inv(0, 2) = (b * f - c * e) / det;
        inv(1, 0) = (f * g - d * i) / det;
        inv(1, 1) = (a * i - c * g) / det;
        inv(1, 2) = (c * d - a * f) / det;
        inv(2, 0) = (d * h - e * g) / det;
        inv(2, 1) = (b * g - a * h) / det;
        inv(2, 2) = (a * e - b * d) / det;

        return inv;
    }

    return DenseMatrix(rows, cols);
}

void DenseMatrix::print() const{
    for(std::size_t i=0; i<rows; ++i){
        for(std::size_t j=0; j<cols; ++j){
            std::cout << (*this)(i,j) << " ";
        }
        std::cout << std::endl;
    }
}

DenseMatrix DenseMatrix::transpose() const{
    DenseMatrix result(cols, rows);
    for(std::size_t i=0; i<rows; ++i){
        for(std::size_t j=0; j<cols; ++j){
            result(j,i) = (*this)(i,j);
        }
    }
    return result;
}

const double& DenseMatrix::operator()(std::size_t i, std::size_t j) const{
    return matrix[i*cols+j];
}

double& DenseMatrix::operator()(std::size_t i, std::size_t j){
    return matrix[i*cols+j];
}

DenseMatrix DenseMatrix::operator+(const DenseMatrix& other) const{
    DenseMatrix result(rows, cols);
    for(std::size_t i=0; i<rows; ++i){
        for(std::size_t j=0; j<cols; ++j){
            result(i,j) = matrix[i*cols+j] + other(i,j);
        }
    }
    return result;
}

DenseMatrix DenseMatrix::operator*(double scalar) const{
    DenseMatrix result(rows, cols);
    for(std::size_t i=0; i<rows; ++i){
        for(std::size_t j=0; j<cols; ++j){
            result(i,j) = matrix[i*cols+j] * scalar;
        }
    }
    return result;
}

MyVector DenseMatrix::operator*(const MyVector &v) const{
    std::vector<double> result(rows);
    for(std::size_t i=0; i<rows; ++i){
        double sum = 0.;
        for(std::size_t j=0; j<cols; ++j){
            sum += matrix[i*cols+j] * v(j);
        }
        result[i] = sum;
    }
    return MyVector(result);
}