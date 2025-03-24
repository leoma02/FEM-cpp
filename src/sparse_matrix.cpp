#include "sparse_matrix.hpp"
#include "dense_matrix.hpp"
#include <iostream>

void SparseMatrix::build(const std::map<int, std::map<int, double>> &m, const std::set<std::pair<int,int>> &i){
    std::size_t dim = i.size();
    matrix.reserve(dim);
    row_index.reserve(dim);
    col_index.reserve(dim);

    for(auto a=i.cbegin(); a!=i.cend(); ++a){
        row_index.push_back(a->first);
        col_index.push_back(a->second);
        matrix.push_back(m.at(a->first).at(a->second));
    }
}

void SparseMatrix::print() const{
    int count = 0;
    for(int i=0; i<row_index.size() && count < 1000000; ++i){
        std::cout << "(" << row_index[i] << ", " << col_index[i] << "): " << matrix[i] << std::endl;
        ++count;
    }
}

MyVector SparseMatrix::operator*(const MyVector &x) const{
    MyVector result(x.getSize());
    for(int i=0; i<matrix.size(); ++i){
        result(row_index[i]) += matrix[i]*x(col_index[i]);
    }

    return result;
}

double& SparseMatrix::operator()(int i, int j){
    for(int i=0; i<matrix.size(); ++i){
        if(row_index[i]==i && col_index[i]==j)
            return matrix[i];
    }

    return buffer;
}

void SparseMatrix::set_boundary(std::set<int> b){
    for(int i=0; i<matrix.size(); ++i){
        bool row = b.find(row_index[i])!=b.end();
        bool col = b.find(col_index[i])!=b.end();
        if(row_index[i]==col_index[i] && col && row){
            matrix[i] = 1;
        }
        else if(row || col){
            matrix[i] = 0;
        }
    }
}