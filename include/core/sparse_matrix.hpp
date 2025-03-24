#ifndef HEADER_SPARSE_MATRIX_HPP
#define HEADER_SPARSE_MATRIX_HPP

#include <map>
#include <set>
#include <vector>
#include "MyVector.hpp"

class SparseMatrix {
private:
    std::vector<double> matrix;
    std::vector<int> row_index;
    std::vector<int> col_index;
    double buffer=0.;
    std::size_t rows;
    std::size_t cols;

public:
    SparseMatrix() = default;

    void build(const std::map<int, std::map<int, double>> &m, const std::set<std::pair<int,int>> &i);
    void print() const;
    void set_boundary(std::set<int> b);
    MyVector operator*(const MyVector &x) const;
    double& operator()(int i, int j);
};

#endif