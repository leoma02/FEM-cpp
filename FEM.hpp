#ifndef HEADER_FEM_HPP
#define HEADER_FEM_HPP

#include <vector>
#include <forward_list>
#include <fstream>
#include <set>
#include "MyVector.hpp"
#include "sparse_matrix.hpp"
#include <functional>

class FEM {
private:
    int dimension;
    int dof=0;
    int nln=0;
    int degree=1;
    double alpha = 0;
    MyVector beta;
    double gamma = 0;

    std::vector<std::vector<double>> nodes;
    std::forward_list<std::vector<int>> elements;
    std::set<int> boundary_nodes;

    SparseMatrix stiffness_matrix;
    std::function<double(MyVector)> f;
    std::function<double(MyVector)> bf;
    MyVector load_vector;
    MyVector lifting;
    MyVector solution;

    void compute_dof();
    void compute_nln();
    std::vector<MyVector> compute_gradients() const;
    void compute_nodes_weights(std::vector<std::vector<double>> &n, std::vector<double> &w) const;
    double compute_basis(std::vector<double> x, int i) const;

public:
    FEM(std::ifstream& nodes_file, int dimension, double a, std::vector<double> b, double g, std::function<double(MyVector)> f, std::function<double(MyVector)> bf);

    void assemble_matrix();
    void assemble_load_vector();
    void assemble_lifting();
    void set_boundary_conditions();
    void solve(double tol, int maxit);
    const MyVector& get_solution() const;
    const std::vector<std::vector<double>>& get_nodes() const;
    const std::forward_list<std::vector<int>>& get_elements() const;
    void print_nodes() const;
    void print_elements() const;
};

#endif