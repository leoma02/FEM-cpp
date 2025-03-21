#include "FEM.hpp"
#include "dense_matrix.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

FEM::FEM(std::ifstream& nodes_file, int dimension, double a, std::vector<double> b, double g, std::function<double(MyVector)> f, std::function<double(MyVector)> bf) 
        : dimension(dimension), alpha(a), beta(b), gamma(g), f(f), bf(bf) {
    std::string line;
    bool flag = true;

    while(getline(nodes_file, line)){
        std::vector<std::string> words;
        std::istringstream reader(line);
        std::string word;
        while(reader >> word){
            words.push_back(word);
        }
        if(flag){
            nodes.reserve(std::stoi(words[0]));
            flag = false;
        }
        if(words.size() == dimension+1){
            std::vector<double> node(dimension);
            for(int i=1; i<=dimension; ++i){
                node[i-1] = std::stod(words[i]);
            }
            nodes.push_back(node);
        }
        if(words.size() == 5+dimension+1){
            std::vector<int> element(dimension+1);
            for(int i=5; i<=5+dimension; ++i){
                element[i-5] = std::stoi(words[i])-1;
            }
            elements.push_front(element);
        }
        if(words.size() == 5+dimension){
            for(int i=5; i<5+dimension; ++i){
                boundary_nodes.insert(std::stoi(words[i])-1);
            }
        }
    }
    
    std::cout << "Mesh generated" << std::endl;
}

void FEM::compute_dof(){
    if(degree == 1){
        dof = nodes.size();
    }
}

void FEM::compute_nln(){
    if(degree == 1){
        nln = dimension+1;
    }
}

std::vector<MyVector> FEM::compute_gradients() const{
    std::vector<MyVector> gradients;
    gradients.resize(nln);
    if(degree == 1){
        gradients[0] = MyVector(std::vector<double>(dimension, -1.));
        for(int i=1; i<nln; ++i){
            std::vector<double> gradient(dimension, 0.);
            gradient[i-1] = 1.;
            gradients[i] = MyVector(gradient);
        }
    }

    return gradients;
}

void FEM::compute_nodes_weights(std::vector<std::vector<double>> &n, std::vector<double> &w) const{
    if(degree == 1){
        if(dimension == 2){
            n = {{1./6., 1./6.}, {2./3., 1./6.}, {1./6., 2./3.}};
            w = {1./6., 1./6., 1./6.};
        }
        if(dimension == 3){
            n = {{0.58541020,0.13819660,0.13819660,0.13819660}, {0.13819660,0.58541020,0.13819660,0.13819660},
                 {0.13819660,0.13819660,0.58541020,0.13819660}, {0.13819660,0.13819660,0.13819660,0.58541020}};
            w = {1./24., 1./24., 1./24., 1./24.};
        }
    } 

}

double FEM::compute_basis(std::vector<double> x, int i) const{
    if(dimension==2){
        if(i==0)
            return 1 - x[0] - x[1];
        else if(i==1)
            return x[0];
        else if(i==2)
            return x[1];
    }
    if(dimension==3){
        if(i==0)
            return 1 - x[0] - x[1] - x[2];
        if(i==1)
            return x[0];
        if(i==2)
            return x[1];
        if(i==3)
            return x[2];
    }
    
    return 0;
}

void FEM::assemble_matrix(){
    compute_dof();
    compute_nln();
    std::map<int, std::map<int, double>> matrix;
    std::set<std::pair<int,int>> indexes;
    //std::vector<Polynomial> basis = compute_basis();

    for(auto a=elements.cbegin(); a!=elements.cend(); ++a){
        std::vector<int> connectivity = *a;
        DenseMatrix bk(nodes[connectivity[0]]);
        DenseMatrix B(dimension, dimension);
        for(int i=1; i<=dimension; ++i){
            std::vector<double> coordinates = nodes[connectivity[i]];
            for(int j=0; j<dimension; ++j){
                B(j,i-1) = coordinates[j] - bk(j,0); //segmentation fault!!
            }
        }
        double det = B.determinant();
        B = B.inverse();
        B = B.transpose();
        std::vector<MyVector> gradients = compute_gradients();
        std::vector<std::vector<double>> n;
        std::vector<double> w;
        compute_nodes_weights(n, w);
        const int nqn = w.size();
        for(int i=0; i<nln; ++i){
            for(int j=0; j<nln; ++j){
                double sum_A = 0;
                double sum_M = 0;
                double sum_B = 0;
                for(int k=0; k<nqn; ++k){
                    MyVector v1 = B*gradients[j];
                    MyVector v2 = B*gradients[i];
                    sum_A += (v1 * v2) * det * w[k];
                    sum_M += compute_basis(n[k],j) * compute_basis(n[k],i) * det * w[k];
                    //sum_M += basis[j](n[k]) * basis[i](n[k]) * det * w[k];
                    sum_B += beta * v1 * compute_basis(n[k],i) * det * w[k];
                }
                matrix[connectivity[i]][connectivity[j]] += alpha*sum_A + gamma*sum_M + sum_B;
                indexes.insert(std::make_pair(connectivity[i], connectivity[j]));
            }
        }
    }

    stiffness_matrix.build(matrix, indexes);
    std::cout << "Stiffness matrix assembled" << std::endl;
    //stiffness_matrix.print();
}

void FEM::assemble_load_vector(){
    load_vector = MyVector(dof);

    for(auto a=elements.cbegin(); a!=elements.cend(); ++a){
        std::vector<int> connectivity = *a;
        MyVector bk(nodes[connectivity[0]]);
        DenseMatrix B(dimension, dimension);
        for(int i=1; i<=dimension; ++i){
            std::vector<double> coordinates = nodes[connectivity[i]];
            for(int j=0; j<dimension; ++j){
                B(j,i-1) = coordinates[j] - bk(j);
            }
        }
        double det = B.determinant();
        std::vector<std::vector<double>> n;
        std::vector<double> w;
        compute_nodes_weights(n, w);
        const int nqn = w.size();
        for(int i=0; i<nln; ++i){
            double sum = 0;
            for(int k=0; k<nqn; ++k){
                MyVector x(n[k]);
                sum += f(B*x + bk) * compute_basis(n[k],i) * det * w[k];
            }
            load_vector(connectivity[i]) += sum;
        }
    }

    std::cout << "Load vector assembled" << std::endl;
    //load_vector.print();
}

void FEM::assemble_lifting(){
    lifting = MyVector(dof);

    for(auto a=boundary_nodes.cbegin(); a!=boundary_nodes.cend(); ++a){
        int index = *a;
        lifting(index) = bf(MyVector(nodes[index]));
    }

    std::cout << "Lifting assembled" << std::endl;
    //lifting.print();
}

void FEM::set_boundary_conditions(){
    MyVector temp = stiffness_matrix*lifting;
    load_vector = load_vector - temp;
    stiffness_matrix.set_boundary(boundary_nodes);
    for(auto i=boundary_nodes.cbegin(); i!=boundary_nodes.cend(); ++i){
        int index = *i;
        load_vector(index) = 0;
    }
    std::cout << "Boundary conditions set" << std::endl;
}

void FEM::solve(double tol, int maxit) {
    double res = tol + 1;
    double alpha = 0;
    double beta  = 0;
    int it=0;
    MyVector xv(dof);
    MyVector xn(dof);
    MyVector rv = load_vector - stiffness_matrix*xv;
    MyVector rn(dof);
    MyVector p = rv;
    if(rv.norm() < tol){
        std::cout << "Iterations: " << it << std::endl;
        solution = xv + lifting;
        return;
    }

    while(rv.norm() > tol && it < maxit){
        alpha = rv*rv / (p*(stiffness_matrix*p));
        xn = xv + alpha*p;
        rn = rv - alpha*(stiffness_matrix*p);
        beta = rn*rn / (rv*rv);
        p = rn + beta*p;
        rv = rn;
        xv = xn;
        ++it;
    }

    std::cout << "Iterations: " << it << std::endl;
    solution = xn + lifting;
    return;
}

const MyVector& FEM::get_solution() const{
    return solution;
}

const std::vector<std::vector<double>>& FEM::get_nodes() const{
    return nodes;
}

const std::forward_list<std::vector<int>>& FEM::get_elements() const{
    return elements;
}

void FEM::print_nodes() const{
    
    auto a = nodes.cbegin();
    auto b = nodes.cend()-1;
    for(int i=0; i<dimension; ++i){
        std::cout << (*a)[i] << " ";
    }
    std::cout << std::endl; 
    for(int i=0; i<dimension; ++i){
        std::cout << (*b)[i] << " ";
    }
    std::cout << std::endl;
    std::cout << nodes.size() << std::endl;
    /*
    for(auto a=nodes.cbegin(); a!=nodes.cend(); ++a){
        for(int i=0; i<dimension; ++i){
            std::cout << (*a)[i] << " ";
        }
        std::cout << std::endl;
    }
    */
}

void FEM::print_elements() const{
    auto a = elements.cbegin();
    for(auto c=a->cbegin(); c!=a->cend(); ++c){
        std::cout << *c << " ";
    }
    std::cout << std::endl;
    auto b = a;
    a++;
    while(a!=elements.cend()){
        ++a;
        ++b;
    }
    for(auto c=b->cbegin(); c!=b->cend(); ++c){
        std::cout << *c << " ";
    }
    std::cout << std::endl;
}