#include "FEM.hpp"
#include <fstream>
#include <iostream>
#include "dense_matrix.hpp"
#include "MyVector.hpp"
#include <sstream>
#include "Dati.hpp"
#include <functional>
#include "Exporting.hpp"
#include <chrono>

int main() {

    std::ifstream mesh_file("Gear_alt.txt");
    double alpha = 1.;
    std::vector<double> beta = {0.,0.,0.};
    double gamma = 0.;
    std::function<double(MyVector)> f  = force;
    std::function<double(MyVector)> bf = boundary_function;

    auto start = std::chrono::high_resolution_clock::now();
    FEM fem(mesh_file, 3, alpha, beta, gamma, f, bf);
    //fem.print_nodes();
    //fem.print_elements();
    fem.assemble_matrix();
    fem.assemble_load_vector();
    fem.assemble_lifting();
    fem.set_boundary_conditions();
    fem.solve(1e-6,1000);
    MyVector u = fem.get_solution();
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

    std::vector<std::vector<double>> nodes = fem.get_nodes();
    std::forward_list<std::vector<int>> elements = fem.get_elements();
    int size = 0;
    for(auto a=elements.cbegin(); a!=elements.cend(); ++a)
        ++size;
    std::vector<int> cellTypes(size, VTK_TETRA);

    writeVTK("risultati.vtk", fem.get_nodes(), u.get_vector(), fem.get_elements(), cellTypes);

    return 0;
    }