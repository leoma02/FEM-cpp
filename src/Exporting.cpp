#include "Exporting.hpp"

#include <fstream>
#include <vector>
#include <array>
#include <stdexcept>
#include <string>

void writeVTK(const std::string &filename,
              const std::vector<std::vector<double>> &nodes,
              const std::vector<double> &solution,
              const std::forward_list<std::vector<int>> &elements,
              const std::vector<int> &cellTypes) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Impossibile aprire il file: " + filename);
    }

    // Intestazione del file VTK
    file << "# vtk DataFile Version 3.0\n";
    file << "Risultati FEM\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // Scrittura dei punti (vertici)
    file << "POINTS " << nodes.size() << " float\n";
    for (const auto &v : nodes) {
        file << v[0] << " " << v[1] << " " << v[2] << "\n";
    }

    // Calcolo del numero totale interi per la sezione CELLS:
    // Per ogni cella: 1 (numero di vertici) + numero di vertici
    int totalIntCount = 0;
    for (const auto &cell : elements) {
        totalIntCount += (1 + cell.size());
    }

    int size = 0;
    for(auto a=elements.cbegin(); a!=elements.cend(); ++a)
        ++size;


    // Scrittura delle celle
    file << "CELLS " << size << " " << totalIntCount << "\n";
    for (const auto &cell : elements) {
        file << cell.size();
        for (auto idx : cell)
            file << " " << idx;
        file << "\n";
    }

    // Scrittura del tipo di ogni cella
    file << "CELL_TYPES " << size << "\n";
    for (auto type : cellTypes) {
        file << type << "\n";
    }

    // Scrittura dei dati scalari associati ai punti (soluzione)
    file << "POINT_DATA " << solution.size() << "\n";
    file << "SCALARS solution float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (auto val : solution) {
        file << val << "\n";
    }

    file.close();
}
