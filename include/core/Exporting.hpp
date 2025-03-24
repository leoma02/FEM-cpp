#ifndef HEADER_EXP_HPP
#define HEADER_EXP_HPP

#include <fstream>
#include <vector>
#include <forward_list>

// Codici VTK per tipi di celle
constexpr int VTK_TRIANGLE = 5;
constexpr int VTK_TETRA    = 10;

void writeVTK(const std::string &filename,
              const std::vector<std::vector<double>> &nodes,
              const std::vector<double> &solution,
              const std::forward_list<std::vector<int>> &elements,
              const std::vector<int> &cellTypes);

#endif