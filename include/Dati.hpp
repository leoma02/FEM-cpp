#ifndef HEADER_DATI
#define HEADER_DATI

#include <cmath>
#include "MyVector.hpp"

inline double force(MyVector x){
    return 1000;
}

inline double boundary_function(MyVector x){
    if(std::pow(x(0),2) + std::pow(x(1),2) < 0.7)
        return 0;
    else
        return 0;
}

#endif