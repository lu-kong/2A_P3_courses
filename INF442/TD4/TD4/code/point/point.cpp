#include "point.hpp"

#include <iostream>
#include <cmath>
#include <cfloat>

int point::d = INT16_MAX;
point::point()
{
    coords = new double[d];
    for (int m = 0; m < d; m++)
        coords[m] = 0.0;
    label = 0;
}

point::~point()
{
    delete[] coords;
}

void point::print()
{
    std::cout << coords[0];

    for (int j = 1; j < d; j++)
        std::cout << '\t' << coords[j];

    std::cout << std::endl;
}

double point::dist(point &q)
{
    double sqd = 0.0;

    for (int m = 0; m < d; m++)
        sqd += (coords[m] - q.coords[m]) * (coords[m] - q.coords[m]);

    return std::sqrt(sqd);
}

int point::get_dim(){
    return d;
}

bool point::set_dim(int _d)
{
    
    if(point::d==INT16_MAX) return false;
    else{
        d = _d;
        return true;
    }
}

