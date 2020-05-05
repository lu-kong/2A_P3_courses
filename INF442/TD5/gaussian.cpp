#include <cmath>
#include <iostream>

#include "point.hpp"
#include "cloud.hpp"
#include "gaussian.hpp"

// TODO
gaussian::gaussian(cloud *data_, double bandwidth_) : radial(data_, bandwidth_){};

double gaussian::volume(){
    int dim = data->get_point(0).get_dim();
    return pow(M_PI*2.0,dim/2.0);
}

double gaussian::profile(double t){
    return pow(M_E,t/-2.0);
}

void gaussian::guess_bandwidth(){
    double sgm = data->standard_deviation();
    bandwidth = 1.06*sgm/pow(data->get_n(),0.2);
}