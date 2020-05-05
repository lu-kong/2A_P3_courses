#include <cmath>

#include "point.hpp"
#include "flat.hpp"

// TODO
flat::flat(cloud * data_, double bandwidth_):radial(data_,bandwidth_){};

double flat::volume(){
    int dim = point::get_dim();
    double res = pow(M_PI,dim/2.0);
    double d = dim/2.0;
    if(dim%2!=0){
        res /= pow(M_PI,0.5);
    }
        // for (int d = dim / 2; d > 0; d--)
        // {
        //     res /= d;
        // }
        
        // d--;
    for(;d>0;d--){res /= d;}
    return res;
}

double flat::profile(double t){
    return (t<=1)? 1:0;
}