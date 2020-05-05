#include <cmath>

#include "point.hpp"
#include "cloud.hpp"
#include "radial.hpp"

// TODO

radial::radial(cloud *data_, double bandwidth_):kernel(data_){
    bandwidth = bandwidth_;
}

 
double radial::density(point &p){
        //for point q 
        // p->dist(&q) give us the distance
        double res = 0;
        int n = data->get_n();
        int dim = p.get_dim();
        while(n-->0){
            double dis = p.dist(data->get_point(n));
            res += profile(dis*dis / bandwidth / bandwidth);
        }
        return res/data->get_n()/pow(bandwidth,dim)/volume();
    }
