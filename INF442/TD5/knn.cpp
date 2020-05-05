#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"
#include "knn.hpp"
#include <algorithm>
// TODO
knn::knn(cloud * data_, int k_, double V_):kernel(data_){
    k = k_;
    V = V_;
}

double knn::density(point &p){
    int n = data->get_n();
    double * dis = new double[n];
    for (int i =0; i<n;i++){
        dis[i] = p.dist(data->get_point(i));
    }
    std::sort(dis,dis+n);


    return k/V/2/n/dis[k-1];
}