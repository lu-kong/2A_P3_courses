#pragma once

#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"

// TODO
class radial : public kernel{
    
    // double bandwidth;


    public:
        double bandwidth;
        
        virtual double volume() = 0;
        virtual double profile(double t) = 0;   

        radial(cloud *data_, double bandwidth_);
        double density(point &p);
};