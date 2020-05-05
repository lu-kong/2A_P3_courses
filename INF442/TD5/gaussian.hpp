#pragma once

#include "cloud.hpp"
#include "radial.hpp"

// TODO
class gaussian: public radial{
    public:
    double volume();
    double profile(double t);
    gaussian(cloud *data_, double bandwidth_);

    void guess_bandwidth();
};