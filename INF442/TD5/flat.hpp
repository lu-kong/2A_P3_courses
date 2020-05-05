#pragma once

#include "cloud.hpp"
#include "radial.hpp"

// TODO
class flat : virtual public radial{
public:
    double volume();
    double profile(double t);
    flat(cloud * data_, double bandwidth_);
};