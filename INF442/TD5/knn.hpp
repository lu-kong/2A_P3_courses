#pragma once

#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"

// TODO
class knn : public kernel{
        int k;
        double V;//volume
    public:

        double density(point &p);
        knn(cloud *data_, int k_, double V_);
};