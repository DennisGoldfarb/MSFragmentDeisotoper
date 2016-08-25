//
// Created by Dennis Goldfarb on 2/10/16.
//

#ifndef MSFRAGMENTDEISOTOPER_AVERAGINEMODEL_H
#define MSFRAGMENTDEISOTOPER_AVERAGINEMODEL_H

#include <vector>
#include "libmercury++.h"
#include "Element.h"

class AveragineModel {

private:

public:

    AveragineModel() {};

    std::vector<unsigned int> estimate_composition(double mass);
    std::vector<unsigned int> estimate_sulfur_corrected_composition(unsigned int num_sulfur, double mass);

    double AVERAGINE_C = 0.044439885;
    double AVERAGINE_H = 0.069815722;
    double AVERAGINE_N = 0.012217729;
    double AVERAGINE_O = 0.01329399;
    double AVERAGINE_S = 0.000375252;

    /*double AVERAGINE_NO_S_C = 0.04497953;
    double AVERAGINE_NO_S_H = 0.070663513;
    double AVERAGINE_NO_S_N = 0.012366092;
    double AVERAGINE_NO_S_O = 0.013455423;*/

    double AVERAGINE_NO_S_C = 0.04475297;
    double AVERAGINE_NO_S_H = 0.06966804;
    double AVERAGINE_NO_S_N = 0.0125457;
    double AVERAGINE_NO_S_O = 0.01357137;

    double AVERAGINE_S_C = 3.96230073;
    double AVERAGINE_S_H = 6.92460146;
    double AVERAGINE_S_N = 1;
    double AVERAGINE_S_O = 1;
    double AVERAGINE_S_S = 1;
    double AVERAGINE_S_MASS = 116.496452;
};


#endif //MSFRAGMENTDEISOTOPER_AVERAGINEMODEL_H
