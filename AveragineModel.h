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

    std::vector<unsigned int> estimate_composition(double mass);

    double AVERAGINE_MASS = 111.1254;
    double AVERAGINE_C = 4.9384;
    double AVERAGINE_H = 7.7583;
    double AVERAGINE_N = 1.3577;
    double AVERAGINE_O = 1.4773;
    double AVERAGINE_S = 0.0417;
};


#endif //MSFRAGMENTDEISOTOPER_AVERAGINEMODEL_H
