//
// Created by Dennis Goldfarb on 2/10/16.
//

#include "AveragineModel.h"

std::vector<unsigned int> AveragineModel::estimate_composition(double mass) {
    double num_averagines = mass / AVERAGINE_MASS;
    std::vector<unsigned int> composition(mercury::MAX_ELEMENTS);
    composition[elements::ELEMENTS::C] = (unsigned int) std::round(num_averagines * AVERAGINE_C);
    composition[elements::ELEMENTS::N] = (unsigned int) std::round(num_averagines * AVERAGINE_N);
    composition[elements::ELEMENTS::O] = (unsigned int) std::round(num_averagines * AVERAGINE_O);
    composition[elements::ELEMENTS::S] = (unsigned int) std::round(num_averagines * AVERAGINE_S);
    double current_avg_mass = mercury::get_average_mass(composition);
    composition[elements::ELEMENTS::H] = (unsigned int) std::round((mass-current_avg_mass) / mercury::elemAtomicMasses[elements::ELEMENTS::H]);

    return composition;
}
