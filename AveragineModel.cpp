//
// Created by Dennis Goldfarb on 2/10/16.
//

#include "AveragineModel.h"

std::vector<unsigned int> AveragineModel::estimate_composition(double mass) {
    std::vector<unsigned int> composition(mercury::MAX_ELEMENTS);
    composition[elements::ELEMENTS::C] = std::round(mass * AVERAGINE_C);
    composition[elements::ELEMENTS::N] = std::round(mass * AVERAGINE_N);
    composition[elements::ELEMENTS::O] = std::round(mass * AVERAGINE_O);
    composition[elements::ELEMENTS::S] = std::round(mass * AVERAGINE_S);

    //double current_avg_mass = mercury::get_average_mass(composition);
    //composition[elements::ELEMENTS::H] = std::max(0.0, std::round((mass-current_avg_mass) / mercury::elemAtomicMasses[elements::ELEMENTS::H]));
    double current_mass = mercury::get_monoisotopic_mass(composition);

    composition[elements::ELEMENTS::H] = std::max(0.0, std::round((mass-current_mass) / mercury::elemMasses[elements::ELEMENTS::H][0]));


    return composition;
}

std::vector<unsigned int> AveragineModel::estimate_sulfur_corrected_composition(unsigned int num_sulfur, double mass) {
    std::vector<unsigned int> composition(mercury::MAX_ELEMENTS);

    double current_mass = mass - num_sulfur*AVERAGINE_S_MASS;

    composition[elements::ELEMENTS::C] = std::round(num_sulfur*AVERAGINE_S_C + current_mass * AVERAGINE_NO_S_C);
    composition[elements::ELEMENTS::N] = std::round(num_sulfur*AVERAGINE_S_N + current_mass * AVERAGINE_NO_S_N);
    composition[elements::ELEMENTS::O] = std::round(num_sulfur*AVERAGINE_S_O + current_mass * AVERAGINE_NO_S_O);
    composition[elements::ELEMENTS::S] = std::round(num_sulfur*AVERAGINE_S_S);

    current_mass = mercury::get_monoisotopic_mass(composition);

    composition[elements::ELEMENTS::H] = std::max(0.0, std::round((mass-current_mass) / mercury::elemMasses[elements::ELEMENTS::H][0]));


    return composition;
}
