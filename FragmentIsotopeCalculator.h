//
// Created by Dennis Goldfarb on 12/6/15.
//

#ifndef MSFRAGMENTDEISOTOPER_FRAGMENTISOTOPECALCULATOR_H
#define MSFRAGMENTDEISOTOPER_FRAGMENTISOTOPECALCULATOR_H

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include "libmercury++.h"

namespace FragmentIsotopeCalculator {

    /* Need isolated precursor distribution
     * Need total composition
     * Need fragment composition
     */
    int fragment_isotopic_distribution(std::vector<double>& msa_mz, /* return value */
                                       std::vector<double>& msa_abundance,  /* return value */
                                       const std::vector<unsigned int>& precursor_isotopes,
                                       std::vector<double>& precursor_abundances,
                                       const std::vector<unsigned int>& total_composition,
                                       const std::vector<unsigned int>& fragment_composition,
                                       const int charge,
                                       const double limit);

    int fragment_isotopic_distribution(std::vector<double> &msa_mz,
                                       std::vector<double> &msa_abundance,
                                       std::vector<double> &tot_mz,
                                       std::vector<double> &tot_abundance,
                                       const std::vector<unsigned int>& precursor_isotopes,
                                       std::vector<double>& precursor_abundances,
                                       const std::vector<unsigned int> &total_composition,
                                       const std::vector<unsigned int> &fragment_composition,
                                       const int charge,
                                       const double limit);
}

#endif //MSFRAGMENTDEISOTOPER_FRAGMENTISOTOPECALCULATOR_H
