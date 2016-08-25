//
// Created by Dennis Goldfarb on 8/17/16.
//

#ifndef MSFRAGMENTDEISOTOPER_FRAGMENTISOTOPEAPPROXIMATOR_H
#define MSFRAGMENTDEISOTOPER_FRAGMENTISOTOPEAPPROXIMATOR_H

#include <map>
#include <algorithm>

#include <xercesc/sax/HandlerBase.hpp>

#include "TensorSplineModel.h"
#include "AveragineModel.h"
#include "libmercury++.h"

class FragmentIsotopeApproximator {
private:

public:
    FragmentIsotopeApproximator() {};

    FragmentIsotopeApproximator(char* infile, xercesc::XercesDOMParser* parser);

    ~FragmentIsotopeApproximator() {
        for (auto itr = models.begin(); itr != models.end(); ) {
            itr = models.erase(itr);
        }
    }

    float calc_expected_probability(std::vector<unsigned int> precursor_composition, std::vector<unsigned int> fragment_composition,
                                    unsigned int precursor_isotope, unsigned int fragment_isotope);

    float calc_probability_spline(unsigned int num_sulfur, unsigned int num_comp_sulfur, unsigned int num_selenium,
                           unsigned int num_comp_selenium, unsigned int precursor_isotope, unsigned int fragment_isotope,
                           float precusor_mass, float fragment_mass);

    float calc_probability_averagine(unsigned int precursor_isotope, unsigned int fragment_isotope,
                                     float precursor_mass, float fragment_mass);

    float calc_probability_sulfur_corrected_averagine(unsigned int num_sulfur, unsigned int num_comp_sulfur,
                                                      unsigned int precursor_isotope, unsigned int fragment_isotope,
                                                      float precursor_mass, float fragment_mass);


    float calc_probability_sulfur_corrected_averagine(unsigned int num_sulfur, unsigned int num_comp_sulfur,
                                                      unsigned int precursor_isotope, unsigned int fragment_isotope,
                                                      float precursor_mass, float fragment_mass, float precursor_probability);


    std::map<ModelAttributes, TensorSplineModel*> models;
};


#endif //MSFRAGMENTDEISOTOPER_FRAGMENTISOTOPEAPPROXIMATOR_H
