//
// Created by Dennis Goldfarb on 4/19/16.
//

#include "y_ion.h"


std::vector<unsigned int> y_ion::get_composition() const {
    std::vector<unsigned int> composition(mercury::MAX_ELEMENTS);

    for (char c : sequence) {
        const Residue* r = residues::name2residue.at(c);
        for (std::pair<int, int> pair : r->molecular_formula.element2count) {
            composition[pair.first] += pair.second;
        }
    }

    composition[elements::ELEMENTS::H] += 1;
    composition[elements::ELEMENTS::O] += 1;

    return composition;
}

std::vector<unsigned int> y_ion::get_composition(int charge) const {
    std::vector<unsigned int> composition = get_composition();
    composition[elements::ELEMENTS::H] += charge;
    return composition;
}

double y_ion::calc_monoisotopic_mass() {
    if (monoisotopic_mass == 0) {
        std::vector<unsigned int> element_counts = get_composition();
        for (int i = 0; i < element_counts.size(); i++) {
            monoisotopic_mass += element_counts[i] * mercury::elemMasses[i][0];
        }
    }
    return monoisotopic_mass;
}

int y_ion::get_most_abundant_isotope() {
    return std::distance(std::max_element(isotope_abundance.begin(), isotope_abundance.end()), isotope_abundance.begin());
}