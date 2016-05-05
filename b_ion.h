//
// Created by Dennis Goldfarb on 4/19/16.
//

#ifndef MSFRAGMENTDEISOTOPER_B_ION_H
#define MSFRAGMENTDEISOTOPER_B_ION_H


#include "Peptide.h"

class b_ion : public Peptide {

public:

    b_ion() {}

    b_ion(std::string sequence, int charge) : Peptide(sequence, charge, true) {
        mercury::mercury(isotope_mz, isotope_abundance, get_composition(), charge, 1e-30);
    }

    std::vector<unsigned int> get_composition() const;
    std::vector<unsigned int> get_composition(int charge) const;
    double calc_monoisotopic_mass();
};


#endif //MSFRAGMENTDEISOTOPER_B_ION_H
