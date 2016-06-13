//
// Created by Dennis Goldfarb on 4/19/16.
//

#ifndef MSFRAGMENTDEISOTOPER_Y_ION_H
#define MSFRAGMENTDEISOTOPER_Y_ION_H


#include "Peptide.h"

class y_ion : public Peptide {

public:

    y_ion() {}

    y_ion(std::string sequence, int charge) : Peptide(sequence, charge, true) {
        mercury::mercury(isotope_mz, isotope_abundance, get_composition(), charge, 1e-30);
    }

    y_ion(std::string sequence, int charge, MolecularFormula loss) : Peptide(sequence, charge, true), loss(loss) {
            mercury::mercury(isotope_mz, isotope_abundance, get_composition(), charge, 1e-30);
    }

    std::vector<unsigned int> get_composition() const;
    std::vector<unsigned int> get_composition(int charge) const;
    double calc_monoisotopic_mass();
    int get_most_abundant_isotope();
    int get_max_isotope();
    MolecularFormula loss;
};



#endif //MSFRAGMENTDEISOTOPER_Y_ION_H
