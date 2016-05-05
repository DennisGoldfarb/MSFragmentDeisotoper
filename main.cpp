#include <iostream>
#include <vector>

#include "FragmentIsotopeCalculator.h"

int main(int argc, const char ** argv) {

    std::vector<double> mz, abundance;
    std::vector<unsigned int> precursor_isotopes = {6};
    std::vector<double> precursor_abundances = {1};
    // 64366168, 61518072, 31141578, 9734178, 2696946, 654956, 566000
    std::vector<unsigned int> total_composition(mercury::MAX_ELEMENTS);
    std::vector<unsigned int> fragment_composition(mercury::MAX_ELEMENTS);
    int charge = 1;
    double limit = 1e-30;

    // H
    total_composition[0] = 121;
    fragment_composition[0] = 105;
    // C
    total_composition[1] = 78;
    fragment_composition[1] = 67;
    // N
    total_composition[2] = 21;
    fragment_composition[2] = 19;
    // O
    total_composition[3] = 20;
    fragment_composition[3] = 17;

    FragmentIsotopeCalculator::fragment_isotopic_distribution(mz, abundance, precursor_isotopes, precursor_abundances,
                                                              total_composition, fragment_composition, charge, limit);

    for (int i = 0; i < mz.size(); i++) {
        std::cout << mz[i] << " " << abundance[i] << std::endl;
    }

    return 0;
}