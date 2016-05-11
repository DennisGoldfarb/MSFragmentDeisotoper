//
// Created by Dennis Goldfarb on 7/12/15.
//

#ifndef MSACQUISITIONSIMULATOR_PEPTIDE_H
#define MSACQUISITIONSIMULATOR_PEPTIDE_H

#include <vector>
#include <map>
#include <algorithm>
#include "Residue.h"
#include "libmercury++.h"
#include <iostream>

class Peptide {

protected:
	void remove_low_probability_isotopes();
	const double MIN_ISOTOPE_ABUNDANCE = 1e-4;

public:
	Peptide() : sequence(""), charge(0)  {}

	Peptide(std::string sequence, int charge, bool no_isotopes) : sequence(sequence), charge(charge) {

	}

	Peptide(std::string sequence, int charge) : sequence(sequence), charge(charge) {
		mercury::mercury(isotope_mz, isotope_abundance, get_composition(), charge, 1e-30);
		//remove_low_probability_isotopes();
	}

	const char& operator[](std::size_t i);
	const int length();
	const int num_isotopes();

	std::vector<unsigned int> get_composition() const;
	virtual std::vector<unsigned int> get_composition(int charge) const;
	std::vector<unsigned int> get_b_ion_composition(int index, int charge) const;
	std::vector<unsigned int> get_y_ion_composition(int index, int charge) const;
	int get_most_abundant_isotope();


	double calc_monoisotopic_mass();
	double monoisotopic_mass=0;
	int charge;
	std::string sequence;
	std::vector<double> isotope_mz;
	std::vector<double> isotope_abundance;

};


#endif //MSACQUISITIONSIMULATOR_PEPTIDE_H
