//
// Created by Dennis Goldfarb on 7/12/15.
//

#include "Peptide.h"

const char& Peptide::operator[](std::size_t i) {
	return sequence[i];
}

std::vector<unsigned int> Peptide::get_composition() const {
	std::vector<unsigned int> composition(mercury::MAX_ELEMENTS);

	for (char c : sequence) {
		const Residue* r = residues::name2residue.at(c);
		for (std::pair<int, int> pair : r->molecular_formula.element2count) {
			composition[pair.first] += pair.second;
		}
	}

	composition[elements::ELEMENTS::H] += 2;
	composition[elements::ELEMENTS::O] += 1;

	return composition;
}

std::vector<unsigned int> Peptide::get_composition(int charge) const {
	std::vector<unsigned int> composition = get_composition();
	composition[elements::ELEMENTS::H] += charge;
	return composition;
}

double Peptide::calc_monoisotopic_mass() {
	if (monoisotopic_mass == 0) {
		std::vector<unsigned int> element_counts = get_composition();
		for (int i = 0; i < element_counts.size(); i++) {
			monoisotopic_mass += element_counts[i] * mercury::elemMasses[i][0];
		}
	}
	return monoisotopic_mass;
}

void Peptide::remove_low_probability_isotopes() {
	auto itr_mz = isotope_mz.begin();
	auto itr_abundance = isotope_abundance.begin();
	for (; itr_mz != isotope_mz.end() && itr_abundance != isotope_abundance.end();) {
		if (*itr_abundance < MIN_ISOTOPE_ABUNDANCE) {
			itr_mz = isotope_mz.erase(itr_mz);
			itr_abundance = isotope_abundance.erase(itr_abundance);
		} else {
			++itr_mz;
			++itr_abundance;
		}
	}
}

const int Peptide::length() {
	return sequence.length();
}

const int Peptide::num_isotopes() {
	return isotope_mz.size();
}

std::vector<unsigned int> Peptide::get_b_ion_composition(int index, int charge) const {
	std::vector<unsigned int> composition(mercury::MAX_ELEMENTS);

	for (int i = 0; i < index && i < sequence.length(); ++i) {
		const Residue* r = residues::name2residue.at(sequence[i]);
		for (std::pair<int, int> pair : r->molecular_formula.element2count) {
			composition[pair.first] += pair.second;
		}
	}

	composition[elements::ELEMENTS::H] += 1;

	composition[elements::ELEMENTS::H] += charge;
	return composition;
}

std::vector<unsigned int> Peptide::get_y_ion_composition(int index, int charge) const {
	std::vector<unsigned int> composition(mercury::MAX_ELEMENTS);

	for (int i = index; i < sequence.length(); ++i) {
		const Residue* r = residues::name2residue.at(sequence[i]);
		for (std::pair<int, int> pair : r->molecular_formula.element2count) {
			composition[pair.first] += pair.second;
		}
	}

	composition[elements::ELEMENTS::H] += 1;
	composition[elements::ELEMENTS::O] += 1;

	composition[elements::ELEMENTS::H] += charge;
	return composition;
}