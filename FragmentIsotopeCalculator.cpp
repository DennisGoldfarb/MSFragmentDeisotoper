//
// Created by Dennis Goldfarb on 12/6/15.
//

#include "FragmentIsotopeCalculator.h"

namespace FragmentIsotopeCalculator {

    int fragment_isotopic_distribution(std::vector<double> &msa_mz,
                                       std::vector<double> &msa_abundance,
                                       const std::vector<unsigned int>& precursor_isotopes,
                                       std::vector<double>& precursor_abundances,
                                       const std::vector<unsigned int> &total_composition,
                                       const std::vector<unsigned int> &fragment_composition,
                                       const int charge,
                                       const double limit) {

        double total_precursor_abundance = std::accumulate(precursor_abundances.begin(), precursor_abundances.end(), 0.0);
        for (int i = 0; i < precursor_abundances.size(); i++) {
            precursor_abundances[i] /= total_precursor_abundance;
        }


        std::vector<unsigned int> complementary_fragment_composition;
        for (int i = 0; i < total_composition.size(); i++) {
            complementary_fragment_composition.push_back(total_composition[i] - fragment_composition[i]);
        }

        std::vector<double> frag_mz, frag_abundance, comp_frag_mz, comp_frag_abundance, tot_mz, tot_abundance;

        mercury::mercury(frag_mz, frag_abundance, fragment_composition, charge, limit);
        mercury::mercury(comp_frag_mz, comp_frag_abundance, complementary_fragment_composition, charge, limit);
        mercury::mercury(tot_mz, tot_abundance, total_composition, charge, limit);

        auto max_isotope = std::max_element(precursor_isotopes.begin(), precursor_isotopes.end());

        // for each possible fragment isotope
        for (unsigned int fragment_isotope = 0; fragment_isotope <= *max_isotope; ++fragment_isotope) {

            double fragment_abundance = 0;

            for (int i = 0; i < precursor_abundances.size(); ++i) {

                unsigned int precursor_isotope = precursor_isotopes[i];
                double precursor_abundance = precursor_abundances[i];

                if (precursor_isotope >= fragment_isotope) {
                    unsigned int comp_isotope = precursor_isotope - fragment_isotope;

                    fragment_abundance += (frag_abundance[fragment_isotope] * comp_frag_abundance[comp_isotope]) * (precursor_abundance / tot_abundance[i]);
               }
            }

            msa_mz.push_back(fragment_isotope/charge);
            msa_abundance.push_back(fragment_abundance);
        }

        double total_abundance = *std::max_element(msa_abundance.begin(), msa_abundance.end());
        for (int i = 0; i < msa_abundance.size(); i++) {
            msa_abundance[i] /= .01*total_abundance;
        }

        return 0;
    }

    int fragment_isotopic_distribution(std::vector<double> &msa_mz,
                                       std::vector<double> &msa_abundance,
                                       std::vector<double> &tot_mz,
                                       std::vector<double> &tot_abundance,
                                       const std::vector<unsigned int>& precursor_isotopes,
                                       std::vector<double>& precursor_abundances,
                                       const std::vector<unsigned int> &total_composition,
                                       const std::vector<unsigned int> &fragment_composition,
                                       const int charge,
                                       const double limit) {

        double total_precursor_abundance = std::accumulate(precursor_abundances.begin(), precursor_abundances.end(), 0.0);
        for (int i = 0; i < precursor_abundances.size(); i++) {
            precursor_abundances[i] /= total_precursor_abundance;
        }


        std::vector<unsigned int> complementary_fragment_composition;
        for (int i = 0; i < total_composition.size(); i++) {
            if (total_composition[i] < fragment_composition[i]) {
                std::cout << "OH NO!!" << std::endl;
            }
            complementary_fragment_composition.push_back(total_composition[i] - fragment_composition[i]);
        }

        std::vector<double> frag_mz, frag_abundance, comp_frag_mz, comp_frag_abundance;

        mercury::mercury(frag_mz, frag_abundance, fragment_composition, charge, limit);
        mercury::mercury(comp_frag_mz, comp_frag_abundance, complementary_fragment_composition, charge, limit);
        auto max_isotope = std::max_element(precursor_isotopes.begin(), precursor_isotopes.end());

        // for each possible fragment isotope
        for (unsigned int fragment_isotope = 0; fragment_isotope <= *max_isotope; ++fragment_isotope) {

            double fragment_abundance = 0;

            for (int i = 0; i < precursor_abundances.size(); ++i) {

                unsigned int precursor_isotope = precursor_isotopes[i];
                double precursor_abundance = precursor_abundances[i];

                if (precursor_isotope >= fragment_isotope) {
                    unsigned int comp_isotope = precursor_isotope - fragment_isotope;

                    fragment_abundance += (frag_abundance[fragment_isotope] * comp_frag_abundance[comp_isotope]) * (precursor_abundance / tot_abundance[i]);
                }
            }

            msa_mz.push_back(fragment_isotope/charge);
            msa_abundance.push_back(fragment_abundance);
        }

        double total_abundance = *std::max_element(msa_abundance.begin(), msa_abundance.end());
        for (int i = 0; i < msa_abundance.size(); i++) {
            msa_abundance[i] /= .01*total_abundance;
        }

        return 0;
    }
}