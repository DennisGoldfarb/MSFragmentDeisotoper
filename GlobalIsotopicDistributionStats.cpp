//
// Created by Dennis Goldfarb on 2/23/16.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string>
#include <random>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "Peptide.h"
#include "b_ion.h"
#include "y_ion.h"


static std::string AMINO_ACIDS = "ADEFGHIKLNPQRSTVWY";
static std::string AMINO_ACIDS_SULFUR = "CM";
static std::string AMINO_ACIDS_SELENIUM = "U";

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> dis_AA(0, AMINO_ACIDS.length()-1);
std::uniform_int_distribution<> dis_S(0, AMINO_ACIDS_SULFUR.length()-1);

MolecularFormula waterloss("H2,O1");
MolecularFormula ammoniumloss("N1,H3");

double dynamic_range = 5000;
const int max_isotope = 30;

void get_precursor_isotopic_ratios(std::string path) {
    double limit = 1e-30;
    int max_isotopes = 7;
    double max_mass = 4000;

    std::ifstream infile(path);
    std::string line;

    while (infile >> line) {
        Peptide p = Peptide(line, 1);

        for (int i = 1; i < max_isotopes; i++) {
            double ratio = -1;
            if (i < p.num_isotopes() && p.calc_monoisotopic_mass() < max_mass) {
                std::cout << p.isotope_abundance[i]/p.isotope_abundance[i-1] << "\t"
                << i-1 << "\t"
                << p.calc_monoisotopic_mass() << "\t"
                << p.get_composition()[elements::ELEMENTS::S] << "\t"
                << std::count(p.sequence.begin(), p.sequence.end(), 'W') << "\t"
                << std::count(p.sequence.begin(), p.sequence.end(), 'Y') << "\t"
                << std::count(p.sequence.begin(), p.sequence.end(), 'F') << "\t"
                << std::count(p.sequence.begin(), p.sequence.end(), 'U') << "\t"
                << std::endl;

            }
        }
    }

    infile.close();
}

void get_fragment_isotopic_ratios(std::string path) {
    double limit = 1e-30;
    int max_isotopes = 7;

    std::cout << "ratio" << "\t" << "ratio_index" << "\t" << "frag.mass" << "\t" << "precursor.isotope" << "\t"
    << "precursor.mass" << "\t" << "ion.type" << "\t" << "sulfurs" << "\t"
    //<< "seleniums" << "\t"
    << "peptide" << "\t" << "fragment" << "\t" << "compliment.sulfurs" << std::endl;

    std::ifstream infile(path);
    std::string line;

    while (infile >> line) {
        Peptide p = Peptide(line, 0);
        if (p.get_composition()[elements::ELEMENTS::Se] > 0) continue;
        //p.sequence[p.sequence.length()-1] = 'A';
        int num_fragments = p.length()-1;

        std::vector<b_ion> b_ions;
        std::vector<y_ion> y_ions;
        for (int i = 1; i <= num_fragments; ++i) {   // for each cleavage site
            b_ions.push_back(b_ion(p.sequence.substr(0,i),0));
            y_ions.push_back(y_ion(p.sequence.substr(i,p.length()),0));
        }

        for (int precursor_isotope = 0; precursor_isotope < p.num_isotopes(); ++precursor_isotope) {
            for (int index = 0; index < num_fragments; ++index) {
                std::vector<double> b_ion_isotope_abundances(precursor_isotope+1);
                std::vector<double> y_ion_isotope_abundances(precursor_isotope+1);

                for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) { // for each fragment isotope
                    unsigned int comp_isotope = precursor_isotope - fragment_isotope;
                    double fragment_isotope_abundance = b_ions[index].isotope_abundance[fragment_isotope] * y_ions[index].isotope_abundance[comp_isotope];

                    b_ion_isotope_abundances[fragment_isotope] = fragment_isotope_abundance;
                    y_ion_isotope_abundances[comp_isotope] = fragment_isotope_abundance;
                }

                for (int fragment_isotope = 1; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
                    double ratio = b_ion_isotope_abundances[fragment_isotope]/b_ion_isotope_abundances[fragment_isotope-1];
                    std::cout << ratio << "\t" << fragment_isotope << "\t" << b_ions[index].calc_monoisotopic_mass() << "\t"
                    << precursor_isotope << "\t" << p.calc_monoisotopic_mass() << "\t" << "b" << "\t"
                    << b_ions[index].get_composition()[elements::ELEMENTS::S] << "\t"
                    //<< b_ions[index].get_composition()[elements::ELEMENTS::Se] << "\t"
                    << p.sequence << "\t" << b_ions[index].sequence << "\t"
                    << p.get_composition()[elements::ELEMENTS::S] - b_ions[index].get_composition()[elements::ELEMENTS::S] << std::endl;

                    ratio = y_ion_isotope_abundances[fragment_isotope]/y_ion_isotope_abundances[fragment_isotope-1];
                    std::cout << ratio << "\t" << fragment_isotope << "\t" << y_ions[index].calc_monoisotopic_mass() << "\t"
                    << precursor_isotope << "\t" << p.calc_monoisotopic_mass()  << "\t" << "y" << "\t"
                    << y_ions[index].get_composition()[elements::ELEMENTS::S] << "\t"
                    //<< y_ions[index].get_composition()[elements::ELEMENTS::Se] << "\t"
                    << p.sequence << "\t" << y_ions[index].sequence << "\t"
                    << p.get_composition()[elements::ELEMENTS::S] - y_ions[index].get_composition()[elements::ELEMENTS::S] << std::endl;

                }
            }
        }
    }

}


bool is_palindrome(std::string &s) {
    return equal(s.begin(), s.begin() + s.size()/2, s.rbegin());
}

std::string create_random_peptide_sequence(int peptide_length, int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium) {
    std::string random_peptide;

    // for insertion of sulfur containing amino acids in fragment
    for (int i = 0; i < num_sulfurs; ++i) {
        random_peptide.push_back(AMINO_ACIDS_SULFUR[dis_S(gen)]);
    }

    // for insertion of selenocysteines in fragment
    for (int i = 0; i < num_selenium; ++i) {
        random_peptide.push_back(AMINO_ACIDS_SELENIUM[0]);
    }

    // random amino acid insertion (non Sulfur and Selenium amino acids)
    for (int aa_index = 0; aa_index < peptide_length; ++aa_index) {
        random_peptide.push_back(AMINO_ACIDS[dis_AA(gen)]);
    }

    // for insertion of sulfur containing amino acids in fragment
    for (int i = 0; i < num_c_sulfurs; ++i) {
        random_peptide.push_back(AMINO_ACIDS_SULFUR[dis_S(gen)]);
    }

    // for insertion of selenocysteines in fragment
    for (int i = 0; i < num_c_selenium; ++i) {
        random_peptide.push_back(AMINO_ACIDS_SELENIUM[0]);
    }


    return random_peptide;
}

void create_fragments(Peptide &p, std::ofstream outfiles[max_isotope][max_isotope-1], int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium, std::string border) {
    int num_fragments = p.length()-1;
    int most_abundant_isotope_index = p.get_most_abundant_isotope();
    double most_abundant_isotope_abundance = p.isotope_abundance[most_abundant_isotope_index];
    double min_abundance = most_abundant_isotope_abundance/dynamic_range;

    int tot_left_SSe = num_sulfurs + num_selenium;
    int tot_right_SSe = num_c_sulfurs + num_c_selenium;

    std::vector<b_ion> b_ions;
    std::vector<y_ion> y_ions;

    std::vector<b_ion> b_ions_waterloss;
    std::vector<y_ion> y_ions_waterloss;

    std::vector<b_ion> b_ions_ammoniumloss;
    std::vector<y_ion> y_ions_ammoniumloss;


    std::vector<int> b_s;
    std::vector<int> b_se;
    std::vector<int> y_s;
    std::vector<int> y_se;

    for (int i = 1; i <= num_fragments; ++i) {   // for each cleavage site
        b_ions.push_back(b_ion(p.sequence.substr(0,i),0));
        y_ions.push_back(y_ion(p.sequence.substr(i,p.length()),0));

        b_ions_waterloss.push_back(b_ion(p.sequence.substr(0,i),0,waterloss));
        y_ions_waterloss.push_back(y_ion(p.sequence.substr(i,p.length()),0,waterloss));

        b_ions_ammoniumloss.push_back(b_ion(p.sequence.substr(0,i),0,ammoniumloss));
        y_ions_ammoniumloss.push_back(y_ion(p.sequence.substr(i,p.length()),0,ammoniumloss));

        b_s.push_back(b_ions[i-1].get_composition()[elements::ELEMENTS::S]);
        b_se.push_back(b_ions[i-1].get_composition()[elements::ELEMENTS::Se]);

        y_s.push_back(y_ions[i-1].get_composition()[elements::ELEMENTS::S]);
        y_se.push_back(y_ions[i-1].get_composition()[elements::ELEMENTS::Se]);
    }

    for (int precursor_isotope = 0; precursor_isotope < max_isotope && precursor_isotope < p.isotope_abundance.size(); ++precursor_isotope) {

        if (p.isotope_abundance[precursor_isotope] >= min_abundance) {

            for (int index = tot_left_SSe; index < num_fragments-tot_right_SSe; ++index) {
                std::vector<double> b_ion_isotope_abundances(precursor_isotope + 1);
                std::vector<double> y_ion_isotope_abundances(precursor_isotope + 1);

                std::vector<double> b_ion_waterloss_isotope_abundances(precursor_isotope + 1);
                std::vector<double> y_ion_waterloss_isotope_abundances(precursor_isotope + 1);

                std::vector<double> b_ion_ammoniumloss_isotope_abundances(precursor_isotope + 1);
                std::vector<double> y_ion_ammoniumloss_isotope_abundances(precursor_isotope + 1);

                for (int fragment_isotope = 0;
                     fragment_isotope <= precursor_isotope; ++fragment_isotope) { // for each fragment isotope
                    unsigned int comp_isotope = precursor_isotope - fragment_isotope;

                    double fragment_isotope_abundance = 0;

                    if (b_ions[index].isotope_abundance.size() > fragment_isotope && y_ions[index].isotope_abundance.size() > comp_isotope) {
                        fragment_isotope_abundance = b_ions[index].isotope_abundance[fragment_isotope] *
                                                     y_ions[index].isotope_abundance[comp_isotope];
                    }

                    b_ion_isotope_abundances[fragment_isotope] = fragment_isotope_abundance;
                    y_ion_isotope_abundances[comp_isotope] = fragment_isotope_abundance;

                    // b-ion water loss
                    fragment_isotope_abundance = 0;
                    if (b_ions_waterloss[index].isotope_abundance.size() > fragment_isotope && y_ions[index].isotope_abundance.size() > comp_isotope) {
                        fragment_isotope_abundance = b_ions_waterloss[index].isotope_abundance[fragment_isotope] *
                                                     y_ions[index].isotope_abundance[comp_isotope];
                    }

                    b_ion_waterloss_isotope_abundances[fragment_isotope] = fragment_isotope_abundance;

                    // b-ion ammonium loss
                    fragment_isotope_abundance = 0;
                    if (b_ions_ammoniumloss[index].isotope_abundance.size() > fragment_isotope && y_ions[index].isotope_abundance.size() > comp_isotope) {
                        fragment_isotope_abundance = b_ions_ammoniumloss[index].isotope_abundance[fragment_isotope] *
                                                     y_ions[index].isotope_abundance[comp_isotope];
                    }

                    b_ion_ammoniumloss_isotope_abundances[fragment_isotope] = fragment_isotope_abundance;

                    // y-ion water loss
                    fragment_isotope_abundance = 0;
                    if (b_ions[index].isotope_abundance.size() > fragment_isotope && y_ions_waterloss[index].isotope_abundance.size() > comp_isotope) {
                        fragment_isotope_abundance = b_ions[index].isotope_abundance[fragment_isotope] *
                                                     y_ions_waterloss[index].isotope_abundance[comp_isotope];
                    }

                    y_ion_waterloss_isotope_abundances[comp_isotope] = fragment_isotope_abundance;

                    // y-ion ammonium loss
                    fragment_isotope_abundance = 0;
                    if (b_ions[index].isotope_abundance.size() > fragment_isotope && y_ions_ammoniumloss[index].isotope_abundance.size() > comp_isotope) {
                        fragment_isotope_abundance = b_ions[index].isotope_abundance[fragment_isotope] *
                                                 y_ions_ammoniumloss[index].isotope_abundance[comp_isotope];
                    }

                    y_ion_ammoniumloss_isotope_abundances[comp_isotope] = fragment_isotope_abundance;

                }

                // calculate minimum abundance for this fragment
                double min_abundance_b = *std::max_element(b_ion_isotope_abundances.begin(), b_ion_isotope_abundances.end())/dynamic_range;
                double min_abundance_y = *std::max_element(y_ion_isotope_abundances.begin(), y_ion_isotope_abundances.end())/dynamic_range;
                double min_abundance_b_waterloss = *std::max_element(b_ion_waterloss_isotope_abundances.begin(), b_ion_waterloss_isotope_abundances.end())/dynamic_range;
                double min_abundance_y_waterloss = *std::max_element(y_ion_waterloss_isotope_abundances.begin(), y_ion_waterloss_isotope_abundances.end())/dynamic_range;
                double min_abundance_b_ammoniumloss = *std::max_element(b_ion_ammoniumloss_isotope_abundances.begin(), b_ion_ammoniumloss_isotope_abundances.end())/dynamic_range;
                double min_abundance_y_ammoniumloss = *std::max_element(y_ion_ammoniumloss_isotope_abundances.begin(), y_ion_ammoniumloss_isotope_abundances.end())/dynamic_range;

                for (int fragment_isotope = 1; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
                    // write to appropriate file: precursor isotope, fragment isotope
                    // check if appropriate given constraints: #S in fragment, #S in complement, #Se in fragment, #Se in complement
                    int ratio_index = fragment_isotope - 1;

                    if (b_s[index] == num_sulfurs && b_se[index] == num_selenium && y_s[index] == num_c_sulfurs &&
                        y_se[index] == num_c_selenium) {

                        // check that both isotopes are greater than the minimum abundance
                        // check that both can have the number of isotopes
                        if (b_ion_isotope_abundances[fragment_isotope] > min_abundance_b &&
                            b_ion_isotope_abundances[fragment_isotope - 1] > min_abundance_b &&
                            b_ions[index].get_max_isotope() >= fragment_isotope &&
                            y_ions[index].get_max_isotope() >= precursor_isotope - fragment_isotope) {

                            double ratio = std::log2(b_ion_isotope_abundances[fragment_isotope]) -
                                           std::log2(b_ion_isotope_abundances[fragment_isotope - 1]);
                            if (!isnan(ratio) && !isinf(ratio)) {
                                outfiles[precursor_isotope][ratio_index] << ratio << "\t" <<
                                                                         b_ions[index].calc_monoisotopic_mass() << "\t"
                                                                         << p.calc_monoisotopic_mass() <<
                                                                         std::endl;
                            }
                        }

                        if (b_ion_waterloss_isotope_abundances[fragment_isotope] > min_abundance_b &&
                            b_ion_waterloss_isotope_abundances[fragment_isotope - 1] > min_abundance_b &&
                            b_ions_waterloss[index].get_max_isotope() >= fragment_isotope &&
                            y_ions[index].get_max_isotope() >= precursor_isotope - fragment_isotope) {

                            double ratio = std::log2(b_ion_waterloss_isotope_abundances[fragment_isotope]) -
                                           std::log2(b_ion_waterloss_isotope_abundances[fragment_isotope - 1]);
                            if (!isnan(ratio) && !isinf(ratio)) {
                                outfiles[precursor_isotope][ratio_index] << ratio << "\t" <<
                                                                         b_ions_waterloss[index].calc_monoisotopic_mass()
                                                                         << "\t" <<
                                                                         p.calc_monoisotopic_mass() <<
                                                                         std::endl;
                            }

                        }

                        if (b_ion_ammoniumloss_isotope_abundances[fragment_isotope] > min_abundance_b &&
                            b_ion_ammoniumloss_isotope_abundances[fragment_isotope - 1] > min_abundance_b &&
                            b_ions_ammoniumloss[index].get_max_isotope() >= fragment_isotope &&
                            y_ions[index].get_max_isotope() >= precursor_isotope - fragment_isotope) {

                            double ratio = std::log2(b_ion_ammoniumloss_isotope_abundances[fragment_isotope]) -
                                           std::log2(b_ion_ammoniumloss_isotope_abundances[fragment_isotope - 1]);
                            if (!isnan(ratio) && !isinf(ratio)) {
                                outfiles[precursor_isotope][ratio_index] << ratio << "\t" <<
                                                                         b_ions_ammoniumloss[index].calc_monoisotopic_mass()
                                                                         << "\t" <<
                                                                         p.calc_monoisotopic_mass() <<
                                                                         std::endl;
                            }
                        }

                    }

                    if (y_s[index] == num_sulfurs && y_se[index] == num_selenium && b_s[index] == num_c_sulfurs &&
                        b_se[index] == num_c_selenium) {

                        if (y_ion_isotope_abundances[fragment_isotope] > min_abundance_b &&
                            y_ion_isotope_abundances[fragment_isotope - 1] > min_abundance_b &&
                            y_ions[index].get_max_isotope() >= fragment_isotope &&
                            b_ions[index].get_max_isotope() >= precursor_isotope - fragment_isotope) {

                            double ratio = std::log2(y_ion_isotope_abundances[fragment_isotope]) -
                                           std::log2(y_ion_isotope_abundances[fragment_isotope - 1]);
                            if (!isnan(ratio) && !isinf(ratio)) {
                                outfiles[precursor_isotope][ratio_index] << ratio << "\t" <<
                                                                         y_ions[index].calc_monoisotopic_mass() << "\t"
                                                                         << p.calc_monoisotopic_mass() << std::endl;
                            }
                        }

                        if (y_ion_waterloss_isotope_abundances[fragment_isotope] > min_abundance_b &&
                            y_ion_waterloss_isotope_abundances[fragment_isotope - 1] > min_abundance_b &&
                            y_ions_waterloss[index].get_max_isotope() >= fragment_isotope &&
                            b_ions[index].get_max_isotope() >= precursor_isotope - fragment_isotope) {

                            double ratio = std::log2(y_ion_waterloss_isotope_abundances[fragment_isotope]) -
                                           std::log2(y_ion_waterloss_isotope_abundances[fragment_isotope - 1]);
                            if (!isnan(ratio) && !isinf(ratio)) {
                                outfiles[precursor_isotope][ratio_index] << ratio << "\t" <<
                                                                         y_ions_waterloss[index].calc_monoisotopic_mass()
                                                                         << "\t" << p.calc_monoisotopic_mass() <<
                                                                         std::endl;
                            }
                        }

                        if (y_ion_ammoniumloss_isotope_abundances[fragment_isotope] > min_abundance_b &&
                            y_ion_ammoniumloss_isotope_abundances[fragment_isotope - 1] > min_abundance_b &&
                            y_ions_ammoniumloss[index].get_max_isotope() >= fragment_isotope &&
                            b_ions[index].get_max_isotope() >= precursor_isotope - fragment_isotope) {

                            double ratio = std::log2(y_ion_ammoniumloss_isotope_abundances[fragment_isotope]) -
                                           std::log2(y_ion_ammoniumloss_isotope_abundances[fragment_isotope - 1]);
                            if (!isnan(ratio) && !isinf(ratio)) {
                                outfiles[precursor_isotope][ratio_index] << ratio << "\t" <<
                                                                         y_ions_ammoniumloss[index].calc_monoisotopic_mass()
                                                                         << "\t" << p.calc_monoisotopic_mass() <<
                                                                         std::endl;
                            }
                        }
                    }

                }
            }

            // create smallest carbon-only fragment
            if (p.length() > 0 && border == "T") {
                for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope-1; ++fragment_isotope) { // for each fragment isotope
                    // the smallest carbon-only fragment has fragment_isotope carbons
                    int num_carbons = 1 + fragment_isotope - 2 * num_sulfurs - 5 * num_selenium;
                    std::vector<unsigned int> composition_smallest(mercury::MAX_ELEMENTS);
                    composition_smallest[elements::ELEMENTS::Se] = num_selenium;
                    composition_smallest[elements::ELEMENTS::S] = num_sulfurs;
                    composition_smallest[elements::ELEMENTS::C] = num_carbons;
                    MolecularFormula loss = MolecularFormula(
                            "Se" + std::to_string(num_selenium) + ",S" + std::to_string(num_sulfurs) + ",C" +
                            std::to_string(num_carbons));

                    std::vector<double> isotope_mz;
                    std::vector<double> isotope_abundance;
                    mercury::mercury(isotope_mz, isotope_abundance, composition_smallest, 0, 1e-30);

                    b_ion c_b = b_ion(p.sequence, 0, loss);
                    y_ion c_y = y_ion(p.sequence, 0, loss);

                    unsigned int comp_isotope = precursor_isotope - fragment_isotope;

                    double c_abundance1 = c_b.isotope_abundance[comp_isotope];
                    double c_abundance2 = c_b.isotope_abundance[comp_isotope - 1];
                    double abundance1 = isotope_abundance[fragment_isotope];
                    double abundance2 = isotope_abundance[fragment_isotope + 1];


                    double fragment_isotope_abundance1 = std::log2(abundance1) + std::log2(c_abundance1);
                    double fragment_isotope_abundance2 = std::log2(abundance2) + std::log2(c_abundance2);

                    double ratio = fragment_isotope_abundance2 - fragment_isotope_abundance1;
                    if (!isnan(ratio) && !isinf(ratio)) {
                        outfiles[precursor_isotope][fragment_isotope + 1] << ratio << "\t" <<
                                                                          loss.get_monoisotopic_mass() << "\t"
                                                                          << p.calc_monoisotopic_mass() <<
                                                                          std::endl;
                        std::cout << ratio << "\t" <<
                                   loss.get_monoisotopic_mass() << "\t"
                                   << p.calc_monoisotopic_mass() <<
                                   std::endl;
                    }
                }
            }
        }
    }
}

void sample_fragment_isotopic_ratios(std::string base_path, float max_mass, int num_samples, int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium, std::string border) {

    if (border == "T") {
        std::cout << "T!" << std::endl;
    } else {
        std::cout << "F!" << std::endl;
    }

    std::ofstream outfiles[max_isotope][max_isotope-1];

    for (int precursor_isotope = 1; precursor_isotope < max_isotope; ++precursor_isotope) {
        for (int ratio_index = 0; ratio_index < precursor_isotope; ++ratio_index) {
            std::string filename = "Precursor" + std::to_string(precursor_isotope) + "_" +
                                   "Ratio" + std::to_string(ratio_index) + ".tab";
            outfiles[precursor_isotope][ratio_index].open(base_path + filename);

            outfiles[precursor_isotope][ratio_index] << "ratio" << "\tfrag.mass" << "\tprecursor.mass" << std::endl; //"\tfrag.a.mass" << "\tprecursor.a.mass" << std::endl;
        }
    }

    int max_length = max_mass/100;

    for (int peptide_length = 0; peptide_length <= max_length; ++peptide_length) {

        for (int sample = 0; sample < num_samples; ++sample) {


            std::string random_sequence = create_random_peptide_sequence(peptide_length, num_sulfurs, num_c_sulfurs,
                                                                         num_selenium, num_c_selenium);

            Peptide p = Peptide(random_sequence, 0);
            if (p.calc_monoisotopic_mass() <= max_mass) {

                create_fragments(p, outfiles, num_sulfurs, num_c_sulfurs, num_selenium, num_c_selenium, border);

                if (!is_palindrome(random_sequence)) {
                    std::reverse(random_sequence.begin(), random_sequence.end());
                    Peptide reverse_p = Peptide(random_sequence, 0);
                    create_fragments(reverse_p, outfiles, num_sulfurs, num_c_sulfurs, num_selenium, num_c_selenium, border);
                }
            }
        }

    }

    for (int precursor_isotope = 1; precursor_isotope < max_isotope; ++precursor_isotope) {
        for (int ratio_index = 0; ratio_index < precursor_isotope; ++ratio_index) {
            outfiles[precursor_isotope][ratio_index].close();
        }
    }
}

int main(int argc, const char ** argv) {
    //get_precursor_isotopic_ratios(argv[1]);
    //get_fragment_isotopic_ratios(argv[1]);
    struct rlimit rlp;
    rlp.rlim_cur = 600;
    setrlimit(RLIMIT_NOFILE, &rlp);

    sample_fragment_isotopic_ratios(argv[1], atof(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), argv[8]);

    return 0;
}
