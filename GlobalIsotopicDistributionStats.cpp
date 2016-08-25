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
#include "AveragineModel.h"


static std::string AMINO_ACIDS = "ADEFGHIKLNPQRSTVWY";
static std::string AMINO_ACIDS_SULFUR = "CM";
static std::string AMINO_ACIDS_SELENIUM = "U";

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> dis_AA(0, AMINO_ACIDS.length()-1);
std::uniform_int_distribution<> dis_S(0, AMINO_ACIDS_SULFUR.length()-1);

MolecularFormula waterloss("H2,O1");
MolecularFormula ammoniumloss("N1,H3");

double dynamic_range = 20000;
double dynamic_range_frag = 20000;
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

std::vector<std::vector<double>> calc_probabilities_averagine(double precursor_mass, double fragment_mass) {
    AveragineModel averagine_model;
    std::vector<std::vector<double>> precursor2fragment_probabilities;

    std::vector<unsigned int> precursor_composition = averagine_model.estimate_composition(precursor_mass);
    std::vector<unsigned int> fragment_composition = averagine_model.estimate_composition(fragment_mass);
    std::vector<unsigned int> complement_composition = averagine_model.estimate_composition(precursor_mass-fragment_mass);

    std::vector<double> frag_mz, frag_abundance, comp_frag_mz, comp_frag_abundance, precursor_mz, precursor_abundance;

    mercury::mercury(frag_mz, frag_abundance, fragment_composition, 1, 1e-30);
    mercury::mercury(comp_frag_mz, comp_frag_abundance, complement_composition, 1, 1e-30);
    mercury::mercury(precursor_mz, precursor_abundance, precursor_composition, 1, 1e-30);

    for (int precursor_isotope = 1; precursor_isotope < 18; ++precursor_isotope) {
        std::vector<double> fragment_probabilities(precursor_isotope+1);
        for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
            int comp_isotope = precursor_isotope - fragment_isotope;
            double probability = std::pow(2,std::log2(frag_abundance[fragment_isotope]) + std::log2(comp_frag_abundance[comp_isotope]) - std::log2(precursor_abundance[precursor_isotope]));
            fragment_probabilities.push_back(probability);
        }
        precursor2fragment_probabilities.push_back(fragment_probabilities);
    }

    return precursor2fragment_probabilities;
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

void create_fragments(Peptide &p, std::ofstream outfiles[max_isotope][max_isotope], int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium, std::string border) {
    int num_fragments = p.length()-1;
    int most_abundant_isotope_index = p.get_most_abundant_isotope();
    double most_abundant_isotope_abundance = p.isotope_abundance[most_abundant_isotope_index];
    double min_abundance = most_abundant_isotope_abundance/dynamic_range;

    int tot_left_SSe = num_sulfurs + num_selenium;
    int tot_right_SSe = num_c_sulfurs + num_c_selenium;

    std::vector<b_ion> b_ions;
    std::vector<y_ion> y_ions;

    std::vector<int> b_s;
    std::vector<int> b_se;
    std::vector<int> y_s;
    std::vector<int> y_se;

    for (int i = 1; i <= num_fragments; ++i) {   // for each cleavage site
        b_ions.push_back(b_ion(p.sequence.substr(0,i),0));
        y_ions.push_back(y_ion(p.sequence.substr(i,p.length()),0));

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


                for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) { // for each fragment isotope
                    unsigned int comp_isotope = precursor_isotope - fragment_isotope;

                    double fragment_isotope_abundance = std::log2(0);

                    if (b_ions[index].isotope_abundance.size() > fragment_isotope && y_ions[index].isotope_abundance.size() > comp_isotope) {
                        fragment_isotope_abundance = std::log2(b_ions[index].isotope_abundance[fragment_isotope]) +
                                                     std::log2(y_ions[index].isotope_abundance[comp_isotope]) -
                                                    std::log2(p.isotope_abundance[precursor_isotope]);
                    }

                    b_ion_isotope_abundances[fragment_isotope] = fragment_isotope_abundance;
                    y_ion_isotope_abundances[comp_isotope] = fragment_isotope_abundance;

                }

                // calculate minimum abundance for this fragment
                double min_abundance_b = *std::max_element(b_ion_isotope_abundances.begin(), b_ion_isotope_abundances.end()) - std::log2(dynamic_range_frag);
                double min_abundance_y = *std::max_element(y_ion_isotope_abundances.begin(), y_ion_isotope_abundances.end()) - std::log2(dynamic_range_frag);

                for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {

                    // check if appropriate given constraints: #S in fragment, #S in complement, #Se in fragment, #Se in complement
                    if (b_s[index] == num_sulfurs && b_se[index] == num_selenium && y_s[index] == num_c_sulfurs &&
                        y_se[index] == num_c_selenium) {

                        // check that the isotope abundance is greater than the minimum abundance
                        // check that both fragment and complementary fragment can have the number of isotopes
                        if (b_ion_isotope_abundances[fragment_isotope] > min_abundance_b &&
                            b_ions[index].get_max_isotope() >= fragment_isotope &&
                            y_ions[index].get_max_isotope() >= precursor_isotope - fragment_isotope) {


                            double abundance = std::pow(2,b_ion_isotope_abundances[fragment_isotope]);

                            if (!isnan(abundance) && !isinf(abundance)) {
                                // write to appropriate file: precursor isotope, fragment isotope
                                outfiles[precursor_isotope][fragment_isotope] << abundance << "\t" <<
                                                                         b_ions[index].calc_monoisotopic_mass() << "\t"
                                                                         << p.calc_monoisotopic_mass() <<
                                                                         std::endl;
                            }
                        }
                    }

                    // check if appropriate given constraints: #S in fragment, #S in complement, #Se in fragment, #Se in complement
                    if (y_s[index] == num_sulfurs && y_se[index] == num_selenium && b_s[index] == num_c_sulfurs &&
                        b_se[index] == num_c_selenium) {

                        // check that the isotope abundance is greater than the minimum abundance
                        // check that both fragment and complementary fragment can have the number of isotopes
                        if (y_ion_isotope_abundances[fragment_isotope] > min_abundance_b &&
                            y_ions[index].get_max_isotope() >= fragment_isotope &&
                            b_ions[index].get_max_isotope() >= precursor_isotope - fragment_isotope) {

                            double abundance = std::pow(2,y_ion_isotope_abundances[fragment_isotope]);

                            if (!isnan(abundance) && !isinf(abundance)) {
                                // write to appropriate file: precursor isotope, fragment isotope
                                outfiles[precursor_isotope][fragment_isotope] << abundance << "\t" <<
                                                                         y_ions[index].calc_monoisotopic_mass() << "\t"
                                                                         << p.calc_monoisotopic_mass() << std::endl;
                            }
                        }

                    }

                }
            }

        }
    }
}

void sample_fragment_isotopic_ratios(std::string base_path, float max_mass, int num_samples, int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium, std::string border) {

    // create all output files and write header to each
    std::ofstream outfiles[max_isotope][max_isotope];
    for (int precursor_isotope = 1; precursor_isotope < max_isotope; ++precursor_isotope) {
        for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
            std::string filename = "Precursor" + std::to_string(precursor_isotope) + "_" +
                                   "Fragment" + std::to_string(fragment_isotope) + ".tab";
            outfiles[precursor_isotope][fragment_isotope].open(base_path + filename);

            outfiles[precursor_isotope][fragment_isotope] << "ratio" << "\tfrag.mass" << "\tprecursor.mass" << std::endl; //"\tfrag.a.mass" << "\tprecursor.a.mass" << std::endl;
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

    // close all output files
    for (int precursor_isotope = 1; precursor_isotope < max_isotope; ++precursor_isotope) {
        for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
            outfiles[precursor_isotope][fragment_isotope].close();
        }
    }
}

int main(int argc, const char ** argv) {
    // Increase the maximum number of open files for this process. Was necessary for me.
    struct rlimit rlp;
    rlp.rlim_cur = 600;
    setrlimit(RLIMIT_NOFILE, &rlp);

    sample_fragment_isotopic_ratios(argv[1], atof(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), argv[8]);

    return 0;
}
