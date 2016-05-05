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

#include "Peptide.h"
#include "b_ion.h"
#include "y_ion.h"


static std::string AMINO_ACIDS = "ADEFGHIKLNPQRSTVWY";
static std::string AMINO_ACIDS_SULFUR = "CM";
static std::string AMINO_ACIDS_SELENIUM = "U";

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

void sample_fragment_isotopic_ratios(std::string base_path, int max_length, int num_sulfurs, int num_c_sulfurs, int num_selenium, int num_c_selenium) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis_AA(0, AMINO_ACIDS.length()-1);
    std::uniform_int_distribution<> dis_S(0, AMINO_ACIDS_SULFUR.length()-1);

    const int max_isotope = 9;
    int tot_sulfur = num_sulfurs + num_c_sulfurs;
    int tot_selenium = num_selenium + num_c_selenium;

    std::ofstream outfiles[max_isotope][max_isotope-1];

    for (int precursor_isotope = 1; precursor_isotope < max_isotope; ++precursor_isotope) {
        for (int ratio_index = 0; ratio_index < precursor_isotope; ++ratio_index) {
            std::string filename = "Precursor" + std::to_string(precursor_isotope) + "_" +
                                   "Ratio" + std::to_string(ratio_index) + ".tab";
                                   //"Ratio" + std::to_string(ratio_index) + "_" +
                                   //"S" + std::to_string(num_sulfurs) + "_" +
                                   //"CS" + std::to_string(num_c_sulfurs) + "_" +
                                   //"Se" + std::to_string(num_selenium) + "_" +
                                   //"CSe" + std::to_string(num_c_selenium) + ".tab";
            outfiles[precursor_isotope][ratio_index].open(base_path + filename);

            outfiles[precursor_isotope][ratio_index] << "ratio" << "\t" << "frag.mass" << "\t" << "precursor.mass" << std::endl;

        }
    }



    for (int peptide_length = tot_selenium + tot_sulfur; peptide_length <= max_length; ++peptide_length) {
        double SSe_factor = (tot_selenium+tot_sulfur)/5.0;
        int num_samples = std::pow(20,std::pow(peptide_length, std::min(1.0,(3.0+SSe_factor)/peptide_length) ));

        for (int sample = 0; sample < num_samples; ++sample) {
            std::string random_peptide;
            int positions[peptide_length];

            for (int aa_index = 0; aa_index < peptide_length; ++aa_index) {
                random_peptide.push_back(AMINO_ACIDS[dis_AA(gen)]);
                positions[aa_index] = aa_index;
            }

            // for random insertion of sulfur and selenium containing amino acids
            if (tot_sulfur + tot_selenium + num_c_sulfurs> 0) {
                std::random_shuffle(positions, positions + peptide_length);

                for (int pos_index = 0; pos_index < tot_sulfur; ++pos_index) {
                    random_peptide[positions[pos_index]] = AMINO_ACIDS_SULFUR[dis_S(gen)];
                }

                for (int pos_index = tot_sulfur; pos_index < tot_sulfur + tot_selenium; ++pos_index) {
                    random_peptide[positions[pos_index]] = AMINO_ACIDS_SELENIUM[0];
                }
            }



            Peptide p = Peptide(random_peptide, 0);
            int num_fragments = p.length()-1;

            std::vector<b_ion> b_ions;
            std::vector<y_ion> y_ions;
            for (int i = 1; i <= num_fragments; ++i) {   // for each cleavage site
                b_ions.push_back(b_ion(p.sequence.substr(0,i),0));
                y_ions.push_back(y_ion(p.sequence.substr(i,p.length()),0));
            }




            for (int precursor_isotope = 0; precursor_isotope < max_isotope; ++precursor_isotope) {
                for (int index = 0; index < num_fragments; ++index) {
                    std::vector<double> b_ion_isotope_abundances(precursor_isotope + 1);
                    std::vector<double> y_ion_isotope_abundances(precursor_isotope + 1);

                    for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) { // for each fragment isotope
                        unsigned int comp_isotope = precursor_isotope - fragment_isotope;
                        double fragment_isotope_abundance = b_ions[index].isotope_abundance[fragment_isotope] *
                                                            y_ions[index].isotope_abundance[comp_isotope];

                        b_ion_isotope_abundances[fragment_isotope] = fragment_isotope_abundance;
                        y_ion_isotope_abundances[comp_isotope] = fragment_isotope_abundance;
                    }


                    for (int fragment_isotope = 1; fragment_isotope <= precursor_isotope; ++fragment_isotope) {


                        // write to appropriate file: precursor isotope, fragment isotope
                        // check if appropriate given constraints: #S in fragment, #S in complement, #Se in fragment, #Se in complement
                        int b_s = b_ions[index].get_composition()[elements::ELEMENTS::S];
                        int b_se = b_ions[index].get_composition()[elements::ELEMENTS::Se];
                        int y_s = y_ions[index].get_composition()[elements::ELEMENTS::S];
                        int y_se = y_ions[index].get_composition()[elements::ELEMENTS::Se];

                        int ratio_index = fragment_isotope-1;

                        if (b_s == num_sulfurs && b_se == num_selenium && y_s == num_c_sulfurs && y_se == num_c_selenium) {

                            double ratio = b_ion_isotope_abundances[fragment_isotope] / b_ion_isotope_abundances[fragment_isotope-1];
                            outfiles[precursor_isotope][ratio_index] << ratio << "\t" << b_ions[index].calc_monoisotopic_mass() << "\t" <<
                            p.calc_monoisotopic_mass() << std::endl;

                        }

                        if (y_s == num_sulfurs && y_se == num_selenium && b_s == num_c_sulfurs && b_se == num_c_selenium) {

                            double ratio = y_ion_isotope_abundances[fragment_isotope] / y_ion_isotope_abundances[fragment_isotope-1];
                            outfiles[precursor_isotope][ratio_index] << ratio << "\t" << y_ions[index].calc_monoisotopic_mass() << "\t" <<
                            p.calc_monoisotopic_mass() << std::endl;
                        }

                    }
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

    sample_fragment_isotopic_ratios(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));

    return 0;
}
