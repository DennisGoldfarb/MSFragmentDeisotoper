//
// Created by Dennis Goldfarb on 2/10/16.
//

#include <iostream>
#include <vector>
#include <numeric>

#include "FragmentIsotopeCalculator.h"
#include "Peptide.h"
#include "AveragineModel.h"


void normalize_vector_L1(std::vector<double> &v) {
    double L1_norm = std::accumulate(v.begin(), v.end(), 0.0);

    for (int i = 0; i < v.size(); i++) {
        v[i]/=L1_norm;
    }
}

void normalize_vector_L2(std::vector<double> &v) {
    double L2_norm = std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));

    for (int i = 0; i < v.size(); i++) {
        v[i]/=L2_norm;
    }
}

void normalize_vector_max(std::vector<double> &v1) {
    double max = *std::max_element(v1.begin(), v1.end());

    for (int i = 0; i < v1.size(); i++) {
        v1[i]/=max;
    }
}

std::vector<double> calc_cross_correlation(std::vector<double> &v1, std::vector<double> &v2) {
    std::vector<double> cross_correlation_result;

    normalize_vector_L2(v1);
    normalize_vector_L2(v2);

    for (int offset = 1-(int)v1.size(); offset < (int)v1.size(); offset++) {
        double sum = 0;
        for (int i = std::max(0,-(int)offset), j = std::max(0,(int)offset); i < v1.size() && j < v2.size(); i++,j++) {
            sum += v1[i]*v2[j];
        }
        cross_correlation_result.push_back(sum);
    }

    return cross_correlation_result;
}

double calculate_xcorr(std::vector<double> &v1, std::vector<double> &v2) {
    std::vector<double> cross_correlation = calc_cross_correlation(v1,v2);

    double max_xcorr = *std::max_element(cross_correlation.begin(), cross_correlation.end());
    //double average_xcorr = std::accumulate(cross_correlation.begin(), cross_correlation.end(), 0.0)/cross_correlation.size();

    //std::cout << max_xcorr << " " << average_xcorr << " " << max_xcorr - average_xcorr << std::endl;

    return max_xcorr;// - average_xcorr;
}

int main(int argc, const char ** argv) {
    double limit = 1e-30;

    AveragineModel averagineModel;

    Peptide peptide1 = Peptide(std::string(argv[1]), atoi(argv[2]));
    Peptide peptide2 = Peptide(std::string(argv[3]), atoi(argv[4]));

    std::vector<unsigned int> estimated_composition1 = averagineModel.estimate_composition(peptide1.calc_monoisotopic_mass());
    std::vector<unsigned int> estimated_composition2 = averagineModel.estimate_composition(peptide2.calc_monoisotopic_mass());

    double min_mass = std::min(peptide1.calc_monoisotopic_mass(), peptide2.calc_monoisotopic_mass());
    int num_fragments = floor(min_mass / averagineModel.AVERAGINE_MASS);
    std::cout << peptide1.calc_monoisotopic_mass() << "\t" << peptide2.calc_monoisotopic_mass() << std::endl;

    std::vector<std::vector<double>> index2frag_mz(num_fragments), index2frag_abundance(num_fragments);
    std::vector<std::vector<double>> index2comp_frag_mz1(num_fragments), index2comp_frag_abundance1(num_fragments);
    std::vector<std::vector<double>> index2comp_frag_mz2(num_fragments), index2comp_frag_abundance2(num_fragments);

    for (int i = 0; i < num_fragments; i++) {
        double mass = (i+1)*averagineModel.AVERAGINE_MASS;
        // estimate fragment elemental composition via averagine model
        std::vector<unsigned int> fragment_composition = averagineModel.estimate_composition(mass);
        std::vector<unsigned int> comp_fragment_composition1 = averagineModel.estimate_composition(peptide1.calc_monoisotopic_mass()-mass);
        std::vector<unsigned int> comp_fragment_composition2 = averagineModel.estimate_composition(peptide2.calc_monoisotopic_mass()-mass);

        mercury::mercury(index2frag_mz[i], index2frag_abundance[i], fragment_composition, 1, limit);
        mercury::mercury(index2comp_frag_mz1[i], index2comp_frag_abundance1[i], comp_fragment_composition1, 1, limit);
        mercury::mercury(index2comp_frag_mz2[i], index2comp_frag_abundance2[i], comp_fragment_composition2, 1, limit);

    }

    for (int start1 = 0; start1 <= peptide1.num_isotopes(); ++start1) {
        for (int end1 = start1; end1 <= peptide1.num_isotopes(); ++end1) {
            for (int start2 = 0; start2 <= peptide2.num_isotopes(); ++start2) {
                for (int end2 = start2; end2 <= peptide2.num_isotopes(); ++end2) {

                    std::vector<double> xcorrs;

                    unsigned int max_isotope = std::max(end1,end2);
                    for (int fragment_index = 0; fragment_index < index2frag_mz.size(); fragment_index++) {

                        std::vector<double> fragment_isotope_abundances1(max_isotope+1);
                        std::vector<double> fragment_isotope_abundances2(max_isotope+1);

                        for (unsigned int fragment_isotope = 0; fragment_isotope <= max_isotope; ++fragment_isotope) {

                            double fragment_isotope_abundance = 0;
                            for (int precursor_isotope = start1; precursor_isotope <= end1; ++precursor_isotope) {

                                if (precursor_isotope >= fragment_isotope) {
                                    unsigned int comp_isotope = precursor_isotope - fragment_isotope;

                                    fragment_isotope_abundance += index2frag_abundance[fragment_index][fragment_isotope] * index2comp_frag_abundance1[fragment_index][comp_isotope];
                                }
                            }

                            double fragment_isotope_abundance2 = 0;
                            for (int precursor_isotope = start2; precursor_isotope <= end2; ++precursor_isotope) {

                                if (precursor_isotope >= fragment_isotope) {
                                    unsigned int comp_isotope = precursor_isotope - fragment_isotope;

                                    fragment_isotope_abundance2 += index2frag_abundance[fragment_index][fragment_isotope] * index2comp_frag_abundance2[fragment_index][comp_isotope];
                                }
                            }

                            fragment_isotope_abundances1[fragment_isotope] = fragment_isotope_abundance;
                            fragment_isotope_abundances2[fragment_isotope] = fragment_isotope_abundance2;
                        }

                        xcorrs.push_back(calculate_xcorr(fragment_isotope_abundances1, fragment_isotope_abundances2));
                    }

                    double average_xcorr = std::accumulate(xcorrs.begin(), xcorrs.end(), 0.0)/xcorrs.size();
                    std::cout << "Peptide 1 isotopes: " << start1 << " " << end1 << "\tPeptide 2 isotopes: " << start2 << " " << end2 << "\taverage xcorr: " << average_xcorr << std::endl;
                }
            }
        }
    }






    /*std::vector<double> frag_mz1, frag_abundance1;
    std::vector<double> frag_mz2, frag_abundance2;
    std::vector<double> total_mz1, total_abundance1;
    std::vector<double> total_mz2, total_abundance2;

    mercury::mercury(total_mz1, total_abundance1, estimated_composition1, peptide1.charge, limit);
    mercury::mercury(total_mz2, total_abundance2, estimated_composition2, peptide1.charge, limit);






    for (int start1 = 0; start1 < peptide1.num_isotopes(); ++start1) {
        for (int end1 = start1; end1 < peptide1.num_isotopes(); ++end1) {
            for (int start2 = 0; start2 < peptide2.num_isotopes(); ++start2) {
                for (int end2 = start2; end2 < peptide2.num_isotopes(); ++end2) {

                    std::cout << "Peptide 1 isotopes: " << start1 << " " << end1 << " Peptide 2 isotopes: " << start2 << " " << end2 << std::endl;

                    std::vector<unsigned int> precursor_isotopes1(end1-start1+1);
                    std::iota(precursor_isotopes1.begin(), precursor_isotopes1.end(), 0);
                    std::vector<double> precursor_abundances1(&peptide1.isotope_abundance[start1], &peptide1.isotope_abundance[end1]);

                    std::vector<unsigned int> precursor_isotopes2(end2-start2+1);
                    std::iota(precursor_isotopes2.begin(), precursor_isotopes2.end(), 0);
                    std::vector<double> precursor_abundances2(&peptide1.isotope_abundance[start2], &peptide1.isotope_abundance[end2]);


                    FragmentIsotopeCalculator::fragment_isotopic_distribution(mz1, abundance1, total_mz1, total_abundance1, precursor_isotopes1, precursor_abundances1,
                                                                              estimated_composition1, fragment_composition, 1, limit);
                    FragmentIsotopeCalculator::fragment_isotopic_distribution(mz2, abundance2, precursor_isotopes2, precursor_abundances2,

                                                                          estimated_composition2, fragment_composition, 1, limit);

                }
            }
        }
    }*/




/*
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
*/
    return 0;
}