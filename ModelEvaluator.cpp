//
// Created by Dennis Goldfarb on 8/2/16.
//

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <numeric>

#include <zlib.h>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>

#include "kseq.h"

#include "TensorSplineModel.h"
#include "FragmentIsotopeApproximator.h"
#include "Peptide.h"
#include "b_ion.h"
#include "y_ion.h"
#include "Histogram.h"

KSEQ_INIT(gzFile, gzread)


void calc_differences(float prob_spline, float prob_averagineS, float prob_averagine, float abundance, std::vector<std::pair<float,ModelAttributes>>& all_probs,
                      FragmentIsotopeApproximator* FIA, Histogram* spline, Histogram* averagine_s,
                      Histogram* averagine, Histogram* num_better_prob, Histogram* best_prob_diff, float& distance_b) {
    if (prob_spline != -1) {
        double diff = std::abs(prob_spline - abundance);
        spline->add_data(diff);
        distance_b += diff;

        int num_better = -1;
        float best_diff = diff;
        for (std::pair<float,ModelAttributes> &v : all_probs) {
            float other_model_diff = std::abs(v.first-prob_spline);
            num_better += other_model_diff <= diff;
            best_diff = std::min(std::abs(v.first-abundance), best_diff);

        }

        num_better_prob->add_data(num_better);
        best_prob_diff->add_data(std::abs(best_diff - diff));
        averagine_s->add_data(std::abs(prob_averagineS - abundance));
        averagine->add_data(std::abs(prob_averagine - abundance));
    }
}

void test_peptide(std::string peptide_seq, FragmentIsotopeApproximator* FIA, Histogram* spline, Histogram* averagine_s,
                  Histogram* averagine, Histogram* statistical_distance, Histogram* num_better_prob, Histogram* best_prob_diff) {

    Peptide p(peptide_seq, 0);
    int num_fragments = p.length() - 1;
    std::vector<b_ion> b_ions;
    std::vector<y_ion> y_ions;

    std::vector<int> b_s;
    std::vector<int> y_s;

    for (int i = 1; i <= num_fragments; ++i) {   // for each cleavage site
        b_ions.push_back(b_ion(p.sequence.substr(0, i), 0));
        y_ions.push_back(y_ion(p.sequence.substr(i, p.length()), 0));

        b_s.push_back(b_ions[i - 1].get_composition()[elements::ELEMENTS::S]);
        y_s.push_back(y_ions[i - 1].get_composition()[elements::ELEMENTS::S]);

    }

    for (int precursor_isotope = 1; precursor_isotope < 6; ++precursor_isotope) {
        for (int index = 0; index < num_fragments; ++index) {
            std::vector<double> b_ion_isotope_abundances(precursor_isotope + 1);
            std::vector<double> y_ion_isotope_abundances(precursor_isotope + 1);

            for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
                int comp_isotope = precursor_isotope - fragment_isotope;

                double fragment_isotope_abundance = std::log2(0);

                if (b_ions[index].isotope_abundance.size() > fragment_isotope &&
                    y_ions[index].isotope_abundance.size() > comp_isotope) {
                    fragment_isotope_abundance = std::log2(b_ions[index].isotope_abundance[fragment_isotope]) +
                                                 std::log2(y_ions[index].isotope_abundance[comp_isotope]) -
                                                 std::log2(p.isotope_abundance[precursor_isotope]);
                }

                b_ion_isotope_abundances[fragment_isotope] = fragment_isotope_abundance;
                y_ion_isotope_abundances[comp_isotope] = fragment_isotope_abundance;

            }

            float distance_y = 0;
            float distance_b = 0;
            std::vector<double> distances_y;
            std::vector<double> distances_b;

            for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
                double abundance = std::pow(2, b_ion_isotope_abundances[fragment_isotope]);
                if (!isnan(abundance) && !isinf(abundance) && index > 0) {

                    float pm = p.calc_monoisotopic_mass(), fm = b_ions[index].calc_monoisotopic_mass();
                    int pi = precursor_isotope, fi = fragment_isotope;
                    int S = b_s[index], CS = p.get_composition()[elements::ELEMENTS::S]-S;

                    float prob_spline = FIA->calc_probability_spline(S,CS,0,0,pi,fi,pm,fm,false);
                    float prob_averagineS = FIA->calc_probability_sulfur_corrected_averagine(S, CS, pi, fi, pm, fm);
                    float prob_averagine = FIA->calc_probability_averagine(pi, fi, pm, fm);
                    std::vector<std::pair<float,ModelAttributes>> all_probs = FIA->get_all_spline_probabilities(pm,fm);

                    calc_differences(prob_spline, prob_averagineS, prob_averagine, abundance, all_probs,
                                     FIA, spline, averagine_s, averagine, num_better_prob, best_prob_diff, distance_b);
                }
                abundance = std::pow(2, y_ion_isotope_abundances[fragment_isotope]);
                if (!isnan(abundance) && !isinf(abundance) && index > 0) {

                    float pm = p.calc_monoisotopic_mass(), fm = y_ions[index].calc_monoisotopic_mass();
                    int pi = precursor_isotope, fi = fragment_isotope;
                    int S = y_s[index], CS = p.get_composition()[elements::ELEMENTS::S]-S;

                    float prob_spline = FIA->calc_probability_spline(S,CS,0,0,pi,fi,pm,fm,false);
                    float prob_averagineS = FIA->calc_probability_sulfur_corrected_averagine(S, CS, pi, fi, pm, fm);
                    float prob_averagine = FIA->calc_probability_averagine(pi, fi, pm, fm);
                    std::vector<std::pair<float,ModelAttributes>> all_probs = FIA->get_all_spline_probabilities(pm,fm);

                    calc_differences(prob_spline, prob_averagineS, prob_averagine, abundance, all_probs,
                                     FIA, spline, averagine_s, averagine, num_better_prob, best_prob_diff, distance_y);

                }
            }

            if (distance_b != 0 && distance_y != 0) {
                statistical_distance->add_data(distance_b);
                statistical_distance->add_data(distance_y);
            }
        }
    }
}


void test_tryptic_peptides(const char* path_FASTA, FragmentIsotopeApproximator* FIA, int job_id, int num_jobs) {
    Histogram* spline = new Histogram("Spline Residual Distribution", "residual", "log2(count)", 0.01);
    Histogram* averagine_s = new Histogram("Averagine S Residual Distribution", "residual", "log2(count)", 0.01);
    Histogram* averagine = new Histogram("Averagine Distribution", "residual", "log2(count)", 0.01);
    Histogram* statistical_distance = new Histogram("Statistical Distance Distribution", "Statistical Distance", "log2(count)", 0.01);
    Histogram* num_better_prob = new Histogram("Number of Better Probabilities Distribution", "# of models", "log2(count)", 1);
    Histogram* best_prob_diff = new Histogram("Correct - Best Model Probability Distribution", "Difference", "log2(count)", 0.01);

    gzFile fp;
    fp = gzopen(path_FASTA, "r");
    if (!fp) {
        std::cerr << "Error opening FASTA file: " << path_FASTA << std::endl;
        throw std::exception();
    }

    kseq_t *seq;
    int l;

    std::map<char,int> char2count;
    std::unordered_set<std::string> peptides;

    int i = 0;
    seq = kseq_init(fp); // initialize seq
    while ((l = kseq_read(seq)) >= 0) { // read sequence
        std::string prot_seq = seq->seq.s;
        unsigned int start_i = 0, end_i = 0;

        while (start_i < prot_seq.size()) {
            if (prot_seq[end_i] == 'R' || prot_seq[end_i] == 'K' || end_i == prot_seq.size()-1) {
                std::string peptide_seq = prot_seq.substr(start_i, end_i-start_i+1);
                bool valid = true;
                for (int k = 0; k < peptide_seq.size(); ++k) {
                    if (peptide_seq[k] == 'B' || peptide_seq[k] == 'X' || peptide_seq[k] == 'Z' || peptide_seq[k] == 'U') {
                        valid = false;
                    }
                }

                start_i = end_i+1;

                if (valid && peptide_seq.size() >= 4 && peptide_seq.size() <= 85) {
                    peptides.insert(peptide_seq);
                }
            }

            end_i++;
        }
    }
    kseq_destroy(seq);
    gzclose(fp);

    int num_test = std::ceil(peptides.size()/num_jobs);
    int offset = (job_id-1)*num_test;

    for (auto itr = peptides.begin(); itr != peptides.end(); ++itr) {
        if (i >= offset && i < offset+num_test) {
            test_peptide(*itr, FIA, spline, averagine_s, averagine, statistical_distance, num_better_prob, best_prob_diff);

           /*if (i%1000 == 0) {
                std::cout << "Number of peptides processed: " << i << std::endl;
                spline->print_histogram();
                averagine_s->print_histogram();
                averagine->print_histogram();
            }*/
        }

        /*for (int ai = 0; ai < itr->size(); ai++) {
            char c = (*itr)[ai];
            if (char2count.find(c) == char2count.end()) char2count[c] = 0;
            char2count[c]++;
        }*/

        i++;
    }


    std::cout << "Number of peptides processed: " << i << std::endl;
    spline->print_histogram();
    averagine_s->print_histogram();
    averagine->print_histogram();
    statistical_distance->print_histogram();
    num_better_prob->print_histogram();
    best_prob_diff->print_histogram();

}

int main(int argc, char ** argv) {
    try {
        xercesc::XMLPlatformUtils::Initialize();
    } catch (const xercesc::XMLException& e) {
        std::cerr << "Error during Xerces initialization. " << e.getMessage() << std::endl;
        return 1;
    }

    xercesc::XercesDOMParser* parser = new xercesc::XercesDOMParser();
    parser->setLoadSchema(false);

    xercesc::ErrorHandler* errHandler = (xercesc::ErrorHandler*) new xercesc::HandlerBase();
    parser->setErrorHandler(errHandler);

    FragmentIsotopeApproximator* FIA = new FragmentIsotopeApproximator(argv[1], parser);

    //std::cout << FIA->calc_probability_spline(1,0,0,0,1,0,2000,1800,true) << std::endl;
    //std::cout << FIA->calc_probability_sulfur_corrected_averagine(1,0,1,0,2000,1800) << std::endl;

    test_tryptic_peptides(argv[2], FIA, atoi(argv[3]), atoi(argv[4]));

    delete parser;
    delete errHandler;
    delete FIA;
    xercesc::XMLPlatformUtils::Terminate();
    return 0;
}
