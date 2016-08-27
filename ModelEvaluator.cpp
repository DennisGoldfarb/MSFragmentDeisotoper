//
// Created by Dennis Goldfarb on 8/2/16.
//

#include <iostream>
#include <fstream>

#include <zlib.h>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <unordered_set>

#include "kseq.h"

#include "TensorSplineModel.h"
#include "FragmentIsotopeApproximator.h"
#include "Peptide.h"
#include "b_ion.h"
#include "y_ion.h"
#include "Histogram.h"

KSEQ_INIT(gzFile, gzread)

void test_peptide(std::string peptide_seq, FragmentIsotopeApproximator* FIA, Histogram* spline, Histogram* averagine_s, Histogram* averagine) {
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

    for (int precursor_isotope = 1; precursor_isotope < 5; ++precursor_isotope) {
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

            for (int fragment_isotope = 0; fragment_isotope <= precursor_isotope; ++fragment_isotope) {
                double abundance = std::pow(2, b_ion_isotope_abundances[fragment_isotope]);
                if (!isnan(abundance) && !isinf(abundance) && index > 0) {

                    float pm = p.calc_monoisotopic_mass(), fm = b_ions[index].calc_monoisotopic_mass();
                    int pi = precursor_isotope, fi = fragment_isotope;
                    int S = b_s[index], CS = p.get_composition()[elements::ELEMENTS::S]-S;

                    float prob = FIA->calc_probability_spline(S,CS,0,0,pi,fi,pm,fm);
                    if (prob > 0) {
                        spline->add_data(std::abs(prob - abundance));

                        if (std::abs(prob - abundance) > .5) {
                            std::cout << p.sequence << " " << index << " " << precursor_isotope << " " << fragment_isotope << " " << prob-abundance << std::endl;
                            FIA->calc_probability_spline(S,CS,0,0,pi,fi,pm,fm);
                        }

                        /*prob = FIA->calc_probability_sulfur_corrected_averagine(S, CS, pi, fi, pm, fm);

                        double diff = std::abs(prob - abundance);
                        averagine_s->add_data(diff);

                        prob = FIA->calc_probability_averagine(pi, fi, pm, fm);
                        averagine->add_data(std::abs(prob - abundance));*/
                    }
                }
                abundance = std::pow(2, y_ion_isotope_abundances[fragment_isotope]);
                if (!isnan(abundance) && !isinf(abundance) && index > 0) {

                    float pm = p.calc_monoisotopic_mass(), fm = y_ions[index].calc_monoisotopic_mass();
                    int pi = precursor_isotope, fi = fragment_isotope;
                    int S = y_s[index], CS = p.get_composition()[elements::ELEMENTS::S]-S;

                    float prob = FIA->calc_probability_spline(S,CS,0,0,pi,fi,pm,fm);
                    if (prob > 0) {
                        spline->add_data(std::abs(prob - abundance));

                        if (std::abs(prob - abundance) > .5) {
                            std::cout << p.sequence << " " << index << " " << precursor_isotope << " " << fragment_isotope << " " << prob-abundance << std::endl;
                            FIA->calc_probability_spline(S,CS,0,0,pi,fi,pm,fm);
                        }
                        /*prob = FIA->calc_probability_sulfur_corrected_averagine(S, CS, pi, fi, pm, fm);

                        double diff = std::abs(prob - abundance);
                        averagine_s->add_data(diff);

                        prob = FIA->calc_probability_averagine(pi, fi, pm, fm);
                        averagine->add_data(std::abs(prob - abundance));*/
                    }
                }
            }
        }
    }
}


void test_tryptic_peptides(const char* path_FASTA, FragmentIsotopeApproximator* FIA, int job_id, int num_jobs) {
    Histogram* spline = new Histogram("Spline Residual Distribution", "residual", "log2(count)");
    Histogram* averagine_s = new Histogram("Averagine S Residual Distribution", "residual", "log2(count)");
    Histogram* averagine = new Histogram("Averagine Distribution", "residual", "log2(count)");

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

                if (valid && peptide_seq.size() >= 4 && peptide_seq.size() <= 50) {
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
        if (i > offset && i < offset+num_test) {
            test_peptide(*itr, FIA, spline, averagine_s, averagine);

            if (i%10000 == 0) {
                std::cout << "Number of peptides processed: " << i << std::endl;
                spline->print_histogram();
                averagine_s->print_histogram();
                averagine->print_histogram();
            }
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


    for (auto c = char2count.begin(); c != char2count.end(); ++c) {
        std::cout << c->first << " " << c->second << std::endl;
    }
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

    test_tryptic_peptides(argv[2], FIA, atoi(argv[3]), atoi(argv[4]));

    delete parser;
    delete errHandler;
    delete FIA;
    xercesc::XMLPlatformUtils::Terminate();
    return 0;
}
