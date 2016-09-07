//
// Created by Dennis Goldfarb on 8/17/16.
//

#include "FragmentIsotopeApproximator.h"

FragmentIsotopeApproximator::FragmentIsotopeApproximator(char *infile, xercesc::XercesDOMParser *parser) {
    try {
        parser->parse(infile);

        xercesc::DOMElement* docRootNode;
        xercesc::DOMDocument* doc;
        xercesc::DOMNodeIterator * walker;
        doc = parser->getDocument();
        docRootNode = doc->getDocumentElement();

        // Create the node iterator, that will walk through each element.
        walker = doc->createNodeIterator(docRootNode,xercesc::DOMNodeFilter::SHOW_ELEMENT,NULL,true);

        xercesc::DOMNode * current_node = NULL;
        std::string thisNodeName;


        for (current_node = walker->nextNode(); current_node != 0; current_node = walker->nextNode()) {

            thisNodeName = xercesc::XMLString::transcode(current_node->getNodeName());

            if (thisNodeName == "model") {
                TensorSplineModel* tsm = new TensorSplineModel(current_node, walker);
                models.insert(std::make_pair(tsm->get_model_attributes(), tsm));
            }
        }

    }
    catch (const xercesc::XMLException& e) {
        char* message = xercesc::XMLString::transcode(e.getMessage());
        std::cout << "Exception message is: \n" << message << "\n";
        xercesc::XMLString::release(&message);
        throw;
    }
    catch (const xercesc::DOMException& e) {
        char* message = xercesc::XMLString::transcode(e.msg);
        std::cout << "Exception message is: \n" << message << "\n";
        xercesc::XMLString::release(&message);
        throw;
    }
    catch (const xercesc::SAXParseException& e) {
        char* message = xercesc::XMLString::transcode(e.getMessage());

        std::cout << "Exception message is: \n" << message << "\n";
        std::cout << e.getLineNumber() << std::endl;
        std::cout << e.getColumnNumber() << std::endl;
        xercesc::XMLString::release(&message);
        throw;
    }
    catch (...) {
        std::cout << "Unexpected Exception \n";
        throw;
    }
}

float FragmentIsotopeApproximator::calc_probability_spline(unsigned int num_sulfur, unsigned int num_comp_sulfur,
                                                    unsigned int num_selenium, unsigned int num_comp_selenium,
                                                    unsigned int precursor_isotope, unsigned int fragment_isotope,
                                                    float precursor_mass, float fragment_mass, bool verbose) {

    if (precursor_isotope == 0 && fragment_isotope == 0) return 1;

    ModelAttributes att(num_sulfur, num_comp_sulfur, num_selenium, num_comp_selenium, precursor_isotope, fragment_isotope);

    if (num_comp_sulfur > num_sulfur || num_comp_selenium > num_selenium) {
        std::swap(att.num_sulfur, att.num_comp_sulfur);
        std::swap(att.num_selenium, att.num_comp_selenium);
        att.fragment_isotope = precursor_isotope - fragment_isotope;
        fragment_mass = precursor_mass - fragment_mass;
    }

    TensorSplineModel* model = models[att];
    if (model == nullptr) {
        return -1;
    }
    return model->evaluate_model(precursor_mass, fragment_mass, verbose);
}

float FragmentIsotopeApproximator::calc_probability_averagine(unsigned int precursor_isotope, unsigned int fragment_isotope, float precursor_mass, float fragment_mass) {

    AveragineModel averagine_model;
    std::vector<std::vector<double>> precursor2fragment_probabilities;

    std::vector<unsigned int> precursor_composition = averagine_model.estimate_composition(precursor_mass);
    std::vector<unsigned int> fragment_composition = averagine_model.estimate_composition(fragment_mass);
    std::vector<unsigned int> complement_composition = averagine_model.estimate_composition(precursor_mass-fragment_mass);

    std::vector<double> frag_mz, frag_abundance, comp_frag_mz, comp_frag_abundance, precursor_mz, precursor_abundance;

    mercury::mercury(frag_mz, frag_abundance, fragment_composition, 1, 1e-30);
    mercury::mercury(comp_frag_mz, comp_frag_abundance, complement_composition, 1, 1e-30);
    mercury::mercury(precursor_mz, precursor_abundance, precursor_composition, 1, 1e-30);

    int comp_isotope = precursor_isotope - fragment_isotope;

    double probability = std::pow(2,std::log2(frag_abundance[fragment_isotope]) + std::log2(comp_frag_abundance[comp_isotope]) - std::log2(precursor_abundance[precursor_isotope]));

    return probability;
}

float FragmentIsotopeApproximator::calc_probability_sulfur_corrected_averagine(unsigned int num_sulfur,
                                                                               unsigned int num_comp_sulfur,
                                                                               unsigned int precursor_isotope,
                                                                               unsigned int fragment_isotope,
                                                                               float precursor_mass,
                                                                               float fragment_mass) {
    AveragineModel averagine_model;
    std::vector<std::vector<double>> precursor2fragment_probabilities;

    std::vector<unsigned int> precursor_composition = averagine_model.estimate_sulfur_corrected_composition(num_sulfur+num_comp_sulfur,precursor_mass);
    std::vector<unsigned int> fragment_composition = averagine_model.estimate_sulfur_corrected_composition(num_sulfur,fragment_mass);
    std::vector<unsigned int> complement_composition = averagine_model.estimate_sulfur_corrected_composition(num_comp_sulfur,precursor_mass-fragment_mass);

    std::vector<double> frag_mz, frag_abundance, comp_frag_mz, comp_frag_abundance, precursor_mz, precursor_abundance;

    mercury::mercury(frag_mz, frag_abundance, fragment_composition, 1, 1e-30);
    mercury::mercury(comp_frag_mz, comp_frag_abundance, complement_composition, 1, 1e-30);
    mercury::mercury(precursor_mz, precursor_abundance, precursor_composition, 1, 1e-30);

    int comp_isotope = precursor_isotope - fragment_isotope;

    double tot_prob = 0;
    for (int fi = 0; fi <= precursor_isotope; fi++) {
        int ci = precursor_isotope - fi;
        tot_prob += frag_abundance[fi] * comp_frag_abundance[ci];
    }


    double probability = std::pow(2,std::log2(frag_abundance[fragment_isotope]) + std::log2(comp_frag_abundance[comp_isotope]) - std::log2(tot_prob));//- std::log2(precursor_abundance[precursor_isotope]));

    return probability;
}

float FragmentIsotopeApproximator::calc_probability_sulfur_corrected_averagine(unsigned int num_sulfur,
                                                                               unsigned int num_comp_sulfur,
                                                                               unsigned int precursor_isotope,
                                                                               unsigned int fragment_isotope,
                                                                               float precursor_mass,
                                                                               float fragment_mass,
                                                                               float precursor_probability) {
    AveragineModel averagine_model;
    std::vector<std::vector<double>> precursor2fragment_probabilities;

    std::vector<unsigned int> precursor_composition = averagine_model.estimate_sulfur_corrected_composition(num_sulfur+num_comp_sulfur,precursor_mass);
    std::vector<unsigned int> fragment_composition = averagine_model.estimate_sulfur_corrected_composition(num_sulfur,fragment_mass);
    std::vector<unsigned int> complement_composition = averagine_model.estimate_sulfur_corrected_composition(num_comp_sulfur,precursor_mass-fragment_mass);

    std::vector<double> frag_mz, frag_abundance, comp_frag_mz, comp_frag_abundance, precursor_mz, precursor_abundance;

    mercury::mercury(frag_mz, frag_abundance, fragment_composition, 1, 1e-30);
    mercury::mercury(comp_frag_mz, comp_frag_abundance, complement_composition, 1, 1e-30);
    mercury::mercury(precursor_mz, precursor_abundance, precursor_composition, 1, 1e-30);

    int comp_isotope = precursor_isotope - fragment_isotope;

    double probability = std::pow(2,std::log2(frag_abundance[fragment_isotope]) + std::log2(comp_frag_abundance[comp_isotope]) - std::log2(precursor_probability));

    return probability;
}

float FragmentIsotopeApproximator::calc_expected_probability(std::vector<unsigned int> precursor_composition,
                                                             std::vector<unsigned int> fragment_composition,
                                                             unsigned int precursor_isotope, unsigned int fragment_isotope) {

    std::vector<double> frag_mz, frag_abundance, comp_frag_mz, comp_frag_abundance, precursor_mz, precursor_abundance;

    std::vector<unsigned int> complement_composition;
    for (int i = 0; i < precursor_composition.size(); ++i) {
        complement_composition.push_back(precursor_composition[i] - fragment_composition[i]);
    }

    mercury::mercury(frag_mz, frag_abundance, fragment_composition, 1, 1e-30);
    mercury::mercury(comp_frag_mz, comp_frag_abundance, complement_composition, 1, 1e-30);
    mercury::mercury(precursor_mz, precursor_abundance, precursor_composition, 1, 1e-30);

    int comp_isotope = precursor_isotope - fragment_isotope;

    double probability = std::pow(2,std::log2(frag_abundance[fragment_isotope]) + std::log2(comp_frag_abundance[comp_isotope]) - std::log2(precursor_abundance[precursor_isotope]));

    return probability;
}

float FragmentIsotopeApproximator::get_closest_spline_probability(float actual_prob, float precursor_mass, float fragment_mass,
                                                                  bool verbose) {
    TensorSplineModel* bestModel;
    float bestDiff = 1;

    for (auto model = models.begin(); model != models.end(); ++model) {
        if (model->second != nullptr) {
            float prob = model->second->evaluate_model(precursor_mass, fragment_mass, false);
            if (prob != -1) {
                float diff = std::abs(prob - actual_prob);
                if (diff < bestDiff) {
                    bestDiff = diff;
                    bestModel = model->second;
                }

            }
            if (model->second->num_sulfur > model->second->num_comp_sulfur) {
                ModelAttributes att(model->second->num_comp_sulfur, model->second->num_sulfur,
                                    model->second->num_comp_selenium, model->second->num_selenium,
                                    model->second->precursor_isotope, model->second->precursor_isotope-model->second->fragment_isotope);

                float prob = model->second->evaluate_model(precursor_mass, precursor_mass - fragment_mass, false);
                if (prob != -1) {
                    float diff = std::abs(prob - actual_prob);
                    if (diff < bestDiff) {
                        bestDiff = diff;
                        bestModel = model->second;
                    }

                }
            }
        }
    }

    if (verbose) {
        std::cout << bestModel->num_sulfur << " " << bestModel->num_comp_sulfur << " " <<
        bestModel->precursor_isotope << " " << bestModel->fragment_isotope << std::endl;
    }

    return bestDiff;
}
